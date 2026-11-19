use std::collections::{HashMap, HashSet};

use wireframes_core::exact::Expansion;
use wireframes_core::predicates::{Sign, collinear3, kept_axes, orient2d_axis, orient3d};
use wireframes_core::topology::{Radial, radial_pairs};
use wireframes_core::{Error, Mesh, Result, Table, units};
use wireframes_ops::{RepairOptions, RepairReport, repair};

use crate::arrangement::{self, Arrangement};
use crate::classify::{self, Status};
use crate::hp;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Op {
    Union,
    Intersection,
    Difference,
    SymmetricDifference,
}

/// `repair`: repair each input with `RepairOptions::default()` first.
#[derive(Clone, Debug, Default)]
pub struct BooleanOptions {
    pub repair: bool,
}

#[derive(Clone, Debug)]
pub struct BooleanResult {
    pub mesh: Mesh,
    /// Per output face: (0 for a, 1 for b, input face index after repair).
    pub face_source: Vec<(u8, u32)>,
    pub repairs: [Option<RepairReport>; 2],
}

/// Repairs `mesh` when asked, then refuses it unless it is a valid solid or
/// has no faces (the empty solid a boolean returns).
pub(crate) fn prepare(
    mesh: &Mesh,
    context: &str,
    fix: bool,
) -> Result<(Mesh, Option<RepairReport>)> {
    let (m, r) = if fix {
        let (m, r) = repair(mesh, &RepairOptions::default())?;
        (m, Some(r))
    } else {
        (mesh.clone(), None)
    };
    let report = m.validate(0.0)?;
    if m.n_faces() > 0 && !report.summary.is_solid() {
        return Err(Error::InvalidMesh {
            context: context.to_string(),
            report: Box::new(report),
        });
    }
    Ok((m, r))
}

/// Sub-triangle ids with a flip flag.
pub(crate) type Piece = Vec<(u32, bool)>;

/// One output shell: sub-triangle ids with a flip flag and, per corner, the
/// umbrella of faces about that corner's point it belongs to. Shells, and
/// umbrellas within a shell, get separate vertex copies, so shells may touch
/// each other and a shell may touch itself without making the output
/// non-manifold.
pub(crate) type Shell = Vec<(u32, bool, [u32; 3])>;

pub(crate) struct Assembled {
    pub vertices: Vec<[f64; 3]>,
    pub triangles: Vec<[u32; 3]>,
    /// Row of the concatenated input vertex tables, per output vertex.
    pub vertex_rows: Vec<Option<u32>>,
    /// (mesh, face) per output triangle.
    pub face_source: Vec<(u8, u32)>,
    /// Outward normal of the exact plane each triangle lies in, rounded to
    /// `f64`; empty when unknown.
    pub normal: Vec<[f64; 3]>,
}

/// Sign of `n . ((b - a) x (c - a))`, exact on the given `f64` values.
fn facing(c: [[f64; 3]; 3], n: [f64; 3]) -> Sign {
    let u = [0, 1, 2].map(|k| Expansion::diff_f64(c[1][k], c[0][k]));
    let v = [0, 1, 2].map(|k| Expansion::diff_f64(c[2][k], c[0][k]));
    (0..3)
        .fold(Expansion::zero(), |acc, k| {
            let (i, j) = kept_axes(k);
            acc.sum(&u[i].product(&v[j]).diff(&u[j].product(&v[i])).scale(n[k]))
        })
        .sign()
}

pub(crate) fn assemble(arr: &Arrangement, pieces: &[Shell]) -> Assembled {
    let mut out = Assembled {
        vertices: vec![],
        triangles: vec![],
        vertex_rows: vec![],
        face_source: vec![],
        normal: vec![],
    };
    for piece in pieces {
        let mut shell = Assembled {
            vertices: vec![],
            triangles: vec![],
            vertex_rows: vec![],
            face_source: vec![],
            normal: vec![],
        };
        let mut map: HashMap<(u32, u32), u32> = HashMap::new();
        let mut point = vec![];
        for &(s, flip, umbrella) in piece {
            let sub = &arr.subs[s as usize];
            let v = [0, 1, 2].map(|c| {
                let p = sub.v[c];
                *map.entry((p, umbrella[c])).or_insert_with(|| {
                    shell.vertices.push(arr.pts[p as usize].to_f64());
                    shell.vertex_rows.push(arr.source_row(p));
                    point.push(p);
                    shell.vertices.len() as u32 - 1
                })
            });
            shell
                .triangles
                .push(if flip { [v[0], v[2], v[1]] } else { v });
            let f = &arr.faces[sub.face as usize];
            shell.face_source.push((f.mesh, f.index));
            let n = arr.planes[sub.face as usize]
                .normal_exp()
                .clone()
                .map(|x| x.estimate());
            shell.normal.push(if flip { n.map(|x| -x) } else { n });
        }
        let (shell, _) = collapse_rounded(shell, &point);
        let base = out.vertices.len() as u32;
        out.vertices.extend(shell.vertices);
        out.vertex_rows.extend(shell.vertex_rows);
        out.triangles
            .extend(shell.triangles.iter().map(|t| t.map(|i| i + base)));
        out.face_source.extend(shell.face_source);
        out.normal.extend(shell.normal);
    }
    out
}

/// One shell after rounding its exact vertices to `f64`. Distinct exact points
/// can round to the same coordinates, which makes faces with a repeated or
/// collinear corner. Vertices of the shell with equal coordinates become one
/// (an edge between them has zero length), faces that repeat a vertex are
/// dropped, opposite faces on the same three vertices cancel, and a face whose
/// corners are collinear is removed by splitting its neighbour across its
/// longest edge at its middle corner. A sliver whose middle corner rounds to
/// strictly inside that coplanar neighbour (a fold) is removed the same way,
/// and so, when `normal` is known, is a face that rounding turned away from
/// its exact plane (`n . ((b - a) x (c - a)) <= 0`); all three tests are exact
/// on the rounded points. Where that split would duplicate an edge, a flat or
/// turned face loses its shortest edge instead: one end merges into the other
/// when the link condition holds and no moved face turns. When `normal` is
/// known, faces that `validate` still finds degenerate, self-intersecting or
/// overlapping, with their edge neighbours, get the split too where the
/// neighbour across the longest edge comes from the same input face and both
/// new faces face its plane; such rounds repeat while the problem count
/// drops, and a round that does not lower it is undone. Each step keeps every
/// edge paired. `point` is the exact point of each vertex; copies of one exact
/// point, where the shell touches itself, stay apart and never merge. Also
/// returns the exact point of each kept vertex.
pub(crate) fn collapse_rounded(shell: Assembled, point: &[u32]) -> (Assembled, Vec<u32>) {
    let n = shell.vertices.len();
    let key = |i: usize| shell.vertices[i].map(|x| (x + 0.0).to_bits());
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by_key(|&i| (key(i), shell.vertex_rows[i].is_none(), i));
    let mut rep: Vec<u32> = (0..n as u32).collect();
    let mut start = 0;
    for w in 1..=n {
        if w < n && key(order[w]) == key(order[start]) {
            continue;
        }
        let mut ps: Vec<u32> = order[start..w].iter().map(|&i| point[i]).collect();
        ps.sort_unstable();
        if ps.windows(2).all(|p| p[0] != p[1]) {
            for &i in &order[start + 1..w] {
                rep[i] = order[start] as u32;
            }
        }
        start = w;
    }
    let has_normals = !shell.normal.is_empty();
    let mut tris: Vec<[u32; 3]> = vec![];
    let mut src: Vec<(u8, u32)> = vec![];
    let mut nrm: Vec<[f64; 3]> = vec![];
    for (f, (t, &s)) in shell.triangles.iter().zip(&shell.face_source).enumerate() {
        let t = t.map(|i| rep[i as usize]);
        if t[0] != t[1] && t[1] != t[2] && t[0] != t[2] {
            tris.push(t);
            src.push(s);
            if has_normals {
                nrm.push(shell.normal[f]);
            }
        }
    }
    let mut alive = vec![true; tris.len()];
    let mut by_set: HashMap<[u32; 3], Vec<usize>> = HashMap::new();
    for (f, t) in tris.iter().enumerate() {
        let mut s = *t;
        s.sort_unstable();
        by_set.entry(s).or_default().push(f);
    }
    let mut sets: Vec<_> = by_set.into_values().filter(|fs| fs.len() > 1).collect();
    sets.sort_unstable();
    let parity = |t: [u32; 3]| {
        let k = (0..3).min_by_key(|&k| t[k]).expect("three corners");
        t[(k + 1) % 3] < t[(k + 2) % 3]
    };
    for fs in sets {
        let (mut up, mut down): (Vec<usize>, Vec<usize>) =
            fs.into_iter().partition(|&f| parity(tris[f]));
        while let (Some(a), Some(b)) = (up.pop(), down.pop()) {
            alive[a] = false;
            alive[b] = false;
        }
    }
    let v = &shell.vertices;
    let flat = |t: [u32; 3]| collinear3(v[t[0] as usize], v[t[1] as usize], v[t[2] as usize]);
    let corners = |t: [u32; 3]| t.map(|i| v[i as usize]);
    let turned = |t: [u32; 3], n: [f64; 3]| facing(corners(t), n) != Sign::Positive;
    let mut edge: HashMap<(u32, u32), usize> = HashMap::new();
    let mut inc: Vec<Vec<usize>> = vec![vec![]; n];
    for (f, t) in tris.iter().enumerate().filter(|(f, _)| alive[*f]) {
        for k in 0..3 {
            edge.insert((t[k], t[(k + 1) % 3]), f);
            inc[t[k] as usize].push(f);
        }
    }
    let inside = |g: [u32; 3], m: u32| {
        let c = g.map(|i| v[i as usize]);
        let at = v[m as usize];
        let Some(k) = (0..3).find(|&k| orient2d_axis(c[0], c[1], c[2], k) != 0.0) else {
            return false;
        };
        let s = Sign::of(orient2d_axis(c[0], c[1], c[2], k));
        orient3d(c[0], c[1], c[2], at) == 0.0
            && (0..3).all(|e| Sign::of(orient2d_axis(c[e], c[(e + 1) % 3], at, k)) == s)
    };
    let mut work: Vec<usize> = (0..tris.len()).filter(|&f| alive[f]).rev().collect();
    // A split flips edge `pq` to `md`. An edge a flip made is never flipped
    // again, so every flip removes an edge of the shell as it came in and the
    // loop ends; without this, flips on rounded points can cycle forever.
    let mut made: HashSet<(u32, u32)> = HashSet::new();
    let mut forced: HashSet<usize> = HashSet::new();
    let mut last = usize::MAX;
    let mut saved = None;
    loop {
        while let Some(f) = work.pop() {
            if !alive[f] {
                continue;
            }
            let t = tris[f];
            let axis = (0..3)
                .max_by(|&a, &b| {
                    let spread = |k: usize| {
                        let c = t.map(|i| v[i as usize][k]);
                        c.iter().copied().fold(f64::NEG_INFINITY, f64::max)
                            - c.iter().copied().fold(f64::INFINITY, f64::min)
                    };
                    spread(a).total_cmp(&spread(b))
                })
                .expect("three axes");
            let mut by: [usize; 3] = [0, 1, 2];
            by.sort_by(|&a, &b| v[t[a] as usize][axis].total_cmp(&v[t[b] as usize][axis]));
            let i = by[1];
            let (m, p, q) = (t[i], t[(i + 1) % 3], t[(i + 2) % 3]);
            let Some(&g) = edge.get(&(q, p)) else {
                continue;
            };
            if g == f || !alive[g] {
                continue;
            }
            let u = tris[g];
            let d = u[(0..3)
                .find(|&k| u[k] != p && u[k] != q)
                .expect("third corner")];
            let broken = flat(t) || (has_normals && turned(t, nrm[f]));
            let refit = !broken
                && forced.contains(&f)
                && src[f] == src[g]
                && d != m
                && !turned([q, m, d], nrm[g])
                && !turned([m, p, d], nrm[g]);
            if !broken && !refit && !inside(u, m) {
                continue;
            }
            if made.contains(&ekey(p, q))
                || edge.contains_key(&(m, d))
                || edge.contains_key(&(d, m))
            {
                if broken && has_normals {
                    let len2 = |(a, b): (u32, u32)| {
                        let (a, b) = (v[a as usize], v[b as usize]);
                        (0..3).map(|k| (a[k] - b[k]) * (a[k] - b[k])).sum::<f64>()
                    };
                    let short = [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])]
                        .into_iter()
                        .min_by(|x, y| len2(*x).total_cmp(&len2(*y)))
                        .expect("three edges");
                    let mut tries = [(short.0, short.1), (short.1, short.0)];
                    if shell.vertex_rows[short.0 as usize].is_some() {
                        tries.reverse();
                    }
                    for (gone, keep) in tries {
                        let Some(nt) = edge_collapse(
                            &tris, &alive, &edge, &inc, &nrm, gone, keep, point, &turned,
                        ) else {
                            continue;
                        };
                        for &(h, _) in &nt {
                            for k in 0..3 {
                                let e = (tris[h][k], tris[h][(k + 1) % 3]);
                                if edge.get(&e) == Some(&h) {
                                    edge.remove(&e);
                                }
                            }
                        }
                        for (h, w) in nt {
                            match w {
                                None => alive[h] = false,
                                Some(w) => {
                                    tris[h] = w;
                                    for k in 0..3 {
                                        edge.insert((w[k], w[(k + 1) % 3]), h);
                                    }
                                    inc[keep as usize].push(h);
                                    work.push(h);
                                }
                            }
                        }
                        break;
                    }
                }
                continue;
            }
            made.insert(ekey(m, d));
            for (h, w) in [(f, t), (g, u)] {
                alive[h] = false;
                for k in 0..3 {
                    if edge.get(&(w[k], w[(k + 1) % 3])) == Some(&h) {
                        edge.remove(&(w[k], w[(k + 1) % 3]));
                    }
                }
            }
            for nt in [[q, m, d], [m, p, d]] {
                if d == m {
                    continue;
                }
                let h = tris.len();
                tris.push(nt);
                src.push(src[g]);
                if has_normals {
                    nrm.push(nrm[g]);
                }
                alive.push(true);
                for k in 0..3 {
                    edge.insert((nt[k], nt[(k + 1) % 3]), h);
                    inc[nt[k] as usize].push(h);
                }
                work.push(h);
            }
        }
        if !has_normals {
            break;
        }
        let (count, bad) = rounding_problems(v, &tris, &alive);
        if count >= last {
            if let Some(state) = saved.take() {
                (tris, alive, src, nrm) = state;
            }
            break;
        }
        if count == 0 {
            break;
        }
        last = count;
        saved = Some((tris.clone(), alive.clone(), src.clone(), nrm.clone()));
        forced = bad
            .iter()
            .flat_map(|&f| {
                let t = tris[f];
                let edge = &edge;
                std::iter::once(f)
                    .chain((0..3).filter_map(move |k| edge.get(&(t[(k + 1) % 3], t[k])).copied()))
            })
            .collect();
        work.extend(forced.iter().copied());
        work.sort_unstable_by(|a, b| b.cmp(a));
    }
    let mut out = Assembled {
        vertices: vec![],
        triangles: vec![],
        vertex_rows: vec![],
        face_source: vec![],
        normal: vec![],
    };
    let mut new_id: Vec<Option<u32>> = vec![None; n];
    let mut kept_point = vec![];
    for (f, t) in tris.iter().enumerate().filter(|(f, _)| alive[*f]) {
        if has_normals {
            out.normal.push(nrm[f]);
        }
        out.triangles.push(t.map(|i| {
            *new_id[i as usize].get_or_insert_with(|| {
                out.vertices.push(shell.vertices[i as usize]);
                out.vertex_rows.push(shell.vertex_rows[i as usize]);
                kept_point.push(point[i as usize]);
                out.vertices.len() as u32 - 1
            })
        }));
        out.face_source.push(src[f]);
    }
    (out, kept_point)
}

/// Collapse of vertex `gone` onto `keep` along their edge, as the faces it
/// changes: `None` for the two faces on the edge, the moved face otherwise.
/// Refused when both are copies of one exact point, when it breaks the link
/// condition (the shell would stop being manifold), or when a moved face would
/// turn away from its exact plane's normal or become flat (exact on the
/// rounded points).
#[allow(clippy::too_many_arguments)]
fn edge_collapse(
    tris: &[[u32; 3]],
    alive: &[bool],
    edge: &HashMap<(u32, u32), usize>,
    inc: &[Vec<usize>],
    nrm: &[[f64; 3]],
    gone: u32,
    keep: u32,
    point: &[u32],
    turned: &dyn Fn([u32; 3], [f64; 3]) -> bool,
) -> Option<Vec<(usize, Option<[u32; 3]>)>> {
    if point[gone as usize] == point[keep as usize] {
        return None;
    }
    let (&f1, &f2) = (edge.get(&(keep, gone))?, edge.get(&(gone, keep))?);
    let third = |f: usize| {
        tris[f]
            .into_iter()
            .find(|&x| x != keep && x != gone)
            .expect("third corner")
    };
    let (o1, o2) = (third(f1), third(f2));
    if f1 == f2 || o1 == o2 {
        return None;
    }
    let around = |x: u32| {
        let mut fs: Vec<usize> = inc[x as usize]
            .iter()
            .copied()
            .filter(|&h| alive[h] && tris[h].contains(&x))
            .collect();
        fs.sort_unstable();
        fs.dedup();
        fs
    };
    let ring = |x: u32| -> HashSet<u32> {
        around(x)
            .iter()
            .flat_map(|&h| tris[h])
            .filter(|&y| y != x)
            .collect()
    };
    let (rk, rg) = (ring(keep), ring(gone));
    if rk.intersection(&rg).any(|&y| y != o1 && y != o2) {
        return None;
    }
    let mut out = vec![(f1, None), (f2, None)];
    for h in around(gone) {
        if h == f1 || h == f2 {
            continue;
        }
        let w = tris[h].map(|x| if x == gone { keep } else { x });
        if turned(w, nrm[h]) {
            return None;
        }
        out.push((h, Some(w)));
    }
    Some(out)
}

/// Problems `validate` finds in the live faces of a shell, counting every kind
/// but unused and duplicate vertices, and the faces of degenerate,
/// self-intersecting and coplanar-overlap problems, in increasing order.
fn rounding_problems(v: &[[f64; 3]], tris: &[[u32; 3]], alive: &[bool]) -> (usize, Vec<usize>) {
    use wireframes_core::MeshProblemKind as K;
    let ids: Vec<usize> = (0..tris.len()).filter(|&f| alive[f]).collect();
    let Ok(report) =
        Mesh::new(v.to_vec(), ids.iter().map(|&f| tris[f]).collect()).and_then(|m| m.validate(0.0))
    else {
        return (usize::MAX, vec![]);
    };
    let mut bad = vec![];
    let mut count = 0;
    for p in &report.problems {
        if matches!(p.kind, K::UnusedVertex | K::DuplicateVertex) {
            continue;
        }
        count += 1;
        if matches!(
            p.kind,
            K::DegenerateFace | K::SelfIntersection | K::CoplanarOverlap
        ) {
            bad.extend(
                [p.face, p.other]
                    .into_iter()
                    .flatten()
                    .map(|f| ids[f as usize]),
            );
        }
    }
    bad.sort_unstable();
    bad.dedup();
    (count, bad)
}

fn ekey(a: u32, b: u32) -> (u32, u32) {
    if a < b { (a, b) } else { (b, a) }
}

/// The shells of a closed piece, each ordered as in `piece` and the shells
/// ordered by first sub-triangle. Sub-triangles join across edges with two of
/// them, and across edges with more in the pairs that bound wedges of solid,
/// found by exact radial order about the edge. Corners of joined
/// sub-triangles at the ends of the edge share an umbrella.
pub(crate) fn shells(arr: &Arrangement, piece: &Piece) -> Result<Vec<Shell>> {
    let tri = |k: usize| {
        let (s, flip) = piece[k];
        let v = arr.subs[s as usize].v;
        if flip { [v[0], v[2], v[1]] } else { v }
    };
    let mut at: HashMap<(u32, u32), Vec<u32>> = HashMap::new();
    for k in 0..piece.len() {
        let t = tri(k);
        for e in 0..3 {
            at.entry(ekey(t[e], t[(e + 1) % 3]))
                .or_default()
                .push(k as u32);
        }
    }
    let mut parent: Vec<u32> = (0..piece.len() as u32).collect();
    fn root(parent: &mut [u32], mut i: u32) -> u32 {
        while parent[i as usize] != i {
            parent[i as usize] = parent[parent[i as usize] as usize];
            i = parent[i as usize];
        }
        i
    }
    let join = |parent: &mut [u32], a: u32, b: u32| {
        let (a, b) = (root(parent, a), root(parent, b));
        parent[a.max(b) as usize] = a.min(b);
    };
    let mut corners: Vec<u32> = (0..3 * piece.len() as u32).collect();
    let corner = |k: u32, p: u32| {
        let v = arr.subs[piece[k as usize].0 as usize].v;
        3 * k + v.iter().position(|&x| x == p).expect("corner on the edge") as u32
    };
    let link = |parent: &mut [u32], corners: &mut [u32], f: u32, g: u32, p: u32, q: u32| {
        join(parent, f, g);
        join(corners, corner(f, p), corner(g, p));
        join(corners, corner(f, q), corner(g, q));
    };
    let mut keys: Vec<(u32, u32)> = at.keys().copied().collect();
    keys.sort_unstable();
    for (p, q) in keys {
        let ks = &at[&(p, q)];
        if ks.len() == 2 {
            link(&mut parent, &mut corners, ks[0], ks[1], p, q);
            continue;
        }
        if ks.len() < 2 {
            continue;
        }
        let d = hp::direction(&arr.pts[p as usize], &arr.pts[q as usize]);
        let forward: Vec<bool> = ks
            .iter()
            .map(|&k| {
                let t = tri(k as usize);
                (0..3).any(|e| t[e] == p && t[(e + 1) % 3] == q)
            })
            .collect();
        let m: Vec<[wireframes_core::exact::Expansion; 3]> = ks
            .iter()
            .zip(&forward)
            .map(|(&k, &fwd)| {
                let (s, flip) = piece[k as usize];
                let n = arr.planes[arr.subs[s as usize].face as usize].normal_exp();
                if flip == fwd {
                    n.clone().map(|x| x.neg())
                } else {
                    n.clone()
                }
            })
            .collect();
        let turn = |a: usize, b: usize| hp::det3(&d, &m[a], &m[b]);
        let same = |a: usize, b: usize| hp::dot3(&m[a], &m[b]) == Sign::Positive;
        match radial_pairs(&forward, turn, same) {
            Radial::Paired(pairs) => {
                for [a, b] in pairs {
                    link(&mut parent, &mut corners, ks[a], ks[b], p, q);
                }
            }
            _ => {
                return Err(Error::kernel(
                    "the faces around an edge of the result do not bound wedges of solid",
                ));
            }
        }
    }
    let mut index: HashMap<u32, usize> = HashMap::new();
    let mut out: Vec<Shell> = vec![];
    for (k, &(s, flip)) in piece.iter().enumerate() {
        let r = root(&mut parent, k as u32);
        let next = out.len();
        let i = *index.entry(r).or_insert(next);
        if i == out.len() {
            out.push(vec![]);
        }
        let umbrella = [0, 1, 2].map(|c| root(&mut corners, 3 * k as u32 + c));
        out[i].push((s, flip, umbrella));
    }
    Ok(out)
}

/// Every edge used once in each direction and every vertex with a single fan.
pub(crate) fn is_closed_manifold(triangles: &[[u32; 3]]) -> bool {
    let mut dir: HashMap<(u32, u32), u32> = HashMap::new();
    for t in triangles {
        for k in 0..3 {
            *dir.entry((t[k], t[(k + 1) % 3])).or_default() += 1;
        }
    }
    if dir
        .iter()
        .any(|(&(a, b), &n)| n != 1 || dir.get(&(b, a)) != Some(&1))
    {
        return false;
    }
    let mut link: HashMap<u32, Vec<(u32, u32)>> = HashMap::new();
    for t in triangles {
        for k in 0..3 {
            link.entry(t[k])
                .or_default()
                .push((t[(k + 1) % 3], t[(k + 2) % 3]));
        }
    }
    link.values().all(|arcs| {
        let next: HashMap<u32, u32> = arcs.iter().copied().collect();
        if next.len() != arcs.len() {
            return false;
        }
        let start = arcs[0].0;
        let mut c = arcs[0].1;
        let mut n = 1;
        while c != start {
            match next.get(&c) {
                Some(&d) => c = d,
                None => return false,
            }
            n += 1;
            if n > arcs.len() {
                return false;
            }
        }
        n == arcs.len()
    })
}

pub(crate) fn build_mesh(
    asm: Assembled,
    vertex_tables: (&Table, &Table),
    face_tables: (&Table, &Table),
    labels: (Option<String>, Option<String>),
) -> Result<Mesh> {
    let vt = Table::concat(vertex_tables.0, vertex_tables.1)?;
    let ft = Table::concat(face_tables.0, face_tables.1)?;
    let n0 = face_tables.0.rows() as u32;
    let rows: Vec<u32> = asm
        .face_source
        .iter()
        .map(|&(m, f)| if m == 0 { f } else { n0 + f })
        .collect();
    let mesh = if asm.triangles.is_empty() {
        Mesh::empty().with_attributes(vt.take(&[]), ft.take(&[]))?
    } else {
        Mesh::new(asm.vertices, asm.triangles)?
            .with_attributes(vt.take_or_null(&asm.vertex_rows), ft.take(&rows))?
    };
    Ok(mesh.with_crs(labels.0).with_length_unit(labels.1))
}

/// `mesh` when `validate` passes it as a solid. Exact results of solids are
/// solids; rounding their new vertices to `f64` can fold slivers thinner than
/// the float spacing into self-intersections or flatten a face, and such a
/// result is an error that carries the report and the mesh, never a return.
/// Whether a split of a closed mesh that is not a solid is still closed and
/// manifold: its edges pair into wedges of solid and its vertices are manifold.
pub(crate) fn is_closed_result(mesh: &Mesh) -> Result<bool> {
    if mesh.is_empty() || is_closed_manifold(mesh.triangles()) {
        return Ok(true);
    }
    let s = mesh.validate(0.0)?.summary;
    Ok(s.is_closed
        && s.boundary_edges == 0
        && s.non_manifold_edges == 0
        && s.non_manifold_vertices == 0
        && s.inconsistent_edges == 0)
}

pub(crate) fn solid_result(mesh: Mesh, what: &str) -> Result<Mesh> {
    if mesh.is_empty() || mesh.validate(0.0)?.summary.is_solid() {
        return Ok(mesh);
    }
    Err(Error::kernel_result(
        format!(
            "rounding {what} to f64 left a mesh that is not a solid (snap rounding is not \
             implemented)"
        ),
        mesh,
    )?)
}

fn keep_rule(op: Op, mesh: usize, st: Status) -> Option<bool> {
    use Status::*;
    match (op, mesh, st) {
        (Op::Union, _, Out) | (Op::Union, 0, OnSame) => Some(false),
        (Op::Intersection, _, In) | (Op::Intersection, 0, OnSame) => Some(false),
        (Op::Difference, 0, Out) | (Op::Difference, 0, OnOpp) => Some(false),
        (Op::Difference, 1, In) => Some(true),
        (Op::SymmetricDifference, _, Out) => Some(false),
        (Op::SymmetricDifference, _, In) => Some(true),
        _ => None,
    }
}

pub fn boolean(a: &Mesh, b: &Mesh, op: Op, opts: &BooleanOptions) -> Result<BooleanResult> {
    let labels = units::check_same("b", (&a.crs, &a.length_unit), (&b.crs, &b.length_unit))?;
    let (a, ra) = prepare(a, "a", opts.repair)?;
    let (b, rb) = prepare(b, "b", opts.repair)?;
    let arr = arrangement::build(&a, &b)?;
    let p = classify::patches(&arr);
    let st = classify::solid_status(&arr, &p)?;
    let kept: Piece = arr
        .subs
        .iter()
        .enumerate()
        .filter_map(|(i, s)| {
            keep_rule(op, classify::mesh_of(&arr, s), st[i]).map(|flip| (i as u32, flip))
        })
        .collect();
    let asm = assemble(&arr, &shells(&arr, &kept)?);
    let face_source = asm.face_source.clone();
    let mesh = build_mesh(
        asm,
        (a.vertex_attributes(), b.vertex_attributes()),
        (a.face_attributes(), b.face_attributes()),
        labels,
    )?;
    let mesh = solid_result(mesh, "the exact result")?;
    Ok(BooleanResult {
        mesh,
        face_source,
        repairs: [ra, rb],
    })
}

pub fn union(a: &Mesh, b: &Mesh, opts: &BooleanOptions) -> Result<BooleanResult> {
    boolean(a, b, Op::Union, opts)
}

pub fn intersection(a: &Mesh, b: &Mesh, opts: &BooleanOptions) -> Result<BooleanResult> {
    boolean(a, b, Op::Intersection, opts)
}

pub fn difference(a: &Mesh, b: &Mesh, opts: &BooleanOptions) -> Result<BooleanResult> {
    boolean(a, b, Op::Difference, opts)
}

/// The parts of `a` outside `b` and of `b` outside `a`. Where they meet along
/// the intersection curve they are separate shells in contact.
pub fn symmetric_difference(a: &Mesh, b: &Mesh, opts: &BooleanOptions) -> Result<BooleanResult> {
    boolean(a, b, Op::SymmetricDifference, opts)
}

#[cfg(test)]
mod tests {
    use super::*;

    const A: [f64; 3] = [0.0, 0.0, 0.0];
    const B: [f64; 3] = [2.0, 0.0, 0.0];
    const C: [f64; 3] = [0.0, 2.0, 0.0];
    const D: [f64; 3] = [0.0, 0.0, 2.0];

    fn shell(vertices: Vec<[f64; 3]>, triangles: Vec<[u32; 3]>) -> Assembled {
        let n = vertices.len();
        Assembled {
            face_source: (0..triangles.len() as u32).map(|f| (0, f)).collect(),
            vertices,
            triangles,
            vertex_rows: (0..n as u32).map(|i| (i < 4).then_some(i)).collect(),
            normal: vec![],
        }
    }

    fn check(out: &Assembled) {
        let m = Mesh::new(out.vertices.clone(), out.triangles.clone()).unwrap();
        assert!(is_closed_manifold(&out.triangles));
        assert!(m.validate(0.0).unwrap().summary.is_solid());
        assert_eq!(out.face_source.len(), out.triangles.len());
        assert_eq!(out.vertex_rows.len(), out.vertices.len());
    }

    #[test]
    fn vertices_that_round_together_become_one() {
        let s = shell(
            vec![A, B, C, D, D],
            vec![
                [0, 2, 1],
                [0, 1, 3],
                [0, 3, 2],
                [1, 2, 4],
                [1, 4, 3],
                [2, 3, 4],
            ],
        );
        let (out, _) = collapse_rounded(s, &[0, 1, 2, 3, 4]);
        check(&out);
        assert_eq!(out.triangles.len(), 4);
        assert_eq!(out.vertex_rows, [Some(0), Some(2), Some(1), Some(3)]);
    }

    #[test]
    fn a_sliver_folded_onto_its_coplanar_neighbour_is_split_away() {
        let m = [1.0, 0.01, 0.0];
        let s = shell(
            vec![A, B, C, D, m],
            vec![
                [0, 2, 1],
                [0, 1, 4],
                [0, 4, 3],
                [4, 1, 3],
                [0, 3, 2],
                [1, 2, 3],
            ],
        );
        let before = Mesh::new(s.vertices.clone(), s.triangles.clone()).unwrap();
        assert!(!before.validate(0.0).unwrap().summary.is_solid());
        let (out, _) = collapse_rounded(s, &[0, 1, 2, 3, 4]);
        check(&out);
        assert_eq!(out.triangles.len(), 6);
    }

    #[test]
    fn a_face_with_collinear_rounded_corners_leaves_a_result_open() {
        let v = vec![
            [60.0, 16.0, 0.0],
            [60.72727272741031, 15.41818181797862, 0.0],
            [60.0, 15.200000000186265, 0.0],
            [62.66666666662786, 16.0, 0.0],
            [62.66666666732635, 16.0, 0.0],
            [62.66666666709352, 16.0, 3.7255176721373573e-10],
            [60.0, 16.0, 4.705882354479968],
            [60.0, 15.362549800425768, 0.9561752990377386],
        ];
        let t = vec![
            [0, 1, 2],
            [0, 3, 1],
            [0, 4, 3],
            [4, 0, 5],
            [5, 0, 6],
            [7, 0, 2],
            [7, 6, 0],
            [4, 5, 3],
            [3, 5, 1],
            [5, 6, 7],
            [2, 1, 7],
            [7, 1, 5],
        ];
        let m = Mesh::new(v, t).unwrap();
        assert!(is_closed_manifold(m.triangles()));
        assert!(!m.validate(0.0).unwrap().summary.is_closed);
        let Err(Error::Kernel { report, mesh, .. }) = solid_result(m.clone(), "it") else {
            panic!("an open result must be a kernel error");
        };
        assert!(!report.unwrap().summary.is_closed);
        assert_eq!(*mesh.unwrap(), m);
    }

    #[test]
    fn copies_of_one_exact_point_stay_apart() {
        let s = shell(vec![A, B, C, D, D], vec![[0, 2, 1], [0, 1, 3], [1, 2, 4]]);
        let (out, _) = collapse_rounded(s, &[0, 1, 2, 3, 3]);
        assert_eq!((out.vertices.len(), out.triangles.len()), (5, 3));
    }

    #[test]
    fn faces_with_collinear_corners_are_split_away() {
        let m = [1.0, 0.0, 0.0];
        let s = shell(
            vec![A, B, C, D, m],
            vec![
                [0, 2, 1],
                [0, 1, 4],
                [0, 4, 3],
                [4, 1, 3],
                [0, 3, 2],
                [1, 2, 3],
            ],
        );
        let (out, _) = collapse_rounded(s, &[0, 1, 2, 3, 4]);
        check(&out);
        assert_eq!(out.triangles.len(), 6);
        let area = |t: &[u32; 3]| {
            let p = t.map(|i| out.vertices[i as usize]);
            !collinear3(p[0], p[1], p[2])
        };
        assert!(out.triangles.iter().all(area));
        assert!(out.face_source.contains(&(0, 0)));
    }

    #[test]
    fn flips_on_rounded_points_end_and_keep_the_shell_closed() {
        let v: Vec<[f64; 3]> = vec![
            [700078.4443665423, 7499882.444366542, 327.3486700859249],
            [700078.2071804033, 7499882.207180403, 327.261322566267],
            [700078.7188038859, 7499882.718803885, 327.44973594655903],
            [700078.7188038859, 7499882.718803885, 327.44973594659086],
            [700079.3243227224, 7499883.324322722, 327.5160718666926],
            [700080.9600533535, 7499884.960053354, 327.6952697543524],
            [700080.9600533535, 7499884.960053354, 327.6952697543533],
            [700077.787231901, 7499881.787231901, 327.1547302266305],
            [700080.9600533535, 7499884.960053354, 327.69526975435275],
            [700075.0, 7499879.0, 326.31824503651836],
            [700075.0, 7499879.0, 326.1635151547103],
            [700076.9973626635, 7499880.997362663, 326.7991542135016],
            [700077.6132896044, 7499881.613289604, 326.99516629889933],
            [700077.8357285517, 7499881.835728551, 327.05523173366186],
            [700077.8815485869, 7499881.881548587, 327.0733674685544],
            [700078.3691364272, 7499882.3691364275, 327.10634971351146],
            [700078.3796673589, 7499882.379667359, 327.1022548467648],
            [700078.4260729539, 7499882.426072954, 327.12142251730654],
            [700079.972448509, 7499883.972448509, 327.4942443004887],
            [700075.3512054217, 7499879.351205422, 327.06499324347294],
            [700075.0, 7499879.0, 326.31824503656975],
            [700075.3512054217, 7499879.351205422, 327.0649932434694],
            [700077.7872319008, 7499881.787231901, 327.15473022662326],
            [700077.7872319008, 7499881.787231901, 327.1547302266264],
            [700077.7872319008, 7499881.787231901, 327.154730226627],
            [700075.3512054224, 7499879.351205423, 327.0649932435096],
            [700075.3512054217, 7499879.351205422, 327.06499324347385],
            [700077.787231901, 7499881.787231901, 327.1547302266468],
            [700077.787231901, 7499881.787231901, 327.1547302266412],
            [700077.6395250404, 7499881.639525041, 327.3048681992224],
            [700077.618429519, 7499881.618429519, 327.32631093261244],
            [700075.3512054217, 7499879.351205422, 327.06499324347294],
            [700080.9600533532, 7499884.960053354, 327.6952697543287],
            [700079.3243227208, 7499883.324322721, 327.51607186652836],
            [700078.2071804032, 7499882.207180403, 327.2613225663012],
            [700078.2071804032, 7499882.207180403, 327.26132256629705],
            [700078.4443665423, 7499882.444366542, 327.3486700860085],
            [700077.6395250404, 7499881.639525041, 327.30486819922487],
            [700077.787231901, 7499881.787231901, 327.15473022669346],
            [700077.787231901, 7499881.787231901, 327.15473022667715],
            [700078.2071804032, 7499882.207180403, 327.2613225662874],
            [700078.2071804032, 7499882.207180403, 327.2613225662879],
            [700078.2071804032, 7499882.207180403, 327.26132256628745],
            [700078.7188038857, 7499882.718803885, 327.44973594659166],
            [700077.7872319008, 7499881.787231901, 327.1547302266699],
        ];
        let t: Vec<[u32; 3]> = vec![
            [0, 1, 2],
            [2, 1, 3],
            [3, 1, 4],
            [4, 1, 5],
            [4, 5, 6],
            [7, 5, 1],
            [7, 8, 5],
            [9, 10, 7],
            [10, 11, 7],
            [11, 12, 7],
            [7, 12, 13],
            [7, 13, 14],
            [7, 14, 15],
            [15, 16, 17],
            [7, 15, 17],
            [7, 17, 18],
            [7, 18, 8],
            [9, 19, 20],
            [20, 19, 21],
            [9, 22, 19],
            [22, 23, 19],
            [9, 7, 22],
            [22, 7, 24],
            [19, 25, 26],
            [19, 23, 25],
            [7, 27, 28],
            [25, 29, 30],
            [28, 24, 7],
            [29, 25, 23],
            [31, 19, 26],
            [32, 4, 6],
            [32, 33, 4],
            [34, 35, 1],
            [34, 1, 0],
            [36, 34, 0],
            [31, 21, 19],
            [30, 29, 37],
            [38, 7, 39],
            [38, 27, 7],
            [39, 7, 40],
            [40, 7, 41],
            [42, 40, 41],
            [41, 7, 35],
            [35, 7, 1],
            [33, 32, 34],
            [43, 33, 36],
            [36, 33, 34],
            [21, 31, 20],
            [31, 26, 20],
            [26, 25, 20],
            [25, 44, 20],
            [44, 25, 37],
            [25, 30, 37],
            [12, 11, 44],
            [10, 9, 11],
            [9, 20, 11],
            [11, 20, 44],
            [13, 12, 44],
            [14, 13, 44],
            [38, 39, 44],
            [39, 14, 44],
            [40, 15, 14],
            [40, 14, 39],
            [42, 15, 40],
            [42, 16, 15],
            [35, 17, 41],
            [34, 17, 35],
            [34, 32, 18],
            [17, 34, 18],
            [18, 5, 8],
            [5, 18, 32],
            [5, 32, 6],
            [17, 16, 42],
            [17, 42, 41],
            [43, 3, 33],
            [33, 3, 4],
            [43, 2, 3],
            [36, 0, 43],
            [43, 0, 2],
            [37, 23, 22],
            [37, 22, 44],
            [29, 23, 37],
            [44, 22, 24],
            [44, 24, 28],
            [44, 28, 38],
            [38, 28, 27],
        ];
        let point: Vec<u32> = vec![
            291, 1, 420, 422, 441, 446, 449, 17, 447, 18, 5, 286, 287, 290, 304, 312, 314, 316,
            362, 19, 89, 259, 272, 273, 275, 267, 265, 281, 279, 277, 269, 263, 443, 442, 355, 358,
            292, 270, 285, 301, 306, 310, 308, 86, 92,
        ];
        let n = v.len();
        let s = Assembled {
            face_source: (0..t.len() as u32).map(|f| (0, f)).collect(),
            vertices: v,
            triangles: t,
            vertex_rows: vec![None; n],
            normal: vec![],
        };
        assert!(is_closed_manifold(&s.triangles));
        let (out, kept) = collapse_rounded(s, &point);
        assert!(is_closed_manifold(&out.triangles));
        assert_eq!(kept.len(), out.vertices.len());
    }
}
