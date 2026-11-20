//! Materialization of exact results as `f64` meshes: method C of
//! `docs/design/rounding.md`. Every output coordinate is the rounding `r(p)`
//! of some exact vertex `p` of the result; the operations below decide with
//! exact tests on those roundings which features rounding cannot place and
//! collapse them, and the result is certified by exact `validate`.

use std::collections::{HashMap, HashSet};

use wireframes_core::exact::{Expansion, shell_volume_sign};
use wireframes_core::predicates::{Sign, collinear3, kept_axes, orient2d_axis, orient3d};
use wireframes_core::{Error, Mesh, Result};

/// What rounding the exact result to `f64` did, and the bounds it proves.
///
/// Lengths are in the mesh's length unit and rounded up. A vertex counts as
/// moved when it ends at a coordinate other than its own rounding, so it may
/// lie farther than half an ulp from its exact position; `max_move` is the
/// largest distance per axis of any output vertex from the exact position of
/// any exact vertex that ended there, and every vertex lies within
/// `move_bound` of its exact position. `volume_bound` bounds
/// `|vol(output) - vol(exact)|`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct RoundingReport {
    /// Output vertices that are rounded intersection points, not input vertices.
    pub rounded_vertices: usize,
    /// Exact vertices that end at a coordinate other than their own rounding.
    pub moved_vertices: usize,
    /// The input vertices among `moved_vertices`.
    pub moved_input_vertices: usize,
    pub max_move: [f64; 3],
    pub move_bound: [f64; 3],
    /// Vertices merged into another of equal coordinates (C1).
    pub merged_vertices: usize,
    /// Edges collapsed (C2).
    pub collapsed_edges: usize,
    /// Rounded vertices inserted into faces they lie within rounding of (C5),
    /// including input vertices inserted into the other input's faces.
    pub inserted_points: usize,
    /// Edges of an input flipped to match the other input's faces where the
    /// two agree to rounding.
    pub flipped_edges: usize,
    /// Shells of zero or wrong-sign rounded volume removed (C3).
    pub removed_shells: usize,
    /// Faces removed in opposite pairs on the same three vertices (C4).
    pub cancelled_faces: usize,
    pub volume_bound: f64,
    /// Operations applied, and the bound on them proven before the first.
    pub steps: usize,
    pub step_bound: usize,
}

/// The exact result before rounding: per vertex copy its exact point's
/// identity (copies of one exact point, where shells touch, share it), its
/// rounding, and a bound per axis on the distance from the rounding to the
/// exact point (zero for input vertices); per face its corners, source and
/// the outward normal of its exact plane, rounded; and the faces of each
/// shell.
#[derive(Clone)]
pub(crate) struct Draft {
    pub point: Vec<u32>,
    pub home: Vec<[f64; 3]>,
    pub err: Vec<[f64; 3]>,
    pub rows: Vec<Option<u32>>,
    pub triangles: Vec<[u32; 3]>,
    pub face_source: Vec<(u8, u32)>,
    pub normal: Vec<[f64; 3]>,
    pub shells: Vec<Vec<u32>>,
    /// The exact result is closed, so the operations over all shells apply.
    pub closed: bool,
}

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

/// Half the spacing of `f64` at `x`: the largest distance from a real number
/// to `x` when `x` is its nearest `f64`.
pub(crate) fn half_ulp(x: f64) -> f64 {
    let a = x.abs();
    ((a.next_up() - a) / 2.0).max(f64::from_bits(1))
}

pub(crate) fn up(x: f64) -> f64 {
    (x * (1.0 + 64.0 * f64::EPSILON)).next_up()
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn norm(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

/// Upper bound on the volume of `t + B(r)` for a triangle `t`, by Steiner's
/// formula `2 A r + (pi / 2) P r^2 + (4 / 3) pi r^3`.
fn tube(t: [[f64; 3]; 3], r: f64) -> f64 {
    if r == 0.0 {
        return 0.0;
    }
    let (u, v) = (sub(t[1], t[0]), sub(t[2], t[0]));
    let c = [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ];
    let area = up(norm(c) / 2.0);
    let perim = up(norm(u) + norm(v) + norm(sub(t[2], t[1])));
    let pi = std::f64::consts::PI;
    up(2.0 * area * r + pi / 2.0 * perim * r * r + 4.0 / 3.0 * pi * r * r * r)
}

/// `a . (b x c)` exactly: six times the signed volume of the tetrahedron on
/// the origin and the face `a b c`.
fn det(t: [[f64; 3]; 3]) -> Expansion {
    let [a, b, c] = t;
    (0..3).fold(Expansion::zero(), |acc, k| {
        let (i, j) = kept_axes(k);
        acc.sum(
            &Expansion::from_f64(b[i])
                .scale(c[j])
                .diff(&Expansion::from_f64(b[j]).scale(c[i]))
                .scale(a[k]),
        )
    })
}

/// `|vol(after) - vol(before)|`, rounded up, from the faces that differ.
fn volume_change(before: &[[[f64; 3]; 3]], after: &[[[f64; 3]; 3]]) -> f64 {
    let key = |t: &[[f64; 3]; 3]| {
        let b = t.map(|p| p.map(|x| (x + 0.0).to_bits()));
        let k = (0..3).min_by_key(|&k| b[k]).expect("three corners");
        [b[k], b[(k + 1) % 3], b[(k + 2) % 3]]
    };
    let mut count: HashMap<[[u64; 3]; 3], (i64, [[f64; 3]; 3])> = HashMap::new();
    for (t, s) in before
        .iter()
        .map(|t| (t, -1))
        .chain(after.iter().map(|t| (t, 1)))
    {
        count.entry(key(t)).or_insert((0, *t)).0 += s;
    }
    let mut diff: Vec<_> = count.into_values().filter(|(n, _)| *n != 0).collect();
    diff.sort_by(|a, b| {
        let bits = |t: &[[f64; 3]; 3]| t.map(|p| p.map(f64::to_bits));
        bits(&a.1).cmp(&bits(&b.1)).then(a.0.cmp(&b.0))
    });
    let total = diff.iter().fold(Expansion::zero(), |acc, (n, t)| {
        acc.sum(&det(*t).scale(*n as f64))
    });
    up(total.estimate().abs() / 6.0)
}

/// Sign of `n . ((b - a) x (c - a))`, exact on the given `f64` values.
pub(crate) fn facing(c: [[f64; 3]; 3], n: [f64; 3]) -> Sign {
    let u = [0, 1, 2].map(|k| Expansion::diff_f64(c[1][k], c[0][k]));
    let v = [0, 1, 2].map(|k| Expansion::diff_f64(c[2][k], c[0][k]));
    (0..3)
        .fold(Expansion::zero(), |acc, k| {
            let (i, j) = kept_axes(k);
            acc.sum(&u[i].product(&v[j]).diff(&u[j].product(&v[i])).scale(n[k]))
        })
        .sign()
}

/// Counts of what a cleanup did.
#[derive(Default)]
pub(crate) struct Counts {
    pub merged: usize,
    pub collapsed: usize,
    pub cancelled: usize,
    pub flips: usize,
    pub rounds: usize,
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
/// point, where the shell touches itself, stay apart and never merge. With
/// `cap`, a merge that would leave a vertex past [`move_cap`] from the exact
/// position of a vertex that ends at it (`err` bounds the rounding of each) is
/// refused. Also returns the exact point of each kept vertex, the output vertex each input
/// vertex ended at (`None` when its faces are gone), and the counts.
pub(crate) fn collapse_rounded(
    shell: Assembled,
    point: &[u32],
    err: &[[f64; 3]],
    cap: bool,
) -> (Assembled, Vec<u32>, Vec<Option<u32>>, Counts) {
    let mut counts = Counts::default();
    let n = shell.vertices.len();
    let mut span: Vec<[f64; 3]> = err.to_vec();
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
                counts.merged += 1;
                for k in 0..3 {
                    span[order[start]][k] = span[order[start]][k].max(span[i][k]);
                }
            }
        }
        start = w;
    }
    let mut to = rep.clone();
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
    counts.cancelled += cancel_pairs(&tris, &mut alive);
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
    // again, so every flip removes an edge of the shell as it came in or one
    // an edge collapse made, and the loop ends; without this, flips on
    // rounded points can cycle forever.
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
                        let (g, k) = (v[gone as usize], v[keep as usize]);
                        let grown: [f64; 3] = [0, 1, 2].map(|a| {
                            up(span[keep as usize][a]
                                .max(span[gone as usize][a] + (g[a] - k[a]).abs()))
                        });
                        if cap && (0..3).any(|a| grown[a] > move_cap(k)) {
                            continue;
                        }
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
                        to[gone as usize] = keep;
                        span[keep as usize] = grown;
                        counts.collapsed += 1;
                        break;
                    }
                }
                continue;
            }
            made.insert(ekey(m, d));
            counts.flips += 1;
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
        counts.rounds += 1;
        let (count, bad) = rounding_problems(v, &tris, &alive);
        if count >= last {
            if let Some((t, a, s, nr, tt)) = saved.take() {
                (tris, alive, src, nrm, to) = (t, a, s, nr, tt);
            }
            break;
        }
        if count == 0 {
            break;
        }
        last = count;
        saved = Some((
            tris.clone(),
            alive.clone(),
            src.clone(),
            nrm.clone(),
            to.clone(),
        ));
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
    let fwd = (0..n)
        .map(|i| {
            let mut j = i;
            while to[j] as usize != j {
                j = to[j] as usize;
            }
            new_id[j]
        })
        .collect();
    (out, kept_point, fwd, counts)
}

/// Kills opposite faces on the same three vertices in pairs; returns how many
/// faces it killed.
fn cancel_pairs(tris: &[[u32; 3]], alive: &mut [bool]) -> usize {
    let mut by_set: HashMap<[u32; 3], Vec<usize>> = HashMap::new();
    for (f, t) in tris.iter().enumerate().filter(|(f, _)| alive[*f]) {
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
    let mut killed = 0;
    for fs in sets {
        let (mut up, mut down): (Vec<usize>, Vec<usize>) =
            fs.into_iter().partition(|&f| parity(tris[f]));
        while let (Some(a), Some(b)) = (up.pop(), down.pop()) {
            alive[a] = false;
            alive[b] = false;
            killed += 2;
        }
    }
    killed
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

/// The whole output after the per-shell cleanup, for the operations that see
/// every shell at once.
#[derive(Clone)]
struct Global {
    v: Vec<[f64; 3]>,
    point: Vec<u32>,
    rows: Vec<Option<u32>>,
    tris: Vec<[u32; 3]>,
    src: Vec<(u8, u32)>,
    nrm: Vec<[f64; 3]>,
    alive: Vec<bool>,
    /// Vertex each vertex was merged or collapsed into; itself when kept.
    to: Vec<u32>,
    /// Per face, the face of the output after the per-shell cleanup it came from.
    orig: Vec<u32>,
    /// Per vertex, the largest distance per axis from its coordinates to the
    /// exact position of any exact vertex that ended at it, rounded up.
    span: Vec<[f64; 3]>,
}

/// A local change [`Global::descend`] tries.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
enum Change {
    Flip(usize, usize),
    Join(usize, usize),
    Collapse(u32, u32),
    Insert(u32, usize, Option<usize>),
}

type Ev = [Expansion; 3];

fn ev(a: [f64; 3], b: [f64; 3]) -> Ev {
    [0, 1, 2].map(|k| Expansion::diff_f64(a[k], b[k]))
}

fn cross_e(u: &Ev, v: &Ev) -> Ev {
    [0, 1, 2].map(|k| {
        let (i, j) = kept_axes(k);
        u[i].product(&v[j]).diff(&u[j].product(&v[i]))
    })
}

fn dot_f(a: &Ev, p: [f64; 3]) -> Expansion {
    (0..3).fold(Expansion::zero(), |acc, k| acc.sum(&a[k].scale(p[k])))
}

/// Whether the projections of `pts` and of the box `[lo, hi]` on `axis`
/// overlap, exactly; a zero axis separates nothing.
fn overlap_on(axis: &Ev, pts: &[[f64; 3]], lo: [f64; 3], hi: [f64; 3]) -> bool {
    if axis.iter().all(|x| x.sign() == Sign::Zero) {
        return true;
    }
    let (mut bmin, mut bmax) = (Expansion::zero(), Expansion::zero());
    for k in 0..3 {
        let (a, b) = (axis[k].scale(lo[k]), axis[k].scale(hi[k]));
        if axis[k].sign() == Sign::Positive {
            (bmin, bmax) = (bmin.sum(&a), bmax.sum(&b));
        } else {
            (bmin, bmax) = (bmin.sum(&b), bmax.sum(&a));
        }
    }
    let proj: Vec<Expansion> = pts.iter().map(|&p| dot_f(axis, p)).collect();
    let below = proj.iter().all(|x| x.diff(&bmin).sign() == Sign::Negative);
    let above = proj.iter().all(|x| x.diff(&bmax).sign() == Sign::Positive);
    !below && !above
}

/// Whether the convex hull of `pts` (a segment or a triangle) meets the
/// closed box `[lo, hi]`, by the separating axis theorem with exact signs.
pub(crate) fn box_meets(pts: &[[f64; 3]], lo: [f64; 3], hi: [f64; 3]) -> bool {
    for k in 0..3 {
        if pts.iter().all(|p| p[k] < lo[k]) || pts.iter().all(|p| p[k] > hi[k]) {
            return false;
        }
    }
    let unit = |k: usize| {
        let mut u = [0.0; 3];
        u[k] = 1.0;
        ev(u, [0.0; 3])
    };
    let edges: Vec<Ev> = (0..pts.len())
        .map(|i| ev(pts[(i + 1) % pts.len()], pts[i]))
        .collect();
    if pts.len() == 3 {
        let n = cross_e(&edges[0], &ev(pts[2], pts[0]));
        if !overlap_on(&n, pts, lo, hi) {
            return false;
        }
    }
    edges
        .iter()
        .all(|e| (0..3).all(|k| overlap_on(&cross_e(e, &unit(k)), pts, lo, hi)))
}

/// The box about `h` within which a face with corners `c` can reach `h` once
/// both are rounded, rounded outward: an ulp of `h` on each axis, or, when
/// `wide`, on every axis the sum of the rounding radii of `h` and of the
/// face's farthest-rounding corner (the length of their half-ulp vectors),
/// which holds every point the separation lemma lets a rounded face reach.
pub(crate) fn hot_box(h: [f64; 3], c: &[[f64; 3]], wide: bool) -> ([f64; 3], [f64; 3]) {
    let radius = |p: [f64; 3]| up(norm(p.map(half_ulp)));
    let w: [f64; 3] = if wide {
        let reach = up(2.0 * (radius(h) + c.iter().map(|&p| radius(p)).fold(0.0, f64::max)));
        [reach; 3]
    } else {
        h.map(|x| 2.0 * half_ulp(x))
    };
    (
        [0, 1, 2].map(|k| (h[k] - w[k]).next_down()),
        [0, 1, 2].map(|k| (h[k] + w[k]).next_up()),
    )
}

/// How far, on every axis, a vertex at `p` may end from the exact position of
/// any exact vertex that ends at it: four ulps of its largest coordinate.
fn move_cap(p: [f64; 3]) -> f64 {
    8.0 * p.map(half_ulp).into_iter().fold(0.0, f64::max)
}

impl Global {
    fn live(&self) -> impl Iterator<Item = usize> + '_ {
        (0..self.tris.len()).filter(|&f| self.alive[f])
    }

    /// The live faces as a mesh of the vertices they use, with the face
    /// index of each.
    fn mesh(&self) -> Result<(Mesh, Vec<usize>)> {
        let ids: Vec<usize> = self.live().collect();
        let mut id: Vec<Option<u32>> = vec![None; self.v.len()];
        let mut v = vec![];
        let tris = ids
            .iter()
            .map(|&f| {
                self.tris[f].map(|i| {
                    *id[i as usize].get_or_insert_with(|| {
                        v.push(self.v[i as usize]);
                        v.len() as u32 - 1
                    })
                })
            })
            .collect();
        Ok((Mesh::new(v, tris)?, ids))
    }

    /// C1: vertices with equal coordinates become one, unless two of them are
    /// copies of one exact point; then faces that repeat a vertex go.
    fn merge_equal(&mut self) -> usize {
        let used: HashSet<u32> = self.live().flat_map(|f| self.tris[f]).collect();
        let mut ids: Vec<u32> = used.into_iter().collect();
        let key = |i: u32| self.v[i as usize].map(|x| (x + 0.0).to_bits());
        ids.sort_unstable_by_key(|&i| (key(i), self.rows[i as usize].is_none(), i));
        let mut merged = 0;
        let mut start = 0;
        for w in 1..=ids.len() {
            if w < ids.len() && key(ids[w]) == key(ids[start]) {
                continue;
            }
            let mut ps: Vec<u32> = ids[start..w]
                .iter()
                .map(|&i| self.point[i as usize])
                .collect();
            ps.sort_unstable();
            if ps.windows(2).all(|p| p[0] != p[1]) {
                for &i in &ids[start + 1..w] {
                    self.to[i as usize] = ids[start];
                    merged += 1;
                }
            }
            start = w;
        }
        if merged > 0 {
            for f in 0..self.tris.len() {
                if !self.alive[f] {
                    continue;
                }
                let t = self.tris[f].map(|i| self.to[i as usize]);
                self.tris[f] = t;
                if t[0] == t[1] || t[1] == t[2] || t[0] == t[2] {
                    self.alive[f] = false;
                }
            }
        }
        merged
    }

    /// C3: removes each edge-connected component whose rounded volume is zero,
    /// or negative while `validate` reports it as an inward shell, or, with
    /// `crossing`, that holds a face of a self-intersecting, overlapping or
    /// flat problem, when that volume is within the bound (V) gives for
    /// rounding it, so its exact volume is too.
    fn remove_flat_shells(&mut self, err: &[f64], crossing: bool) -> Result<usize> {
        let (m, ids) = self.mesh()?;
        let report = m.validate(0.0)?;
        let mut parent: Vec<u32> = (0..ids.len() as u32).collect();
        fn root(parent: &mut [u32], mut i: u32) -> u32 {
            while parent[i as usize] != i {
                parent[i as usize] = parent[parent[i as usize] as usize];
                i = parent[i as usize];
            }
            i
        }
        let mut at: HashMap<(u32, u32), u32> = HashMap::new();
        for (k, &f) in ids.iter().enumerate() {
            let t = self.tris[f];
            for e in 0..3 {
                let key = ekey(t[e], t[(e + 1) % 3]);
                if let Some(&o) = at.get(&key) {
                    let (a, b) = (root(&mut parent, o), root(&mut parent, k as u32));
                    parent[a.max(b) as usize] = a.min(b);
                } else {
                    at.insert(key, k as u32);
                }
            }
        }
        let mut groups: HashMap<u32, Vec<usize>> = HashMap::new();
        for k in 0..ids.len() {
            groups
                .entry(root(&mut parent, k as u32))
                .or_default()
                .push(k);
        }
        let inward: HashSet<usize> = report
            .problems
            .iter()
            .filter(|p| p.kind == wireframes_core::MeshProblemKind::InwardShell)
            .filter_map(|p| p.face.map(|f| f as usize))
            .collect();
        use wireframes_core::MeshProblemKind as K;
        let troubled: HashSet<usize> = report
            .problems
            .iter()
            .filter(|p| {
                matches!(
                    p.kind,
                    K::SelfIntersection | K::CoplanarOverlap | K::DegenerateFace
                )
            })
            .flat_map(|p| [p.face, p.other])
            .flatten()
            .map(|f| f as usize)
            .collect();
        let mut removed = 0;
        let mut keys: Vec<u32> = groups.keys().copied().collect();
        keys.sort_unstable();
        for g in keys {
            let ks = &groups[&g];
            let faces: Vec<[u32; 3]> = ks.iter().map(|&k| self.tris[ids[k]]).collect();
            let sign = shell_volume_sign(&self.v, &faces);
            let flagged = ks.iter().any(|k| inward.contains(k));
            let crossed = crossing && ks.iter().any(|k| troubled.contains(k));
            if !(sign == Sign::Zero || (sign == Sign::Negative && flagged) || crossed) {
                continue;
            }
            let corners: Vec<[[f64; 3]; 3]> = faces
                .iter()
                .map(|t| t.map(|i| self.v[i as usize]))
                .collect();
            let vol = volume_change(&corners, &[]);
            let allowed: f64 = ks
                .iter()
                .zip(&corners)
                .map(|(&k, c)| {
                    let r = self.tris[ids[k]]
                        .iter()
                        .map(|&i| err[i as usize])
                        .fold(0.0, f64::max);
                    tube(*c, 2.0 * r)
                })
                .sum();
            if vol <= allowed {
                for &k in ks {
                    self.alive[ids[k]] = false;
                }
                removed += 1;
            }
        }
        Ok(removed)
    }

    /// C5 for each self-intersecting or overlapping pair `validate` reports:
    /// a corner `h` of one face whose hot box meets the other face `g`, not at
    /// a corner, is inserted into `g`, on the edge of `g` the box meets when
    /// there is one (splitting both faces on that edge), else inside `g`.
    /// Each (exact point, original face) pair is inserted at most once.
    /// Returns how many points it inserted.
    fn insert_hot_points(
        &mut self,
        report: &wireframes_core::MeshReport,
        ids: &[usize],
        done: &mut HashSet<(u32, u32)>,
        wide: bool,
    ) -> (usize, usize) {
        use wireframes_core::MeshProblemKind as K;
        let mut inserted = 0;
        let mut collapsed = 0;
        for p in &report.problems {
            if p.kind == K::DegenerateFace
                && let Some(f) = p.face
            {
                inserted += usize::from(self.unflatten(ids[f as usize], done));
                continue;
            }
            if !matches!(p.kind, K::SelfIntersection | K::CoplanarOverlap) {
                continue;
            }
            let (Some(f), Some(o)) = (p.face, p.other) else {
                continue;
            };
            let (f, o) = (ids[f as usize], ids[o as usize]);
            if self.join_opposite(f, o) {
                collapsed += 1;
                continue;
            }
            'pair: for (a, b) in [(f, o), (o, f)] {
                if !self.alive[a] || !self.alive[b] {
                    break;
                }
                let tb = self.tris[b];
                let cb = tb.map(|i| self.v[i as usize]);
                for h in self.tris[a] {
                    let at = self.v[h as usize];
                    if tb.contains(&h) || cb.contains(&at) {
                        continue;
                    }
                    let key = (self.point[h as usize], self.orig[b]);
                    if done.contains(&key) {
                        continue;
                    }
                    let (lo, hi) = hot_box(at, &cb, wide);
                    if !box_meets(&cb, lo, hi) {
                        continue;
                    }
                    if let Some(k) = (0..3).find(|&k| box_meets(&[cb[k]], lo, hi)) {
                        if self.collapse(h, tb[k]) {
                            collapsed += 1;
                            break 'pair;
                        }
                        continue;
                    }
                    let on_edge = (0..3).find(|&e| box_meets(&[cb[e], cb[(e + 1) % 3]], lo, hi));
                    if self.insert(h, b, on_edge) {
                        done.insert(key);
                        inserted += 1;
                        break 'pair;
                    }
                }
            }
        }
        (inserted, collapsed)
    }

    /// C4 for faces `f` and `o` on the same three coordinates with opposite
    /// orientation whose corners differ by copies of one exact point: the
    /// copies in `o` become those of `f`, so the faces cancel. Returns whether
    /// it joined any.
    fn join_opposite(&mut self, f: usize, o: usize) -> bool {
        if !self.alive[f] || !self.alive[o] {
            return false;
        }
        let (tf, to) = (self.tris[f], self.tris[o]);
        let at = |i: u32| self.v[i as usize].map(|x| (x + 0.0).to_bits());
        let Some(k) = (0..3).find(|&k| at(to[k]) == at(tf[0])) else {
            return false;
        };
        let opposite = at(to[(k + 2) % 3]) == at(tf[1]) && at(to[(k + 1) % 3]) == at(tf[2]);
        if !opposite {
            return false;
        }
        let pairs = [
            (tf[0], to[k]),
            (tf[1], to[(k + 2) % 3]),
            (tf[2], to[(k + 1) % 3]),
        ];
        if pairs
            .iter()
            .any(|&(x, y)| x != y && self.point[x as usize] != self.point[y as usize])
        {
            return false;
        }
        let mut joined = false;
        for (keep, gone) in pairs {
            if keep == gone {
                continue;
            }
            joined = true;
            self.to[gone as usize] = keep;
            for g in 0..self.tris.len() {
                if self.alive[g] && self.tris[g].contains(&gone) {
                    self.tris[g] = self.tris[g].map(|i| if i == gone { keep } else { i });
                }
            }
        }
        joined
    }

    /// Removes each edge-connected run of flat faces and splits the edges of
    /// the faces around it at the run's vertices, ordered along the axis the
    /// run spreads most on, so the faces on either side meet edge for edge. A
    /// face is flat when its rounded corners are collinear, or when the hot
    /// box of its middle corner meets the segment of the other two; a run is
    /// removed when one of its faces is collinear or in `seeds` and every
    /// vertex's wide hot box meets the segment between its two extreme
    /// vertices.
    /// A run is put back when the edges then do not pair. Returns how many faces it
    /// removed.
    fn zip_flat_runs(&mut self, seeds: &HashSet<usize>) -> usize {
        let flat = |g: &Global, f: usize| {
            let c = g.tris[f].map(|i| g.v[i as usize]);
            if collinear3(c[0], c[1], c[2]) {
                return true;
            }
            let axis = (0..3)
                .max_by(|&a, &b| {
                    let spread = |k: usize| {
                        c.iter().map(|p| p[k]).fold(f64::NEG_INFINITY, f64::max)
                            - c.iter().map(|p| p[k]).fold(f64::INFINITY, f64::min)
                    };
                    spread(a).total_cmp(&spread(b))
                })
                .expect("three axes");
            let mut by = [0usize, 1, 2];
            by.sort_by(|&a, &b| c[a][axis].total_cmp(&c[b][axis]));
            let (lo, hi) = hot_box(c[by[1]], &[], false);
            box_meets(&[c[by[0]], c[by[2]]], lo, hi)
        };
        let flats: Vec<usize> = self.live().filter(|&f| flat(self, f)).collect();
        if flats.is_empty() {
            return 0;
        }
        let mut parent: HashMap<usize, usize> = flats.iter().map(|&f| (f, f)).collect();
        fn root(parent: &mut HashMap<usize, usize>, mut f: usize) -> usize {
            while parent[&f] != f {
                let up = parent[&parent[&f]];
                parent.insert(f, up);
                f = up;
            }
            f
        }
        let mut by_edge: HashMap<(u32, u32), usize> = HashMap::new();
        for &f in &flats {
            let t = self.tris[f];
            for k in 0..3 {
                let e = ekey(t[k], t[(k + 1) % 3]);
                if let Some(&g) = by_edge.get(&e) {
                    let (a, b) = (root(&mut parent, f), root(&mut parent, g));
                    parent.insert(a.max(b), a.min(b));
                } else {
                    by_edge.insert(e, f);
                }
            }
        }
        let mut runs: HashMap<usize, Vec<usize>> = HashMap::new();
        for &f in &flats {
            runs.entry(root(&mut parent, f)).or_default().push(f);
        }
        let mut keys: Vec<usize> = runs.keys().copied().collect();
        keys.sort_unstable();
        let mut removed = 0;
        for key in keys {
            let run = &runs[&key];
            if run.iter().any(|&f| !self.alive[f]) {
                continue;
            }
            let exact = |g: &Global, f: usize| {
                let c = g.tris[f].map(|i| g.v[i as usize]);
                collinear3(c[0], c[1], c[2])
            };
            if !run.iter().any(|&f| seeds.contains(&f) || exact(self, f)) {
                continue;
            }
            let mut on: Vec<u32> = run.iter().flat_map(|&f| self.tris[f]).collect();
            on.sort_unstable();
            on.dedup();
            let axis = (0..3)
                .max_by(|&a, &b| {
                    let spread = |k: usize| {
                        let c = on.iter().map(|&i| self.v[i as usize][k]);
                        c.clone().fold(f64::NEG_INFINITY, f64::max)
                            - c.fold(f64::INFINITY, f64::min)
                    };
                    spread(a).total_cmp(&spread(b))
                })
                .expect("three axes");
            let at = |g: &Global, i: u32| g.v[i as usize][axis];
            let ends = [
                *on.iter()
                    .min_by(|&&p, &&q| at(self, p).total_cmp(&at(self, q)))
                    .expect("corners"),
                *on.iter()
                    .max_by(|&&p, &&q| at(self, p).total_cmp(&at(self, q)))
                    .expect("corners"),
            ]
            .map(|i| self.v[i as usize]);
            let straight = on.iter().all(|&i| {
                let (lo, hi) = hot_box(self.v[i as usize], &ends, true);
                box_meets(&ends, lo, hi)
            });
            if !straight {
                continue;
            }
            for &f in run {
                self.alive[f] = false;
            }
            let inside: HashSet<(u32, u32)> = run
                .iter()
                .flat_map(|&f| {
                    let t = self.tris[f];
                    (0..3).map(move |k| (t[k], t[(k + 1) % 3]))
                })
                .collect();
            let around: Vec<usize> = self
                .live()
                .filter(|&h| {
                    let t = self.tris[h];
                    (0..3).any(|k| inside.contains(&(t[(k + 1) % 3], t[k])))
                })
                .collect();
            let (mut made, mut killed) = (vec![], vec![]);
            for h in around {
                let t = self.tris[h];
                let Some(k) = (0..3).find(|&k| inside.contains(&(t[(k + 1) % 3], t[k]))) else {
                    continue;
                };
                let (x, y, c) = (t[k], t[(k + 1) % 3], t[(k + 2) % 3]);
                let (lo, hi) = (at(self, x).min(at(self, y)), at(self, x).max(at(self, y)));
                let mut mid: Vec<u32> = on
                    .iter()
                    .copied()
                    .filter(|&w| w != x && w != y && lo < at(self, w) && at(self, w) < hi)
                    .collect();
                if mid.is_empty() {
                    continue;
                }
                mid.sort_by(|&p, &q| at(self, p).total_cmp(&at(self, q)));
                if at(self, x) > at(self, y) {
                    mid.reverse();
                }
                self.alive[h] = false;
                killed.push(h);
                let mut prev = x;
                for w in mid.into_iter().chain(std::iter::once(y)) {
                    self.push([prev, w, c], h);
                    made.push(self.tris.len() - 1);
                    prev = w;
                }
            }
            let line: HashSet<u32> = on.iter().copied().collect();
            let mut count: HashMap<(u32, u32), i32> = HashMap::new();
            for f in self.live() {
                let t = self.tris[f];
                for k in 0..3 {
                    let (a, b) = (t[k], t[(k + 1) % 3]);
                    if line.contains(&a) && line.contains(&b) {
                        *count.entry(ekey(a, b)).or_default() += if a < b { 1 } else { -1 };
                    }
                }
            }
            let still_flat = made
                .iter()
                .filter(|&&f| self.alive[f] && flat(self, f))
                .count();
            if count.values().any(|&n| n != 0) || still_flat >= run.len() {
                for f in made {
                    self.alive[f] = false;
                }
                for f in run.iter().copied().chain(killed) {
                    self.alive[f] = true;
                }
                continue;
            }
            removed += run.len();
        }
        removed
    }

    /// C5 for a face `f` whose corners are collinear: its middle corner, on
    /// the segment of the other two, is inserted on that edge of the face
    /// across it, which leaves `f` with a repeated corner, so it goes.
    fn unflatten(&mut self, f: usize, done: &mut HashSet<(u32, u32)>) -> bool {
        if !self.alive[f] {
            return false;
        }
        let t = self.tris[f];
        let c = t.map(|i| self.v[i as usize]);
        if !collinear3(c[0], c[1], c[2]) {
            return false;
        }
        let axis = (0..3)
            .max_by(|&a, &b| {
                let spread = |k: usize| {
                    c.iter().map(|p| p[k]).fold(f64::NEG_INFINITY, f64::max)
                        - c.iter().map(|p| p[k]).fold(f64::INFINITY, f64::min)
                };
                spread(a).total_cmp(&spread(b))
            })
            .expect("three axes");
        let mut by = [0usize, 1, 2];
        by.sort_by(|&a, &b| c[a][axis].total_cmp(&c[b][axis]));
        let i = by[1];
        let (m, e) = (t[i], (i + 1) % 3);
        let (a, b) = (t[e], t[(e + 1) % 3]);
        let Some(g) = self.live().find(|&x| {
            let u = self.tris[x];
            (0..3).any(|k| u[k] == b && u[(k + 1) % 3] == a)
        }) else {
            return false;
        };
        let key = (self.point[m as usize], self.orig[g]);
        if done.contains(&key) || !self.insert(m, f, Some(e)) {
            return false;
        }
        done.insert(key);
        true
    }

    /// C2: vertex `h` and vertex `c`, within each other's hot box, become one:
    /// the input vertex is kept when one of them is, else `c`. Refused for
    /// copies of one exact point, when the kept coordinates would lie past
    /// [`move_cap`] from the exact position of a vertex that ends there, and
    /// when a face would hold two copies of one exact point.
    fn collapse(&mut self, h: u32, c: u32) -> bool {
        if self.point[h as usize] == self.point[c as usize] {
            return false;
        }
        let (gone, keep) = if self.rows[h as usize].is_some() && self.rows[c as usize].is_none() {
            (c, h)
        } else {
            (h, c)
        };
        let (g, k) = (self.v[gone as usize], self.v[keep as usize]);
        let span: [f64; 3] = [0, 1, 2].map(|a| {
            up(self.span[keep as usize][a].max(self.span[gone as usize][a] + (g[a] - k[a]).abs()))
        });
        if (0..3).any(|a| span[a] > move_cap(k)) {
            return false;
        }
        let pk = self.point[keep as usize];
        let joins_a_copy = self.live().any(|f| {
            let t = self.tris[f];
            t.contains(&gone)
                && t.iter()
                    .any(|&x| x != keep && x != gone && self.point[x as usize] == pk)
        });
        if joins_a_copy {
            return false;
        }
        self.span[keep as usize] = span;
        self.to[gone as usize] = keep;
        for f in 0..self.tris.len() {
            if !self.alive[f] || !self.tris[f].contains(&gone) {
                continue;
            }
            let t = self.tris[f].map(|i| if i == gone { keep } else { i });
            self.tris[f] = t;
            if t[0] == t[1] || t[1] == t[2] || t[0] == t[2] {
                self.alive[f] = false;
            }
        }
        true
    }

    /// The problems `validate` finds among the live faces whose bounding
    /// boxes meet `[lo, hi]`, of the kinds a change of faces inside the box
    /// can make or remove.
    fn local_problems(&self, lo: [f64; 3], hi: [f64; 3]) -> usize {
        use wireframes_core::MeshProblemKind as K;
        let mut id: HashMap<u32, u32> = HashMap::new();
        let mut v = vec![];
        let mut tris = vec![];
        for f in self.live() {
            let c = self.tris[f].map(|i| self.v[i as usize]);
            if !(0..3).all(|k| c.iter().any(|p| p[k] >= lo[k]) && c.iter().any(|p| p[k] <= hi[k])) {
                continue;
            }
            tris.push(self.tris[f].map(|i| {
                *id.entry(i).or_insert_with(|| {
                    v.push(self.v[i as usize]);
                    v.len() as u32 - 1
                })
            }));
        }
        let Ok(r) = Mesh::new(v, tris).and_then(|m| m.validate(0.0)) else {
            return usize::MAX;
        };
        r.problems
            .iter()
            .filter(|p| {
                matches!(
                    p.kind,
                    K::DegenerateFace
                        | K::DuplicateFace
                        | K::BoundaryEdge
                        | K::NonManifoldEdge
                        | K::InconsistentWinding
                        | K::SelfIntersection
                        | K::CoplanarOverlap
                )
            })
            .count()
    }

    /// The box of the corners of every face that differs between `self` and
    /// `next`, a copy of `self` after a local change.
    fn changed_box(&self, next: &Global) -> Option<([f64; 3], [f64; 3])> {
        let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
        let mut any = false;
        let mut take = |g: &Global, f: usize| {
            any = true;
            for i in g.tris[f] {
                for k in 0..3 {
                    lo[k] = lo[k].min(g.v[i as usize][k]);
                    hi[k] = hi[k].max(g.v[i as usize][k]);
                }
            }
        };
        for f in 0..next.tris.len() {
            let old = f < self.tris.len();
            let same = old && self.alive[f] == next.alive[f] && self.tris[f] == next.tris[f];
            if same {
                continue;
            }
            if old && self.alive[f] {
                take(self, f);
            }
            if next.alive[f] {
                take(next, f);
            }
        }
        any.then_some((lo, hi))
    }

    /// Flips edge `e` of face `f` (from corner `e` to the next) with the one
    /// face across it.
    fn flip(&mut self, f: usize, e: usize) -> bool {
        let t = self.tris[f];
        let (p, q, m) = (t[e], t[(e + 1) % 3], t[(e + 2) % 3]);
        let opp: Vec<usize> = self
            .live()
            .filter(|&x| {
                let u = self.tris[x];
                (0..3).any(|k| u[k] == q && u[(k + 1) % 3] == p)
            })
            .collect();
        let [g] = opp[..] else {
            return false;
        };
        let u = self.tris[g];
        let k = (0..3).find(|&k| u[k] == q).expect("runs the edge");
        let d = u[(k + 2) % 3];
        if d == m {
            return false;
        }
        self.alive[f] = false;
        self.alive[g] = false;
        self.push([p, d, m], f);
        self.push([d, q, m], g);
        true
    }

    /// The local changes tried for the problem pairs and flat faces of
    /// `report`: joins, collapses and insertions of a corner of one face
    /// into the other within its wide hot box, and flips of their edges.
    fn candidates(&self, report: &wireframes_core::MeshReport, ids: &[usize]) -> Vec<Change> {
        use wireframes_core::MeshProblemKind as K;
        let mut out = vec![];
        for p in &report.problems {
            let faces: Vec<usize> = match p.kind {
                K::SelfIntersection | K::CoplanarOverlap | K::DegenerateFace => [p.face, p.other]
                    .into_iter()
                    .flatten()
                    .map(|f| ids[f as usize])
                    .collect(),
                _ => continue,
            };
            for &f in &faces {
                for e in 0..3 {
                    out.push(Change::Flip(f, e));
                }
            }
            let [f, o] = faces[..] else {
                continue;
            };
            out.push(Change::Join(f, o));
            for (a, b) in [(f, o), (o, f)] {
                let tb = self.tris[b];
                let cb = tb.map(|i| self.v[i as usize]);
                for h in self.tris[a] {
                    if tb.contains(&h) {
                        continue;
                    }
                    let (lo, hi) = hot_box(self.v[h as usize], &cb, true);
                    if !box_meets(&cb, lo, hi) {
                        continue;
                    }
                    for k in 0..3 {
                        if box_meets(&[cb[k]], lo, hi) {
                            out.push(Change::Collapse(h, tb[k]));
                        }
                        if box_meets(&[cb[k], cb[(k + 1) % 3]], lo, hi) {
                            out.push(Change::Insert(h, b, Some(k)));
                        }
                    }
                    out.push(Change::Insert(h, b, None));
                }
            }
        }
        let mut seen = HashSet::new();
        out.retain(|c| seen.insert(*c));
        out
    }

    fn apply(&mut self, c: Change) -> bool {
        match c {
            Change::Flip(f, e) => self.alive[f] && self.flip(f, e),
            Change::Join(f, o) => self.join_opposite(f, o),
            Change::Collapse(h, c) => self.collapse(h, c),
            Change::Insert(h, g, e) => self.alive[g] && self.insert(h, g, e),
        }
    }

    /// Applies the change among [`Global::candidates`] that, with the
    /// cancellations it allows, lowers the count of problems most; the first
    /// such in candidate order on ties, and returns it; `None`, changing
    /// nothing, when none lowers it. Each applied change lowers the count of
    /// problems of the whole mesh by at least one, since faces outside the
    /// changed box and their pairs are the same before and after.
    fn descend(&mut self, report: &wireframes_core::MeshReport, ids: &[usize]) -> Option<Change> {
        let mut best: Option<(usize, Change, Global)> = None;
        for c in self.candidates(report, ids) {
            let mut next = self.clone();
            if !next.apply(c) {
                continue;
            }
            cancel_pairs(&next.tris, &mut next.alive);
            let Some((lo, hi)) = self.changed_box(&next) else {
                continue;
            };
            let (before, after) = (self.local_problems(lo, hi), next.local_problems(lo, hi));
            if after < before && best.as_ref().is_none_or(|(d, ..)| before - after > *d) {
                best = Some((before - after, c, next));
            }
        }
        let (_, c, next) = best?;
        *self = next;
        Some(c)
    }

    fn push(&mut self, t: [u32; 3], like: usize) {
        self.tris.push(t);
        self.alive
            .push(t[0] != t[1] && t[1] != t[2] && t[0] != t[2]);
        self.src.push(self.src[like]);
        if !self.nrm.is_empty() {
            self.nrm.push(self.nrm[like]);
        }
        self.orig.push(self.orig[like]);
    }

    /// Inserts vertex `h` into face `g`: on its edge `e` (from corner `e` to
    /// the next), splitting the one other face that runs that edge the other
    /// way too, or inside `g` when `e` is `None`. Faces that then repeat a
    /// vertex go. False when the edge has no single opposite face.
    fn insert(&mut self, h: u32, g: usize, e: Option<usize>) -> bool {
        let t = self.tris[g];
        let ph = self.point[h as usize];
        let copy = |x: u32| x != h && self.point[x as usize] == ph;
        if t.iter().any(|&x| copy(x)) {
            return false;
        }
        let twin = |g2: &Global, n: [u32; 3]| {
            g2.live().any(|x| {
                let u = g2.tris[x];
                (0..3).any(|k| u == [n[k], n[(k + 1) % 3], n[(k + 2) % 3]])
            })
        };
        match e {
            None => {
                if (0..3).any(|k| twin(self, [t[k], t[(k + 1) % 3], h])) {
                    return false;
                }
                self.alive[g] = false;
                for k in 0..3 {
                    self.push([t[k], t[(k + 1) % 3], h], g);
                }
            }
            Some(e) => {
                let (a, b, c) = (t[e], t[(e + 1) % 3], t[(e + 2) % 3]);
                let opp: Vec<usize> = self
                    .live()
                    .filter(|&x| {
                        let u = self.tris[x];
                        (0..3).any(|k| u[k] == b && u[(k + 1) % 3] == a)
                    })
                    .collect();
                if opp.len() != 1 {
                    return false;
                }
                let g2 = opp[0];
                let u = self.tris[g2];
                let k = (0..3).find(|&k| u[k] == b).expect("runs the edge");
                let d = u[(k + 2) % 3];
                if copy(d)
                    || [[a, h, c], [h, b, c], [b, h, d], [h, a, d]]
                        .into_iter()
                        .any(|n| n[0] != n[2] && n[1] != n[2] && twin(self, n))
                {
                    return false;
                }
                self.alive[g] = false;
                self.alive[g2] = false;
                self.push([a, h, c], g);
                self.push([h, b, c], g);
                self.push([b, h, d], g2);
                self.push([h, a, d], g2);
            }
        }
        true
    }
}

/// The step bound of `docs/design/rounding.md` for `v0` vertices and `f0`
/// faces: the per-shell cleanup flips each edge it did not make at most once
/// and makes at most `3 f0 / 2 + v0 f0` such edges per round, with at most
/// `f0^2 + 1` rounds; merges, collapses and cancellations each remove a
/// vertex or a face.
fn step_bound(v0: usize, f0: usize) -> usize {
    let flips = (3 * f0 / 2 + 1).saturating_add(v0.saturating_mul(f0));
    let rounds = f0.saturating_mul(f0).saturating_add(1);
    flips
        .saturating_mul(rounds)
        .saturating_add(v0)
        .saturating_add(f0)
        .saturating_add(rounds)
}

/// A rounded result, its report, and the output vertex each draft vertex
/// ended at (`None` when its faces are gone).
pub(crate) struct Rounded {
    pub mesh: Assembled,
    pub report: RoundingReport,
    /// The mesh is a solid reached within the move cap (or the draft is not
    /// closed); otherwise it is the closed mesh of the cleanup without the cap,
    /// for the error to carry.
    pub certified: bool,
    #[cfg_attr(not(test), allow(dead_code))]
    pub ends: Vec<Option<u32>>,
}

/// The removals of flat runs over all shells for `v` vertices and `f` faces
/// that the loop allows: it counts them and takes none past this.
fn zip_bound(v: usize, f: usize) -> usize {
    f.saturating_add(v.saturating_mul(f).saturating_mul(4))
}

/// The descents over all shells for `v` vertices and `f` faces that the loop
/// allows: it counts them and takes none past this.
fn descent_bound(v: usize, f: usize) -> usize {
    v.saturating_add(f).saturating_mul(4)
}

/// The bound on the operations over all shells for `v` vertices and `f`
/// faces: at most `i = v f` insertions (one per exact point and original
/// face), each adding at most two faces; merges, collapses, cancellations and
/// removals each remove a vertex or a face of the at most `v + f + 2 i` there
/// ever are; at most [`zip_bound`] removals of flat runs, a count the loop
/// enforces; one round per insertion, collapse or removal plus the last; and
/// per descent (at most [`descent_bound`], a count the loop enforces) its
/// change, its round, and the removal of the at most two faces it adds.
fn global_step_bound(v: usize, f: usize) -> usize {
    let i = v.saturating_mul(f);
    i.saturating_mul(9)
        .saturating_add(v.saturating_mul(3))
        .saturating_add(f.saturating_mul(2))
        .saturating_add(descent_bound(v, f).saturating_mul(4))
        .saturating_add(1)
}

/// Rounds a draft to an `f64` mesh with the report of what it did. When the
/// operations cannot certify a closed draft, the result is the per-shell
/// cleanup without the move cap, a closed mesh for the error to carry, with
/// `certified` false and its report's `move_bound` raised to its `max_move`.
pub(crate) fn materialize(d: Draft, descend: bool) -> Result<Rounded> {
    if d.home.iter().flatten().any(|x| !x.is_finite()) {
        return Err(Error::kernel(
            "non-finite: an exact coordinate of the result rounds outside the f64 range",
        ));
    }
    let (r, solid) = run(&d, true, descend)?;
    if solid || !d.closed {
        return Ok(Rounded {
            certified: true,
            ..r
        });
    }
    let (mut r, _) = run(&d, false, false)?;
    for k in 0..3 {
        r.report.move_bound[k] = r.report.move_bound[k].max(r.report.max_move[k]);
    }
    Ok(Rounded {
        certified: false,
        ..r
    })
}

/// One rounding of `d`: the per-shell cleanup, with the move cap when `cap`,
/// then, with `cap` and a closed draft that is not yet a solid, the operations
/// over all shells. Also returns whether the mesh is a solid.
fn run(d: &Draft, cap: bool, descend: bool) -> Result<(Rounded, bool)> {
    let nv = d.home.len();
    let rho: Vec<f64> = d.err.iter().map(|e| up(norm(*e))).collect();
    let raw: Vec<[[f64; 3]; 3]> = d
        .triangles
        .iter()
        .map(|t| t.map(|i| d.home[i as usize]))
        .collect();
    let rounding_bound: f64 = d
        .triangles
        .iter()
        .zip(&raw)
        .map(|(t, c)| {
            tube(
                *c,
                2.0 * t.iter().map(|&i| rho[i as usize]).fold(0.0, f64::max),
            )
        })
        .sum();
    let mut report = RoundingReport {
        step_bound: step_bound(nv, d.triangles.len()),
        ..RoundingReport::default()
    };
    let mut g = Global {
        v: vec![],
        point: vec![],
        rows: vec![],
        tris: vec![],
        src: vec![],
        nrm: vec![],
        alive: vec![],
        to: vec![],
        orig: vec![],
        span: vec![],
    };
    let mut origin: Vec<Option<u32>> = vec![None; nv];
    let mut err_g: Vec<f64> = vec![];
    for faces in &d.shells {
        let mut local: HashMap<u32, u32> = HashMap::new();
        let mut verts: Vec<u32> = vec![];
        let tris: Vec<[u32; 3]> = faces
            .iter()
            .map(|&f| {
                d.triangles[f as usize].map(|i| {
                    *local.entry(i).or_insert_with(|| {
                        verts.push(i);
                        verts.len() as u32 - 1
                    })
                })
            })
            .collect();
        let shell = Assembled {
            vertices: verts.iter().map(|&i| d.home[i as usize]).collect(),
            triangles: tris,
            vertex_rows: verts.iter().map(|&i| d.rows[i as usize]).collect(),
            face_source: faces.iter().map(|&f| d.face_source[f as usize]).collect(),
            normal: if d.normal.is_empty() {
                vec![]
            } else {
                faces.iter().map(|&f| d.normal[f as usize]).collect()
            },
        };
        let point: Vec<u32> = verts.iter().map(|&i| d.point[i as usize]).collect();
        let err: Vec<[f64; 3]> = verts.iter().map(|&i| d.err[i as usize]).collect();
        let (out, kept, fwd, c) = collapse_rounded(shell, &point, &err, cap);
        report.merged_vertices += c.merged;
        report.collapsed_edges += c.collapsed;
        report.cancelled_faces += c.cancelled;
        report.steps += c.merged + c.collapsed + c.cancelled / 2 + c.flips + c.rounds;
        let base = g.v.len() as u32;
        for (k, &i) in verts.iter().enumerate() {
            origin[i as usize] = fwd[k].map(|o| o + base);
        }
        let mut best: Vec<f64> = vec![0.0; out.vertices.len()];
        for (k, &i) in verts.iter().enumerate() {
            if let Some(o) = fwd[k] {
                best[o as usize] = best[o as usize].max(rho[i as usize]);
            }
        }
        err_g.extend(best);
        g.v.extend(out.vertices);
        g.point.extend(kept);
        g.rows.extend(out.vertex_rows);
        g.tris
            .extend(out.triangles.iter().map(|t| t.map(|i| i + base)));
        g.alive.extend(out.triangles.iter().map(|_| true));
        g.src.extend(out.face_source);
        g.nrm.extend(out.normal);
    }
    g.to = (0..g.v.len() as u32).collect();
    g.orig = (0..g.tris.len() as u32).collect();
    g.span = vec![[0.0; 3]; g.v.len()];
    for (i, o) in origin.iter().enumerate() {
        if let Some(o) = *o {
            let at = g.v[o as usize];
            for k in 0..3 {
                let s = up((at[k] - d.home[i][k]).abs() + d.err[i][k]);
                g.span[o as usize][k] = g.span[o as usize][k].max(s);
            }
        }
    }
    report.step_bound = report
        .step_bound
        .saturating_add(global_step_bound(g.v.len(), g.tris.len()));
    if cap && d.closed && !is_solid(&g)? {
        let before = (g.clone(), report.clone());
        let mut done: HashSet<(u32, u32)> = HashSet::new();
        let (mut zips, zip_limit) = (0usize, zip_bound(g.v.len(), g.tris.len()));
        let (mut descents, descent_limit) = (0usize, descent_bound(g.v.len(), g.tris.len()));
        loop {
            let merged = g.merge_equal();
            report.merged_vertices += merged;
            let cancelled = cancel_pairs(&g.tris, &mut g.alive);
            report.cancelled_faces += cancelled;
            let removed = g.remove_flat_shells(&err_g, descend)?;
            report.removed_shells += removed;
            report.steps += merged + cancelled / 2 + removed + 1;
            let (m, ids) = g.mesh()?;
            if ids.is_empty() {
                break;
            }
            let found = m.validate(0.0)?;
            if found.summary.is_solid() {
                break;
            }
            if descents == 0 && zips < zip_limit {
                let seeds: HashSet<usize> = found
                    .problems
                    .iter()
                    .flat_map(|p| [p.face, p.other])
                    .flatten()
                    .map(|f| ids[f as usize])
                    .collect();
                let zipped = g.zip_flat_runs(&seeds);
                if zipped > 0 {
                    zips += 1;
                    report.cancelled_faces += zipped;
                    report.steps += 1;
                    continue;
                }
            }
            let (mut inserted, mut collapsed) = (0, 0);
            if descents == 0 {
                (inserted, collapsed) = g.insert_hot_points(&found, &ids, &mut done, false);
                if inserted + collapsed == 0 {
                    (inserted, collapsed) = g.insert_hot_points(&found, &ids, &mut done, true);
                }
            }
            report.inserted_points += inserted;
            report.collapsed_edges += collapsed;
            report.steps += inserted + collapsed;
            if inserted + collapsed == 0 {
                if descend
                    && descents < descent_limit
                    && let Some(c) = g.descend(&found, &ids)
                {
                    descents += 1;
                    report.steps += 1;
                    match c {
                        Change::Flip(..) => report.flipped_edges += 1,
                        Change::Join(..) | Change::Collapse(..) => report.collapsed_edges += 1,
                        Change::Insert(..) => report.inserted_points += 1,
                    }
                    continue;
                }
                break;
            }
        }
        if !is_solid(&g)? {
            (g, report) = before;
        }
    }
    assert!(
        report.steps <= report.step_bound,
        "rounding took {} steps, past its bound {}",
        report.steps,
        report.step_bound
    );
    let resolve = |mut i: u32| {
        while g.to[i as usize] != i {
            i = g.to[i as usize];
        }
        i
    };
    let mut out = Assembled {
        vertices: vec![],
        triangles: vec![],
        vertex_rows: vec![],
        face_source: vec![],
        normal: vec![],
    };
    let mut new_id: Vec<Option<u32>> = vec![None; g.v.len()];
    for f in g.live().collect::<Vec<_>>() {
        out.triangles.push(g.tris[f].map(|i| {
            *new_id[i as usize].get_or_insert_with(|| {
                out.vertices.push(g.v[i as usize]);
                out.vertex_rows.push(g.rows[i as usize]);
                out.vertices.len() as u32 - 1
            })
        }));
        out.face_source.push(g.src[f]);
        if !g.nrm.is_empty() {
            out.normal.push(g.nrm[f]);
        }
    }
    report.rounded_vertices = out.vertex_rows.iter().filter(|r| r.is_none()).count();
    let mut moved: HashSet<u32> = HashSet::new();
    let mut moved_input: HashSet<u32> = HashSet::new();
    let ends: Vec<Option<u32>> = (0..nv)
        .map(|i| origin[i].map(resolve).and_then(|o| new_id[o as usize]))
        .collect();
    for i in 0..nv {
        let Some(o) = ends[i] else {
            continue;
        };
        let at = out.vertices[o as usize];
        for k in 0..3 {
            let dist = up((at[k] - d.home[i][k]).abs() + d.err[i][k]);
            report.max_move[k] = report.max_move[k].max(dist);
        }
        if at != d.home[i] {
            moved.insert(d.point[i]);
            if d.rows[i].is_some() {
                moved_input.insert(d.point[i]);
            }
        }
    }
    report.moved_vertices = moved.len();
    report.moved_input_vertices = moved_input.len();
    for p in &out.vertices {
        for k in 0..3 {
            report.move_bound[k] = report.move_bound[k].max(move_cap(*p));
        }
    }
    let after: Vec<[[f64; 3]; 3]> = out
        .triangles
        .iter()
        .map(|t| t.map(|i| out.vertices[i as usize]))
        .collect();
    report.volume_bound = up(rounding_bound + volume_change(&raw, &after));
    let solid = is_solid(&g)?;
    Ok((
        Rounded {
            mesh: out,
            report,
            certified: solid,
            ends,
        },
        solid,
    ))
}

fn is_solid(g: &Global) -> Result<bool> {
    let (m, ids) = g.mesh()?;
    Ok(ids.is_empty() || m.validate(0.0)?.summary.is_solid())
}

/// `mesh` when the rounding certified it and `validate` passes it as a solid,
/// else the `uncertified` kernel error carrying the report and the mesh.
pub(crate) fn certify(mesh: Mesh, certified: bool, what: &str) -> Result<Mesh> {
    if certified && (mesh.is_empty() || mesh.validate(0.0)?.summary.is_solid()) {
        return Ok(mesh);
    }
    Err(Error::kernel_result(
        format!(
            "uncertified: rounding {what} to f64 left a mesh that is not a solid, and the \
             rounding operations could not certify it"
        ),
        mesh,
    )?)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ops::{is_closed_manifold, solid_result};

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
        let (out, _, _, _) = collapse_rounded(s, &[0, 1, 2, 3, 4], &[[0.0; 3]; 5], true);
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
        let (out, _, _, _) = collapse_rounded(s, &[0, 1, 2, 3, 4], &[[0.0; 3]; 5], true);
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
        let Err(Error::Kernel { report, mesh, .. }) = solid_result(m.clone(), true, "it") else {
            panic!("an open result must be a kernel error");
        };
        assert!(!report.unwrap().summary.is_closed);
        assert_eq!(*mesh.unwrap(), m);
    }

    #[test]
    fn copies_of_one_exact_point_stay_apart() {
        let s = shell(vec![A, B, C, D, D], vec![[0, 2, 1], [0, 1, 3], [1, 2, 4]]);
        let (out, _, _, _) = collapse_rounded(s, &[0, 1, 2, 3, 3], &[[0.0; 3]; 5], true);
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
        let (out, _, _, _) = collapse_rounded(s, &[0, 1, 2, 3, 4], &[[0.0; 3]; 5], true);
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
        let (out, kept, _, _) = collapse_rounded(s, &point, &vec![[0.0; 3]; point.len()], true);
        assert!(is_closed_manifold(&out.triangles));
        assert_eq!(kept.len(), out.vertices.len());
    }
}

/// Exact checks of a rounding's report against the exact result, with
/// rationals.
#[cfg(test)]
pub(crate) mod checks {
    use num_rational::BigRational as Q;

    use super::*;
    use crate::arrangement::Arrangement;
    use crate::hp;
    use crate::ops::{self, Shell};

    pub(crate) fn pair(name: &str) -> (Mesh, Mesh) {
        let path = format!(
            "{}/tests/data/rounding/{name}.txt",
            env!("CARGO_MANIFEST_DIR")
        );
        let text = std::fs::read_to_string(&path).expect("reduced case");
        let mut lines = text.lines();
        let mut read = || {
            let head: Vec<usize> = lines
                .next()
                .expect("header")
                .split(' ')
                .skip(1)
                .map(|x| x.parse().expect("count"))
                .collect();
            let row = |l: &str| -> Vec<f64> {
                l.split(' ').map(|x| x.parse().expect("number")).collect()
            };
            let v: Vec<[f64; 3]> = (0..head[0])
                .map(|_| {
                    let p = row(lines.next().expect("vertex"));
                    [p[0], p[1], p[2]]
                })
                .collect();
            let t: Vec<[u32; 3]> = (0..head[1])
                .map(|_| {
                    let p = row(lines.next().expect("face"));
                    [p[0] as u32, p[1] as u32, p[2] as u32]
                })
                .collect();
            Mesh::new(v, t).expect("mesh")
        };
        let a = read();
        (a, read())
    }

    fn q(x: f64) -> Q {
        Q::from_float(x).expect("finite")
    }

    fn det(a: &[Q; 3], b: &[Q; 3], c: &[Q; 3]) -> Q {
        &a[0] * (&b[1] * &c[2] - &b[2] * &c[1])
            + &a[1] * (&b[2] * &c[0] - &b[0] * &c[2])
            + &a[2] * (&b[0] * &c[1] - &b[1] * &c[0])
    }

    fn abs(x: Q) -> Q {
        if x < q(0.0) { -x } else { x }
    }

    /// Rounds the exact result of `pieces` and checks, in exact arithmetic,
    /// that every exact vertex ends within `max_move <= move_bound` of its
    /// exact position, within half an ulp unless counted as moved, and that
    /// `|vol(out) - vol(exact)| <= volume_bound`.
    pub(crate) fn rounded_exactly(arr: &Arrangement, pieces: &[Shell]) -> Rounded {
        let d = ops::draft(arr, pieces);
        let exact: Vec<[Q; 3]> = d
            .point
            .iter()
            .map(|&p| hp::rational(&arr.pts[p as usize]))
            .collect();
        let six = q(6.0);
        let vol_exact = d.triangles.iter().fold(q(0.0), |acc, t| {
            acc + det(
                &exact[t[0] as usize],
                &exact[t[1] as usize],
                &exact[t[2] as usize],
            )
        }) / &six;
        let (home, rows, point) = (d.home.clone(), d.rows.clone(), d.point.clone());
        let r = materialize(d, false).expect("finite");
        let out = &r.mesh.vertices;
        let vol_out = r.mesh.triangles.iter().fold(q(0.0), |acc, t| {
            let c = t.map(|i| out[i as usize].map(q));
            acc + det(&c[0], &c[1], &c[2])
        }) / &six;
        assert!(
            abs(vol_out - vol_exact) <= q(r.report.volume_bound),
            "volume bound {}",
            r.report.volume_bound
        );
        let mut moved = HashSet::new();
        let mut moved_input = HashSet::new();
        for (i, e) in r.ends.iter().enumerate() {
            let Some(o) = e else { continue };
            let at = out[*o as usize];
            for k in 0..3 {
                let dist = abs(q(at[k]) - &exact[i][k]);
                assert!(dist <= q(r.report.max_move[k]), "max_move axis {k}");
                assert!(r.report.max_move[k] <= r.report.move_bound[k]);
                if at == home[i] {
                    assert!(
                        dist <= q(half_ulp(at[k])),
                        "unmoved vertex past half an ulp"
                    );
                }
            }
            if at != home[i] {
                moved.insert(point[i]);
                if rows[i].is_some() {
                    moved_input.insert(point[i]);
                }
            }
        }
        assert_eq!(r.report.moved_vertices, moved.len());
        assert_eq!(r.report.moved_input_vertices, moved_input.len());
        assert!(r.report.steps <= r.report.step_bound);
        r
    }

    pub(crate) const REDUCED: [&str; 9] = [
        "clipped_boxes",
        "pit_pushback",
        "touching_slivers",
        "runaway_union",
        "runaway_intersection",
        "runaway_difference",
        "runaway_difference_reversed",
        "rotated_stopes",
        "crosscut_chain",
    ];

    fn ulp_prisms() -> (Mesh, Mesh) {
        let faces = vec![
            [0, 2, 1],
            [1, 2, 3],
            [4, 5, 6],
            [5, 7, 6],
            [0, 1, 4],
            [1, 5, 4],
            [2, 6, 3],
            [3, 6, 7],
            [0, 4, 2],
            [2, 4, 6],
            [1, 3, 5],
            [3, 7, 5],
        ];
        let prism = |half: f64, bottom: f64| {
            let mut v = vec![];
            for top in [false, true] {
                for y in [-half, half] {
                    for x in [-half, half] {
                        v.push([x, y, if top { 5.0 + x / 4.0 + y / 2.0 } else { bottom }]);
                    }
                }
            }
            v
        };
        let a = Mesh::new(prism(10.0, -10.0), faces.clone()).expect("prism");
        let mut v = prism(6.0, -5.0);
        v[4][2] = v[4][2].next_up();
        (a, Mesh::new(v, faces).expect("prism"))
    }

    fn booleans_round_within_their_bounds(a: &Mesh, b: &Mesh) {
        let arr = crate::arrangement::build(a, b).expect("arrangement");
        for op in [ops::Op::Union, ops::Op::Intersection, ops::Op::Difference] {
            let kept = ops::kept(&arr, op).expect("classified");
            rounded_exactly(&arr, &ops::shells(&arr, &kept).expect("shells"));
        }
    }

    #[test]
    fn prisms_whose_caps_differ_by_one_ulp_round_within_their_bounds() {
        let (a, b) = ulp_prisms();
        booleans_round_within_their_bounds(&a, &b);
        booleans_round_within_their_bounds(&b, &a);
    }

    #[test]
    fn reduced_cases_round_within_their_bounds() {
        for name in REDUCED {
            let (a, b) = pair(name);
            booleans_round_within_their_bounds(&a, &b);
        }
    }

    /// Six times the exact volume of the exact result of `pieces`.
    fn exact_volume6(arr: &Arrangement, pieces: &[Shell]) -> Q {
        let d = ops::draft(arr, pieces);
        let exact: Vec<[Q; 3]> = d
            .point
            .iter()
            .map(|&p| hp::rational(&arr.pts[p as usize]))
            .collect();
        d.triangles.iter().fold(q(0.0), |acc, t| {
            acc + det(
                &exact[t[0] as usize],
                &exact[t[1] as usize],
                &exact[t[2] as usize],
            )
        })
    }

    fn volume6(m: &Mesh) -> Q {
        m.triangles().iter().fold(q(0.0), |acc, t| {
            let c = t.map(|i| m.vertices()[i as usize].map(q));
            acc + det(&c[0], &c[1], &c[2])
        })
    }

    /// The public booleans of every reduced case, conforming included, are
    /// solids within `volume_bound` of the exact result of the inputs as
    /// given, and leave the inputs as they were.
    #[test]
    fn public_booleans_of_reduced_cases_stay_within_their_bounds() {
        let mut cases: Vec<(Mesh, Mesh)> = REDUCED.iter().map(|n| pair(n)).collect();
        cases.push(ulp_prisms());
        for (a, b) in cases {
            let (a0, b0) = (a.clone(), b.clone());
            let arr = crate::arrangement::build(&a, &b).expect("arrangement");
            for op in [ops::Op::Union, ops::Op::Intersection, ops::Op::Difference] {
                let kept = ops::kept(&arr, op).expect("classified");
                let exact = exact_volume6(&arr, &ops::shells(&arr, &kept).expect("shells"));
                let r = ops::boolean(&a, &b, op, &ops::BooleanOptions::default()).expect("solid");
                assert!(
                    r.mesh.is_empty() || r.mesh.validate(0.0).expect("report").summary.is_solid()
                );
                let gap = abs(volume6(&r.mesh) - exact);
                assert!(gap <= q(6.0 * r.rounding.volume_bound), "{op:?}");
                assert!(r.rounding.steps <= r.rounding.step_bound);
                for k in 0..3 {
                    assert!(r.rounding.max_move[k] <= r.rounding.move_bound[k]);
                }
            }
            assert_eq!((a, b), (a0, b0));
        }
    }

    #[test]
    fn reduced_clip_rounds_within_its_bounds() {
        let (pit, ground) = pair("square_pit_tilted_ground");
        let arr = crate::arrangement::build(&pit, &ground).expect("arrangement");
        for piece in crate::clip::pieces(&arr).expect("pieces") {
            let r = rounded_exactly(&arr, &ops::shells(&arr, &piece).expect("shells"));
            let m = Mesh::new(r.mesh.vertices, r.mesh.triangles).expect("mesh");
            assert!(m.validate(0.0).expect("report").summary.is_solid());
        }
    }
}

#[cfg(test)]
mod operations {
    use super::*;
    use crate::ops::is_closed_manifold;

    const A: [f64; 3] = [0.0, 0.0, 0.0];
    const B: [f64; 3] = [2.0, 0.0, 0.0];
    const C: [f64; 3] = [0.0, 2.0, 0.0];
    const D: [f64; 3] = [0.0, 0.0, 2.0];
    const TETRA: [[u32; 3]; 4] = [[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]];

    fn global(v: Vec<[f64; 3]>, tris: Vec<[u32; 3]>, point: Vec<u32>) -> Global {
        let n = v.len();
        let f = tris.len();
        Global {
            span: vec![[0.0; 3]; n],
            rows: (0..n as u32).map(|i| (i < 4).then_some(i)).collect(),
            to: (0..n as u32).collect(),
            point,
            v,
            src: vec![(0, 0); f],
            nrm: vec![],
            alive: vec![true; f],
            orig: (0..f as u32).collect(),
            tris,
        }
    }

    fn live(g: &Global) -> Vec<[u32; 3]> {
        g.live().map(|f| g.tris[f]).collect()
    }

    /// A tetrahedron whose bottom face has a rounded point `E` an ulp across
    /// edge `AB`, joined to the three bottom corners: the sliver `B A E`
    /// folds over the other two faces. Each descent lowers the problem count
    /// and the last leaves a solid.
    #[test]
    fn descents_lower_the_problem_count_until_none_does() {
        let e = [1.0, -1e-16, 0.0];
        let mut g = global(
            vec![A, B, C, D, e],
            vec![
                [0, 1, 3],
                [0, 3, 2],
                [1, 2, 3],
                [0, 2, 4],
                [2, 1, 4],
                [1, 0, 4],
            ],
            vec![0, 1, 2, 3, 4],
        );
        let problems = |g: &Global| {
            let (m, ids) = g.mesh().unwrap();
            (m.validate(0.0).unwrap(), ids)
        };
        let (mut found, mut ids) = problems(&g);
        assert!(!found.summary.is_solid());
        let mut steps = 0;
        let mut count = found.problems.len();
        while g.descend(&found, &ids).is_some() {
            cancel_pairs(&g.tris, &mut g.alive);
            (found, ids) = problems(&g);
            assert!(found.problems.len() < count);
            count = found.problems.len();
            steps += 1;
            assert!(steps <= descent_bound(5, 6));
        }
        assert!(found.summary.is_solid(), "{:?}", found.summary);
        assert!(is_closed_manifold(&live(&g)));
        assert!(g.descend(&found, &ids).is_none());
    }

    #[test]
    fn boxes_meet_triangles_and_segments_exactly() {
        let t = [A, B, C];
        assert!(box_meets(&t, [0.5, 0.5, -1.0], [0.6, 0.6, 1.0]));
        assert!(!box_meets(&t, [0.5, 0.5, 0.1], [0.6, 0.6, 1.0]));
        assert!(box_meets(&t, [1.0, 1.0, 0.0], [2.0, 2.0, 0.0]));
        assert!(!box_meets(
            &t,
            [1.0, 1.0, 0.0].map(f64::next_up),
            [2.0, 2.0, 0.0]
        ));
        assert!(box_meets(&[A, D], [-1.0, -1.0, 1.0], [0.0, 0.0, 1.0]));
        assert!(!box_meets(
            &[A, D],
            [0.0f64.next_up(), -1.0, 0.0],
            [1.0, 1.0, 1.0]
        ));
    }

    #[test]
    fn a_point_inserted_inside_a_face_keeps_the_shell_closed() {
        let e = [0.5, 0.5, 0.0];
        let mut g = global(vec![A, B, C, D, e], TETRA.to_vec(), vec![0, 1, 2, 3, 4]);
        assert!(g.insert(4, 0, None));
        assert!(is_closed_manifold(&live(&g)));
        assert_eq!(live(&g).len(), 6);
        assert!(g.orig[4..].iter().all(|&o| o == 0));
    }

    #[test]
    fn a_point_inserted_on_an_edge_splits_both_faces() {
        let e = [1.0, 0.0, 0.0];
        let mut g = global(vec![A, B, C, D, e], TETRA.to_vec(), vec![0, 1, 2, 3, 4]);
        assert!(g.insert(4, 1, Some(0)));
        assert!(is_closed_manifold(&live(&g)));
        assert_eq!(live(&g).len(), 6);
        assert!(live(&g).iter().all(|t| !(t.contains(&0) && t.contains(&1))));
    }

    #[test]
    fn inserting_the_far_corner_on_an_edge_flips_it() {
        let m = [1.0, 1.0, 0.0];
        let tris = vec![
            [0, 2, 1],
            [1, 2, 4],
            [0, 1, 3],
            [0, 3, 2],
            [1, 4, 3],
            [4, 2, 3],
        ];
        let mut g = global(vec![A, B, C, D, m], tris, vec![0, 1, 2, 3, 4]);
        assert!(g.insert(4, 0, Some(1)));
        assert!(is_closed_manifold(&live(&g)));
        assert_eq!(live(&g).len(), 6);
    }

    #[test]
    fn collapses_keep_input_vertices_and_refuse_copies_and_long_moves() {
        let near = [B[0].next_up(), B[1], B[2]];
        let tris = vec![[0, 2, 4], [0, 4, 3], [0, 3, 2], [4, 2, 3]];
        let mut g = global(vec![A, B, C, D, near], tris.clone(), vec![0, 1, 2, 3, 4]);
        assert!(g.collapse(4, 1));
        assert_eq!(g.to[4], 1);
        assert!(is_closed_manifold(&live(&g)));
        assert!(live(&g).iter().all(|t| !t.contains(&4)));
        let mut g = global(vec![A, B, C, D, near], tris.clone(), vec![0, 1, 2, 3, 1]);
        assert!(!g.collapse(4, 1));
        let far = [2.001, 0.0, 0.0];
        let mut g = global(vec![A, B, C, D, far], tris, vec![0, 1, 2, 3, 4]);
        assert!(!g.collapse(4, 1));
    }

    #[test]
    fn the_step_bounds_grow_polynomially_and_saturate() {
        assert!(global_step_bound(10, 20) >= 5 * 10 * 20);
        assert_eq!(global_step_bound(usize::MAX, 2), usize::MAX);
        assert!(step_bound(10, 20) > 20 * 20);
    }
}
