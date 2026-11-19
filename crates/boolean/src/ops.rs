use std::collections::HashMap;

use wireframes_core::predicates::{Sign, collinear3};
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

/// Repairs `mesh` when asked, then refuses it unless it is a valid solid.
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
    if !report.summary.is_solid() {
        return Err(Error::InvalidMesh {
            context: context.to_string(),
            report: Box::new(report),
        });
    }
    Ok((m, r))
}

/// One output piece: sub-triangle ids with a flip flag. Pieces get separate
/// vertex copies, so two pieces may touch without making the output
/// non-manifold.
pub(crate) type Piece = Vec<(u32, bool)>;

pub(crate) struct Assembled {
    pub vertices: Vec<[f64; 3]>,
    pub triangles: Vec<[u32; 3]>,
    /// Row of the concatenated input vertex tables, per output vertex.
    pub vertex_rows: Vec<Option<u32>>,
    /// (mesh, face) per output triangle.
    pub face_source: Vec<(u8, u32)>,
}

pub(crate) fn assemble(arr: &Arrangement, pieces: &[Piece]) -> Assembled {
    let mut out = Assembled {
        vertices: vec![],
        triangles: vec![],
        vertex_rows: vec![],
        face_source: vec![],
    };
    for piece in pieces {
        let mut shell = Assembled {
            vertices: vec![],
            triangles: vec![],
            vertex_rows: vec![],
            face_source: vec![],
        };
        let mut map: HashMap<u32, u32> = HashMap::new();
        for &(s, flip) in piece {
            let sub = &arr.subs[s as usize];
            let v = sub.v.map(|p| {
                *map.entry(p).or_insert_with(|| {
                    shell.vertices.push(arr.pts[p as usize].to_f64());
                    shell.vertex_rows.push(arr.source_row(p));
                    shell.vertices.len() as u32 - 1
                })
            });
            shell
                .triangles
                .push(if flip { [v[0], v[2], v[1]] } else { v });
            let f = &arr.faces[sub.face as usize];
            shell.face_source.push((f.mesh, f.index));
        }
        let shell = collapse_rounded(shell);
        let base = out.vertices.len() as u32;
        out.vertices.extend(shell.vertices);
        out.vertex_rows.extend(shell.vertex_rows);
        out.triangles
            .extend(shell.triangles.iter().map(|t| t.map(|i| i + base)));
        out.face_source.extend(shell.face_source);
    }
    out
}

/// One shell after rounding its exact vertices to `f64`. Distinct exact points
/// can round to the same coordinates, which makes faces with a repeated or
/// collinear corner. Vertices of the shell with equal coordinates become one
/// (an edge between them has zero length), faces that repeat a vertex are
/// dropped, opposite faces on the same three vertices cancel, and a face whose
/// corners are collinear is removed by splitting its neighbour across its
/// longest edge at its middle corner. Each step keeps every edge paired, and
/// none moves a point, so the surface is the same set of points.
fn collapse_rounded(shell: Assembled) -> Assembled {
    let n = shell.vertices.len();
    let key = |i: usize| shell.vertices[i].map(|x| (x + 0.0).to_bits());
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by_key(|&i| (key(i), shell.vertex_rows[i].is_none(), i));
    let mut rep: Vec<u32> = (0..n as u32).collect();
    for w in 1..n {
        if key(order[w]) == key(order[w - 1]) {
            rep[order[w]] = rep[order[w - 1]];
        }
    }
    let mut tris: Vec<[u32; 3]> = vec![];
    let mut src: Vec<(u8, u32)> = vec![];
    for (t, &s) in shell.triangles.iter().zip(&shell.face_source) {
        let t = t.map(|i| rep[i as usize]);
        if t[0] != t[1] && t[1] != t[2] && t[0] != t[2] {
            tris.push(t);
            src.push(s);
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
    let mut edge: HashMap<(u32, u32), usize> = HashMap::new();
    for (f, t) in tris.iter().enumerate().filter(|(f, _)| alive[*f]) {
        for k in 0..3 {
            edge.insert((t[k], t[(k + 1) % 3]), f);
        }
    }
    let mut work: Vec<usize> = (0..tris.len())
        .filter(|&f| alive[f] && flat(tris[f]))
        .rev()
        .collect();
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
            alive.push(true);
            for k in 0..3 {
                edge.insert((nt[k], nt[(k + 1) % 3]), h);
            }
            if flat(nt) {
                work.push(h);
            }
        }
    }
    let mut out = Assembled {
        vertices: vec![],
        triangles: vec![],
        vertex_rows: vec![],
        face_source: vec![],
    };
    let mut new_id: Vec<Option<u32>> = vec![None; n];
    for (f, t) in tris.iter().enumerate().filter(|(f, _)| alive[*f]) {
        out.triangles.push(t.map(|i| {
            *new_id[i as usize].get_or_insert_with(|| {
                out.vertices.push(shell.vertices[i as usize]);
                out.vertex_rows.push(shell.vertex_rows[i as usize]);
                out.vertices.len() as u32 - 1
            })
        }));
        out.face_source.push(src[f]);
    }
    out
}

fn ekey(a: u32, b: u32) -> (u32, u32) {
    if a < b { (a, b) } else { (b, a) }
}

/// The shells of a closed piece, each a piece ordered as in `piece` and the
/// shells ordered by first sub-triangle. Sub-triangles join across edges with
/// two of them, and across edges with more in the pairs that bound wedges of
/// solid, found by exact radial order about the edge.
pub(crate) fn shells(arr: &Arrangement, piece: &Piece) -> Result<Vec<Piece>> {
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
    let mut keys: Vec<(u32, u32)> = at.keys().copied().collect();
    keys.sort_unstable();
    for (p, q) in keys {
        let ks = &at[&(p, q)];
        if ks.len() == 2 {
            join(&mut parent, ks[0], ks[1]);
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
                    join(&mut parent, ks[a], ks[b]);
                }
            }
            _ => {
                return Err(Error::invalid(
                    "internal: the faces around an edge of the result do not bound wedges of solid",
                ));
            }
        }
    }
    let mut index: HashMap<u32, usize> = HashMap::new();
    let mut out: Vec<Piece> = vec![];
    for (k, &entry) in piece.iter().enumerate() {
        let r = root(&mut parent, k as u32);
        let next = out.len();
        let i = *index.entry(r).or_insert(next);
        if i == out.len() {
            out.push(vec![]);
        }
        out[i].push(entry);
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

pub(crate) const SELF_TOUCH: &str =
    "result is not manifold: a shell touches itself along an edge or vertex";

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
    if !asm.triangles.is_empty() && !is_closed_manifold(&asm.triangles) {
        return Err(Error::invalid(SELF_TOUCH));
    }
    let face_source = asm.face_source.clone();
    let mesh = build_mesh(
        asm,
        (a.vertex_attributes(), b.vertex_attributes()),
        (a.face_attributes(), b.face_attributes()),
        labels,
    )?;
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
        let out = collapse_rounded(s);
        check(&out);
        assert_eq!(out.triangles.len(), 4);
        assert_eq!(out.vertex_rows, [Some(0), Some(2), Some(1), Some(3)]);
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
        let out = collapse_rounded(s);
        check(&out);
        assert_eq!(out.triangles.len(), 6);
        let area = |t: &[u32; 3]| {
            let p = t.map(|i| out.vertices[i as usize]);
            !collinear3(p[0], p[1], p[2])
        };
        assert!(out.triangles.iter().all(area));
        assert!(out.face_source.contains(&(0, 0)));
    }
}
