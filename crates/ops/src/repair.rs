use std::collections::{HashMap, HashSet};

use rayon::prelude::*;
use wireframes_core::predicates::collinear3;
use wireframes_core::topology::{Shells, locate_point};
use wireframes_core::validate::{degenerate_mask, earliest_within};
use wireframes_core::{
    Aabb, EdgeMap, Mesh, MeshProblemKind, MeshSummary, Result, tolerance, validate,
};

use crate::measure::{cross, dot, norm, sub, sum_by};

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum HoleFill {
    Fan,
    MinArea,
}

#[derive(Clone, Debug, PartialEq)]
pub struct RepairOptions {
    /// None: no weld; `Some(0.0)`: exactly equal coordinates.
    pub weld: Option<f64>,
    pub remove_degenerate: bool,
    pub remove_duplicate_faces: bool,
    pub remove_unused_vertices: bool,
    /// Consistent winding, outward by nesting parity.
    pub orient: bool,
    pub fill_holes: Option<HoleFill>,
    /// Weld also merges the coincident vertices where closed shells touch,
    /// which makes their contact share vertex indices. By default each closed
    /// shell keeps its own vertices there.
    pub weld_touching_shells: bool,
}

impl Default for RepairOptions {
    fn default() -> RepairOptions {
        RepairOptions {
            weld: Some(0.0),
            remove_degenerate: true,
            remove_duplicate_faces: true,
            remove_unused_vertices: true,
            orient: true,
            fill_holes: None,
            weld_touching_shells: false,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum RepairAction {
    VertexMerged {
        vertex: u32,
        into: u32,
    },
    VertexRemoved {
        vertex: u32,
    },
    /// Input `vertex` was shared by closed shells that touch there; output
    /// vertex `copy` is the copy it gets for one more of them.
    VertexSplit {
        vertex: u32,
        copy: u32,
    },
    /// `reason` is `DegenerateFace` or `DuplicateFace`; `face` indexes the input.
    FaceRemoved {
        face: u32,
        reason: MeshProblemKind,
    },
    /// `face` indexes the input.
    FaceFlipped {
        face: u32,
    },
    /// `face` indexes the output.
    FaceAdded {
        face: u32,
        hole: u32,
    },
}

#[derive(Clone, Debug, PartialEq)]
pub struct RepairReport {
    /// In step order: weld (merges, then splits), degenerate, duplicate,
    /// unused, orient, fill.
    pub actions: Vec<RepairAction>,
    /// Output vertex of each input vertex; a split vertex maps to its first copy.
    pub vertex_map: Vec<Option<u32>>,
    pub face_map: Vec<Option<u32>>,
    pub before: MeshSummary,
    pub after: MeshSummary,
    /// Loops not filled: pinched at a vertex used by two boundary edges, or
    /// whose fill would add a degenerate face, a duplicate face, or a diagonal
    /// on an edge the mesh already has.
    pub holes_left: usize,
}

fn root(parent: &mut [u32], mut i: u32) -> u32 {
    while parent[i as usize] != i {
        parent[i as usize] = parent[parent[i as usize] as usize];
        i = parent[i as usize];
    }
    i
}

/// Smallest index of each vertex's cluster, the transitive closure of pairs
/// within `tol`.
fn clusters(v: &[[f64; 3]], tol: f64) -> Vec<u32> {
    if tol == 0.0 {
        return earliest_within(v, 0.0);
    }
    let cell = |p: [f64; 3]| p.map(|x| (x / tol).floor() as i64);
    let mut grid: HashMap<[i64; 3], Vec<u32>> = HashMap::new();
    let mut parent: Vec<u32> = (0..v.len() as u32).collect();
    let t2 = tol * tol;
    for (i, &p) in v.iter().enumerate() {
        let c = cell(p);
        for dx in -1..=1i64 {
            for dy in -1..=1i64 {
                for dz in -1..=1i64 {
                    let k = [
                        c[0].saturating_add(dx),
                        c[1].saturating_add(dy),
                        c[2].saturating_add(dz),
                    ];
                    for &j in grid.get(&k).into_iter().flatten() {
                        let d = sub(p, v[j as usize]);
                        if dot(d, d) <= t2 {
                            let (a, b) = (root(&mut parent, i as u32), root(&mut parent, j));
                            parent[a.max(b) as usize] = a.min(b);
                        }
                    }
                }
            }
        }
        grid.entry(c).or_default().push(i as u32);
    }
    (0..v.len() as u32).map(|i| root(&mut parent, i)).collect()
}

struct Component {
    faces: Vec<u32>,
    parity: Vec<bool>,
    orientable: bool,
    closed: bool,
}

/// Edge-connected components through edges with exactly two faces, each face's
/// winding parity relative to the component's first face.
fn components(mesh: &Mesh, skip: &[bool], edges: &EdgeMap) -> Vec<Component> {
    let tris = mesh.triangles();
    let forward = |f: u32, e: u32| {
        let t = tris[f as usize];
        let k = edges.face_edges[f as usize]
            .iter()
            .position(|&x| x == e)
            .expect("face holds edge");
        t[k] < t[(k + 1) % 3]
    };
    let mut parity: Vec<Option<bool>> = vec![None; tris.len()];
    let mut out = vec![];
    for start in (0..tris.len() as u32).filter(|&f| !skip[f as usize]) {
        if parity[start as usize].is_some() {
            continue;
        }
        parity[start as usize] = Some(false);
        let mut c = Component {
            faces: vec![],
            parity: vec![],
            orientable: true,
            closed: true,
        };
        let mut stack = vec![start];
        while let Some(f) = stack.pop() {
            c.faces.push(f);
            let pf = parity[f as usize].expect("set before push");
            for &e in &edges.face_edges[f as usize] {
                let [a, b] = match *edges.faces_of_edge(e) {
                    [a, b] => [a, b],
                    _ => {
                        c.closed = false;
                        continue;
                    }
                };
                let g = if a == f { b } else { a };
                let want = pf ^ (forward(f, e) == forward(g, e));
                match parity[g as usize] {
                    None => {
                        parity[g as usize] = Some(want);
                        stack.push(g);
                    }
                    Some(p) if p != want => c.orientable = false,
                    Some(_) => {}
                }
            }
        }
        c.faces.sort_unstable();
        c.parity = c
            .faces
            .iter()
            .map(|&f| parity[f as usize].expect("visited"))
            .collect();
        out.push(c);
    }
    out
}

/// Per face, whether `orient` flips it.
fn orient(mesh: &Mesh) -> Vec<bool> {
    let tris = mesh.triangles();
    let v = mesh.vertices();
    let skip = degenerate_mask(mesh);
    let edges = EdgeMap::new(mesh, Some(&skip));
    let comps = components(mesh, &skip, &edges);
    let o = mesh.bounds().map_or([0.0; 3], |b| b.min);
    let closed: Vec<usize> = (0..comps.len())
        .filter(|&c| comps[c].orientable && comps[c].closed)
        .collect();
    let face_lists: Vec<Vec<[u32; 3]>> = closed
        .iter()
        .map(|&c| comps[c].faces.iter().map(|&f| tris[f as usize]).collect())
        .collect();
    let boxes: Vec<Aabb> = face_lists
        .iter()
        .map(|fs| {
            let pts: Vec<[f64; 3]> = fs.iter().flatten().map(|&i| v[i as usize]).collect();
            Aabb::from_points(&pts)
        })
        .collect();
    let flip_closed: Vec<bool> = (0..closed.len())
        .into_par_iter()
        .map(|s| {
            let c = &comps[closed[s]];
            let vol = sum_by(c.faces.len(), |k| {
                let p = tris[c.faces[k] as usize].map(|i| sub(v[i as usize], o));
                let d = dot(p[0], cross(p[1], p[2]));
                [if c.parity[k] { -d } else { d }]
            })[0];
            if vol == 0.0 {
                return false;
            }
            let depth = (0..closed.len())
                .filter(|&t| t != s)
                .filter(|&t| {
                    let (a, b) = (&boxes[t], &boxes[s]);
                    (0..3).all(|k| a.min[k] <= b.min[k] && b.max[k] <= a.max[k])
                        && face_lists[s]
                            .iter()
                            .flatten()
                            .map(|&i| locate_point(v, &face_lists[t], v[i as usize]))
                            .find(|&l| l != 0)
                            == Some(-1)
                })
                .count();
            (depth % 2 == 0) != (vol > 0.0)
        })
        .collect();
    let mut flip = vec![false; tris.len()];
    for (k, c) in comps.iter().enumerate() {
        if !c.orientable {
            continue;
        }
        let all = match closed.binary_search(&k) {
            Ok(s) => flip_closed[s],
            Err(_) => 2 * c.parity.iter().filter(|&&p| p).count() > c.parity.len(),
        };
        for (&f, &p) in c.faces.iter().zip(&c.parity) {
            flip[f as usize] = p != all;
        }
    }
    flip
}

/// Boundary loops in face direction and whether each can be filled: closed,
/// at least 3 vertices, and no vertex with two outgoing boundary half-edges.
fn boundary_loops(mesh: &Mesh) -> Vec<(Vec<u32>, bool)> {
    let tris = mesh.triangles();
    let skip = degenerate_mask(mesh);
    let edges = EdgeMap::new(mesh, Some(&skip));
    let mut half: Vec<(u32, u32)> = vec![];
    for e in 0..edges.edges.len() as u32 {
        if let [f] = *edges.faces_of_edge(e) {
            let k = edges.face_edges[f as usize]
                .iter()
                .position(|&x| x == e)
                .expect("face holds edge");
            half.push((3 * f + k as u32, tris[f as usize][k]));
        }
    }
    half.sort_unstable();
    let target = |h: u32| tris[(h / 3) as usize][((h % 3 + 1) % 3) as usize];
    let mut leaving: Vec<(u32, u32)> = half.iter().map(|&(h, o)| (o, h)).collect();
    leaving.sort_unstable();
    let out_degree = |v: u32| {
        let from = leaving.partition_point(|&(o, _)| o < v);
        leaving[from..].iter().take_while(|&&(o, _)| o == v).count()
    };
    let mut used: HashMap<u32, bool> = half.iter().map(|&(h, _)| (h, false)).collect();
    let mut loops = vec![];
    for &(start, first) in &half {
        if used[&start] {
            continue;
        }
        let mut ring = vec![];
        let mut h = start;
        let closed = loop {
            used.insert(h, true);
            ring.push(tris[(h / 3) as usize][(h % 3) as usize]);
            let v = target(h);
            let from = leaving.partition_point(|&(o, _)| o < v);
            match leaving[from..]
                .iter()
                .take_while(|&&(o, _)| o == v)
                .find(|&&(_, g)| !used[&g])
            {
                Some(&(_, g)) => h = g,
                None => break v == first,
            }
        };
        let ok = closed && ring.len() >= 3 && ring.iter().all(|&v| out_degree(v) == 1);
        loops.push((ring, ok));
    }
    loops
}

fn tri_area2(v: &[[f64; 3]], t: [u32; 3]) -> f64 {
    let [a, b, c] = t.map(|i| v[i as usize]);
    norm(cross(sub(b, a), sub(c, a)))
}

/// Triangles closing polygon `p` (already in the winding the new faces take).
fn fill(v: &[[f64; 3]], p: &[u32], how: HoleFill) -> Vec<[u32; 3]> {
    let n = p.len();
    if how == HoleFill::Fan || n > 256 {
        return (1..n - 1).map(|i| [p[0], p[i], p[i + 1]]).collect();
    }
    let o = Aabb::from_points(&p.iter().map(|&i| v[i as usize]).collect::<Vec<_>>()).min;
    let local: Vec<[f64; 3]> = p.iter().map(|&i| sub(v[i as usize], o)).collect();
    let cost = |i: usize, m: usize, j: usize| {
        let degenerate =
            wireframes_core::predicates::collinear3(local[i], local[m], local[j]) as usize;
        (
            degenerate,
            tri_area2(&local, [i as u32, m as u32, j as u32]),
        )
    };
    let mut w = vec![(0usize, 0.0f64); n * n];
    let mut choice = vec![0usize; n * n];
    for len in 2..n {
        for i in 0..n - len {
            let j = i + len;
            let mut best: Option<((usize, f64), usize)> = None;
            for m in i + 1..j {
                let (a, b, c) = (w[i * n + m], w[m * n + j], cost(i, m, j));
                let total = (a.0 + b.0 + c.0, a.1 + b.1 + c.1);
                if best.is_none_or(|(b, _)| total < b) {
                    best = Some((total, m));
                }
            }
            let (total, m) = best.expect("len >= 2");
            w[i * n + j] = total;
            choice[i * n + j] = m;
        }
    }
    let mut out = vec![];
    let mut stack = vec![(0, n - 1)];
    while let Some((i, j)) = stack.pop() {
        if j - i < 2 {
            continue;
        }
        let m = choice[i * n + j];
        out.push([p[i], p[m], p[j]]);
        stack.push((m, j));
        stack.push((i, m));
    }
    out
}

const OPEN: u32 = u32::MAX - 1;
const UNUSED: u32 = u32::MAX;

/// Per face, the closed shell of `welded` it lies in when `split`, so each
/// closed shell keeps its own vertices where it touches another; `OPEN` for
/// every other face.
fn shell_buckets(welded: &Mesh, keep: &[bool], degenerate: &[bool], split: bool) -> Vec<u32> {
    let nf = welded.n_faces();
    if !split {
        return vec![OPEN; nf];
    }
    let skip: Vec<bool> = (0..nf).map(|f| !keep[f] || degenerate[f]).collect();
    let edges = EdgeMap::new(welded, Some(&skip));
    let shells = Shells::new(welded, &edges);
    let tris = welded.triangles();
    let mut closed = vec![true; shells.count];
    for id in 0..edges.edges.len() as u32 {
        let fs = edges.faces_of_edge(id);
        let open = match *fs {
            [_] => true,
            [a, b] => {
                let [p, q] = edges.edges[id as usize];
                let dir = |f: u32| {
                    let t = tris[f as usize];
                    (0..3).any(|k| t[k] == p && t[(k + 1) % 3] == q)
                };
                dir(a) == dir(b)
            }
            _ => shells.non_manifold.binary_search(&id).is_ok(),
        };
        if open {
            for &f in fs {
                closed[shells.of[f as usize] as usize] = false;
            }
        }
    }
    (0..nf)
        .map(|f| match shells.of[f] {
            s if s != u32::MAX && closed[s as usize] => s,
            _ => OPEN,
        })
        .collect()
}

pub fn repair(mesh: &Mesh, opts: &RepairOptions) -> Result<(Mesh, RepairReport)> {
    use MeshProblemKind::{DegenerateFace, DuplicateFace};
    let weld = opts.weld.map(|t| tolerance::check("weld", t)).transpose()?;
    let tol = weld.unwrap_or(0.0);
    let before = validate(mesh, tol)?.summary;
    let (nv, nf) = (mesh.n_vertices(), mesh.n_faces());
    let mut actions = vec![];

    let rep: Vec<u32> = match weld {
        Some(t) => clusters(mesh.vertices(), t),
        None => (0..nv as u32).collect(),
    };
    let welded = Mesh::new(
        mesh.vertices().to_vec(),
        mesh.triangles()
            .iter()
            .map(|t| t.map(|i| rep[i as usize]))
            .collect(),
    )?;
    let tris = welded.triangles();

    let degenerate = degenerate_mask(&welded);
    let mut keep = vec![true; nf];
    if opts.remove_degenerate {
        for f in (0..nf).filter(|&f| degenerate[f]) {
            keep[f] = false;
            actions.push(RepairAction::FaceRemoved {
                face: f as u32,
                reason: DegenerateFace,
            });
        }
    }
    if opts.remove_duplicate_faces {
        let mut seen: HashMap<[u32; 3], u32> = HashMap::new();
        for f in 0..nf {
            if !keep[f] || degenerate[f] {
                continue;
            }
            let mut s = tris[f];
            s.sort_unstable();
            if *seen.entry(s).or_insert(f as u32) != f as u32 {
                keep[f] = false;
                actions.push(RepairAction::FaceRemoved {
                    face: f as u32,
                    reason: DuplicateFace,
                });
            }
        }
    }

    let bucket = shell_buckets(
        &welded,
        &keep,
        &degenerate,
        weld.is_some() && !opts.weld_touching_shells,
    );
    let orig = mesh.triangles();
    let kept_faces: Vec<u32> = (0..nf as u32).filter(|&f| keep[f as usize]).collect();
    let mut first_row: HashMap<(u32, u32), u32> = HashMap::new();
    for &f in &kept_faces {
        for &j in &orig[f as usize] {
            let row = first_row
                .entry((rep[j as usize], bucket[f as usize]))
                .or_insert(j);
            *row = (*row).min(j);
        }
    }
    let mut by_cluster: HashMap<u32, Vec<(u32, u32)>> = HashMap::new();
    for (&(c, b), &row) in &first_row {
        by_cluster.entry(c).or_default().push((row, b));
    }
    let mut idents: Vec<(u32, u32, u32)> = vec![];
    for (&c, list) in by_cluster.iter_mut() {
        list.sort_unstable();
        for (k, &(row, b)) in list.iter().enumerate() {
            idents.push((if k == 0 { c } else { row }, b, c));
        }
    }
    for c in (0..nv as u32).filter(|&c| rep[c as usize] == c && !by_cluster.contains_key(&c)) {
        if opts.remove_unused_vertices {
            continue;
        }
        idents.push((c, UNUSED, c));
    }
    idents.sort_unstable();
    let mut out_of: HashMap<(u32, u32), u32> = HashMap::new();
    let mut first_of_cluster: HashMap<u32, u32> = HashMap::new();
    let mut kept_vertices = vec![];
    let mut splits = vec![];
    let mut seen_row: HashMap<u32, u32> = HashMap::new();
    for (k, &(row, b, c)) in idents.iter().enumerate() {
        out_of.insert((c, b), k as u32);
        first_of_cluster.entry(c).or_insert(k as u32);
        kept_vertices.push(row);
        if seen_row.insert(row, k as u32).is_some() {
            splits.push(RepairAction::VertexSplit {
                vertex: row,
                copy: k as u32,
            });
        }
    }
    let mut first_face: Vec<Option<u32>> = vec![None; nv];
    for &f in &kept_faces {
        for &j in &orig[f as usize] {
            first_face[j as usize].get_or_insert(f);
        }
    }
    let vertex_map: Vec<Option<u32>> = (0..nv)
        .map(|i| {
            let c = rep[i];
            match first_face[i] {
                Some(f) => Some(out_of[&(c, bucket[f as usize])]),
                None => first_of_cluster.get(&c).copied(),
            }
        })
        .collect();
    let mut weld_actions = vec![];
    for i in 0..nv {
        let into = match vertex_map[i] {
            Some(k) if first_face[i].is_some() => kept_vertices[k as usize],
            _ => rep[i],
        };
        if into != i as u32 {
            weld_actions.push(RepairAction::VertexMerged {
                vertex: i as u32,
                into,
            });
        }
    }
    weld_actions.extend(splits);
    actions.splice(0..0, weld_actions);
    for c in (0..nv as u32).filter(|&c| rep[c as usize] == c && !first_of_cluster.contains_key(&c))
    {
        actions.push(RepairAction::VertexRemoved { vertex: c });
    }
    let mut face_map = vec![None; nf];
    for (k, &f) in kept_faces.iter().enumerate() {
        face_map[f as usize] = Some(k as u32);
    }
    let vertices: Vec<[f64; 3]> = kept_vertices
        .iter()
        .map(|&i| mesh.vertices()[i as usize])
        .collect();
    let mut out_tris: Vec<[u32; 3]> = kept_faces
        .iter()
        .map(|&f| orig[f as usize].map(|j| out_of[&(rep[j as usize], bucket[f as usize])]))
        .collect();

    let reverse = |t: &mut [u32; 3]| *t = [t[0], t[2], t[1]];
    let mut flipped = vec![false; out_tris.len()];
    if opts.orient {
        flipped = orient(&Mesh::new(vertices.clone(), out_tris.clone())?);
        for (t, _) in out_tris.iter_mut().zip(&flipped).filter(|(_, f)| **f) {
            reverse(t);
        }
    }

    let loops = boundary_loops(&Mesh::new(vertices.clone(), out_tris.clone())?);
    let edge = |a: u32, b: u32| [a.min(b), a.max(b)];
    let sorted = |mut t: [u32; 3]| {
        t.sort_unstable();
        t
    };
    let mut edges: HashSet<[u32; 2]> = out_tris
        .iter()
        .flat_map(|t| (0..3).map(move |k| edge(t[k], t[(k + 1) % 3])))
        .collect();
    let mut faces: HashSet<[u32; 3]> = out_tris.iter().map(|&t| sorted(t)).collect();
    let mut holes_left = 0;
    let mut added = vec![];
    for (hole, (ring, ok)) in loops.iter().enumerate() {
        let Some(how) = opts.fill_holes.filter(|_| *ok) else {
            holes_left += 1;
            continue;
        };
        let mut p = ring.clone();
        p[1..].reverse();
        let new = fill(&vertices, &p, how);
        let rim: HashSet<[u32; 2]> = (0..p.len())
            .map(|k| edge(p[k], p[(k + 1) % p.len()]))
            .collect();
        let valid = new.iter().all(|&t| {
            let [a, b, c] = t.map(|i| vertices[i as usize]);
            !collinear3(a, b, c)
                && !faces.contains(&sorted(t))
                && (0..3)
                    .map(|k| edge(t[k], t[(k + 1) % 3]))
                    .all(|e| rim.contains(&e) || !edges.contains(&e))
        });
        if !valid {
            holes_left += 1;
            continue;
        }
        for t in new {
            edges.extend((0..3).map(|k| edge(t[k], t[(k + 1) % 3])));
            faces.insert(sorted(t));
            added.push(RepairAction::FaceAdded {
                face: out_tris.len() as u32,
                hole: hole as u32,
            });
            out_tris.push(t);
        }
    }
    // A fill can close a shell or change nesting, so orient again to keep
    // repair idempotent; a face flipped twice nets to no action.
    if opts.orient && !added.is_empty() {
        let again = orient(&Mesh::new(vertices.clone(), out_tris.clone())?);
        for (k, t) in out_tris.iter_mut().enumerate().filter(|&(k, _)| again[k]) {
            reverse(t);
            if let Some(f) = flipped.get_mut(k) {
                *f = !*f;
            }
        }
    }
    actions.extend(
        (0..flipped.len())
            .filter(|&k| flipped[k])
            .map(|k| RepairAction::FaceFlipped {
                face: kept_faces[k],
            }),
    );
    actions.extend(added);

    let mut face_rows: Vec<Option<u32>> = kept_faces.iter().map(|&f| Some(f)).collect();
    face_rows.resize(out_tris.len(), None);
    let out = Mesh::new(vertices, out_tris)?
        .with_attributes(
            mesh.vertex_attributes().take(&kept_vertices),
            mesh.face_attributes().take_or_null(&face_rows),
        )?
        .with_crs(mesh.crs.clone())
        .with_length_unit(mesh.length_unit.clone());
    let after = validate(&out, tol)?.summary;
    Ok((
        out,
        RepairReport {
            actions,
            vertex_map,
            face_map,
            before,
            after,
            holes_left,
        },
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::measure::volume;
    use proptest::prelude::*;
    use wireframes_core::fixtures::{self, Rng, broken};
    use wireframes_core::{Column, Error};

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= 1e-9 * 1f64.max(a.abs()).max(b.abs())
    }

    fn fix(m: &Mesh) -> (Mesh, RepairReport) {
        repair(m, &RepairOptions::default()).unwrap()
    }

    fn check_maps(input: &Mesh, out: &Mesh, r: &RepairReport) {
        assert_eq!(r.vertex_map.len(), input.n_vertices());
        assert_eq!(r.face_map.len(), input.n_faces());
        for (f, m) in r.face_map.iter().enumerate() {
            let Some(g) = m else { continue };
            let mut want: Vec<u32> = input.triangles()[f]
                .iter()
                .map(|&i| r.vertex_map[i as usize].expect("corner of a kept face"))
                .collect();
            let mut got = out.triangles()[*g as usize].to_vec();
            want.sort_unstable();
            got.sort_unstable();
            assert_eq!(want, got, "face {f}");
        }
        for (i, m) in r.vertex_map.iter().enumerate() {
            if let Some(j) = m {
                let (a, b) = (input.vertices()[i], out.vertices()[*j as usize]);
                let tol = 1e-12;
                assert!((0..3).all(|k| (a[k] - b[k]).abs() <= tol) || r.actions.iter().any(
                    |x| matches!(x, RepairAction::VertexMerged { vertex, .. } if *vertex == i as u32)
                ));
            }
        }
    }

    #[test]
    fn valid_meshes_are_unchanged() {
        for m in [
            fixtures::cube(),
            fixtures::tetrahedron(),
            fixtures::sphere(2),
            fixtures::torus(12, 8),
            fixtures::nested_shells(),
            fixtures::open_cube(),
            fixtures::pit_surface(12, 3),
        ] {
            let (out, r) = fix(&m);
            assert!(r.actions.is_empty(), "{:?}", r.actions);
            assert_eq!(out, m);
            assert_eq!(r.before, r.after);
            assert!(
                r.vertex_map
                    .iter()
                    .enumerate()
                    .all(|(i, x)| *x == Some(i as u32))
            );
        }
    }

    #[test]
    fn broken_fixtures_become_solids_and_report_every_change() {
        use RepairAction::*;
        let cases: Vec<(Mesh, Vec<RepairAction>)> = vec![
            (
                broken::degenerate_face(),
                vec![
                    FaceRemoved {
                        face: 12,
                        reason: MeshProblemKind::DegenerateFace,
                    },
                    VertexRemoved { vertex: 8 },
                ],
            ),
            (
                broken::duplicate_face(),
                vec![FaceRemoved {
                    face: 12,
                    reason: MeshProblemKind::DuplicateFace,
                }],
            ),
            (
                broken::duplicate_vertex(),
                vec![VertexMerged { vertex: 8, into: 0 }],
            ),
            (broken::unused_vertex(), vec![VertexRemoved { vertex: 8 }]),
            (broken::flipped_face(), vec![FaceFlipped { face: 0 }]),
            (
                broken::inward_cube(),
                (0..12).map(|face| FaceFlipped { face }).collect(),
            ),
        ];
        for (m, want) in cases {
            let s = m.validate(0.0).unwrap().summary;
            assert!(!s.is_solid() || s.unused_vertices > 0);
            let (out, r) = fix(&m);
            assert_eq!(out.validate(0.0).unwrap().summary.unused_vertices, 0);
            assert_eq!(r.actions, want);
            assert!(r.after.is_solid(), "{:?}", r.after);
            assert!(close(volume(&out).unwrap(), 1.0));
            check_maps(&m, &out, &r);
        }
    }

    #[test]
    fn unrepairable_defects_stay_reported() {
        for m in [
            broken::non_manifold_fin(),
            broken::bowtie(),
            broken::edge_overlapping_shells(),
            broken::overlapping_cubes(),
            broken::touching_cubes(),
        ] {
            let (out, r) = fix(&m);
            assert!(!r.after.is_solid());
            let (_, again) = fix(&out);
            assert!(again.actions.is_empty(), "{:?}", again.actions);
        }
    }

    #[test]
    fn touching_shells_keep_their_own_vertices() {
        let apart = fixtures::edge_touching_cubes(false);
        let (out, r) = fix(&apart);
        assert!(r.actions.is_empty(), "{:?}", r.actions);
        assert_eq!(out, apart);

        let welded = fixtures::edge_touching_cubes(true);
        let (out, r) = fix(&welded);
        let splits: Vec<&RepairAction> = r
            .actions
            .iter()
            .filter(|a| matches!(a, RepairAction::VertexSplit { .. }))
            .collect();
        assert_eq!(splits.len(), 2, "{:?}", r.actions);
        assert_eq!(r.actions.len(), 2, "{:?}", r.actions);
        assert_eq!(out.n_vertices(), 16);
        assert!(r.after.is_solid() && r.before.is_solid(), "{:?}", r.after);
        assert_eq!(
            (
                r.after.touching_shell_edges,
                r.after.touching_shell_vertices
            ),
            (1, 2)
        );
        let (again, r2) = fix(&out);
        assert!(r2.actions.is_empty(), "{:?}", r2.actions);
        assert_eq!(again, out);

        let soup = fixtures::soup(&apart);
        let (out, r) = fix(&soup);
        assert!(r.after.is_solid(), "{:?}", r.after);
        assert_eq!((out.n_vertices(), r.after.shells), (16, 2));
        assert!(close(volume(&out).unwrap(), 2.0));
        assert!(fix(&out).1.actions.is_empty());

        let opts = RepairOptions {
            weld_touching_shells: true,
            ..RepairOptions::default()
        };
        let (out, r) = repair(&apart, &opts).unwrap();
        assert_eq!(out.n_vertices(), 14);
        assert_eq!(
            r.actions
                .iter()
                .filter(|a| matches!(a, RepairAction::VertexMerged { .. }))
                .count(),
            2
        );
        assert!(r.after.is_solid(), "{:?}", r.after);
        check_maps(&apart, &out, &r);
    }

    #[test]
    fn soup_welds_back_into_the_solid() {
        let m = fixtures::soup(&fixtures::sphere(2));
        let (out, r) = fix(&m);
        assert!(r.after.is_solid());
        assert_eq!(out.n_vertices(), fixtures::sphere(2).n_vertices());
        assert!(close(
            volume(&out).unwrap(),
            volume(&fixtures::sphere(2)).unwrap()
        ));
        check_maps(&m, &out, &r);
    }

    #[test]
    fn tolerance_weld_is_transitive() {
        let m = Mesh::new(
            vec![[0.0; 3], [0.6, 0.0, 0.0], [1.2, 0.0, 0.0], [5.0, 0.0, 0.0]],
            vec![],
        )
        .unwrap();
        let opts = RepairOptions {
            weld: Some(1.0),
            remove_unused_vertices: false,
            ..RepairOptions::default()
        };
        let (out, r) = repair(&m, &opts).unwrap();
        assert_eq!(out.n_vertices(), 2);
        assert_eq!(r.vertex_map, vec![Some(0), Some(0), Some(0), Some(1)]);
        let bad = RepairOptions {
            weld: Some(-1.0),
            ..RepairOptions::default()
        };
        assert!(matches!(repair(&m, &bad), Err(Error::InvalidInput(_))));
    }

    #[test]
    fn holes_are_filled_with_outward_faces() {
        for how in [HoleFill::Fan, HoleFill::MinArea] {
            let m = fixtures::open_cube()
                .with_face_column("id", Column::i64((0..10).collect()))
                .unwrap();
            let opts = RepairOptions {
                fill_holes: Some(how),
                ..RepairOptions::default()
            };
            let (out, r) = repair(&m, &opts).unwrap();
            assert_eq!(r.holes_left, 0);
            assert_eq!(
                r.actions,
                vec![
                    RepairAction::FaceAdded { face: 10, hole: 0 },
                    RepairAction::FaceAdded { face: 11, hole: 0 },
                ]
            );
            assert!(r.after.is_solid(), "{:?}", r.after);
            assert!(close(volume(&out).unwrap(), 1.0));
            let id = out.face_attributes().column("id").unwrap();
            assert!(id.is_valid(9) && !id.is_valid(10) && !id.is_valid(11));
            let (_, again) = repair(&out, &opts).unwrap();
            assert!(again.actions.is_empty());
        }
        let (_, r) = fix(&fixtures::open_cube());
        assert_eq!(r.holes_left, 1);
    }

    fn fan_fill_twice(v: Vec<[f64; 3]>, t: Vec<[u32; 3]>) -> (RepairReport, RepairReport) {
        let opts = RepairOptions {
            fill_holes: Some(HoleFill::Fan),
            ..RepairOptions::default()
        };
        let (out, r) = repair(&Mesh::new(v, t).unwrap(), &opts).unwrap();
        let (_, again) = repair(&out, &opts).unwrap();
        (r, again)
    }

    #[test]
    fn fan_fill_of_a_flat_surface_is_idempotent() {
        let square = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let (r, again) = fan_fill_twice(square.clone(), vec![[0, 1, 2], [0, 2, 3]]);
        assert!(r.actions.is_empty(), "{:?}", r.actions);
        assert_eq!(r.holes_left, 1);
        assert_eq!(r.after.duplicate_faces, 0);
        assert!(again.actions.is_empty(), "{:?}", again.actions);
        let (r, again) = fan_fill_twice(square, vec![[0, 1, 2]]);
        assert_eq!(r.holes_left, 1);
        assert!(again.actions.is_empty(), "{:?}", again.actions);
    }

    #[test]
    fn fan_fill_does_not_add_degenerate_faces() {
        let (r, again) = fan_fill_twice(
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
            ],
            vec![[0, 1, 3], [1, 2, 3]],
        );
        assert_eq!(r.after.degenerate_faces, 0);
        assert_eq!(r.holes_left, 1);
        assert!(again.actions.is_empty(), "{:?}", again.actions);
    }

    #[test]
    fn min_area_avoids_collinear_triangles() {
        let ring = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let p: Vec<u32> = (0..5).collect();
        let t = fill(&ring, &p, HoleFill::MinArea);
        assert_eq!(t.len(), 3);
        for f in t {
            let c = f.map(|i| ring[i as usize]);
            assert!(!wireframes_core::predicates::collinear3(c[0], c[1], c[2]));
        }
    }

    fn scrambled(seed: u64, base: &Mesh) -> Mesh {
        let mut rng = Rng::new(seed);
        let mut faces = base.triangles().to_vec();
        let extra = rng.below(4);
        for _ in 0..extra {
            let f = faces[rng.below(faces.len())];
            faces.push(f);
        }
        rng.shuffle(&mut faces);
        for t in faces.iter_mut() {
            if rng.below(3) == 0 {
                *t = [t[0], t[2], t[1]];
            }
        }
        let mut perm: Vec<u32> = (0..base.n_vertices() as u32).collect();
        rng.shuffle(&mut perm);
        let mut v = vec![[0.0; 3]; base.n_vertices()];
        for (i, &p) in perm.iter().enumerate() {
            v[p as usize] = base.vertices()[i];
        }
        let t = faces.iter().map(|t| t.map(|i| perm[i as usize])).collect();
        let m = Mesh::new(v, t).unwrap();
        fixtures::soup(&m)
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(48))]
        #[test]
        fn repair_restores_scrambled_solids(seed in any::<u64>(), which in 0usize..4) {
            let base = match which {
                0 => fixtures::cube(),
                1 => fixtures::sphere(1),
                2 => fixtures::torus(8, 6),
                _ => fixtures::nested_shells(),
            };
            let m = scrambled(seed, &base);
            let (out, r) = fix(&m);
            prop_assert!(r.after.is_solid(), "{:?}", r.after);
            prop_assert!(close(volume(&out).unwrap(), volume(&base).unwrap()));
            check_maps(&m, &out, &r);
            let (again, r2) = fix(&out);
            prop_assert!(r2.actions.is_empty(), "{:?}", r2.actions);
            prop_assert_eq!(again, out);
        }
    }
}
