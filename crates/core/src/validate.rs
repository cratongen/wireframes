use std::collections::HashMap;

use rayon::prelude::*;

use crate::exact::shell_volume_sign;
use crate::intersect::{SegTri, TriTri, seg_tri, tri_tri};
use crate::predicates::{Sign, collinear3, orient3d};
use crate::topology::{Radial, Shells, radial_pairs_at, shell_location};
use crate::{Aabb, Bvh, EdgeMap, Mesh, Result, tolerance};

/// Declaration order is the sort order of problems.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum MeshProblemKind {
    /// Repeated corner index, or three exactly collinear corners.
    DegenerateFace,
    /// Same corner set as an earlier face.
    DuplicateFace,
    /// Within the tolerance of an earlier vertex.
    DuplicateVertex,
    /// Edge of exactly one face.
    BoundaryEdge,
    /// Edge of three or more faces.
    NonManifoldEdge,
    /// Its faces form more than one edge-connected fan.
    NonManifoldVertex,
    /// Edge whose two faces run it the same way.
    InconsistentWinding,
    /// Closed shell wound against its nesting parity.
    InwardShell,
    /// Two faces meet where they should not (crossing or touching).
    SelfIntersection,
    /// Two coplanar faces overlap with positive area.
    CoplanarOverlap,
    /// No face uses the vertex.
    UnusedVertex,
}

impl MeshProblemKind {
    pub fn name(&self) -> &'static str {
        match self {
            MeshProblemKind::DegenerateFace => "degenerate_face",
            MeshProblemKind::DuplicateFace => "duplicate_face",
            MeshProblemKind::DuplicateVertex => "duplicate_vertex",
            MeshProblemKind::BoundaryEdge => "boundary_edge",
            MeshProblemKind::NonManifoldEdge => "non_manifold_edge",
            MeshProblemKind::NonManifoldVertex => "non_manifold_vertex",
            MeshProblemKind::InconsistentWinding => "inconsistent_winding",
            MeshProblemKind::InwardShell => "inward_shell",
            MeshProblemKind::SelfIntersection => "self_intersection",
            MeshProblemKind::CoplanarOverlap => "coplanar_overlap",
            MeshProblemKind::UnusedVertex => "unused_vertex",
        }
    }
}

/// One defect. Edge problems give the `face` holding the edge, its local `edge`
/// (corner `edge` to `(edge + 1) % 3`) and that edge's first `vertex`.
#[derive(Clone, Debug, PartialEq)]
pub struct MeshProblem {
    pub kind: MeshProblemKind,
    pub face: Option<u32>,
    pub vertex: Option<u32>,
    pub edge: Option<u8>,
    /// Earlier duplicate, other face of an edge, or earlier face of a pair.
    pub other: Option<u32>,
}

#[derive(Clone, Debug, PartialEq, Default)]
pub struct MeshSummary {
    pub degenerate_faces: usize,
    pub duplicate_faces: usize,
    pub duplicate_vertices: usize,
    pub boundary_edges: usize,
    pub non_manifold_edges: usize,
    pub non_manifold_vertices: usize,
    pub inconsistent_edges: usize,
    pub shells: usize,
    pub inward_shells: usize,
    pub self_intersections: usize,
    pub coplanar_overlaps: usize,
    pub unused_vertices: usize,
    /// Edge locations where closed shells touch each other or themselves:
    /// exactly coincident edges with more than two faces that pair into
    /// wedges of solid. Information, not a defect.
    pub touching_shell_edges: usize,
    /// Point locations where more than one vertex umbrella of closed shells
    /// meets, ends of touching edges included. Information, not a defect.
    pub touching_shell_vertices: usize,
    pub is_closed: bool,
}

impl MeshSummary {
    /// Closed, and every count except `shells`, `unused_vertices` and the
    /// touching counts is zero.
    pub fn is_solid(&self) -> bool {
        self.is_closed
            && self.degenerate_faces == 0
            && self.duplicate_faces == 0
            && self.duplicate_vertices == 0
            && self.boundary_edges == 0
            && self.non_manifold_edges == 0
            && self.non_manifold_vertices == 0
            && self.inconsistent_edges == 0
            && self.inward_shells == 0
            && self.self_intersections == 0
            && self.coplanar_overlaps == 0
    }

    /// The nonzero defect counts as `name=count`, comma separated, with
    /// `not closed` when open and nothing else is counted; `no defects` for a
    /// solid.
    pub fn defects(&self) -> String {
        let info = [
            "shells",
            "unused_vertices",
            "touching_shell_edges",
            "touching_shell_vertices",
            "is_closed",
        ];
        let mut out: Vec<String> = self
            .to_pairs()
            .into_iter()
            .filter(|(n, v)| *v != 0 && !info.contains(n))
            .map(|(n, v)| format!("{n}={v}"))
            .collect();
        if !self.is_closed && out.is_empty() {
            out.push("not closed".into());
        }
        if out.is_empty() {
            "no defects".into()
        } else {
            out.join(", ")
        }
    }

    /// Counts in a stable order, `is_closed` last as 0 or 1.
    pub fn to_pairs(&self) -> Vec<(&'static str, i64)> {
        vec![
            ("degenerate_faces", self.degenerate_faces as i64),
            ("duplicate_faces", self.duplicate_faces as i64),
            ("duplicate_vertices", self.duplicate_vertices as i64),
            ("boundary_edges", self.boundary_edges as i64),
            ("non_manifold_edges", self.non_manifold_edges as i64),
            ("non_manifold_vertices", self.non_manifold_vertices as i64),
            ("inconsistent_edges", self.inconsistent_edges as i64),
            ("shells", self.shells as i64),
            ("inward_shells", self.inward_shells as i64),
            ("self_intersections", self.self_intersections as i64),
            ("coplanar_overlaps", self.coplanar_overlaps as i64),
            ("unused_vertices", self.unused_vertices as i64),
            ("touching_shell_edges", self.touching_shell_edges as i64),
            (
                "touching_shell_vertices",
                self.touching_shell_vertices as i64,
            ),
            ("is_closed", self.is_closed as i64),
        ]
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct MeshReport {
    pub summary: MeshSummary,
    /// Sorted by kind, then face, vertex, other and edge.
    pub problems: Vec<MeshProblem>,
}

/// True for faces with a repeated corner index or exactly collinear corners.
pub fn degenerate_mask(mesh: &Mesh) -> Vec<bool> {
    mesh.triangles()
        .iter()
        .map(|&[a, b, c]| {
            if a == b || b == c || a == c {
                return true;
            }
            let v = mesh.vertices();
            collinear3(v[a as usize], v[b as usize], v[c as usize])
        })
        .collect()
}

fn root(parent: &mut [u32], mut i: u32) -> u32 {
    while parent[i as usize] != i {
        parent[i as usize] = parent[parent[i as usize] as usize];
        i = parent[i as usize];
    }
    i
}

fn join(parent: &mut [u32], a: u32, b: u32) {
    let (a, b) = (root(parent, a), root(parent, b));
    parent[a.max(b) as usize] = a.min(b);
}

fn key(p: [f64; 3]) -> [u64; 3] {
    p.map(|x| tolerance::canonical(x).to_bits())
}

/// For each vertex, the earliest vertex within `tolerance` of it (itself if none).
pub fn earliest_within(vertices: &[[f64; 3]], tolerance: f64) -> Vec<u32> {
    earliest_within_by(vertices, tolerance, |_, _| true)
}

/// For each vertex `i`, the earliest vertex `j` within `tolerance` of it with
/// `keep(i, j)` (itself if none).
fn earliest_within_by(
    vertices: &[[f64; 3]],
    tolerance: f64,
    keep: impl Fn(u32, u32) -> bool,
) -> Vec<u32> {
    if tolerance == 0.0 {
        let mut seen: HashMap<[u64; 3], Vec<u32>> = HashMap::with_capacity(vertices.len());
        return vertices
            .iter()
            .enumerate()
            .map(|(i, &p)| {
                let group = seen.entry(key(p)).or_default();
                let w = group
                    .iter()
                    .copied()
                    .find(|&j| keep(i as u32, j))
                    .unwrap_or(i as u32);
                group.push(i as u32);
                w
            })
            .collect();
    }
    let cell = |x: f64| (x / tolerance).floor() as i64;
    let mut grid: HashMap<[i64; 3], Vec<u32>> = HashMap::new();
    let t2 = tolerance * tolerance;
    let mut out = Vec::with_capacity(vertices.len());
    for (i, &p) in vertices.iter().enumerate() {
        let c = p.map(cell);
        let mut best = i as u32;
        for dx in -1..=1i64 {
            for dy in -1..=1i64 {
                for dz in -1..=1i64 {
                    let k = [
                        c[0].saturating_add(dx),
                        c[1].saturating_add(dy),
                        c[2].saturating_add(dz),
                    ];
                    for &j in grid.get(&k).into_iter().flatten() {
                        let q = vertices[j as usize];
                        let d2: f64 = (0..3).map(|a| (p[a] - q[a]) * (p[a] - q[a])).sum();
                        if d2 <= t2 && j < best && keep(i as u32, j) {
                            best = j;
                        }
                    }
                }
            }
        }
        grid.entry(c).or_default().push(i as u32);
        out.push(best);
    }
    out
}

fn problem(kind: MeshProblemKind) -> MeshProblem {
    MeshProblem {
        kind,
        face: None,
        vertex: None,
        edge: None,
        other: None,
    }
}

/// Degenerate, duplicate, edge, vertex, orientation and intersection problems,
/// and the contacts between closed shells. `tolerance` is the distance within
/// which a vertex duplicates an earlier one of its own shell.
///
/// Closed shells may touch each other, and a closed shell may touch itself, at
/// exactly coincident vertices and edges, by shared or by separate vertex
/// indices, when their interiors do not overlap there; such contacts are
/// counted in `touching_shell_vertices` and `touching_shell_edges` and are not
/// problems.
pub fn validate(mesh: &Mesh, tolerance: f64) -> Result<MeshReport> {
    use MeshProblemKind::*;
    let tolerance = tolerance::check("tolerance", tolerance)?;
    let (v, tris) = (mesh.vertices(), mesh.triangles());
    let degenerate = degenerate_mask(mesh);
    let mut problems = vec![];
    for f in (0..tris.len()).filter(|&f| degenerate[f]) {
        problems.push(MeshProblem {
            face: Some(f as u32),
            ..problem(DegenerateFace)
        });
    }
    let valid = || (0..tris.len()).filter(|&f| !degenerate[f]);

    let mut seen: HashMap<[u32; 3], u32> = HashMap::new();
    for f in valid() {
        let mut s = tris[f];
        s.sort_unstable();
        let first = *seen.entry(s).or_insert(f as u32);
        if first != f as u32 {
            problems.push(MeshProblem {
                face: Some(f as u32),
                other: Some(first),
                ..problem(DuplicateFace)
            });
        }
    }
    let mut used = vec![false; v.len()];
    for t in tris {
        for &i in t {
            used[i as usize] = true;
        }
    }
    for i in (0..v.len()).filter(|&i| !used[i]) {
        problems.push(MeshProblem {
            vertex: Some(i as u32),
            ..problem(UnusedVertex)
        });
    }

    let edges = EdgeMap::new(mesh, Some(&degenerate));
    let shells = Shells::new(mesh, &edges);
    let local = |f: u32, id: u32| {
        edges.face_edges[f as usize]
            .iter()
            .position(|&e| e == id)
            .expect("face holds edge") as u8
    };
    let mut inconsistent = vec![false; tris.len()];
    let mut open_or_pinched = vec![false; tris.len()];
    for id in 0..edges.edges.len() as u32 {
        let fs = edges.faces_of_edge(id);
        let at = |f: u32| {
            let e = local(f, id);
            MeshProblem {
                face: Some(f),
                edge: Some(e),
                vertex: Some(tris[f as usize][e as usize]),
                ..problem(BoundaryEdge)
            }
        };
        match *fs {
            [one] => {
                open_or_pinched[one as usize] = true;
                problems.push(at(one));
            }
            [a, b] => {
                if tris[a as usize][local(a, id) as usize]
                    == tris[b as usize][local(b, id) as usize]
                {
                    inconsistent[a as usize] = true;
                    inconsistent[b as usize] = true;
                    problems.push(MeshProblem {
                        kind: InconsistentWinding,
                        other: Some(a),
                        ..at(b)
                    });
                }
            }
            _ if shells.non_manifold.binary_search(&id).is_ok() => {
                for &f in fs {
                    open_or_pinched[f as usize] = true;
                }
                problems.push(MeshProblem {
                    kind: NonManifoldEdge,
                    ..at(fs[0])
                });
            }
            _ => {}
        }
    }
    let mut closed_shell = vec![true; shells.count];
    for f in valid() {
        if inconsistent[f] || open_or_pinched[f] {
            closed_shell[shells.of[f] as usize] = false;
        }
    }

    let mut umbrellas: Vec<usize> = (0..v.len()).map(|i| used[i] as usize).collect();
    for vert in 0..v.len() as u32 {
        let fan = edges.faces_of_vertex(vert);
        if fan.len() < 2 {
            continue;
        }
        let mut ids: Vec<u32> = fan
            .iter()
            .flat_map(|&f| {
                let t = tris[f as usize];
                let ids = edges.face_edges[f as usize];
                (0..3)
                    .filter(move |&e| t[e] == vert || t[(e + 1) % 3] == vert)
                    .map(move |e| ids[e])
            })
            .collect();
        ids.sort_unstable();
        ids.dedup();
        let mut parent: Vec<u32> = (0..fan.len() as u32).collect();
        for id in ids {
            for [a, b] in shells.links(&edges, id) {
                let [i, j] =
                    [a, b].map(|g| fan.binary_search(&g).expect("g holds the vertex") as u32);
                join(&mut parent, i, j);
            }
        }
        let fans: Vec<u32> = (0..fan.len() as u32)
            .filter(|&i| root(&mut parent, i) == i)
            .map(|i| shells.of[fan[i as usize] as usize])
            .collect();
        umbrellas[vert as usize] = fans.len();
        if fans.len() < 2 {
            continue;
        }
        if fans.iter().any(|&s| !closed_shell[s as usize]) {
            problems.push(MeshProblem {
                vertex: Some(vert),
                ..problem(NonManifoldVertex)
            });
        }
    }

    let shells_of = |i: usize| -> Vec<u32> {
        let mut s: Vec<u32> = edges
            .faces_of_vertex(i as u32)
            .iter()
            .map(|&f| shells.of[f as usize])
            .collect();
        s.sort_unstable();
        s.dedup();
        s
    };
    let vertex_shells: Vec<Vec<u32>> = (0..v.len()).map(shells_of).collect();
    let all_closed = |s: &[u32]| !s.is_empty() && s.iter().all(|&k| closed_shell[k as usize]);
    let apart = |i: u32, j: u32| {
        let (a, b) = (&vertex_shells[i as usize], &vertex_shells[j as usize]);
        key(v[i as usize]) == key(v[j as usize]) && all_closed(a) && all_closed(b)
    };
    for (i, &w) in earliest_within_by(v, tolerance, |i, j| !apart(i, j))
        .iter()
        .enumerate()
    {
        if w != i as u32 {
            problems.push(MeshProblem {
                vertex: Some(i as u32),
                other: Some(w),
                ..problem(DuplicateVertex)
            });
        }
    }

    let canon = earliest_within(v, 0.0);
    let mut at_point: HashMap<u32, (usize, Vec<u32>)> = HashMap::new();
    for i in 0..v.len() {
        let e = at_point.entry(canon[i]).or_default();
        e.0 += umbrellas[i];
        e.1.extend_from_slice(&vertex_shells[i]);
    }
    let mut contact_points: Vec<u32> = at_point
        .into_iter()
        .filter_map(|(c, (n, mut s))| {
            s.sort_unstable();
            s.dedup();
            (n > 1 && all_closed(&s)).then_some(c)
        })
        .collect();
    contact_points.sort_unstable();
    let is_contact = |i: u32| contact_points.binary_search(&canon[i as usize]).is_ok();
    let mut at_edge: HashMap<[u32; 2], Vec<u32>> = HashMap::new();
    for f in valid() {
        let t = tris[f];
        for e in 0..3 {
            let (a, b) = (t[e], t[(e + 1) % 3]);
            let id = edges.face_edges[f][e];
            if is_contact(a) || is_contact(b) || shells.contacts.binary_search(&id).is_ok() {
                let (a, b) = (canon[a as usize], canon[b as usize]);
                at_edge
                    .entry([a.min(b), a.max(b)])
                    .or_default()
                    .push(f as u32);
            }
        }
    }
    let mut contact_edges: Vec<([u32; 2], Vec<u32>)> = at_edge
        .into_iter()
        .filter(|(_, fs)| {
            let mut s: Vec<u32> = fs.iter().map(|&f| shells.of[f as usize]).collect();
            s.sort_unstable();
            s.dedup();
            fs.len() > 2 && all_closed(&s)
        })
        .collect();
    contact_edges.sort_unstable();
    let mut touching_shell_edges = 0;
    let mut edge_ends = contact_points.clone();
    for ([a, b], fs) in &contact_edges {
        match radial_pairs_at(mesh, fs, v[*a as usize], v[*b as usize]) {
            Radial::Paired(_) => {
                touching_shell_edges += 1;
                edge_ends.extend([*a, *b]);
            }
            Radial::Coincident => {}
            Radial::Unbalanced => problems.push(MeshProblem {
                face: Some(fs[0].max(fs[1])),
                other: Some(fs[0].min(fs[1])),
                ..problem(SelfIntersection)
            }),
        }
    }

    edge_ends.sort_unstable();
    edge_ends.dedup();

    let mut members: Vec<Vec<u32>> = vec![vec![]; shells.count];
    for f in valid() {
        members[shells.of[f] as usize].push(f as u32);
    }
    let closed_ids: Vec<usize> = (0..shells.count).filter(|&s| closed_shell[s]).collect();
    let closed: Vec<(Vec<[u32; 3]>, Aabb, Sign)> = closed_ids
        .iter()
        .map(|&s| {
            let faces: Vec<[u32; 3]> = members[s].iter().map(|&f| tris[f as usize]).collect();
            let pts: Vec<[f64; 3]> = faces.iter().flatten().map(|&i| v[i as usize]).collect();
            let vol = shell_volume_sign(v, &faces);
            (faces, Aabb::from_points(&pts), vol)
        })
        .collect();
    let wrong: Vec<bool> = (0..closed.len())
        .into_par_iter()
        .map(|s| {
            let (faces, bounds, vol) = &closed[s];
            let depth = (0..closed.len())
                .filter(|&t| t != s)
                .filter(|&t| {
                    let (tf, tb, _) = &closed[t];
                    let fits =
                        (0..3).all(|k| tb.min[k] <= bounds.min[k] && bounds.max[k] <= tb.max[k]);
                    fits && shell_location(v, faces, tf) == Some(-1)
                })
                .count();
            (depth % 2 == 0) != (*vol == Sign::Positive)
        })
        .collect();
    for (s, w) in wrong.into_iter().enumerate() {
        if w {
            problems.push(MeshProblem {
                face: Some(members[closed_ids[s]][0]),
                ..problem(InwardShell)
            });
        }
    }

    let boxes: Vec<Aabb> = (0..tris.len() as u32)
        .map(|f| {
            if degenerate[f as usize] {
                Aabb::empty()
            } else {
                Aabb::from_points(&mesh.corners(f))
            }
        })
        .collect();
    let mut pairs = vec![];
    Bvh::new(&boxes).self_pairs(|f, g| pairs.push((f, g)));
    pairs.sort_unstable();
    let found: Vec<Option<MeshProblem>> = pairs
        .par_iter()
        .map(|&(f, g)| {
            pair_problem(mesh, f, g).map(|kind| MeshProblem {
                face: Some(g),
                other: Some(f),
                ..problem(kind)
            })
        })
        .collect();
    problems.extend(found.into_iter().flatten());

    problems.sort_by_key(|p| (p.kind, p.face, p.vertex, p.other, p.edge));
    problems.dedup();
    let count = |kind| problems.iter().filter(|p| p.kind == kind).count();
    let summary = MeshSummary {
        degenerate_faces: count(DegenerateFace),
        duplicate_faces: count(DuplicateFace),
        duplicate_vertices: count(DuplicateVertex),
        boundary_edges: count(BoundaryEdge),
        non_manifold_edges: count(NonManifoldEdge),
        non_manifold_vertices: count(NonManifoldVertex),
        inconsistent_edges: count(InconsistentWinding),
        shells: shells.count,
        inward_shells: count(InwardShell),
        self_intersections: count(SelfIntersection),
        coplanar_overlaps: count(CoplanarOverlap),
        unused_vertices: count(UnusedVertex),
        touching_shell_edges,
        touching_shell_vertices: edge_ends.len(),
        is_closed: degenerate.contains(&false)
            && count(BoundaryEdge) == 0
            && count(NonManifoldEdge) == 0,
    };
    Ok(MeshReport { summary, problems })
}

/// The problem a pair of non-degenerate faces `f < g` forms, if any. Corners
/// are shared when they have the same index or exactly equal coordinates.
fn pair_problem(mesh: &Mesh, f: u32, g: u32) -> Option<MeshProblemKind> {
    let (a, b) = (mesh.corners(f), mesh.corners(g));
    let (ta, tb) = (mesh.triangles()[f as usize], mesh.triangles()[g as usize]);
    let same = |i: usize, j: usize| ta[i] == tb[j] || key(a[i]) == key(b[j]);
    let shared: Vec<(usize, usize)> = (0..3)
        .flat_map(|i| (0..3).map(move |j| (i, j)))
        .filter(|&(i, j)| same(i, j))
        .collect();
    let coplanar = b.iter().all(|&p| orient3d(a[0], a[1], a[2], p) == 0.0);
    let k = coplanar.then(|| tri_tri(a, b));
    if k == Some(TriTri::CoplanarOverlap) {
        let same_corners = ta.iter().all(|i| tb.contains(i));
        return (!same_corners).then_some(MeshProblemKind::CoplanarOverlap);
    }
    let hit = match shared[..] {
        [] => k.unwrap_or_else(|| tri_tri(a, b)) != TriTri::Disjoint,
        [(i, j)] => {
            let opposite = |t: [[f64; 3]; 3], c: usize| (t[(c + 1) % 3], t[(c + 2) % 3]);
            let (p, q) = opposite(a, i);
            let (r, s) = opposite(b, j);
            seg_tri(p, q, b).0 != SegTri::Miss || seg_tri(r, s, a).0 != SegTri::Miss
        }
        _ => false,
    };
    hit.then_some(MeshProblemKind::SelfIntersection)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fixtures::{self, broken};
    use proptest::prelude::*;

    fn summary(m: &Mesh) -> MeshSummary {
        validate(m, 0.0).unwrap().summary
    }

    #[test]
    fn defects_list_only_nonzero_defect_counts() {
        assert_eq!(summary(&fixtures::cube()).defects(), "no defects");
        let open = MeshSummary {
            is_closed: false,
            shells: 1,
            unused_vertices: 2,
            ..MeshSummary::default()
        };
        assert_eq!(open.defects(), "not closed");
        let bad = MeshSummary {
            boundary_edges: 18,
            self_intersections: 29,
            touching_shell_edges: 1,
            ..open
        };
        assert_eq!(bad.defects(), "boundary_edges=18, self_intersections=29");
    }

    #[test]
    fn clean_solids_are_solid() {
        for m in [
            fixtures::cube(),
            fixtures::tetrahedron(),
            fixtures::sphere(2),
            fixtures::torus(12, 8),
        ] {
            let r = validate(&m, 0.0).unwrap();
            assert!(r.problems.is_empty(), "{:?}", r.problems);
            assert!(r.summary.is_solid());
            assert_eq!(r.summary.shells, 1);
        }
        let nested = fixtures::nested_shells();
        let s = summary(&nested);
        assert!(s.is_solid(), "{s:?}");
        assert_eq!(s.shells, 2);
    }

    #[test]
    fn removing_one_face_opens_three_edges() {
        for m in [fixtures::cube(), fixtures::sphere(1)] {
            let keep: Vec<bool> = (0..m.n_faces()).map(|f| f != 0).collect();
            let s = summary(&m.filter_faces(&keep).unwrap());
            assert_eq!(s.boundary_edges, 3);
            assert!(!s.is_closed && !s.is_solid());
        }
    }

    #[test]
    fn broken_meshes_report_their_problem() {
        use MeshProblemKind::*;
        let cases = [
            (broken::degenerate_face(), DegenerateFace),
            (broken::duplicate_face(), DuplicateFace),
            (broken::duplicate_vertex(), DuplicateVertex),
            (fixtures::open_cube(), BoundaryEdge),
            (broken::non_manifold_fin(), NonManifoldEdge),
            (broken::bowtie(), NonManifoldVertex),
            (broken::flipped_face(), InconsistentWinding),
            (broken::inward_cube(), InwardShell),
            (broken::overlapping_cubes(), SelfIntersection),
            (broken::folded_face(), CoplanarOverlap),
            (broken::unused_vertex(), UnusedVertex),
        ];
        for (m, kind) in cases {
            let r = validate(&m, 0.0).unwrap();
            assert!(
                r.problems.iter().any(|p| p.kind == kind),
                "{kind:?} not in {:?}",
                r.problems
            );
            assert!(!r.summary.is_solid() || kind == UnusedVertex, "{kind:?}");
        }
        assert!(summary(&broken::unused_vertex()).is_solid());
    }

    #[test]
    fn counts_of_the_broken_cases() {
        let s = summary(&broken::flipped_face());
        assert_eq!(s.inconsistent_edges, 3);
        let s = summary(&broken::bowtie());
        assert_eq!((s.non_manifold_vertices, s.shells), (1, 2), "{s:?}");
        let s = summary(&broken::overlapping_cubes());
        assert!(s.self_intersections > 0 && s.is_closed);
        let s = summary(&broken::touching_cubes());
        assert!(
            s.self_intersections > 0 && s.coplanar_overlaps > 0 && !s.is_solid(),
            "{s:?}"
        );
        for (name, a, b) in fixtures::coplanar_pairs() {
            let s = summary(&Mesh::concat(&[&a, &b]).unwrap());
            assert!(
                s.self_intersections + s.coplanar_overlaps > 0,
                "{name}: {s:?}"
            );
        }
        let s = summary(&broken::non_manifold_fin());
        assert_eq!(s.non_manifold_edges, 1);
    }

    #[test]
    fn closed_shells_touching_at_edges_and_vertices_are_solid() {
        for (name, m, edges, points) in [
            ("tetrahedra by index", fixtures::touching_tetrahedra(), 0, 1),
            ("cubes at a vertex", fixtures::vertex_touching_cubes(), 0, 1),
            (
                "cubes at an edge",
                fixtures::edge_touching_cubes(false),
                1,
                2,
            ),
            (
                "welded cubes at an edge",
                fixtures::edge_touching_cubes(true),
                1,
                2,
            ),
        ] {
            let r = validate(&m, 0.0).unwrap();
            let s = &r.summary;
            assert!(s.is_solid(), "{name}: {:?}", r.problems);
            assert!(r.problems.is_empty(), "{name}: {:?}", r.problems);
            assert_eq!(s.shells, 2, "{name}");
            assert_eq!(
                (s.touching_shell_edges, s.touching_shell_vertices),
                (edges, points),
                "{name}"
            );
        }
        let s = summary(&fixtures::cube());
        assert_eq!((s.touching_shell_edges, s.touching_shell_vertices), (0, 0));
    }

    /// `m` with a separate copy of each vertex per umbrella of faces around it.
    fn split_umbrellas(m: &Mesh) -> Mesh {
        let edges = EdgeMap::new(m, None);
        let shells = Shells::new(m, &edges);
        let tris = m.triangles();
        let mut parent: Vec<u32> = (0..3 * tris.len() as u32).collect();
        let corner =
            |f: u32, i: u32| 3 * f + tris[f as usize].iter().position(|&x| x == i).unwrap() as u32;
        for id in 0..edges.edges.len() as u32 {
            for [f, g] in shells.links(&edges, id) {
                for i in edges.edges[id as usize] {
                    join(&mut parent, corner(f, i), corner(g, i));
                }
            }
        }
        let mut v = vec![];
        let mut ids: HashMap<u32, u32> = HashMap::new();
        let t = (0..tris.len() as u32)
            .map(|f| {
                tris[f as usize].map(|i| {
                    let r = root(&mut parent, corner(f, i));
                    *ids.entry(r).or_insert_with(|| {
                        v.push(m.vertices()[i as usize]);
                        v.len() as u32 - 1
                    })
                })
            })
            .collect();
        Mesh::new(v, t).unwrap()
    }

    #[test]
    fn a_closed_shell_touching_itself_is_solid() {
        for (name, cells, edges, points) in [
            ("edge", &fixtures::EDGE_LOOP[..], 1, 2),
            ("vertex", &fixtures::VERTEX_CHAIN[..], 0, 1),
        ] {
            let shared = fixtures::voxels(cells);
            let apart = split_umbrellas(&shared);
            assert!(apart.n_vertices() > shared.n_vertices(), "{name}");
            for m in [shared, apart] {
                let r = validate(&m, 0.0).unwrap();
                let s = &r.summary;
                assert!(
                    s.is_solid() && r.problems.is_empty(),
                    "{name}: {:?}",
                    r.problems
                );
                assert_eq!(s.shells, 1, "{name}");
                assert_eq!(
                    (s.touching_shell_edges, s.touching_shell_vertices),
                    (edges, points),
                    "{name}"
                );
            }
        }
    }

    #[test]
    fn a_shell_touching_itself_along_an_edge_with_single_umbrella_ends_is_solid() {
        let mut cells = fixtures::EDGE_LOOP.to_vec();
        cells.extend([[0, 0, -1], [1, 0, -1], [1, 1, -1]]);
        let r = validate(&fixtures::voxels(&cells), 0.0).unwrap();
        let s = &r.summary;
        assert!(s.is_solid() && r.problems.is_empty(), "{:?}", r.problems);
        assert_eq!(
            (s.shells, s.touching_shell_edges, s.touching_shell_vertices),
            (1, 1, 2)
        );
    }

    #[test]
    fn a_shell_touching_itself_with_overlapping_interiors_stays_a_defect() {
        let mut cells = fixtures::EDGE_LOOP.to_vec();
        cells.push([1, 0, 0]);
        let whole = fixtures::voxels(&cells);
        assert!(summary(&whole).is_solid());
        let loop_ = fixtures::voxels(&fixtures::EDGE_LOOP);
        let both = Mesh::concat(&[&loop_, &fixtures::voxels(&[[1, 0, 0]])]).unwrap();
        assert!(!summary(&both).is_solid());
    }

    #[test]
    fn a_cavity_touching_at_every_vertex_nests_by_an_edge_point() {
        let m = fixtures::vertex_touching_cavity();
        let r = validate(&m, 0.0).unwrap();
        assert!(r.summary.is_solid(), "{:?}", r.problems);
        assert_eq!(
            (r.summary.shells, r.summary.touching_shell_vertices),
            (2, 4)
        );
        let outward: Vec<[u32; 3]> = m
            .triangles()
            .iter()
            .enumerate()
            .map(|(f, &[a, b, c])| if f < 12 { [a, b, c] } else { [a, c, b] })
            .collect();
        let s = summary(&Mesh::new(m.vertices().to_vec(), outward).unwrap());
        assert_eq!(s.inward_shells, 1, "{s:?}");
    }

    #[test]
    fn cavities_touching_along_an_edge_are_solid() {
        for (m, shells) in [
            (fixtures::edge_touching_cavity(), 2),
            (fixtures::edge_touching_cavities(), 3),
        ] {
            let r = validate(&m, 0.0).unwrap();
            assert!(
                r.summary.is_solid() && r.problems.is_empty(),
                "{:?}",
                r.problems
            );
            assert_eq!(
                (
                    r.summary.shells,
                    r.summary.touching_shell_edges,
                    r.summary.touching_shell_vertices
                ),
                (shells, 1, 2)
            );
        }
        let outward = Mesh::concat(&[
            &fixtures::box_mesh([0.0; 3], [4.0; 3]),
            &fixtures::box_mesh([1.0; 3], [2.0; 3]),
            &fixtures::flipped(&fixtures::box_mesh([2.0, 2.0, 1.0], [3.0, 3.0, 2.0])),
        ])
        .unwrap();
        assert!(!summary(&outward).is_solid());
    }

    #[test]
    fn contacts_with_overlapping_interiors_stay_defects() {
        let r = validate(&broken::edge_overlapping_shells(), 0.0).unwrap();
        assert!(!r.summary.is_solid(), "{:?}", r.summary);
        assert!(r.summary.self_intersections > 0, "{:?}", r.summary);
        assert_eq!(r.summary.touching_shell_edges, 0);
        let s = summary(&broken::touching_cubes());
        assert!(s.coplanar_overlaps > 0 && !s.is_solid(), "{s:?}");
        let shared = Mesh::concat(&[
            &fixtures::cube(),
            &fixtures::translated(&fixtures::cube(), [1.0, 0.0, 0.0]),
        ])
        .unwrap();
        let s = summary(&shared);
        assert!(s.coplanar_overlaps > 0 && !s.is_solid(), "{s:?}");
    }

    #[test]
    fn duplicates_inside_one_shell_stay_defects() {
        let m = fixtures::edge_touching_cubes(false);
        let mut v = m.vertices().to_vec();
        v.push(v[0]);
        let mut t = m.triangles().to_vec();
        let k = t.iter().position(|t| t.contains(&0)).unwrap();
        t[k] = t[k].map(|i| if i == 0 { v.len() as u32 - 1 } else { i });
        let s = summary(&Mesh::new(v, t).unwrap());
        assert!(s.duplicate_vertices > 0 && !s.is_solid(), "{s:?}");
    }

    #[test]
    fn radial_pairs_follow_the_wedges() {
        use crate::topology::{Radial, radial_pairs_at};
        let m = fixtures::edge_touching_cubes(false);
        let (p, q) = ([1.0, 1.0, 0.0], [1.0, 1.0, 1.0]);
        let on: Vec<u32> = (0..m.n_faces() as u32)
            .filter(|&f| {
                let c = m.corners(f);
                c.contains(&p) && c.contains(&q)
            })
            .collect();
        assert_eq!(on.len(), 4);
        let Radial::Paired(pairs) = radial_pairs_at(&m, &on, p, q) else {
            panic!("not paired")
        };
        for [a, b] in pairs {
            assert_eq!(on[a] < 12, on[b] < 12, "{pairs:?}", pairs = (on[a], on[b]));
        }
        let flipped = fixtures::flipped(&m);
        let Radial::Paired(again) = radial_pairs_at(&flipped, &on, p, q) else {
            panic!("not paired")
        };
        assert!(
            again.iter().all(|&[a, b]| (on[a] < 12) != (on[b] < 12)),
            "{again:?}"
        );
        let one: Vec<u32> = on
            .iter()
            .copied()
            .filter(|&f| f < 12)
            .chain([on[3]])
            .collect();
        assert_eq!(radial_pairs_at(&m, &one, p, q), Radial::Unbalanced);
    }

    #[test]
    fn soup_reports_duplicates_not_crossings() {
        let soup = fixtures::soup(&fixtures::cube());
        let s = summary(&soup);
        assert_eq!(s.self_intersections, 0);
        assert_eq!(s.coplanar_overlaps, 0);
        assert_eq!(s.boundary_edges, 36);
        assert_eq!(s.duplicate_vertices, 36 - 8);
    }

    #[test]
    fn tolerance_merges_near_vertices() {
        let v = [
            [0.0, 0.0, 0.0],
            [0.05, 0.0, 0.0],
            [0.09, 0.0, 0.0],
            [1.0, 0.0, 0.0],
        ];
        assert_eq!(earliest_within(&v, 0.1), vec![0, 0, 0, 3]);
        assert_eq!(earliest_within(&v, 0.0), vec![0, 1, 2, 3]);
        assert_eq!(
            earliest_within(&[[0.0; 3], [-0.0, 0.0, 0.0]], 0.0),
            vec![0, 0]
        );
        assert!(validate(&fixtures::cube(), -1.0).is_err());
    }

    #[test]
    fn pairs_and_names_are_stable() {
        let p = MeshSummary::default().to_pairs();
        assert_eq!(p.len(), 15);
        assert_eq!(p.last().unwrap().0, "is_closed");
        assert_eq!(MeshProblemKind::CoplanarOverlap.name(), "coplanar_overlap");
    }

    fn pair_problem_always_tri_tri(mesh: &Mesh, f: u32, g: u32) -> Option<MeshProblemKind> {
        let (a, b) = (mesh.corners(f), mesh.corners(g));
        let (ta, tb) = (mesh.triangles()[f as usize], mesh.triangles()[g as usize]);
        let same = |i: usize, j: usize| ta[i] == tb[j] || key(a[i]) == key(b[j]);
        let shared: Vec<(usize, usize)> = (0..3)
            .flat_map(|i| (0..3).map(move |j| (i, j)))
            .filter(|&(i, j)| same(i, j))
            .collect();
        let k = tri_tri(a, b);
        if k == TriTri::CoplanarOverlap {
            return (!ta.iter().all(|i| tb.contains(i)))
                .then_some(MeshProblemKind::CoplanarOverlap);
        }
        let hit = match shared[..] {
            [] => k != TriTri::Disjoint,
            [(i, j)] => {
                let opposite = |t: [[f64; 3]; 3], c: usize| (t[(c + 1) % 3], t[(c + 2) % 3]);
                let (p, q) = opposite(a, i);
                let (r, s) = opposite(b, j);
                seg_tri(p, q, b).0 != SegTri::Miss || seg_tri(r, s, a).0 != SegTri::Miss
            }
            _ => false,
        };
        hit.then_some(MeshProblemKind::SelfIntersection)
    }

    #[test]
    fn pair_problem_skips_tri_tri_only_where_it_cannot_change_the_answer() {
        let mut rng = fixtures::Rng::new(7);
        let mut kinds = HashMap::new();
        for _ in 0..20_000 {
            let v: Vec<[f64; 3]> = (0..6)
                .map(|_| [0, 1, 2].map(|_| rng.below(3) as f64))
                .collect();
            let tb = [0, 1, 2].map(|_| rng.below(6) as u32);
            let tris = vec![[0, 1, 2], tb];
            let Ok(m) = Mesh::new(v, tris) else { continue };
            if degenerate_mask(&m).contains(&true) {
                continue;
            }
            let got = pair_problem(&m, 0, 1);
            assert_eq!(got, pair_problem_always_tri_tri(&m, 0, 1), "{:?}", m);
            *kinds.entry(got).or_insert(0) += 1;
        }
        assert!(
            kinds.len() == 3 && kinds.values().all(|&n| n > 100),
            "{kinds:?}"
        );
    }

    #[test]
    fn same_result_on_one_and_four_threads() {
        let m = broken::overlapping_cubes();
        let run = |n| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(n)
                .build()
                .unwrap()
                .install(|| validate(&m, 0.0).unwrap())
        };
        assert_eq!(run(1), run(4));
    }

    fn permuted(m: &Mesh, seed: u64) -> Mesh {
        let mut rng = fixtures::Rng::new(seed);
        let mut fo: Vec<usize> = (0..m.n_faces()).collect();
        rng.shuffle(&mut fo);
        let mut vo: Vec<usize> = (0..m.n_vertices()).collect();
        rng.shuffle(&mut vo);
        let mut inv = vec![0u32; vo.len()];
        for (new, &old) in vo.iter().enumerate() {
            inv[old] = new as u32;
        }
        let verts = vo.iter().map(|&i| m.vertices()[i]).collect();
        let tris = fo
            .iter()
            .map(|&f| m.triangles()[f].map(|i| inv[i as usize]))
            .collect();
        Mesh::new(verts, tris).unwrap()
    }

    proptest! {
        #[test]
        fn counts_survive_relabelling_and_integer_shifts(seed in 0u64..1000, which in 0usize..6, shift in proptest::array::uniform3(-1000i32..1000)) {
            let m = [
                fixtures::cube(),
                broken::flipped_face(),
                broken::overlapping_cubes(),
                broken::non_manifold_fin(),
                fixtures::nested_shells(),
                fixtures::pit_surface(6, 7),
            ][which].clone();
            let base = summary(&m);
            prop_assert_eq!(&summary(&permuted(&m, seed)), &base);
            let moved: Vec<[f64; 3]> = m.vertices().iter().map(|p| [0, 1, 2].map(|k| p[k] + shift[k] as f64)).collect();
            prop_assert_eq!(&summary(&m.clone().with_vertices(moved).unwrap()), &base);
        }
    }
}
