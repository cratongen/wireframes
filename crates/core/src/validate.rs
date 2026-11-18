use std::collections::HashMap;

use rayon::prelude::*;

use crate::intersect::{SegTri, TriTri, seg_tri, tri_tri};
use crate::predicates::collinear3;
use crate::topology::locate_point;
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
    pub is_closed: bool,
}

impl MeshSummary {
    /// Closed, and every count except `shells` and `unused_vertices` is zero.
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
    if tolerance == 0.0 {
        let mut first: HashMap<[u64; 3], u32> = HashMap::with_capacity(vertices.len());
        return vertices
            .iter()
            .enumerate()
            .map(|(i, &p)| *first.entry(key(p)).or_insert(i as u32))
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
                        if d2 <= t2 && j < best {
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

/// Signed volume of `faces` (times 6), relative to `o`.
fn volume6(v: &[[f64; 3]], faces: &[[u32; 3]], o: [f64; 3]) -> f64 {
    faces
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|i| {
                let p = v[i as usize];
                [p[0] - o[0], p[1] - o[1], p[2] - o[2]]
            });
            a[0] * (b[1] * c[2] - b[2] * c[1])
                + a[1] * (b[2] * c[0] - b[0] * c[2])
                + a[2] * (b[0] * c[1] - b[1] * c[0])
        })
        .sum()
}

/// Degenerate, duplicate, edge, vertex, orientation and intersection problems.
/// `tolerance` is the distance within which a vertex duplicates an earlier one.
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
    for (i, &w) in earliest_within(v, tolerance).iter().enumerate() {
        if w != i as u32 {
            problems.push(MeshProblem {
                vertex: Some(i as u32),
                other: Some(w),
                ..problem(DuplicateVertex)
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
    let local = |f: u32, id: u32| {
        edges.face_edges[f as usize]
            .iter()
            .position(|&e| e == id)
            .expect("face holds edge") as u8
    };
    let mut shell_of: Vec<u32> = (0..tris.len() as u32).collect();
    let mut inconsistent = vec![false; tris.len()];
    let mut open_or_pinched = vec![false; tris.len()];
    for id in 0..edges.edges.len() as u32 {
        let fs = edges.faces_of_edge(id);
        for &g in &fs[1..] {
            join(&mut shell_of, fs[0], g);
        }
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
            _ => {
                for &f in fs {
                    open_or_pinched[f as usize] = true;
                }
                problems.push(MeshProblem {
                    kind: NonManifoldEdge,
                    ..at(fs[0])
                });
            }
        }
    }

    for vert in 0..v.len() as u32 {
        let fan = edges.faces_of_vertex(vert);
        if fan.len() < 2 {
            continue;
        }
        let mut parent: Vec<u32> = (0..fan.len() as u32).collect();
        for (i, &f) in fan.iter().enumerate() {
            for (e, &id) in edges.face_edges[f as usize].iter().enumerate() {
                let t = tris[f as usize];
                if t[e] != vert && t[(e + 1) % 3] != vert {
                    continue;
                }
                for g in edges.faces_of_edge(id) {
                    let j = fan.binary_search(g).expect("g holds the vertex");
                    join(&mut parent, i as u32, j as u32);
                }
            }
        }
        if (0..fan.len() as u32)
            .filter(|&i| root(&mut parent, i) == i)
            .count()
            > 1
        {
            problems.push(MeshProblem {
                vertex: Some(vert),
                ..problem(NonManifoldVertex)
            });
        }
    }

    let mut shells: Vec<(u32, Vec<u32>)> = vec![];
    let mut index_of: HashMap<u32, usize> = HashMap::new();
    for f in valid() {
        let r = root(&mut shell_of, f as u32);
        let k = *index_of.entry(r).or_insert_with(|| {
            shells.push((r, vec![]));
            shells.len() - 1
        });
        shells[k].1.push(f as u32);
    }
    let origin = mesh.bounds().map_or([0.0; 3], |b| b.min);
    let closed: Vec<(Vec<[u32; 3]>, Aabb, f64)> = shells
        .iter()
        .filter(|(_, fs)| {
            fs.iter()
                .all(|&f| !inconsistent[f as usize] && !open_or_pinched[f as usize])
        })
        .map(|(_, fs)| {
            let faces: Vec<[u32; 3]> = fs.iter().map(|&f| tris[f as usize]).collect();
            let pts: Vec<[f64; 3]> = faces.iter().flatten().map(|&i| v[i as usize]).collect();
            let vol = volume6(v, &faces, origin);
            (faces, Aabb::from_points(&pts), vol)
        })
        .collect();
    let closed_first: Vec<u32> = shells
        .iter()
        .filter(|(_, fs)| {
            fs.iter()
                .all(|&f| !inconsistent[f as usize] && !open_or_pinched[f as usize])
        })
        .map(|(_, fs)| fs[0])
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
                    fits && faces
                        .iter()
                        .flatten()
                        .map(|&i| locate_point(v, tf, v[i as usize]))
                        .find(|&l| l != 0)
                        == Some(-1)
                })
                .count();
            (depth % 2 == 0) != (*vol > 0.0)
        })
        .collect();
    for (s, w) in wrong.into_iter().enumerate() {
        if w {
            problems.push(MeshProblem {
                face: Some(closed_first[s]),
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
    let count = |kind| problems.iter().filter(|p| p.kind == kind).count();
    let summary = MeshSummary {
        degenerate_faces: count(DegenerateFace),
        duplicate_faces: count(DuplicateFace),
        duplicate_vertices: count(DuplicateVertex),
        boundary_edges: count(BoundaryEdge),
        non_manifold_edges: count(NonManifoldEdge),
        non_manifold_vertices: count(NonManifoldVertex),
        inconsistent_edges: count(InconsistentWinding),
        shells: shells.len(),
        inward_shells: count(InwardShell),
        self_intersections: count(SelfIntersection),
        coplanar_overlaps: count(CoplanarOverlap),
        unused_vertices: count(UnusedVertex),
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
    let k = tri_tri(a, b);
    if k == TriTri::CoplanarOverlap {
        let same_corners = ta.iter().all(|i| tb.contains(i));
        return (!same_corners).then_some(MeshProblemKind::CoplanarOverlap);
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fixtures::{self, broken};
    use proptest::prelude::*;

    fn summary(m: &Mesh) -> MeshSummary {
        validate(m, 0.0).unwrap().summary
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
            (broken::touching_tetrahedra(), NonManifoldVertex),
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
        let s = summary(&broken::touching_tetrahedra());
        assert_eq!(
            (s.non_manifold_vertices, s.self_intersections, s.shells),
            (1, 0, 2)
        );
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
        assert_eq!(p.len(), 13);
        assert_eq!(p.last().unwrap().0, "is_closed");
        assert_eq!(MeshProblemKind::CoplanarOverlap.name(), "coplanar_overlap");
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
