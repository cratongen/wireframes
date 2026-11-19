//! Seeded synthetic meshes for tests in every wireframes crate. Solids wind
//! counter-clockwise seen from outside.

use crate::Mesh;

/// SplitMix64: small, seeded, and identical on every platform.
#[derive(Clone, Debug)]
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Rng {
        Rng(seed)
    }

    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform in `[lo, hi)`.
    pub fn uniform(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * ((self.next_u64() >> 11) as f64 / (1u64 << 53) as f64)
    }

    pub fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }

    pub fn shuffle<T>(&mut self, v: &mut [T]) {
        for i in (1..v.len()).rev() {
            v.swap(i, self.below(i + 1));
        }
    }
}

fn mesh(vertices: Vec<[f64; 3]>, triangles: Vec<[u32; 3]>) -> Mesh {
    Mesh::new(vertices, triangles).expect("fixture is valid")
}

/// Axis-aligned box; faces in the order bottom, top, front (`-y`), back, left (`-x`), right, two each.
pub fn box_mesh(min: [f64; 3], max: [f64; 3]) -> Mesh {
    let v = (0..8)
        .map(|i| {
            [
                if i & 1 == 0 { min[0] } else { max[0] },
                if i & 2 == 0 { min[1] } else { max[1] },
                if i & 4 == 0 { min[2] } else { max[2] },
            ]
        })
        .collect();
    let t = vec![
        [0, 2, 3],
        [0, 3, 1],
        [4, 5, 7],
        [4, 7, 6],
        [0, 1, 5],
        [0, 5, 4],
        [2, 6, 7],
        [2, 7, 3],
        [0, 4, 6],
        [0, 6, 2],
        [1, 3, 7],
        [1, 7, 5],
    ];
    mesh(v, t)
}

/// The unit cube `[0, 1]^3`.
pub fn cube() -> Mesh {
    box_mesh([0.0; 3], [1.0; 3])
}

/// Corners at the origin and the three unit points.
pub fn tetrahedron() -> Mesh {
    mesh(
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
        vec![[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]],
    )
}

/// Unit icosphere at the origin, `subdivisions` times split (20 * 4^s faces).
pub fn sphere(subdivisions: u32) -> Mesh {
    let t = (1.0 + 5f64.sqrt()) / 2.0;
    let mut v: Vec<[f64; 3]> = vec![
        [-1.0, t, 0.0],
        [1.0, t, 0.0],
        [-1.0, -t, 0.0],
        [1.0, -t, 0.0],
        [0.0, -1.0, t],
        [0.0, 1.0, t],
        [0.0, -1.0, -t],
        [0.0, 1.0, -t],
        [t, 0.0, -1.0],
        [t, 0.0, 1.0],
        [-t, 0.0, -1.0],
        [-t, 0.0, 1.0],
    ];
    let mut f: Vec<[u32; 3]> = vec![
        [0, 11, 5],
        [0, 5, 1],
        [0, 1, 7],
        [0, 7, 10],
        [0, 10, 11],
        [1, 5, 9],
        [5, 11, 4],
        [11, 10, 2],
        [10, 7, 6],
        [7, 1, 8],
        [3, 9, 4],
        [3, 4, 2],
        [3, 2, 6],
        [3, 6, 8],
        [3, 8, 9],
        [4, 9, 5],
        [2, 4, 11],
        [6, 2, 10],
        [8, 6, 7],
        [9, 8, 1],
    ];
    let unit = |p: [f64; 3]| {
        let n = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        p.map(|x| x / n)
    };
    v = v.into_iter().map(unit).collect();
    for _ in 0..subdivisions {
        let mut mid = std::collections::HashMap::new();
        let mut split = |a: u32, b: u32, v: &mut Vec<[f64; 3]>| {
            *mid.entry((a.min(b), a.max(b))).or_insert_with(|| {
                let (p, q) = (v[a as usize], v[b as usize]);
                v.push(unit([0, 1, 2].map(|k| 0.5 * (p[k] + q[k]))));
                v.len() as u32 - 1
            })
        };
        let mut next = Vec::with_capacity(4 * f.len());
        for [a, b, c] in f {
            let (ab, bc, ca) = (
                split(a, b, &mut v),
                split(b, c, &mut v),
                split(c, a, &mut v),
            );
            next.extend([[a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]]);
        }
        f = next;
    }
    mesh(v, f)
}

/// Torus around the z axis with radii 2 and 0.5, `major x minor` quads.
pub fn torus(major: u32, minor: u32) -> Mesh {
    let (big, small) = (2.0, 0.5);
    let mut v = vec![];
    for i in 0..major {
        let th = i as f64 / major as f64 * std::f64::consts::TAU;
        for j in 0..minor {
            let ph = j as f64 / minor as f64 * std::f64::consts::TAU;
            let r = big + small * ph.cos();
            v.push([r * th.cos(), r * th.sin(), small * ph.sin()]);
        }
    }
    let id = |i: u32, j: u32| (i % major) * minor + (j % minor);
    let mut f = vec![];
    for i in 0..major {
        for j in 0..minor {
            let (a, b, c, d) = (id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1));
            f.push([a, b, c]);
            f.push([a, c, d]);
        }
    }
    mesh(v, f)
}

/// Every face reversed.
pub fn flipped(m: &Mesh) -> Mesh {
    let t = m.triangles().iter().map(|&[a, b, c]| [a, c, b]).collect();
    mesh(m.vertices().to_vec(), t)
}

pub fn translated(m: &Mesh, d: [f64; 3]) -> Mesh {
    let v = m
        .vertices()
        .iter()
        .map(|p| [p[0] + d[0], p[1] + d[1], p[2] + d[2]])
        .collect();
    m.clone().with_vertices(v).expect("same count")
}

/// The box `[0, 3]^3` with the cavity `[1, 2]^3` wound inward.
pub fn nested_shells() -> Mesh {
    let outer = box_mesh([0.0; 3], [3.0; 3]);
    let inner = flipped(&box_mesh([1.0; 3], [2.0; 3]));
    Mesh::concat(&[&outer, &inner]).expect("no labels")
}

/// Each face with its own three vertices.
pub fn soup(m: &Mesh) -> Mesh {
    let v = (0..m.n_faces() as u32).flat_map(|f| m.corners(f)).collect();
    let t = (0..m.n_faces() as u32)
        .map(|f| [3 * f, 3 * f + 1, 3 * f + 2])
        .collect();
    mesh(v, t)
}

/// The unit cube without its top face (4 boundary edges).
pub fn open_cube() -> Mesh {
    let keep: Vec<bool> = (0..12).map(|f| f != 2 && f != 3).collect();
    cube().filter_faces(&keep).expect("12 faces")
}

/// Open height field over `[0, 100]^2` on an `n x n` vertex grid: a conical
/// pit 30 deep around `(50, 50)` at elevation 100, with seeded noise of 0.5.
/// Faces wind counter-clockwise in plan.
pub fn pit_surface(n: u32, seed: u64) -> Mesh {
    assert!(n >= 2, "a grid needs at least 2 x 2 vertices");
    let mut rng = Rng::new(seed);
    let step = 100.0 / (n - 1) as f64;
    let mut v = vec![];
    for j in 0..n {
        for i in 0..n {
            let (x, y) = (i as f64 * step, j as f64 * step);
            let r = ((x - 50.0).powi(2) + (y - 50.0).powi(2)).sqrt();
            let z = 100.0 - 30.0 * (1.0 - r / 40.0).max(0.0) + rng.uniform(-0.5, 0.5);
            v.push([x, y, z]);
        }
    }
    let id = |i: u32, j: u32| j * n + i;
    let mut f = vec![];
    for j in 0..n - 1 {
        for i in 0..n - 1 {
            let (a, b, c, d) = (id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1));
            f.push([a, b, c]);
            f.push([a, c, d]);
        }
    }
    mesh(v, f)
}

/// `n` points uniform in the box `[min, max)`.
pub fn random_points(seed: u64, n: usize, min: [f64; 3], max: [f64; 3]) -> Vec<[f64; 3]> {
    let mut rng = Rng::new(seed);
    (0..n)
        .map(|_| [0, 1, 2].map(|k| rng.uniform(min[k], max[k])))
        .collect()
}

/// Pairs of solids with coplanar faces, for booleans: `shared_face` (unit
/// cubes side by side), `half_shift` (shifted by half a cell in x and y),
/// `nested` (a smaller box on the bottom face, inside), `opposite_face` (a
/// smaller box standing on the top face).
pub fn coplanar_pairs() -> Vec<(&'static str, Mesh, Mesh)> {
    vec![
        ("shared_face", cube(), translated(&cube(), [1.0, 0.0, 0.0])),
        ("half_shift", cube(), translated(&cube(), [0.5, 0.5, 0.0])),
        (
            "nested",
            cube(),
            box_mesh([0.25, 0.25, 0.0], [0.75, 0.75, 0.5]),
        ),
        (
            "opposite_face",
            cube(),
            box_mesh([0.25, 0.25, 1.0], [0.75, 0.75, 2.0]),
        ),
    ]
}

/// Two tetrahedra sharing only vertex 0, by index: closed shells touching
/// at a vertex.
pub fn touching_tetrahedra() -> Mesh {
    let a = tetrahedron();
    let mut v = a.vertices().to_vec();
    v.extend([[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]);
    let mut t = a.triangles().to_vec();
    t.extend([[0, 4, 5], [0, 6, 4], [0, 5, 6], [4, 6, 5]]);
    mesh(v, t)
}

/// The unit cube and its copy moved by `(1, 1, 0)`: closed shells touching
/// along the edge from `(1, 1, 0)` to `(1, 1, 1)`. `welded`: the two shells
/// share the contact vertices by index, else they hold coincident copies.
pub fn edge_touching_cubes(welded: bool) -> Mesh {
    let m = Mesh::concat(&[&cube(), &translated(&cube(), [1.0, 1.0, 0.0])]).expect("no labels");
    if !welded {
        return m;
    }
    let first = crate::validate::earliest_within(m.vertices(), 0.0);
    let t = m
        .triangles()
        .iter()
        .map(|t| t.map(|i| first[i as usize]))
        .collect();
    let welded = mesh(m.vertices().to_vec(), t);
    welded
        .filter_faces(&vec![true; welded.n_faces()])
        .expect("same faces")
}

/// Boundary of a union of unit cells `[i, i + 1] x [j, j + 1] x [k, k + 1]`,
/// one vertex index per lattice point, faces wound outward.
pub fn voxels(cells: &[[i32; 3]]) -> Mesh {
    let has = |c: [i32; 3]| cells.contains(&c);
    let mut ids: std::collections::BTreeMap<[i32; 3], u32> = Default::default();
    let mut v = vec![];
    let mut t = vec![];
    for &c in cells {
        for axis in 0..3 {
            for side in [0, 1] {
                let mut n = c;
                n[axis] += 2 * side - 1;
                if has(n) {
                    continue;
                }
                let (a, b) = ((axis + 1) % 3, (axis + 2) % 3);
                let q: Vec<u32> = [(0, 0), (1, 0), (1, 1), (0, 1)]
                    .iter()
                    .map(|&(da, db)| {
                        let mut p = c;
                        p[axis] += side;
                        p[a] += da;
                        p[b] += db;
                        *ids.entry(p).or_insert_with(|| {
                            v.push(p.map(f64::from));
                            v.len() as u32 - 1
                        })
                    })
                    .collect();
                if side == 1 {
                    t.extend([[q[0], q[1], q[2]], [q[0], q[2], q[3]]]);
                } else {
                    t.extend([[q[0], q[2], q[1]], [q[0], q[3], q[2]]]);
                }
            }
        }
    }
    mesh(v, t)
}

/// Five cells in a loop whose ends touch along the edge from `(1, 1, 0)` to
/// `(1, 1, 1)`: one closed shell touching itself along an edge.
pub const EDGE_LOOP: [[i32; 3]; 5] = [[0, 0, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [1, 1, 0]];

/// Eight cells in a chain whose ends touch at the point `(1, 1, 1)`: one
/// closed shell touching itself at a vertex.
pub const VERTEX_CHAIN: [[i32; 3]; 8] = [
    [0, 0, 0],
    [-1, 0, 0],
    [-1, 0, 1],
    [-1, 0, 2],
    [0, 0, 2],
    [1, 0, 2],
    [1, 1, 2],
    [1, 1, 1],
];

/// The unit cube and its copy moved by `(1, 1, 1)`: closed shells touching
/// at the vertex `(1, 1, 1)`, with coincident copies of it.
pub fn vertex_touching_cubes() -> Mesh {
    Mesh::concat(&[&cube(), &translated(&cube(), [1.0, 1.0, 1.0])]).expect("no labels")
}

/// The tetrahedron scaled by 3 with each face split at its centroid, and the
/// cavity spanned by the four centroids, wound inward: a cavity touching the
/// outer shell at every vertex it has, with coincident copies. Volume
/// `4.5 - 1 / 6`.
pub fn vertex_touching_cavity() -> Mesh {
    let t = tetrahedron();
    let mut v: Vec<[f64; 3]> = t.vertices().iter().map(|p| p.map(|x| 3.0 * x)).collect();
    let mut f = vec![];
    for &[a, b, c] in t.triangles() {
        let i = v.len() as u32;
        v.push([0, 1, 2].map(|k| (v[a as usize][k] + v[b as usize][k] + v[c as usize][k]) / 3.0));
        f.extend([[a, b, i], [b, c, i], [c, a, i]]);
    }
    let inner: Vec<[f64; 3]> = v[4..8].to_vec();
    let base = v.len() as u32;
    v.extend_from_slice(&inner);
    for skip in 0..4 {
        let [a, b, c] = [0, 1, 2].map(|k| (skip + 1 + k) % 4);
        let p = [a, b, c, skip].map(|i| inner[i]);
        let d = |i: usize| [0, 1, 2].map(|k| p[i][k] - p[0][k]);
        let (u, w, x) = (d(1), d(2), d(3));
        let n = [
            u[1] * w[2] - u[2] * w[1],
            u[2] * w[0] - u[0] * w[2],
            u[0] * w[1] - u[1] * w[0],
        ];
        let inward = n[0] * x[0] + n[1] * x[1] + n[2] * x[2] > 0.0;
        let t = [a, b, c].map(|i| base + i as u32);
        f.push(if inward { t } else { [t[0], t[2], t[1]] });
    }
    mesh(v, f)
}

/// Meshes with one deliberate defect each.
pub mod broken {
    use super::*;

    fn cube_plus(extra_vertices: &[[f64; 3]], extra_faces: &[[u32; 3]]) -> Mesh {
        let c = cube();
        let mut v = c.vertices().to_vec();
        v.extend_from_slice(extra_vertices);
        let mut t = c.triangles().to_vec();
        t.extend_from_slice(extra_faces);
        mesh(v, t)
    }

    /// A cube plus a face whose corners are collinear.
    pub fn degenerate_face() -> Mesh {
        cube_plus(&[[0.5, 0.0, 0.0]], &[[0, 8, 1]])
    }

    pub fn duplicate_face() -> Mesh {
        cube_plus(&[], &[[2, 3, 0]])
    }

    /// A cube plus an unused copy of vertex 0.
    pub fn duplicate_vertex() -> Mesh {
        cube_plus(&[[0.0; 3]], &[])
    }

    pub fn unused_vertex() -> Mesh {
        cube_plus(&[[5.0, 5.0, 5.0]], &[])
    }

    /// A cube plus a triangle hanging off edge 0-1 (three faces on that edge).
    pub fn non_manifold_fin() -> Mesh {
        cube_plus(&[[0.5, -1.0, -1.0]], &[[0, 1, 8]])
    }

    /// Two triangles sharing only vertex 0: an open surface pinched there.
    pub fn bowtie() -> Mesh {
        mesh(
            vec![
                [0.0; 3],
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [-1.0, 0.0, 0.0],
                [-1.0, -1.0, 0.0],
            ],
            vec![[0, 1, 2], [0, 3, 4]],
        )
    }

    /// The unit cube and a tetrahedron inside it that shares the cube's
    /// edge from `(1, 1, 0)` to `(1, 1, 1)`: closed shells in contact along
    /// an edge with overlapping interiors.
    pub fn edge_overlapping_shells() -> Mesh {
        let t = mesh(
            vec![
                [1.0, 1.0, 0.0],
                [1.0, 1.0, 1.0],
                [0.5, 0.75, 0.5],
                [0.75, 0.5, 0.5],
            ],
            vec![[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]],
        );
        let t = if t.validate(0.0).expect("valid").summary.inward_shells > 0 {
            flipped(&t)
        } else {
            t
        };
        Mesh::concat(&[&cube(), &t]).expect("no labels")
    }

    /// A cube with face 0 reversed (3 inconsistent edges).
    pub fn flipped_face() -> Mesh {
        let c = cube();
        let mut t = c.triangles().to_vec();
        t[0] = [t[0][0], t[0][2], t[0][1]];
        mesh(c.vertices().to_vec(), t)
    }

    pub fn inward_cube() -> Mesh {
        flipped(&cube())
    }

    /// Two unit cubes overlapping in volume.
    pub fn overlapping_cubes() -> Mesh {
        Mesh::concat(&[&cube(), &translated(&cube(), [0.5, 0.3, 0.2])]).expect("no labels")
    }

    /// Two unit cubes touching face to face, with separate vertices.
    pub fn touching_cubes() -> Mesh {
        Mesh::concat(&[&cube(), &translated(&cube(), [1.0, 0.5, 0.0])]).expect("no labels")
    }

    /// An open surface: a triangle and a coplanar neighbour folded back over it.
    pub fn folded_face() -> Mesh {
        mesh(
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.25, 0.25, 0.0],
            ],
            vec![[0, 1, 2], [1, 0, 3]],
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rng_is_seeded_and_uniform_in_range() {
        let (mut a, mut b) = (Rng::new(7), Rng::new(7));
        assert_eq!(a.next_u64(), b.next_u64());
        for _ in 0..1000 {
            let x = a.uniform(-2.0, 3.0);
            assert!((-2.0..3.0).contains(&x));
        }
        let mut v: Vec<u32> = (0..50).collect();
        a.shuffle(&mut v);
        let mut s = v.clone();
        s.sort_unstable();
        assert_eq!(s, (0..50).collect::<Vec<_>>());
    }

    #[test]
    fn fixture_sizes() {
        assert_eq!(sphere(2).n_faces(), 320);
        assert_eq!(torus(12, 8).n_faces(), 192);
        assert_eq!(pit_surface(5, 1).n_faces(), 32);
        assert_eq!(pit_surface(5, 1), pit_surface(5, 1));
        assert_ne!(pit_surface(5, 1), pit_surface(5, 2));
        assert_eq!(random_points(3, 10, [0.0; 3], [1.0; 3]).len(), 10);
        assert_eq!(coplanar_pairs().len(), 4);
    }
}
