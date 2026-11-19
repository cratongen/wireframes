use rayon::prelude::*;
use wireframes_core::validate::degenerate_mask;
use wireframes_core::{EdgeMap, Error, Mesh, Result, validate};

const CHUNK: usize = 4096;

fn neumaier<const K: usize>(items: impl Iterator<Item = [f64; K]>) -> [f64; K] {
    let mut s = [0.0; K];
    let mut c = [0.0; K];
    for x in items {
        for k in 0..K {
            let t = s[k] + x[k];
            c[k] += if s[k].abs() >= x[k].abs() {
                (s[k] - t) + x[k]
            } else {
                (x[k] - t) + s[k]
            };
            s[k] = t;
        }
    }
    std::array::from_fn(|k| s[k] + c[k])
}

/// Sum of `f(0..n)` in fixed chunks of 4096 with Neumaier compensation, chunk
/// partials added in chunk order, so the result does not depend on thread count.
pub(crate) fn sum_by<const K: usize>(n: usize, f: impl Fn(usize) -> [f64; K] + Sync) -> [f64; K] {
    let partials: Vec<[f64; K]> = (0..n.div_ceil(CHUNK))
        .into_par_iter()
        .map(|c| neumaier((c * CHUNK..n.min((c + 1) * CHUNK)).map(&f)))
        .collect();
    neumaier(partials.into_iter())
}

pub(crate) fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

pub(crate) fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

pub(crate) fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub(crate) fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

fn origin(mesh: &Mesh) -> [f64; 3] {
    mesh.bounds().map_or([0.0; 3], |b| b.min)
}

fn local_corners(mesh: &Mesh, f: usize, o: [f64; 3]) -> [[f64; 3]; 3] {
    mesh.corners(f as u32).map(|p| sub(p, o))
}

fn area2(c: [[f64; 3]; 3]) -> f64 {
    norm(cross(sub(c[1], c[0]), sub(c[2], c[0])))
}

pub fn area(mesh: &Mesh) -> f64 {
    let o = origin(mesh);
    0.5 * sum_by(mesh.n_faces(), |f| [area2(local_corners(mesh, f, o))])[0]
}

pub fn face_areas(mesh: &Mesh) -> Vec<f64> {
    (0..mesh.n_faces())
        .into_par_iter()
        .map(|f| 0.5 * area2(mesh.corners(f as u32)))
        .collect()
}

/// Unit normals; `[0, 0, 0]` for degenerate faces.
pub fn face_normals(mesh: &Mesh) -> Vec<[f64; 3]> {
    let degenerate = degenerate_mask(mesh);
    (0..mesh.n_faces())
        .into_par_iter()
        .map(|f| {
            let c = mesh.corners(f as u32);
            let n = cross(sub(c[1], c[0]), sub(c[2], c[0]));
            let l = norm(n);
            if degenerate[f] || l == 0.0 {
                [0.0; 3]
            } else {
                n.map(|x| x / l)
            }
        })
        .collect()
}

pub fn face_centroids(mesh: &Mesh) -> Vec<[f64; 3]> {
    (0..mesh.n_faces())
        .into_par_iter()
        .map(|f| {
            let [a, b, c] = mesh.corners(f as u32);
            std::array::from_fn(|k| a[k] + ((b[k] - a[k]) + (c[k] - a[k])) / 3.0)
        })
        .collect()
}

fn angle(u: [f64; 3], v: [f64; 3]) -> f64 {
    norm(cross(u, v)).atan2(dot(u, v))
}

/// Angle-weighted unit normals; `[0, 0, 0]` for isolated vertices.
pub fn vertex_normals(mesh: &Mesh) -> Vec<[f64; 3]> {
    let normals = face_normals(mesh);
    let mut acc = vec![[0.0; 3]; mesh.n_vertices()];
    for (f, t) in mesh.triangles().iter().enumerate() {
        let c = mesh.corners(f as u32);
        for k in 0..3 {
            let w = angle(sub(c[(k + 1) % 3], c[k]), sub(c[(k + 2) % 3], c[k]));
            for a in 0..3 {
                acc[t[k] as usize][a] += w * normals[f][a];
            }
        }
    }
    acc.into_iter()
        .map(|n| {
            let l = norm(n);
            if l == 0.0 { [0.0; 3] } else { n.map(|x| x / l) }
        })
        .collect()
}

/// Per face `(dip, dip direction)` in degrees: dip 0 to 90, direction 0 to 360
/// clockwise from `+y`; NaN direction for horizontal faces.
pub fn dip(mesh: &Mesh) -> Vec<(f64, f64)> {
    face_normals(mesh)
        .into_iter()
        .map(|n| {
            if n == [0.0; 3] {
                return (f64::NAN, f64::NAN);
            }
            let up = if n[2] < 0.0 { n.map(|x| -x) } else { n };
            let dip = up[0].hypot(up[1]).atan2(up[2]).to_degrees();
            let direction = if up[0] == 0.0 && up[1] == 0.0 {
                f64::NAN
            } else {
                let d = up[0].atan2(up[1]).to_degrees().rem_euclid(360.0);
                if d < 360.0 { d } else { 0.0 }
            };
            (dip, direction)
        })
        .collect()
}

/// `InvalidMesh` unless every edge of the non-degenerate faces has exactly two
/// faces that traverse it in opposite directions.
pub(crate) fn check_closed_oriented(mesh: &Mesh, what: &str) -> Result<Vec<bool>> {
    let degenerate = degenerate_mask(mesh);
    let edges = EdgeMap::new(mesh, Some(&degenerate));
    let tris = mesh.triangles();
    let forward = |f: u32, e: u32| {
        let t = tris[f as usize];
        let k = edges.face_edges[f as usize]
            .iter()
            .position(|&x| x == e)
            .expect("face holds edge");
        t[k] < t[(k + 1) % 3]
    };
    let ok = degenerate.contains(&false)
        && (0..edges.edges.len() as u32).all(|e| match *edges.faces_of_edge(e) {
            [a, b] => forward(a, e) != forward(b, e),
            _ => false,
        });
    if ok {
        Ok(degenerate)
    } else {
        Err(Error::InvalidMesh {
            context: format!("{what} needs a closed, consistently wound mesh"),
            report: Box::new(validate(mesh, 0.0)?),
        })
    }
}

fn det(c: [[f64; 3]; 3]) -> f64 {
    dot(c[0], cross(c[1], c[2]))
}

/// Signed volume; `InvalidMesh` unless closed and consistently wound.
pub fn volume(mesh: &Mesh) -> Result<f64> {
    let degenerate = check_closed_oriented(mesh, "volume")?;
    let o = origin(mesh);
    let s = sum_by(mesh.n_faces(), |f| {
        if degenerate[f] {
            [0.0]
        } else {
            [det(local_corners(mesh, f, o))]
        }
    });
    Ok(s[0] / 6.0)
}

/// Volume centroid, same precondition as `volume`; `InvalidInput` when the volume is 0.
pub fn centroid(mesh: &Mesh) -> Result<[f64; 3]> {
    let degenerate = check_closed_oriented(mesh, "centroid")?;
    let o = origin(mesh);
    let s = sum_by(mesh.n_faces(), |f| {
        if degenerate[f] {
            return [0.0; 4];
        }
        let c = local_corners(mesh, f, o);
        let d = det(c);
        [
            d * (c[0][0] + c[1][0] + c[2][0]),
            d * (c[0][1] + c[1][1] + c[2][1]),
            d * (c[0][2] + c[1][2] + c[2][2]),
            d,
        ]
    });
    if s[3] == 0.0 {
        return Err(Error::invalid("centroid needs a non-zero volume, got 0"));
    }
    Ok(std::array::from_fn(|k| o[k] + s[k] / (4.0 * s[3])))
}

/// `None` when the area is 0.
pub fn area_centroid(mesh: &Mesh) -> Option<[f64; 3]> {
    let o = origin(mesh);
    let s = sum_by(mesh.n_faces(), |f| {
        let c = local_corners(mesh, f, o);
        let a = area2(c);
        [
            a * (c[0][0] + c[1][0] + c[2][0]),
            a * (c[0][1] + c[1][1] + c[2][1]),
            a * (c[0][2] + c[1][2] + c[2][2]),
            a,
        ]
    });
    (s[3] != 0.0).then(|| std::array::from_fn(|k| o[k] + s[k] / (3.0 * s[3])))
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;
    use wireframes_core::fixtures::{self, broken};

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= 1e-9 * 1f64.max(a.abs()).max(b.abs())
    }

    #[test]
    fn dip_direction_stays_below_360() {
        let m = Mesh::new(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 1e-17], [0.0, 1.0, -1.0]],
            vec![[0, 1, 2]],
        )
        .unwrap();
        let (d, dir) = dip(&m)[0];
        assert!((d - 45.0).abs() < 1e-9);
        assert!((0.0..360.0).contains(&dir), "dip direction {dir}");
    }

    #[test]
    fn box_measures_match_analytic_values() {
        let m = fixtures::box_mesh([1.0, 2.0, 3.0], [3.0, 5.0, 7.0]);
        assert!(close(volume(&m).unwrap(), 24.0));
        assert!(close(area(&m), 2.0 * (6.0 + 8.0 + 12.0)));
        let c = centroid(&m).unwrap();
        for (a, b) in c.iter().zip([2.0, 3.5, 5.0]) {
            assert!(close(*a, b));
        }
        let c = area_centroid(&m).unwrap();
        for (a, b) in c.iter().zip([2.0, 3.5, 5.0]) {
            assert!(close(*a, b));
        }
    }

    #[test]
    fn divergence_volume_matches_analytic_solids() {
        assert!(close(volume(&fixtures::tetrahedron()).unwrap(), 1.0 / 6.0));
        assert!(close(volume(&fixtures::cube()).unwrap(), 1.0));
        assert!(close(volume(&fixtures::nested_shells()).unwrap(), 26.0));
        let c = centroid(&fixtures::tetrahedron()).unwrap();
        for x in c {
            assert!(close(x, 0.25));
        }
        let s = fixtures::sphere(4);
        let v = volume(&s).unwrap();
        let exact = 4.0 / 3.0 * std::f64::consts::PI;
        assert!(v < exact && (exact - v) / exact < 5e-3, "{v}");
        let t = fixtures::torus(96, 48);
        let v = volume(&t).unwrap();
        let exact = 2.0 * std::f64::consts::PI.powi(2) * 2.0 * 0.25;
        assert!(((exact - v) / exact).abs() < 5e-3, "{v}");
    }

    #[test]
    fn volume_refuses_open_or_inconsistent_meshes() {
        for m in [fixtures::open_cube(), broken::flipped_face()] {
            assert!(matches!(volume(&m), Err(Error::InvalidMesh { .. })));
            assert!(matches!(centroid(&m), Err(Error::InvalidMesh { .. })));
        }
        assert!(close(volume(&broken::inward_cube()).unwrap(), -1.0));
    }

    #[test]
    fn normals_dip_and_centroids() {
        let m = fixtures::cube();
        for (n, d) in face_normals(&m).iter().zip(dip(&m)) {
            assert!((norm(*n) - 1.0).abs() < 1e-12);
            if n[2].abs() == 1.0 {
                assert_eq!(d.0, 0.0);
                assert!(d.1.is_nan());
            } else {
                assert!(close(d.0, 90.0));
            }
        }
        let tilted = Mesh::new(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, -1.0]],
            vec![[0, 1, 2]],
        )
        .unwrap();
        let (d, az) = dip(&tilted)[0];
        assert!(close(d, 45.0) && close(az, 0.0), "{d} {az}");
        let east = Mesh::new(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, -1.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        )
        .unwrap();
        let (d, az) = dip(&east)[0];
        assert!(close(d, 45.0) && close(az, 90.0), "{d} {az}");
        for n in vertex_normals(&m) {
            for x in n {
                assert!(close(x.abs(), 1.0 / 3f64.sqrt()));
            }
        }
        assert_eq!(vertex_normals(&broken::unused_vertex())[8], [0.0; 3]);
        assert_eq!(face_normals(&broken::degenerate_face())[12], [0.0; 3]);
        let fc = face_centroids(&fixtures::tetrahedron());
        assert_eq!(fc.len(), 4);
        assert!(area_centroid(&Mesh::empty()).is_none());
    }

    fn quantized(m: &Mesh) -> Mesh {
        let q = |x: f64| (x * 1048576.0).round() / 1048576.0;
        let v = m.vertices().iter().map(|p| p.map(q)).collect();
        m.clone().with_vertices(v).unwrap()
    }

    #[test]
    fn volume_survives_utm_sized_coordinates() {
        let s = quantized(&fixtures::sphere(3));
        let far = fixtures::translated(&s, [524288.0, 7340032.0, 1024.0]);
        let (a, b) = (volume(&s).unwrap(), volume(&far).unwrap());
        assert!(close(a, b), "{a} {b}");
        let (ca, cb) = (centroid(&s).unwrap(), centroid(&far).unwrap());
        assert!((cb[1] - 7340032.0 - ca[1]).abs() < 1e-9);
    }

    #[test]
    fn sums_do_not_depend_on_thread_count() {
        let s = fixtures::translated(&fixtures::sphere(5), [500000.25, 7500000.5, 800.0]);
        let one = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .unwrap()
            .install(|| (volume(&s).unwrap(), area(&s)));
        let many = rayon::ThreadPoolBuilder::new()
            .num_threads(7)
            .build()
            .unwrap()
            .install(|| (volume(&s).unwrap(), area(&s)));
        assert_eq!(one.0.to_bits(), many.0.to_bits());
        assert_eq!(one.1.to_bits(), many.1.to_bits());
    }

    proptest! {
        #[test]
        fn box_volume_area_centroid_are_analytic(
            min in proptest::array::uniform3(-1000i32..1000),
            size in proptest::array::uniform3(1i32..500),
        ) {
            let lo = min.map(|x| x as f64 / 8.0);
            let d = size.map(|x| x as f64 / 16.0);
            let m = fixtures::box_mesh(lo, [lo[0] + d[0], lo[1] + d[1], lo[2] + d[2]]);
            prop_assert!(close(volume(&m).unwrap(), d[0] * d[1] * d[2]));
            prop_assert!(close(area(&m), 2.0 * (d[0] * d[1] + d[1] * d[2] + d[0] * d[2])));
            let c = centroid(&m).unwrap();
            for k in 0..3 {
                prop_assert!(close(c[k], lo[k] + d[k] / 2.0));
            }
        }

        #[test]
        fn volume_is_translation_invariant_and_flips_sign(
            d in proptest::array::uniform3(-(1i64 << 31)..(1i64 << 31)),
            sub in 0u32..3,
        ) {
            let s = quantized(&fixtures::sphere(sub));
            let shift = d.map(|x| x as f64 / 256.0);
            let far = fixtures::translated(&s, shift);
            let (a, b) = (volume(&s).unwrap(), volume(&far).unwrap());
            prop_assert!(close(a, b), "{} {}", a, b);
            prop_assert!(close(volume(&fixtures::flipped(&far)).unwrap(), -b));
        }
    }
}
