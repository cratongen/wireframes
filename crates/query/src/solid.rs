use rayon::prelude::*;
use wireframes_core::exact::Expansion;
use wireframes_core::intersect::point_in_triangle;
use wireframes_core::predicates::{Sign, orient2d};
use wireframes_core::topology::ray_crosses;
use wireframes_core::{Aabb, Error, Mesh, Result};

use crate::ray::{Line, Param, V, dot, normal};
use crate::{Blocks, Interval, SurfaceIndex, Trace};

/// `ray_crosses` for the exact point `x / w`, `w > 0`, with the same plan
/// perturbation. The point must not lie on the closed face.
fn ray_crosses_h(c: [[f64; 3]; 3], x: &V, w: &Expansion) -> bool {
    let plan = c.map(|v| [v[0], v[1]]);
    let s = Sign::of(orient2d(plan[0], plan[1], plan[2]));
    if s == Sign::Zero {
        return false;
    }
    let side = |a: [f64; 2], b: [f64; 2]| {
        let o = Expansion::diff_f64(b[0], a[0])
            .product(&x[1].diff(&w.scale(a[1])))
            .diff(&Expansion::diff_f64(b[1], a[1]).product(&x[0].diff(&w.scale(a[0]))));
        match o.sign() {
            Sign::Zero => match Sign::of(a[1] - b[1]) {
                Sign::Zero => Sign::of(b[0] - a[0]),
                t => t,
            },
            t => t,
        }
    };
    if (0..3).any(|i| side(plan[i], plan[(i + 1) % 3]) != s) {
        return false;
    }
    let below = [0, 1, 2].map(|k| w.scale(c[0][k]).diff(&x[k]));
    dot(&normal(&c), &below).sign() == s
}

/// Sign agrees with signed distance.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(i8)]
pub enum Location {
    Inside = -1,
    Boundary = 0,
    Outside = 1,
}

pub struct SolidIndex {
    pub(crate) surface: SurfaceIndex,
}

impl SolidIndex {
    /// `InvalidMesh { context: "solid" }` unless `validate(0.0)` says it is a solid.
    pub fn new(solid: &Mesh) -> Result<SolidIndex> {
        let report = solid.validate(0.0)?;
        if !report.summary.is_solid() {
            return Err(Error::InvalidMesh {
                context: "solid".into(),
                report: Box::new(report),
            });
        }
        Ok(SolidIndex {
            surface: SurfaceIndex::new(solid),
        })
    }

    pub fn mesh(&self) -> &Mesh {
        &self.surface.mesh
    }

    /// Exact: `Boundary` on any closed face, else parity of the `+z` ray with
    /// plan ties broken by symbolic perturbation. `p` must be finite.
    pub fn locate(&self, p: [f64; 3]) -> Location {
        let mesh = self.mesh();
        let (mut boundary, mut odd) = (false, false);
        self.surface
            .bvh
            .ray(p, [0.0, 0.0, 1.0], f64::INFINITY, |f| {
                if boundary {
                    return;
                }
                let c = mesh.corners(f);
                if point_in_triangle(p, c) {
                    boundary = true;
                } else if ray_crosses(c, p) {
                    odd = !odd;
                }
            });
        match (boundary, odd) {
            (true, _) => Location::Boundary,
            (false, true) => Location::Inside,
            (false, false) => Location::Outside,
        }
    }

    /// Inside or outside for the exact point of `line` at `t`, which must not
    /// lie on the surface.
    pub(crate) fn locate_on_line(&self, line: &Line, t: &Param) -> Location {
        let mesh = self.mesh();
        let (x, w) = line.point_h(t);
        let est = line.at(t.value());
        let o = line.origin();
        let m = (0..3).fold(0.0f64, |m, k| m.max(o[k].abs() + est[k].abs()));
        let pad = 64.0 * f64::EPSILON * m + f64::MIN_POSITIVE;
        let column = Aabb {
            min: [est[0] - pad, est[1] - pad, est[2] - pad],
            max: [est[0] + pad, est[1] + pad, f64::INFINITY],
        };
        let mut odd = false;
        self.surface.bvh.query(&column, |f| {
            if ray_crosses_h(mesh.corners(f), &x, &w) {
                odd = !odd;
            }
        });
        if odd {
            Location::Inside
        } else {
            Location::Outside
        }
    }

    pub fn locate_all(&self, points: &[[f64; 3]]) -> Vec<Location> {
        points.par_iter().map(|&p| self.locate(p)).collect()
    }

    /// `boundary_inside` treats the solid as a closed set, otherwise as open.
    pub fn contains_all(&self, points: &[[f64; 3]], boundary_inside: bool) -> Vec<bool> {
        points
            .par_iter()
            .map(|&p| match self.locate(p) {
                Location::Inside => true,
                Location::Boundary => boundary_inside,
                Location::Outside => false,
            })
            .collect()
    }

    /// Negative inside, 0 on the boundary.
    pub fn signed_distance(&self, points: &[[f64; 3]]) -> Vec<f64> {
        points
            .par_iter()
            .map(|&p| {
                let d = self.surface.closest_one(p).distance;
                // An off-surface point whose distance rounds to 0 keeps its sign.
                match self.locate(p) {
                    Location::Boundary => 0.0,
                    Location::Inside => -d.max(f64::MIN_POSITIVE),
                    Location::Outside => d.max(f64::MIN_POSITIVE),
                }
            })
            .collect()
    }

    pub fn block_fraction(&self, blocks: &Blocks) -> Result<Vec<f64>> {
        crate::fraction::fractions(&self.surface, blocks)
    }

    pub fn trace_intervals(&self, trace: &Trace) -> Result<Vec<Interval>> {
        crate::trace::intervals(self, trace)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;
    use wireframes_core::fixtures;
    use wireframes_core::topology::locate_point;

    #[test]
    fn locate_on_faces_edges_and_vertices() {
        let s = SolidIndex::new(&fixtures::cube()).unwrap();
        use Location::*;
        for (p, want) in [
            ([0.5, 0.5, 0.5], Inside),
            ([0.5, 0.5, 1.0], Boundary),
            ([1.0, 0.5, 0.5], Boundary),
            ([1.0, 1.0, 0.5], Boundary),
            ([1.0, 1.0, 1.0], Boundary),
            ([0.0, 0.0, 0.0], Boundary),
            ([0.5, 0.5, 0.0], Boundary),
            ([0.5, 0.5, -0.5], Outside),
            ([1.0, 1.0, -0.5], Outside),
            ([0.0, 0.5, -0.5], Outside),
            ([0.5, 0.5, 1.5], Outside),
            ([1.5, 0.5, 0.5], Outside),
            ([0.5, 0.5, 0.5f64.next_up()], Inside),
            ([1f64.next_up(), 0.5, 0.5], Outside),
        ] {
            assert_eq!(s.locate(p), want, "{p:?}");
        }
        let c = s.contains_all(&[[1.0, 0.5, 0.5], [0.5; 3]], true);
        assert_eq!(c, [true, true]);
        let c = s.contains_all(&[[1.0, 0.5, 0.5], [0.5; 3]], false);
        assert_eq!(c, [false, true]);
    }

    #[test]
    fn signed_distance_of_a_cube() {
        let s = SolidIndex::new(&fixtures::cube()).unwrap();
        let d = s.signed_distance(&[[0.5, 0.5, 0.25], [0.5, 0.5, 1.0], [0.5, 0.5, 3.0]]);
        assert_eq!(d, [-0.25, 0.0, 2.0]);
    }

    #[test]
    fn invalid_solids_are_refused() {
        for m in [fixtures::open_cube(), fixtures::broken::inward_cube()] {
            match SolidIndex::new(&m) {
                Err(Error::InvalidMesh { context, .. }) => assert_eq!(context, "solid"),
                _ => panic!("expected InvalidMesh"),
            }
        }
    }

    #[test]
    fn same_result_on_one_and_four_threads() {
        let s = SolidIndex::new(&fixtures::sphere(3)).unwrap();
        let pts = fixtures::random_points(7, 5000, [-1.2; 3], [1.2; 3]);
        let run = |n| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(n)
                .build()
                .unwrap()
                .install(|| (s.locate_all(&pts), s.signed_distance(&pts)))
        };
        let (a, b) = (run(1), run(4));
        assert_eq!(a.0, b.0);
        assert!(
            a.1.iter()
                .zip(&b.1)
                .all(|(x, y)| x.to_bits() == y.to_bits())
        );
    }

    fn grid_points(seed: u64, n: usize) -> Vec<[f64; 3]> {
        let mut rng = fixtures::Rng::new(seed);
        (0..n)
            .map(|_| [0, 1, 2].map(|_| rng.below(9) as f64 * 0.25 - 0.5))
            .collect()
    }

    proptest! {
        #[test]
        fn locate_matches_brute_force_and_sign_of_distance(seed in 0u64..500) {
            for mesh in [fixtures::torus(10, 6), fixtures::cube(), fixtures::nested_shells()] {
                let s = SolidIndex::new(&mesh).unwrap();
                let mut pts = fixtures::random_points(seed, 30, [-1.5; 3], [1.5; 3]);
                pts.extend(grid_points(seed, 30));
                let loc = s.locate_all(&pts);
                let sd = s.signed_distance(&pts);
                for ((p, l), d) in pts.iter().zip(&loc).zip(&sd) {
                    prop_assert_eq!(*l as i8, locate_point(mesh.vertices(), mesh.triangles(), *p));
                    prop_assert_eq!(*l as i8 as f64, d.signum() * (*d != 0.0) as u8 as f64);
                }
            }
        }

        #[test]
        fn vertices_locate_on_the_boundary(sub in 0u32..3) {
            let m = fixtures::sphere(sub);
            let s = SolidIndex::new(&m).unwrap();
            prop_assert!(s.locate_all(m.vertices()).iter().all(|&l| l == Location::Boundary));
        }
    }
}
