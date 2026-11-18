use drusa_core::predicates::Sign;
use drusa_core::validate::degenerate_mask;
use drusa_core::{Aabb, Bvh, Error, Mesh, Result};
use rayon::prelude::*;

use crate::Hit;
use crate::ray::{Line, RawHit, face_hit, hit_point, sort_and_collapse};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Closest {
    pub point: [f64; 3],
    pub face: u32,
    pub distance: f64,
}

pub struct SurfaceIndex {
    pub(crate) mesh: Mesh,
    pub(crate) bvh: Bvh,
    degenerate: Vec<bool>,
    closed: bool,
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn lerp(a: [f64; 3], d: [f64; 3], t: f64) -> [f64; 3] {
    [a[0] + t * d[0], a[1] + t * d[1], a[2] + t * d[2]]
}

/// Closest point of segment `ab` to the origin.
fn closest_on_segment(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    let ab = sub(b, a);
    let l2 = dot(ab, ab);
    if l2 == 0.0 {
        return a;
    }
    lerp(a, ab, (-dot(a, ab) / l2).clamp(0.0, 1.0))
}

/// Closest point of triangle `abc` to the origin (Ericson, Real-Time Collision
/// Detection, 5.1.5). Degenerate triangles use their three edges.
fn closest_on_triangle(t: [[f64; 3]; 3], degenerate: bool) -> [f64; 3] {
    let [a, b, c] = t;
    if degenerate {
        return [
            closest_on_segment(a, b),
            closest_on_segment(b, c),
            closest_on_segment(c, a),
        ]
        .into_iter()
        .min_by(|x, y| dot(*x, *x).total_cmp(&dot(*y, *y)))
        .expect("three edges");
    }
    let (ab, ac) = (sub(b, a), sub(c, a));
    let ap = [-a[0], -a[1], -a[2]];
    let (d1, d2) = (dot(ab, ap), dot(ac, ap));
    if d1 <= 0.0 && d2 <= 0.0 {
        return a;
    }
    let bp = [-b[0], -b[1], -b[2]];
    let (d3, d4) = (dot(ab, bp), dot(ac, bp));
    if d3 >= 0.0 && d4 <= d3 {
        return b;
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        return lerp(a, ab, d1 / (d1 - d3));
    }
    let cp = [-c[0], -c[1], -c[2]];
    let (d5, d6) = (dot(ab, cp), dot(ac, cp));
    if d6 >= 0.0 && d5 <= d6 {
        return c;
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        return lerp(a, ac, d2 / (d2 - d6));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && d4 - d3 >= 0.0 && d5 - d6 >= 0.0 {
        return lerp(b, sub(c, b), (d4 - d3) / ((d4 - d3) + (d5 - d6)));
    }
    let denom = 1.0 / (va + vb + vc);
    let (v, w) = (vb * denom, vc * denom);
    [0, 1, 2].map(|k| a[k] + ab[k] * v + ac[k] * w)
}

pub(crate) fn check_finite(what: &str, points: &[[f64; 3]]) -> Result<()> {
    match points.iter().position(|p| p.iter().any(|x| !x.is_finite())) {
        None => Ok(()),
        Some(i) => Err(Error::invalid(format!(
            "{what} must be finite: {what} {i} is {:?}",
            points[i]
        ))),
    }
}

fn check_pairs(names: [&str; 2], a: &[[f64; 3]], b: &[[f64; 3]]) -> Result<()> {
    if a.len() != b.len() {
        return Err(Error::invalid(format!(
            "{} and {} must have the same length: {} vs {}",
            names[0],
            names[1],
            a.len(),
            b.len()
        )));
    }
    check_finite(names[0], a)?;
    check_finite(names[1], b)
}

impl SurfaceIndex {
    /// Any mesh. Degenerate faces count for `closest` but never for ray hits.
    pub fn new(mesh: &Mesh) -> SurfaceIndex {
        SurfaceIndex {
            mesh: mesh.clone(),
            bvh: Bvh::new(&mesh.face_aabbs()),
            degenerate: degenerate_mask(mesh),
            closed: mesh.is_closed(),
        }
    }

    /// Distances use corners relative to the query point. A mesh without faces
    /// gives a NaN point, face `u32::MAX` and infinite distance.
    pub(crate) fn closest_one(&self, p: [f64; 3]) -> Closest {
        let rel = |f: u32| {
            let c = self.mesh.corners(f).map(|v| sub(v, p));
            closest_on_triangle(c, self.degenerate[f as usize])
        };
        match self.bvh.nearest(p, |f| {
            let q = rel(f);
            dot(q, q)
        }) {
            Some((face, d2)) => {
                let q = rel(face);
                Closest {
                    point: [p[0] + q[0], p[1] + q[1], p[2] + q[2]],
                    face,
                    distance: d2.sqrt(),
                }
            }
            None => Closest {
                point: [f64::NAN; 3],
                face: u32::MAX,
                distance: f64::INFINITY,
            },
        }
    }

    /// Points must be finite.
    pub fn closest(&self, points: &[[f64; 3]]) -> Vec<Closest> {
        points.par_iter().map(|&p| self.closest_one(p)).collect()
    }

    pub fn distance(&self, points: &[[f64; 3]]) -> Vec<f64> {
        points
            .par_iter()
            .map(|&p| self.closest_one(p).distance)
            .collect()
    }

    pub(crate) fn raw_hits(&self, line: &Line) -> Vec<RawHit> {
        let mut faces = vec![];
        match *line {
            Line::Ray(o, d) => self.bvh.ray(o, d, f64::INFINITY, |f| faces.push(f)),
            Line::Segment(p, q) => self
                .bvh
                .query(&Aabb::from_points(&[p, q]), |f| faces.push(f)),
        }
        faces.sort_unstable();
        faces
            .into_iter()
            .filter(|&f| !self.degenerate[f as usize])
            .filter_map(|f| face_hit(&self.mesh.corners(f), f, line))
            .collect()
    }

    fn hits(&self, line: Line) -> Vec<Hit> {
        let mut raw = self.raw_hits(&line);
        sort_and_collapse(&mut raw);
        raw.iter()
            .map(|h| Hit {
                t: h.lo.value(),
                face: h.face,
                point: hit_point(&self.mesh.corners(h.face), &line, h),
                kind: h.kind,
                entering: match (self.closed, h.crossing) {
                    (false, _) | (_, Sign::Zero) => None,
                    (true, s) => Some(s == Sign::Negative),
                },
            })
            .collect()
    }

    /// Hits with `t >= 0`, sorted by `t` then face, one per distinct exact point.
    pub fn ray_hits(&self, origins: &[[f64; 3]], dirs: &[[f64; 3]]) -> Result<Vec<Vec<Hit>>> {
        check_pairs(["origins", "dirs"], origins, dirs)?;
        if let Some(i) = dirs.iter().position(|d| *d == [0.0; 3]) {
            return Err(Error::invalid(format!(
                "dirs must be nonzero: dir {i} is zero"
            )));
        }
        Ok(origins
            .par_iter()
            .zip(dirs.par_iter())
            .map(|(&o, &d)| self.hits(Line::Ray(o, d)))
            .collect())
    }

    /// Hits with `0 <= t <= 1` along `a + t * (b - a)`.
    pub fn segment_hits(&self, a: &[[f64; 3]], b: &[[f64; 3]]) -> Result<Vec<Vec<Hit>>> {
        check_pairs(["a", "b"], a, b)?;
        if let Some(i) = a.iter().zip(b).position(|(p, q)| p == q) {
            return Err(Error::invalid(format!(
                "segments must have nonzero length: segment {i} has equal ends"
            )));
        }
        Ok(a.par_iter()
            .zip(b.par_iter())
            .map(|(&p, &q)| self.hits(Line::Segment(p, q)))
            .collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::HitKind;
    use drusa_core::fixtures;
    use proptest::prelude::*;

    fn brute(mesh: &Mesh, p: [f64; 3]) -> (u32, f64) {
        let deg = degenerate_mask(mesh);
        (0..mesh.n_faces() as u32)
            .map(|f| {
                let q = closest_on_triangle(mesh.corners(f).map(|v| sub(v, p)), deg[f as usize]);
                (f, dot(q, q))
            })
            .fold(
                (u32::MAX, f64::INFINITY),
                |b, x| if x.1 < b.1 { x } else { b },
            )
    }

    #[test]
    fn closest_on_a_cube() {
        let s = SurfaceIndex::new(&fixtures::cube());
        let c = s.closest(&[[0.5, 0.5, 2.0], [2.0, 2.0, 2.0], [0.5, 0.5, 0.25]]);
        assert_eq!(c[0].point, [0.5, 0.5, 1.0]);
        assert_eq!(c[0].distance, 1.0);
        assert_eq!(c[1].point, [1.0, 1.0, 1.0]);
        assert_eq!(c[1].distance, 3f64.sqrt());
        assert_eq!(c[2].distance, 0.25);
        assert_eq!(c[2].point, [0.5, 0.5, 0.0]);
        let empty = SurfaceIndex::new(&Mesh::empty());
        assert_eq!(empty.closest(&[[0.0; 3]])[0].face, u32::MAX);
    }

    #[test]
    fn rays_through_a_cube_collapse_shared_edges() {
        let s = SurfaceIndex::new(&fixtures::cube());
        let h = s
            .ray_hits(
                &[[0.25, 0.75, -1.0], [0.5, 0.5, -1.0]],
                &[[0.0, 0.0, 1.0]; 2],
            )
            .unwrap();
        assert_eq!(h[0].len(), 2);
        assert_eq!(h[0][0].entering, Some(true));
        assert_eq!(h[0][1].entering, Some(false));
        assert_eq!((h[0][0].t, h[0][1].t), (1.0, 2.0));
        assert_eq!(
            h[1].len(),
            2,
            "diagonal edge on bottom and top: one hit each"
        );
        assert_eq!(h[1][0].kind, HitKind::Edge);
        assert_eq!(h[1][0].point, [0.5, 0.5, 0.0]);
        let open = SurfaceIndex::new(&fixtures::open_cube());
        let h = open
            .ray_hits(&[[0.25, 0.75, -1.0]], &[[0.0, 0.0, 1.0]])
            .unwrap();
        assert!(h[0].iter().all(|x| x.entering.is_none()));
        assert!(s.ray_hits(&[[0.0; 3]], &[[0.0; 3]]).is_err());
        assert!(s.ray_hits(&[[0.0; 3]], &[]).is_err());
        assert!(
            s.ray_hits(&[[f64::NAN, 0.0, 0.0]], &[[1.0, 0.0, 0.0]])
                .is_err()
        );
    }

    #[test]
    fn segments_and_coplanar_rays() {
        let s = SurfaceIndex::new(&fixtures::cube());
        let h = s
            .segment_hits(&[[0.25, 0.75, -1.0]], &[[0.25, 0.75, 0.5]])
            .unwrap();
        assert_eq!(h[0].len(), 1);
        assert_eq!(h[0][0].point, [0.25, 0.75, 0.0]);
        assert!((h[0][0].t - 2.0 / 3.0).abs() < 1e-15);
        let h = s
            .ray_hits(&[[-1.0, 0.25, 0.0]], &[[1.0, 0.0, 0.0]])
            .unwrap();
        assert!(
            h[0].iter()
                .any(|x| x.kind == HitKind::Coplanar && x.entering.is_none())
        );
        assert_eq!(h[0][0].t, 1.0);
        assert!(s.segment_hits(&[[0.0; 3]], &[[0.0; 3]]).is_err());
    }

    proptest! {
        #[test]
        fn closest_matches_brute_force(seed in 0u64..1000) {
            let mesh = fixtures::torus(12, 8);
            let s = SurfaceIndex::new(&mesh);
            for p in fixtures::random_points(seed, 20, [-2.0; 3], [2.0; 3]) {
                let c = s.closest(&[p])[0];
                let (_, d2) = brute(&mesh, p);
                prop_assert_eq!(c.distance, d2.sqrt());
                for v in mesh.vertices() {
                    prop_assert!(c.distance <= dot(sub(*v, p), sub(*v, p)).sqrt());
                }
            }
        }
    }
}
