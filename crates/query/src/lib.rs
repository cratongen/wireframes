//! Point-in-solid, distances, ray and drillhole intersections and block fractions.

pub mod distance;
pub mod fraction;
pub mod ray;
pub mod solid;
pub mod trace;

use wireframes_core::{Mesh, Result};

pub use distance::{Closest, SurfaceIndex};
pub use fraction::Blocks;
pub use ray::{Hit, HitKind};
pub use solid::{Location, SolidIndex};
pub use trace::{Interval, Trace};

use distance::check_finite;

pub fn locate(solid: &Mesh, points: &[[f64; 3]]) -> Result<Vec<Location>> {
    check_finite("points", points)?;
    Ok(SolidIndex::new(solid)?.locate_all(points))
}

pub fn signed_distance(solid: &Mesh, points: &[[f64; 3]]) -> Result<Vec<f64>> {
    check_finite("points", points)?;
    Ok(SolidIndex::new(solid)?.signed_distance(points))
}

pub fn block_fraction(solid: &Mesh, blocks: &Blocks) -> Result<Vec<f64>> {
    SolidIndex::new(solid)?.block_fraction(blocks)
}

pub fn trace_intervals(solid: &Mesh, trace: &Trace) -> Result<Vec<Interval>> {
    SolidIndex::new(solid)?.trace_intervals(trace)
}

/// Points must be finite.
pub fn closest(mesh: &Mesh, points: &[[f64; 3]]) -> Vec<Closest> {
    SurfaceIndex::new(mesh).closest(points)
}

pub fn ray_hits(mesh: &Mesh, origins: &[[f64; 3]], dirs: &[[f64; 3]]) -> Result<Vec<Vec<Hit>>> {
    SurfaceIndex::new(mesh).ray_hits(origins, dirs)
}
