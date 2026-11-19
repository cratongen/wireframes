//! Exact mesh arrangements, booleans, surface clipping, plane slicing and cut and fill.

mod arrangement;
mod classify;
pub mod clip;
pub mod cutfill;
mod hp;
pub mod ops;
pub mod slice;
mod triangulate;

pub use clip::{ClipResult, clip_by_surface};
pub use cutfill::{CutFill, cut_fill};
pub use ops::{
    BooleanOptions, BooleanResult, Op, boolean, difference, intersection, symmetric_difference,
    union,
};
pub use slice::{slice_by_plane, split_by_plane};

/// Sub-triangle count of the arrangement of `a` and `b`, for benchmarks.
/// Both must already pass `validate(0.0)` as solids; nothing here checks it.
#[doc(hidden)]
pub fn arrangement_size(
    a: &wireframes_core::Mesh,
    b: &wireframes_core::Mesh,
) -> wireframes_core::Result<usize> {
    arrangement::build(a, b).map(|arr| arr.subs.len())
}
