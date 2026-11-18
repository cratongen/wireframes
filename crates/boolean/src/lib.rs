//! Exact mesh arrangements, booleans, surface clipping, plane slicing and cut and fill.

#![allow(unused_variables)]

mod arrangement;
mod classify;
pub mod clip;
pub mod cutfill;
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
