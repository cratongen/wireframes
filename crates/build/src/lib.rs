//! Surface and solid construction and mesh editing.

pub mod close;
pub mod contour;
pub mod decimate;
pub mod dtm;
pub mod extrude;
mod geom;
pub mod loft;
pub mod offset;
pub mod polygon;
pub mod remesh;
pub mod smooth;

pub use close::close_to_elevation;
pub use contour::contours;
pub use decimate::{DecimateOptions, decimate};
pub use dtm::{DtmOptions, dtm};
pub use extrude::extrude;
pub use loft::{LoftOptions, loft};
pub use offset::offset;
pub use polygon::triangulate_polygons;
pub use remesh::remesh;
pub use smooth::smooth_taubin;
