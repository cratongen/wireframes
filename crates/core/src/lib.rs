//! Core containers, exact predicates, spatial index, topology and validation.

mod aabb;
mod bvh;
mod error;
pub mod exact;
pub mod fixtures;
pub mod intersect;
mod mesh;
#[cfg(test)]
mod oracle;
mod points;
mod polylines;
pub mod predicates;
mod table;
pub mod tolerance;
pub mod topology;
pub mod units;
pub mod validate;

pub use aabb::Aabb;
pub use bvh::Bvh;
pub use error::{Error, Result};
pub use mesh::Mesh;
pub use points::PointSet;
pub use polylines::Polylines;
pub use table::{Column, ColumnData, Table};
pub use topology::{EdgeMap, HalfEdges};
pub use validate::{MeshProblem, MeshProblemKind, MeshReport, MeshSummary, validate};
