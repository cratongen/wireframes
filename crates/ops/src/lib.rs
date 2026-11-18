//! Measures, slab volumes, thickness, repair, components and transforms.

pub mod components;
pub mod measure;
pub mod repair;
pub mod slices;
pub mod thickness;
pub mod transform;

pub use components::{component_labels, components};
pub use measure::{
    area, area_centroid, centroid, dip, face_areas, face_centroids, face_normals, vertex_normals,
    volume,
};
pub use repair::{HoleFill, RepairAction, RepairOptions, RepairReport, repair};
pub use slices::{box_volume, volume_by_elevation};
pub use thickness::{ThicknessMode, surface_z, thickness, vertical_distance};
pub use transform::{flip, to_length_unit, transform, translate};
