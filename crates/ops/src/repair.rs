use drusa_core::{Mesh, MeshProblemKind, MeshSummary, Result};

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum HoleFill {
    Fan,
    MinArea,
}

#[derive(Clone, Debug, PartialEq)]
pub struct RepairOptions {
    /// None: no weld; `Some(0.0)`: exactly equal coordinates.
    pub weld: Option<f64>,
    pub remove_degenerate: bool,
    pub remove_duplicate_faces: bool,
    pub remove_unused_vertices: bool,
    /// Consistent winding, outward by nesting parity.
    pub orient: bool,
    pub fill_holes: Option<HoleFill>,
}

impl Default for RepairOptions {
    fn default() -> RepairOptions {
        RepairOptions {
            weld: Some(0.0),
            remove_degenerate: true,
            remove_duplicate_faces: true,
            remove_unused_vertices: true,
            orient: true,
            fill_holes: None,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum RepairAction {
    VertexMerged {
        vertex: u32,
        into: u32,
    },
    VertexRemoved {
        vertex: u32,
    },
    /// `reason` is `DegenerateFace` or `DuplicateFace`.
    FaceRemoved {
        face: u32,
        reason: MeshProblemKind,
    },
    FaceFlipped {
        face: u32,
    },
    /// `face` indexes the output.
    FaceAdded {
        face: u32,
        hole: u32,
    },
}

#[derive(Clone, Debug, PartialEq)]
pub struct RepairReport {
    /// In step order: weld, degenerate, duplicate, unused, orient, fill.
    pub actions: Vec<RepairAction>,
    pub vertex_map: Vec<Option<u32>>,
    pub face_map: Vec<Option<u32>>,
    pub before: MeshSummary,
    pub after: MeshSummary,
    /// Loops not filled (pinched at a vertex used by two boundary edges).
    pub holes_left: usize,
}

pub fn repair(mesh: &Mesh, opts: &RepairOptions) -> Result<(Mesh, RepairReport)> {
    todo!()
}
