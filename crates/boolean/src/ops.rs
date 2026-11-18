use drusa_core::{Mesh, Result};
use drusa_ops::RepairReport;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Op {
    Union,
    Intersection,
    Difference,
    SymmetricDifference,
}

/// `repair`: repair each input with `RepairOptions::default()` first.
#[derive(Clone, Debug, Default)]
pub struct BooleanOptions {
    pub repair: bool,
}

#[derive(Clone, Debug)]
pub struct BooleanResult {
    pub mesh: Mesh,
    /// Per output face: (0 for a, 1 for b, input face index after repair).
    pub face_source: Vec<(u8, u32)>,
    pub repairs: [Option<RepairReport>; 2],
}

pub fn boolean(a: &Mesh, b: &Mesh, op: Op, opts: &BooleanOptions) -> Result<BooleanResult> {
    todo!()
}

pub fn union(a: &Mesh, b: &Mesh, opts: &BooleanOptions) -> Result<BooleanResult> {
    boolean(a, b, Op::Union, opts)
}

pub fn intersection(a: &Mesh, b: &Mesh, opts: &BooleanOptions) -> Result<BooleanResult> {
    boolean(a, b, Op::Intersection, opts)
}

pub fn difference(a: &Mesh, b: &Mesh, opts: &BooleanOptions) -> Result<BooleanResult> {
    boolean(a, b, Op::Difference, opts)
}

pub fn symmetric_difference(a: &Mesh, b: &Mesh, opts: &BooleanOptions) -> Result<BooleanResult> {
    boolean(a, b, Op::SymmetricDifference, opts)
}
