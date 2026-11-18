#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HitKind {
    Interior,
    Edge,
    Vertex,
    Coplanar,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Hit {
    pub t: f64,
    pub face: u32,
    pub point: [f64; 3],
    pub kind: HitKind,
    /// `Some` for closed meshes.
    pub entering: Option<bool>,
}
