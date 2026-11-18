use drusa_core::{Mesh, Polylines, Result};

#[derive(Clone, Debug, Default)]
pub struct DtmOptions<'a> {
    pub breaklines: Option<&'a Polylines>,
    /// Closed rings, even-odd; triangles outside are dropped.
    pub boundary: Option<&'a Polylines>,
    /// Plan duplicate merge distance; 0 merges exact duplicates only.
    pub min_spacing: f64,
    /// Drops boundary triangles with a plan edge longer than this, repeatedly.
    pub max_edge: Option<f64>,
}

pub fn dtm(points: &[[f64; 3]], opts: &DtmOptions) -> Result<Mesh> {
    todo!()
}
