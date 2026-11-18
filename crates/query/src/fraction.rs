#[derive(Clone, Debug, PartialEq)]
pub struct Blocks {
    pub centroids: Vec<[f64; 3]>,
    /// Positive and finite.
    pub size: [f64; 3],
    /// Columns are the block axes in world coordinates; orthonormal with determinant +1.
    pub rotation: Option<[[f64; 3]; 3]>,
}
