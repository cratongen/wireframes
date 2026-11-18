use drusa_core::{Mesh, Polylines, Result};

#[derive(Clone, Debug)]
pub struct LoftOptions {
    pub cap_ends: bool,
}

impl Default for LoftOptions {
    fn default() -> LoftOptions {
        LoftOptions { cap_ends: true }
    }
}

pub fn loft(sections: &Polylines, opts: &LoftOptions) -> Result<Mesh> {
    todo!()
}
