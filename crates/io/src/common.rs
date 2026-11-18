#![allow(unused_variables)]

use drusa_core::Mesh;

use crate::Format;

/// `mesh` in the canonical form a write then read through `format` returns
/// unchanged (spec section 4.3). Test helper.
#[doc(hidden)]
pub fn canonicalize(mesh: &Mesh, format: Format) -> Mesh {
    todo!()
}
