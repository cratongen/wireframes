#![allow(unused_variables)]

use drusa_core::Result;

use crate::{Data, WriteOptions};

/// Reads a `_pt.csv` and `_tr.csv` pair.
pub fn read(pt: &[u8], tr: &[u8]) -> Result<Data> {
    todo!()
}

/// Writes the `(pt, tr)` pair.
pub fn write(data: &Data, opts: &WriteOptions) -> Result<(Vec<u8>, Vec<u8>)> {
    todo!()
}
