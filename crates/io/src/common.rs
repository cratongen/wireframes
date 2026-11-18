#![allow(unused_variables)]

use std::ops::Range;

use drusa_core::{ColumnData, Error, Mesh, Result, Table};

use crate::Format;

/// `mesh` in the canonical form a write then read through `format` returns
/// unchanged (spec section 4.3). Test helper.
#[doc(hidden)]
pub fn canonicalize(mesh: &Mesh, format: Format) -> Mesh {
    todo!()
}

/// Face ranges that came from one polygon: consecutive faces with the same valid
/// `polygon` id that continue the fan of the first one.
pub(crate) fn polygon_groups(mesh: &Mesh) -> Result<Vec<Range<usize>>> {
    let tris = mesh.triangles();
    let ids = match mesh.face_attributes().column("polygon") {
        None => None,
        Some(c) => match &c.data {
            ColumnData::Int64(v) => Some((v, c)),
            _ => {
                return Err(Error::invalid(format!(
                    "face column \"polygon\" must be int64, found {}",
                    c.type_name()
                )));
            }
        },
    };
    let mut out = Vec::with_capacity(tris.len());
    let mut i = 0;
    while i < tris.len() {
        let mut j = i + 1;
        if let Some((v, c)) = ids
            && c.is_valid(i)
        {
            while j < tris.len()
                && c.is_valid(j)
                && v[j] == v[i]
                && tris[j][0] == tris[i][0]
                && tris[j][1] == tris[j - 1][2]
            {
                j += 1;
            }
        }
        out.push(i..j);
        i = j;
    }
    Ok(out)
}

pub(crate) fn polygon_corners(tris: &[[u32; 3]], g: &Range<usize>, out: &mut Vec<u32>) {
    out.clear();
    out.extend_from_slice(&tris[g.start]);
    out.extend(tris[g.start + 1..g.end].iter().map(|t| t[2]));
}

/// Faces regrouped into the polygons a reader fanned, with the face attributes
/// of each polygon's first triangle; the `polygon` column is consumed.
pub(crate) fn polygons(mesh: &Mesh) -> Result<(Vec<Vec<u32>>, Table)> {
    let groups = polygon_groups(mesh)?;
    let mut polys = Vec::with_capacity(groups.len());
    for g in &groups {
        let mut p = vec![];
        polygon_corners(mesh.triangles(), g, &mut p);
        polys.push(p);
    }
    let first: Vec<u32> = groups.iter().map(|g| g.start as u32).collect();
    let table = mesh.face_attributes().clone().without_column("polygon");
    Ok((polys, table.take(&first)))
}
