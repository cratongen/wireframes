use std::collections::HashMap;
use std::collections::hash_map::Entry;
use std::ops::Range;

use drusa_core::{Column, ColumnData, Error, Mesh, Result, Table};

use crate::ply::{F64_EXACT_INT, regenerated_by_reader};
use crate::{Data, Format, WriteOptions};

/// `mesh` in the canonical form a write then read through `format` returns
/// unchanged (spec section 4.3): columns the format cannot store are dropped,
/// polygons are renumbered or regrouped as the reader produces them, and for STL
/// vertices are rounded to `f32` (out of range values clamp to `f32::MAX`) and
/// exact-welded in first-use order. STL keeps both `solid` (stored by ASCII) and
/// `stl_attribute` (stored by binary); each mode round-trips the mesh without the
/// other column. Formats outside OBJ, PLY, OFF and STL go through their writer
/// (binary, skipping unsupported columns) and reader. Test helper; panics on a
/// mesh the format cannot write.
#[doc(hidden)]
pub fn canonicalize(mesh: &Mesh, format: Format) -> Mesh {
    let mut m = mesh.clone().with_crs(None).with_length_unit(None);
    if m.face_attributes()
        .column("polygon")
        .is_some_and(|c| !matches!(c.data, ColumnData::Int64(_)))
    {
        let (vt, ft) = (
            m.vertex_attributes().clone(),
            m.face_attributes().clone().without_column("polygon"),
        );
        m = m.with_attributes(vt, ft).expect("same row counts");
    }
    let groups = polygon_groups(&m).expect("polygon is int64");
    let firsts: Vec<u32> = groups
        .iter()
        .flat_map(|g| std::iter::repeat_n(g.start as u32, g.len()))
        .collect();
    let per_polygon = m.face_attributes().take(&firsts);
    let renumbered = || {
        groups.iter().any(|g| g.len() > 1).then(|| {
            Column::i64(
                groups
                    .iter()
                    .enumerate()
                    .flat_map(|(k, g)| std::iter::repeat_n(k as i64, g.len()))
                    .collect(),
            )
        })
    };
    let n = m.n_faces();
    let result = match format {
        Format::Obj | Format::Off => {
            let mut ft = Table::new(n);
            if format == Format::Obj && n > 0 {
                for name in ["object", "group"] {
                    if let Some(c) = per_polygon.column(name)
                        && matches!(c.data, ColumnData::Utf8(_))
                    {
                        ft = ft.with_column(name, c.clone()).expect("same row counts");
                    }
                }
            }
            if let Some(c) = renumbered() {
                ft = ft.with_column("polygon", c).expect("same row counts");
            }
            let vt = Table::new(m.n_vertices());
            m.with_attributes(vt, ft)
        }
        Format::Ply => {
            let regenerated = regenerated_by_reader(&m, &groups);
            let vt = ply_table(m.vertex_attributes(), &["x", "y", "z"]);
            let mut ft =
                ply_table(m.face_attributes(), &["vertex_indices", "vertex_index"]).take(&firsts);
            if regenerated {
                ft = ft.without_column("polygon");
            }
            let structural = if regenerated || ft.column("polygon").is_none() {
                renumbered()
            } else {
                None
            };
            let ft = match structural {
                Some(c) => ft.without_column("polygon").with_column("polygon", c),
                None => Ok(ft),
            };
            m.with_attributes(vt, ft.expect("same row counts"))
        }
        Format::Stl => return stl_canonical(&m),
        _ => return round_trip(&m, format),
    };
    result.expect("same row counts")
}

fn ply_table(t: &Table, reserved: &[&str]) -> Table {
    let mut out = Table::new(t.rows());
    for (name, c) in t.names().iter().zip(t.columns()) {
        let bad = name.is_empty()
            || name
                .bytes()
                .any(|b| b.is_ascii_whitespace() || !b.is_ascii())
            || reserved.contains(&name.as_str());
        let valid: Vec<bool> = (0..c.len()).map(|i| c.is_valid(i)).collect();
        let col = match &c.data {
            _ if bad => None,
            ColumnData::Float64(v) => {
                let valid = v
                    .iter()
                    .zip(&valid)
                    .map(|(x, &ok)| ok && !x.is_nan())
                    .collect();
                Column::f64(v.clone()).with_valid(valid).ok()
            }
            ColumnData::Int64(v) => {
                if c.valid.is_none() && v.iter().all(|&x| i32::try_from(x).is_ok()) {
                    Some(Column::i64(v.clone()))
                } else if v.iter().all(|&x| x.unsigned_abs() <= F64_EXACT_INT as u64) {
                    let f = Column::f64(v.iter().map(|&x| x as f64).collect());
                    if name == "polygon" {
                        Some(c.clone())
                    } else {
                        f.with_valid(valid).ok()
                    }
                } else {
                    None
                }
            }
            ColumnData::Bool(v) => Some(Column::i64(
                v.iter()
                    .zip(&valid)
                    .map(|(&b, &ok)| i64::from(ok && b))
                    .collect(),
            )),
            ColumnData::Utf8(_) => None,
        };
        if let Some(col) = col {
            out = out.with_column(name, col).expect("same row count");
        }
    }
    out
}

fn stl_canonical(m: &Mesh) -> Mesh {
    let round = |x: f64| f64::from((x as f32).clamp(-f32::MAX, f32::MAX));
    let mut ids: HashMap<[u64; 3], u32> = HashMap::new();
    let mut verts = Vec::new();
    let tris = m
        .triangles()
        .iter()
        .map(|t| {
            t.map(|i| {
                let p = m.vertices()[i as usize].map(round);
                match ids.entry(p.map(|x| if x == 0.0 { 0 } else { x.to_bits() })) {
                    Entry::Occupied(e) => *e.get(),
                    Entry::Vacant(e) => {
                        verts.push(p);
                        *e.insert(verts.len() as u32 - 1)
                    }
                }
            })
        })
        .collect();
    let mut out = Mesh::new(verts, tris).expect("finite and in range");
    if out.n_faces() == 0 {
        return out;
    }
    let fa = m.face_attributes();
    if let Some(c) = fa.column("solid")
        && matches!(c.data, ColumnData::Utf8(_))
    {
        out = out
            .with_face_column("solid", c.clone())
            .expect("one row per face");
    }
    if let Some(c) = fa.column("stl_attribute")
        && let ColumnData::Int64(v) = &c.data
    {
        let a: Vec<i64> = (0..v.len())
            .map(|i| if c.is_valid(i) { v[i] } else { 0 })
            .collect();
        if a.iter().all(|&x| u16::try_from(x).is_ok()) && a.iter().any(|&x| x != 0) {
            out = out
                .with_face_column("stl_attribute", Column::i64(a))
                .expect("one row per face");
        }
    }
    out
}

fn round_trip(m: &Mesh, format: Format) -> Mesh {
    let data = Data {
        mesh: Some(m.clone()),
        ..Data::default()
    };
    let opts = WriteOptions {
        skip_unsupported_columns: true,
        ..WriteOptions::default()
    };
    let back = match format {
        Format::Dxf => crate::dxf::write(&data, &opts).and_then(|b| crate::dxf::read(&b)),
        Format::Vtk => crate::vtk::write(&data, &opts).and_then(|b| crate::vtk::read(&b)),
        Format::Vtu => crate::vtu::write(&data, &opts).and_then(|b| crate::vtu::read(&b)),
        Format::Vtp => crate::vtu::write_vtp(&data, &opts).and_then(|b| crate::vtu::read(&b)),
        Format::DatamineCsv => crate::datamine_csv::write(&data, &opts)
            .and_then(|(pt, tr)| crate::datamine_csv::read(&pt, &tr)),
        Format::Obj | Format::Ply | Format::Off | Format::Stl => unreachable!("handled above"),
    };
    match back {
        Ok(d) => d.mesh.unwrap_or_else(Mesh::empty),
        Err(e) => panic!("canonicalize: {format:?} cannot write this mesh: {e}"),
    }
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

const NULL_NAME: &[u8] = b"%null";
const EMPTY_NAME: &[u8] = b"%empty";

fn hex(b: u8) -> Option<u8> {
    (b as char).to_digit(16).map(|d| d as u8)
}

/// Inverse of `encode_name`; `None` is a null name.
pub(crate) fn decode_name(raw: &[u8]) -> Option<String> {
    match raw {
        NULL_NAME => None,
        EMPTY_NAME => Some(String::new()),
        _ => {
            let mut out = Vec::with_capacity(raw.len());
            let mut i = 0;
            while i < raw.len() {
                if raw[i] == b'%'
                    && let (Some(h), Some(l)) = (
                        raw.get(i + 1).copied().and_then(hex),
                        raw.get(i + 2).copied().and_then(hex),
                    )
                {
                    out.push(h << 4 | l);
                    i += 3;
                } else {
                    out.push(raw[i]);
                    i += 1;
                }
            }
            Some(String::from_utf8_lossy(&out).into_owned())
        }
    }
}

/// Writes a name that reads back exactly from a line whose outer whitespace is
/// trimmed: `%`, `#`, control bytes, a trailing `\`, outer spaces and, with
/// `single_spaces`, every space next to another space are written as `%XX`.
pub(crate) fn encode_name(name: Option<&str>, single_spaces: bool, out: &mut Vec<u8>) {
    let Some(name) = name else {
        out.extend_from_slice(NULL_NAME);
        return;
    };
    if name.is_empty() {
        out.extend_from_slice(EMPTY_NAME);
        return;
    }
    let b = name.as_bytes();
    let last = b.len() - 1;
    for (i, &c) in b.iter().enumerate() {
        let raw = match c {
            b'%' | b'#' | 0..=0x1f | 0x7f => false,
            b' ' => i > 0 && i < last && !(single_spaces && (b[i - 1] == b' ' || b[i + 1] == b' ')),
            b'\\' => i < last,
            _ => true,
        };
        if raw {
            out.push(c);
        } else {
            out.extend_from_slice(format!("%{c:02X}").as_bytes());
        }
    }
}
