#![allow(dead_code)]

use drusa_core::{Column, Mesh, PointSet, Polylines, Table};
use drusa_io::{Data, WriteOptions};
use proptest::prelude::*;

pub fn ascii() -> WriteOptions {
    WriteOptions {
        ascii: true,
        ..WriteOptions::default()
    }
}

pub fn zlib() -> WriteOptions {
    WriteOptions {
        zlib: true,
        ..WriteOptions::default()
    }
}

pub fn mesh_data(m: Mesh) -> Data {
    Data {
        mesh: Some(m),
        ..Data::default()
    }
}

pub fn coord() -> impl Strategy<Value = f64> {
    prop_oneof![
        any::<f64>().prop_filter("finite", |x| x.is_finite()),
        -1e7..1e7f64,
        (-1000i32..1000).prop_map(f64::from),
        Just(0.0),
        Just(5e-324),
    ]
}

/// Fans `polys` and adds a `polygon` column when some polygon has more than
/// three corners, the convention every reader follows.
pub fn fan(verts: Vec<[f64; 3]>, polys: &[Vec<u32>]) -> (Mesh, Vec<u32>) {
    let mut tris = Vec::new();
    let mut ids = Vec::new();
    let mut first = Vec::new();
    for (k, c) in polys.iter().enumerate() {
        for i in 1..c.len() - 1 {
            tris.push([c[0], c[i], c[i + 1]]);
            ids.push(k as i64);
            first.push(k as u32);
        }
    }
    let mut m = Mesh::new(verts, tris).unwrap();
    if polys.iter().any(|c| c.len() > 3) {
        m = m.with_face_column("polygon", Column::i64(ids)).unwrap();
    }
    (m, first)
}

/// Adds face columns (one value per polygon, repeated over its triangles) ahead
/// of `polygon`, as readers order them.
pub fn with_faces(m: Mesh, per_face: Table) -> Mesh {
    let poly = m.face_attributes().column("polygon").cloned();
    let vt = m.vertex_attributes().clone();
    let mut m = m.with_attributes(vt, per_face).unwrap();
    if let Some(p) = poly {
        m = m.with_face_column("polygon", p).unwrap();
    }
    m
}

/// A mesh with distinct vertices, all used, in first-use order: the canonical
/// form for formats that weld (DXF). Polygons have 3 or 4 distinct corners.
pub fn canonical_mesh(max_corners: usize) -> impl Strategy<Value = Mesh> {
    (3usize..12)
        .prop_flat_map(move |n| {
            let poly = prop::collection::vec(0..n as u32, 3..=max_corners);
            (
                prop::collection::vec([coord(), coord(), coord()], n),
                prop::collection::vec(poly, 1..8),
            )
        })
        .prop_map(|(raw, polys)| {
            let polys: Vec<Vec<u32>> = polys
                .into_iter()
                .map(|mut p| {
                    let mut seen = Vec::new();
                    p.retain(|v| {
                        let new = !seen.contains(v);
                        seen.push(*v);
                        new
                    });
                    p
                })
                .filter(|p| p.len() >= 3)
                .collect();
            let mut order: Vec<u32> = Vec::new();
            for p in &polys {
                for &v in p {
                    if !order.contains(&v) {
                        order.push(v);
                    }
                }
            }
            let verts: Vec<[f64; 3]> = order
                .iter()
                .map(|&v| {
                    let p = raw[v as usize];
                    [v as f64 + 0.5, p[1], p[2]]
                })
                .collect();
            let polys: Vec<Vec<u32>> = polys
                .iter()
                .map(|p| {
                    p.iter()
                        .map(|v| order.iter().position(|o| o == v).unwrap() as u32)
                        .collect()
                })
                .collect();
            fan(verts, &polys).0
        })
}

/// Vertices with Float64 (some null) and Int64 columns, polygons of 3 to 6
/// corners with Float64 and Int64 face columns; unused vertices allowed.
pub fn attributed_mesh() -> impl Strategy<Value = Mesh> {
    (1usize..10)
        .prop_flat_map(|n| {
            let poly = prop::collection::vec(0..n as u32, 3..7);
            (
                prop::collection::vec([coord(), coord(), coord()], n),
                prop::collection::vec(poly, 0..6),
                prop::collection::vec(prop::option::of(coord()), n),
                prop::collection::vec(any::<i64>(), n),
                prop::collection::vec((coord(), any::<i64>()), 6),
            )
        })
        .prop_map(|(verts, polys, vf, vi, ff)| {
            let n = verts.len();
            let (m, first) = fan(verts, &polys);
            let valid: Vec<bool> = vf.iter().map(Option::is_some).collect();
            let mut f = Column::f64(vf.iter().map(|x| x.unwrap_or(f64::NAN)).collect());
            if valid.iter().any(|v| !v) {
                f = f.with_valid(valid).unwrap();
            }
            let vt = Table::new(n)
                .with_column("grade", f)
                .unwrap()
                .with_column("code", Column::i64(vi))
                .unwrap();
            let fcol = Column::f64(first.iter().map(|&k| ff[k as usize].0).collect());
            let icol = Column::i64(first.iter().map(|&k| ff[k as usize].1).collect());
            let ft = Table::new(first.len())
                .with_column("density", fcol)
                .unwrap()
                .with_column("domain", icol)
                .unwrap();
            let m = with_faces(m, ft);
            let ft = m.face_attributes().clone();
            m.with_attributes(vt, ft).unwrap()
        })
}

pub fn polylines(attrs: bool) -> impl Strategy<Value = Polylines> {
    let part = (
        any::<bool>(),
        prop::collection::vec([coord(), coord(), coord()], 3..6),
    );
    prop::collection::vec(part, 1..5).prop_map(move |parts| {
        let n = parts.len();
        let closed = parts.iter().map(|p| p.0).collect();
        let coords = parts.into_iter().map(|p| p.1).collect();
        let pl = Polylines::from_parts(coords, closed, None).unwrap();
        if attrs {
            let ids = (0..n as i64).collect();
            pl.with_column("id", Column::i64(ids)).unwrap()
        } else {
            pl
        }
    })
}

pub fn points() -> impl Strategy<Value = PointSet> {
    prop::collection::vec([coord(), coord(), coord()], 1..6).prop_map(|c| PointSet::new(c).unwrap())
}

/// Every prefix and a set of single-byte corruptions must return, never panic.
pub fn never_panics(bytes: &[u8], read: impl Fn(&[u8])) {
    for end in 0..bytes.len() {
        read(&bytes[..end]);
    }
    for i in (0..bytes.len()).step_by(bytes.len() / 200 + 1) {
        for v in [b'0', b'\n', 0xff, b'-', b'<', b'9'] {
            let mut b = bytes.to_vec();
            b[i] = v;
            read(&b);
        }
    }
}
