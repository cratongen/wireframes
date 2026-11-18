#![allow(dead_code)]

use std::collections::HashSet;
use std::path::PathBuf;

use drusa_core::{Column, Mesh, PointSet, Polylines};
use drusa_io::{Data, WriteOptions};
use proptest::prelude::*;

pub fn fixture(name: &str) -> Vec<u8> {
    let p = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/plain_fixtures")
        .join(name);
    std::fs::read(&p).unwrap_or_else(|e| panic!("{}: {e}", p.display()))
}

pub fn mesh_data(m: Mesh) -> Data {
    Data {
        mesh: Some(m),
        ..Data::default()
    }
}

pub fn ascii() -> WriteOptions {
    WriteOptions {
        ascii: true,
        ..WriteOptions::default()
    }
}

pub fn binary() -> WriteOptions {
    WriteOptions::default()
}

pub fn skipping(ascii: bool) -> WriteOptions {
    WriteOptions {
        ascii,
        skip_unsupported_columns: true,
        ..WriteOptions::default()
    }
}

pub fn coord() -> impl Strategy<Value = f64> {
    prop_oneof![
        any::<f64>().prop_filter("finite", |x| x.is_finite()),
        -1e7..1e7f64,
        (-1000i32..1000).prop_map(f64::from),
        Just(0.0),
        Just(-0.0),
        Just(f64::MIN_POSITIVE),
        Just(5e-324),
    ]
}

pub fn point() -> impl Strategy<Value = [f64; 3]> {
    [coord(), coord(), coord()]
}

/// Vertices and polygons of 3 to 6 corners, fanned with a `polygon` column when
/// some polygon is larger than a triangle, the convention the readers use.
pub fn polygon_mesh(max_v: usize, max_p: usize) -> impl Strategy<Value = (Mesh, Vec<usize>)> {
    prop::collection::vec(point(), 1..max_v).prop_flat_map(move |verts| {
        let n = verts.len() as u32;
        let poly = prop::collection::vec(0..n, 3..7);
        (Just(verts), prop::collection::vec(poly, 0..max_p)).prop_map(|(verts, polys)| {
            let mut tris = Vec::new();
            let mut ids = Vec::new();
            let mut sizes = Vec::new();
            for (k, c) in polys.iter().enumerate() {
                for i in 1..c.len() - 1 {
                    tris.push([c[0], c[i], c[i + 1]]);
                    ids.push(k as i64);
                }
                sizes.push(c.len() - 2);
            }
            let big = polys.iter().any(|c| c.len() > 3);
            let mut m = Mesh::new(verts, tris).unwrap();
            if big {
                m = m.with_face_column("polygon", Column::i64(ids)).unwrap();
            }
            (m, sizes)
        })
    })
}

/// Moves the `polygon` column, if any, to the end, the order readers produce.
pub fn polygon_last(m: Mesh) -> Mesh {
    match m.face_attributes().column("polygon").cloned() {
        None => m,
        Some(c) => {
            let ft = m.face_attributes().clone().without_column("polygon");
            let vt = m.vertex_attributes().clone();
            m.with_attributes(vt, ft)
                .unwrap()
                .with_face_column("polygon", c)
                .unwrap()
        }
    }
}

/// Per-face values repeated over the triangles of each polygon.
pub fn per_polygon<T: Clone>(sizes: &[usize], values: &[T]) -> Vec<T> {
    sizes
        .iter()
        .zip(values)
        .flat_map(|(&n, v)| std::iter::repeat_n(v.clone(), n))
        .collect()
}

pub fn polylines() -> impl Strategy<Value = Polylines> {
    let part = (any::<bool>(), prop::collection::vec(point(), 3..6));
    prop::collection::vec(part, 1..4).prop_map(|parts| {
        let closed = parts.iter().map(|p| p.0).collect();
        let coords = parts.into_iter().map(|p| p.1).collect();
        Polylines::from_parts(coords, closed, None).unwrap()
    })
}

pub fn points() -> impl Strategy<Value = PointSet> {
    prop::collection::vec(point(), 1..6).prop_map(|c| PointSet::new(c).unwrap())
}

/// A mesh canonical for STL: distinct f32-representable vertices, all used, in
/// first-use order over faces.
pub fn stl_mesh(max_f: usize) -> impl Strategy<Value = Mesh> {
    let f = || {
        any::<f32>().prop_filter("finite, not -0", |x| {
            x.is_finite() && x.to_bits() != (-0.0f32).to_bits()
        })
    };
    prop::collection::vec([f(), f(), f()], 1..12).prop_flat_map(move |raw| {
        let mut seen = HashSet::new();
        let verts: Vec<[f64; 3]> = raw
            .into_iter()
            .filter(|p| seen.insert(p.map(f32::to_bits)))
            .map(|p| p.map(f64::from))
            .collect();
        let n = verts.len() as u32;
        (
            Just(verts),
            prop::collection::vec([0..n, 0..n, 0..n], 1..max_f),
        )
            .prop_map(|(verts, tris)| {
                let mut map = vec![u32::MAX; verts.len()];
                let mut out_v = Vec::new();
                let tris = tris
                    .into_iter()
                    .map(|t| {
                        t.map(|i| {
                            if map[i as usize] == u32::MAX {
                                map[i as usize] = out_v.len() as u32;
                                out_v.push(verts[i as usize]);
                            }
                            map[i as usize]
                        })
                    })
                    .collect();
                Mesh::new(out_v, tris).unwrap()
            })
    })
}

pub fn cube() -> Mesh {
    drusa_core::fixtures::cube()
}
