mod plain_support;

use plain_support::*;
use proptest::prelude::*;
use wireframes_core::{Column, ColumnData, Error, Mesh, Table};
use wireframes_io::{Data, ply};

#[test]
fn unknown_element_between_vertex_and_face() {
    let d = ply::read(&fixture("unknown_element.ply")).unwrap();
    let m = d.mesh.as_ref().unwrap();
    assert_eq!(m.n_vertices(), 4);
    assert_eq!(m.triangles(), &[[0, 1, 2], [0, 2, 3], [0, 1, 3]]);
    let va = m.vertex_attributes();
    assert_eq!(
        va.column("red").unwrap(),
        &Column::i64(vec![255, 0, 10, 20])
    );
    let q = va.column("quality").unwrap();
    assert_eq!(q.null_count(), 1);
    assert!(!q.is_valid(1));
    let fa = m.face_attributes();
    assert_eq!(fa.names(), &["label", "polygon"]);
    assert_eq!(fa.column("label").unwrap(), &Column::i64(vec![-3, -3, 9]));
    assert_eq!(fa.column("polygon").unwrap(), &Column::i64(vec![0, 0, 1]));
}

#[test]
fn big_endian_with_tristrips() {
    let d = ply::read(&fixture("tristrips_be.ply")).unwrap();
    let m = d.mesh.as_ref().unwrap();
    assert_eq!(m.vertices()[4], [2.0, 0.0, 0.0]);
    assert_eq!(m.triangles(), &[[0, 1, 2], [0, 2, 3], [1, 4, 2], [0, 1, 3]]);
    let va = m.vertex_attributes();
    assert_eq!(
        va.column("id").unwrap(),
        &Column::i64(vec![-2, -1, 0, 1, 2])
    );
    let w = m.face_attributes().column("weight").unwrap();
    assert_eq!(
        (w.is_valid(0), w.is_valid(1), w.is_valid(2), w.is_valid(3)),
        (true, true, false, false)
    );
    let ColumnData::Float64(v) = &w.data else {
        panic!()
    };
    assert_eq!(v[0], 0.25);
    assert_eq!(
        m.face_attributes().column("polygon").unwrap(),
        &Column::i64(vec![0, 0, 1, 2])
    );
}

#[test]
fn fixtures_read_write_read_identically() {
    for name in ["unknown_element.ply", "tristrips_be.ply"] {
        let d = ply::read(&fixture(name)).unwrap();
        for opts in [ascii(), binary()] {
            let again = ply::read(&ply::write(&d, &opts).unwrap()).unwrap();
            assert_eq!(again, d, "{name}");
        }
    }
}

#[test]
fn polygons_write_back_as_polygons() {
    let d = ply::read(&fixture("unknown_element.ply")).unwrap();
    let text = String::from_utf8(ply::write(&d, &ascii()).unwrap()).unwrap();
    assert!(text.ends_with("4 0 1 2 3 -3\n3 0 1 3 9\n"), "{text}");
}

#[test]
fn wide_polygons_use_int_counts() {
    let n = 300u32;
    let verts: Vec<[f64; 3]> = (0..n)
        .map(|i| [f64::from(i).cos(), f64::from(i).sin(), 0.0])
        .collect();
    let tris: Vec<[u32; 3]> = (1..n - 1).map(|k| [0, k, k + 1]).collect();
    let m = Mesh::new(verts, tris)
        .unwrap()
        .with_face_column("polygon", Column::i64(vec![0; n as usize - 2]))
        .unwrap();
    let d = mesh_data(m);
    for opts in [ascii(), binary()] {
        let bytes = ply::write(&d, &opts).unwrap();
        assert!(String::from_utf8_lossy(&bytes).contains("property list int int vertex_indices"));
        assert_eq!(ply::read(&bytes).unwrap(), d);
    }
}

#[test]
fn column_type_rules() {
    let m = cube()
        .with_vertex_column(
            "flag",
            Column::bool(vec![true, false, true, false, true, false, true, false]),
        )
        .unwrap()
        .with_vertex_column("big", Column::i64(vec![1 << 40, 0, 0, 0, 0, 0, 0, 0]))
        .unwrap()
        .with_vertex_column(
            "gappy",
            Column::i64(vec![1; 8])
                .with_valid(vec![true, false, true, true, true, true, true, true])
                .unwrap(),
        )
        .unwrap();
    for opts in [ascii(), binary()] {
        let back = ply::read(&ply::write(&mesh_data(m.clone()), &opts).unwrap()).unwrap();
        let va = back.mesh.as_ref().unwrap().vertex_attributes().clone();
        assert_eq!(
            va.column("flag").unwrap(),
            &Column::i64(vec![1, 0, 1, 0, 1, 0, 1, 0])
        );
        assert_eq!(va.column("big").unwrap().type_name(), "float64");
        let g = va.column("gappy").unwrap();
        assert_eq!((g.type_name(), g.null_count()), ("float64", 1));
    }
    let m = cube()
        .with_face_column("name", Column::utf8(vec![String::new(); 12]))
        .unwrap();
    assert!(matches!(
        ply::write(&mesh_data(m.clone()), &binary()),
        Err(Error::Unsupported(_))
    ));
    let back = ply::read(&ply::write(&mesh_data(m), &skipping(false)).unwrap()).unwrap();
    assert_eq!(back, mesh_data(cube()));
    let m = cube()
        .with_vertex_column("x", Column::f64(vec![0.0; 8]))
        .unwrap();
    assert!(matches!(
        ply::write(&mesh_data(m), &binary()),
        Err(Error::Unsupported(_))
    ));
}

#[test]
fn malformed_inputs_are_errors() {
    let head = "ply\nformat ascii 1.0\nelement vertex 3\nproperty float x\nproperty float y\nproperty float z\nelement face 1\nproperty list uchar int vertex_indices\nend_header\n";
    let cases: Vec<Vec<u8>> = vec![
        b"plx\n".to_vec(),
        b"ply\nformat ascii 1.0\nelement vertex 1\n".to_vec(),
        b"ply\nformat weird 1.0\nend_header\n".to_vec(),
        b"ply\nelement vertex 0\nend_header\n".to_vec(),
        b"ply\nformat ascii 1.0\nproperty float x\nend_header\n".to_vec(),
        b"ply\nformat ascii 1.0\nelement vertex 1\nproperty float x\nend_header\n0\n".to_vec(),
        b"ply\nformat ascii 1.0\nelement vertex 1\nproperty blob x\nend_header\n".to_vec(),
        b"ply\nformat ascii 1.0\nelement vertex 1\nproperty list float int x\nend_header\n".to_vec(),
        format!("{head}0 0 0\n1 0 0\n0 1 0\n3 0 1 7\n").into_bytes(),
        format!("{head}0 0 0\n1 0 0\n0 1 0\n3 0 1 -1\n").into_bytes(),
        format!("{head}0 0 0\n1 0 0\n0 1 0\n3 0 1\n").into_bytes(),
        format!("{head}0 0 0\n1 0 0\n0 1 nan\n3 0 1 2\n").into_bytes(),
        format!("{head}0 0 0\n1 0 0\n0 1 0\n3 0 1 x\n").into_bytes(),
        format!("{head}0 0 0\n1 0 0\n").into_bytes(),
        b"ply\nformat binary_little_endian 1.0\nelement vertex 2\nproperty float x\nproperty float y\nproperty float z\nend_header\n\0\0\0\0".to_vec(),
        b"ply\nformat binary_little_endian 1.0\nelement face 1\nproperty list uint int vertex_indices\nend_header\n\xff\xff\xff\xff".to_vec(),
        b"ply\nformat ascii 1.0\nelement face 1\nproperty list uchar int vertex_indices\nend_header\n3 0 1 2\n".to_vec(),
    ];
    for c in &cases {
        let r = ply::read(c);
        assert!(
            matches!(r, Err(Error::Format { .. })),
            "{:?}: {r:?}",
            String::from_utf8_lossy(c)
        );
    }
    let r = ply::read(b"ply\nformat ascii 1.0\nelement vertex 1\nproperty int64 x\nend_header\n");
    assert!(matches!(r, Err(Error::Unsupported(_))));
}

#[test]
fn huge_declared_counts_fail_fast() {
    let r = ply::read(b"ply\nformat binary_little_endian 1.0\nelement vertex 18446744073709551615\nproperty float x\nproperty float y\nproperty float z\nend_header\n");
    assert!(r.is_err());
    let r = ply::read(b"ply\nformat ascii 1.0\nelement blob 18446744073709551615\nend_header\n");
    assert_eq!(r.unwrap(), Data::default());
}

fn attr_mesh() -> impl Strategy<Value = Mesh> {
    (polygon_mesh(10, 6), any::<bool>(), any::<bool>()).prop_flat_map(|((m, sizes), fv, iv)| {
        let nv = m.n_vertices();
        let np = sizes.len();
        let fcol = prop::collection::vec(prop::option::of(coord()), nv);
        let icol = prop::collection::vec(any::<i32>(), np);
        (Just((m, sizes)), fcol, icol).prop_map(move |((m, sizes), f, i)| {
            let poly = m.face_attributes().column("polygon").cloned();
            let mut vt = Table::new(m.n_vertices());
            if fv {
                let valid = f.iter().map(Option::is_some).collect();
                let v = f.iter().map(|x| x.unwrap_or(f64::NAN)).collect();
                vt = vt
                    .with_column("grade", Column::f64(v).with_valid(valid).unwrap())
                    .unwrap();
            }
            let mut ft = Table::new(m.n_faces());
            if iv {
                let vals: Vec<i64> = i.iter().map(|&x| i64::from(x)).collect();
                ft = ft
                    .with_column("domain", Column::i64(per_polygon(&sizes, &vals)))
                    .unwrap();
            }
            let mut m = m.with_attributes(vt, ft).unwrap();
            if let Some(c) = poly {
                m = m.with_face_column("polygon", c).unwrap();
            }
            m
        })
    })
}

proptest! {
    #[test]
    fn write_then_read_is_identity(m in attr_mesh(), asc in any::<bool>()) {
        let d = mesh_data(m);
        let opts = if asc { ascii() } else { binary() };
        let back = ply::read(&ply::write(&d, &opts).unwrap()).unwrap();
        prop_assert_eq!(back, d);
    }

    #[test]
    fn arbitrary_bytes_never_panic(b in prop::collection::vec(any::<u8>(), 0..400)) {
        let _ = ply::read(&b);
    }

    #[test]
    fn truncated_files_never_panic(m in attr_mesh(), asc in any::<bool>(), cut in 0.0..1.0f64) {
        let opts = if asc { ascii() } else { binary() };
        let bytes = ply::write(&mesh_data(m), &opts).unwrap();
        let n = (bytes.len() as f64 * cut) as usize;
        let _ = ply::read(&bytes[..n]);
    }

    #[test]
    fn corrupted_binary_never_panics(m in attr_mesh(), at in any::<prop::sample::Index>(), byte in any::<u8>()) {
        let mut bytes = ply::write(&mesh_data(m), &binary()).unwrap();
        let i = at.index(bytes.len());
        bytes[i] = byte;
        let _ = ply::read(&bytes);
    }
}

fn tri() -> wireframes_core::Mesh {
    wireframes_core::Mesh::new(
        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        vec![[0, 1, 2]],
    )
    .unwrap()
}

#[test]
fn face_property_named_polygon_survives_read_write_read() {
    let src = b"ply\nformat ascii 1.0\nelement vertex 4\nproperty float x\nproperty float y\nproperty float z\nelement face 2\nproperty list uchar int vertex_indices\nproperty int polygon\nend_header\n0 0 0\n1 0 0\n1 1 0\n0 1 0\n3 0 1 2 7\n3 0 2 3 7\n";
    let d = ply::read(src).unwrap();
    for opts in [ascii(), binary()] {
        let again = ply::read(&ply::write(&d, &opts).unwrap()).unwrap();
        assert_eq!(again, d);
    }
}

#[test]
fn face_property_named_polygon_is_not_overwritten() {
    let src = b"ply\nformat ascii 1.0\nelement vertex 4\nproperty float x\nproperty float y\nproperty float z\nelement face 1\nproperty list uchar int vertex_indices\nproperty int polygon\nend_header\n0 0 0\n1 0 0\n1 1 0\n0 1 0\n4 0 1 2 3 9\n";
    let d = ply::read(src).unwrap();
    let c = d
        .mesh
        .unwrap()
        .face_attributes()
        .column("polygon")
        .cloned()
        .unwrap();
    assert_eq!(c, Column::i64(vec![9, 9]));
}

#[test]
fn polygon_ids_with_nulls_and_gaps_round_trip() {
    let m = wireframes_core::Mesh::new(
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ],
        vec![[0, 1, 2], [0, 2, 3], [0, 1, 3]],
    )
    .unwrap()
    .with_face_column(
        "polygon",
        Column::i64(vec![40, 40, 0])
            .with_valid(vec![true, true, false])
            .unwrap(),
    )
    .unwrap();
    let d = mesh_data(m);
    for opts in [ascii(), binary()] {
        assert_eq!(ply::read(&ply::write(&d, &opts).unwrap()).unwrap(), d);
    }
}

#[test]
fn int64_beyond_f64_precision_is_refused_not_rounded() {
    let big = (1i64 << 53) + 1;
    let m = tri()
        .with_face_column("id", Column::i64(vec![big]))
        .unwrap();
    assert!(matches!(
        ply::write(&mesh_data(m.clone()), &binary()),
        Err(Error::Unsupported(_))
    ));
    let back = ply::read(&ply::write(&mesh_data(m), &skipping(false)).unwrap()).unwrap();
    assert_eq!(back, mesh_data(tri()));
    let exact = tri()
        .with_face_column("id", Column::i64(vec![1i64 << 53]))
        .unwrap();
    let back = ply::read(&ply::write(&mesh_data(exact), &binary()).unwrap()).unwrap();
    assert_eq!(
        back.mesh.unwrap().face_attributes().column("id").unwrap(),
        &Column::f64(vec![(1i64 << 53) as f64])
    );
}

proptest! {
    #[test]
    fn polygon_ids_anywhere_round_trip(
        m in attr_mesh(),
        off in -1000i64..1000,
        first in any::<bool>(),
        asc in any::<bool>(),
    ) {
        let mut m = m;
        let fa = m.face_attributes().clone();
        if let Some(ColumnData::Int64(v)) = fa.column("polygon").map(|c| &c.data) {
            let shifted = Column::i64(v.iter().map(|x| x + off).collect());
            let rest = fa.clone().without_column("polygon");
            let mut ft = Table::new(m.n_faces());
            if first {
                ft = ft.with_column("polygon", shifted.clone()).unwrap();
            }
            for (name, col) in rest.names().iter().zip(rest.columns()) {
                ft = ft.with_column(name, col.clone()).unwrap();
            }
            if !first {
                ft = ft.with_column("polygon", shifted).unwrap();
            }
            let vt = m.vertex_attributes().clone();
            m = m.with_attributes(vt, ft).unwrap();
        }
        let d = mesh_data(m);
        let opts = if asc { ascii() } else { binary() };
        prop_assert_eq!(ply::read(&ply::write(&d, &opts).unwrap()).unwrap(), d);
    }

    #[test]
    fn int64_columns_are_exact_or_refused(v in prop::collection::vec(any::<i64>(), 1), null in any::<bool>()) {
        let c = Column::i64(v.clone()).with_valid(vec![!null]).unwrap();
        let d = mesh_data(tri().with_face_column("id", c).unwrap());
        match ply::write(&d, &binary()) {
            Err(e) => prop_assert!(matches!(e, Error::Unsupported(_)) && v[0].unsigned_abs() > 1 << 53),
            Ok(bytes) => {
                let back = ply::read(&bytes).unwrap().mesh.unwrap();
                let c = back.face_attributes().column("id").unwrap();
                prop_assert_eq!(c.is_valid(0), !null);
                if !null {
                    let got = match &c.data {
                        ColumnData::Int64(w) => i128::from(w[0]),
                        ColumnData::Float64(w) => w[0] as i128,
                        _ => unreachable!(),
                    };
                    prop_assert_eq!(got, i128::from(v[0]));
                }
            }
        }
    }
}

#[test]
fn binary_error_offsets_count_from_the_start_of_the_file() {
    let head = b"ply\nformat binary_little_endian 1.0\nelement vertex 3\nproperty double x\nproperty double y\nproperty double z\nelement face 1\nproperty list uchar int vertex_indices\nend_header\n";
    let mut full = head.to_vec();
    for p in [[0.0f64, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]] {
        p.iter().for_each(|x| full.extend(x.to_le_bytes()));
    }
    full.push(3);
    [0i32, 1, 2]
        .iter()
        .for_each(|i| full.extend(i.to_le_bytes()));
    assert!(ply::read(&full).is_ok());
    let at = |e: Error| -> usize {
        let s = e.to_string();
        s[s.find("at byte ").expect(&s) + 8..s.len() - 1]
            .parse()
            .unwrap()
    };
    let mut nan = full.clone();
    let x1 = head.len() + 24;
    nan[x1..x1 + 8].copy_from_slice(&f64::NAN.to_le_bytes());
    assert_eq!(at(ply::read(&nan).unwrap_err()), x1);
    let cut = &full[..full.len() - 3];
    let off = at(ply::read(cut).unwrap_err());
    assert!((full.len() - 13..full.len()).contains(&off), "{off}");
}
