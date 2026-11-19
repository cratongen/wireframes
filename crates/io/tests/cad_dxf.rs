mod cad_common;

use cad_common::*;
use proptest::prelude::*;
use wireframes_core::{Column, ColumnData, Error, Mesh, PointSet, Polylines, Table};
use wireframes_io::{Data, DxfEntity, WriteOptions, dxf};

fn polyface() -> WriteOptions {
    WriteOptions {
        dxf_entity: DxfEntity::Polyface,
        ..WriteOptions::default()
    }
}

fn layered(m: Mesh, seed: u64) -> Mesh {
    let ids: Vec<i64> = match m.face_attributes().column("polygon").map(|c| &c.data) {
        Some(ColumnData::Int64(v)) => v.clone(),
        _ => (0..m.n_faces() as i64).collect(),
    };
    let layer = ids
        .iter()
        .map(|&k| ["0", "ORE", "WASTE"][((k as u64 + seed) % 3) as usize].to_string())
        .collect();
    let color: Vec<i64> = ids
        .iter()
        .map(|&k| if k % 2 == 0 { 256 } else { k % 255 + 1 })
        .collect();
    let mut ft = Table::new(ids.len())
        .with_column("layer", Column::utf8(layer))
        .unwrap();
    if color.iter().any(|&c| c != 256) {
        ft = ft.with_column("color", Column::i64(color)).unwrap();
    }
    with_faces(m, ft)
}

fn unit_mesh(m: Mesh) -> Mesh {
    m.with_length_unit(Some("m".into()))
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    #[test]
    fn canonical_meshes_round_trip(m in canonical_mesh(4), seed in 0u64..3, pf in any::<bool>()) {
        prop_assume!(m.n_faces() > 0);
        let m = unit_mesh(layered(m, seed));
        let opts = if pf { polyface() } else { WriteOptions::default() };
        let data = Data { length_unit: Some("m".into()), ..mesh_data(m.clone()) };
        let back = dxf::read(&dxf::write(&data, &opts).unwrap()).unwrap();
        prop_assert_eq!(back.length_unit.as_deref(), Some("m"));
        prop_assert_eq!(back.mesh, Some(m));
    }

    #[test]
    fn polylines_and_points_round_trip(pl in polylines(false, false), ps in points()) {
        let n = pl.len();
        let pl = pl
            .with_column("layer", Column::utf8((0..n).map(|i| format!("S{i}")).collect()))
            .unwrap()
            .with_column("color", Column::i64((0..n as i64).map(|i| i + 1).collect()))
            .unwrap();
        let ps = ps.clone()
            .with_column("layer", Column::utf8(vec!["PTS".into(); ps.len()]))
            .unwrap();
        let data = Data { polylines: Some(pl), points: Some(ps), ..Data::default() };
        let back = dxf::read(&dxf::write(&data, &Default::default()).unwrap()).unwrap();
        prop_assert_eq!(back.polylines, data.polylines);
        prop_assert_eq!(back.points, data.points);
        prop_assert!(back.mesh.is_none());
    }
}

#[test]
fn polyface_splits_large_layers() {
    let n = 40_000u32;
    let mut verts = Vec::new();
    for i in 0..n + 2 {
        verts.push([i as f64, (i % 2) as f64, 0.0]);
    }
    let tris: Vec<[u32; 3]> = (0..n).map(|i| [i, i + 1, i + 2]).collect();
    let m = Mesh::new(verts, tris)
        .unwrap()
        .with_face_column("layer", Column::utf8(vec!["0".into(); n as usize]))
        .unwrap();
    let bytes = dxf::write(&mesh_data(m.clone()), &polyface()).unwrap();
    let text = String::from_utf8(bytes.clone()).unwrap();
    assert_eq!(text.matches("\nPOLYLINE\n").count(), 2);
    assert_eq!(dxf::read(&bytes).unwrap().mesh, Some(m));
}

#[test]
fn unsupported_columns_error_unless_skipped() {
    let m = Mesh::new(
        vec![[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        vec![[0, 1, 2]],
    )
    .unwrap()
    .with_face_column("grade", Column::f64(vec![1.0]))
    .unwrap();
    assert!(matches!(
        dxf::write(&mesh_data(m.clone()), &Default::default()),
        Err(Error::Unsupported(_))
    ));
    let opts = WriteOptions {
        skip_unsupported_columns: true,
        ..WriteOptions::default()
    };
    let back = dxf::read(&dxf::write(&mesh_data(m), &opts).unwrap()).unwrap();
    assert_eq!(back.mesh.unwrap().face_attributes().names(), &["layer"]);
}

fn groups(pairs: &[(i32, &str)]) -> String {
    pairs
        .iter()
        .map(|(c, v)| format!("{c:>3}\n{v}\n"))
        .collect()
}

fn entities(body: &str) -> String {
    groups(&[(0, "SECTION"), (2, "ENTITIES")]) + body + &groups(&[(0, "ENDSEC"), (0, "EOF")])
}

#[test]
fn reads_polyface_with_negative_indices() {
    let mut body = groups(&[
        (0, "POLYLINE"),
        (8, "PF"),
        (66, "1"),
        (70, "64"),
        (71, "4"),
        (72, "2"),
    ]);
    for p in [[0, 0], [1, 0], [1, 1], [0, 1]] {
        body += &groups(&[
            (0, "VERTEX"),
            (10, &p[0].to_string()),
            (20, &p[1].to_string()),
            (30, "0"),
            (70, "192"),
        ]);
    }
    body += &groups(&[
        (0, "VERTEX"),
        (70, "128"),
        (71, "1"),
        (72, "-2"),
        (73, "3"),
        (74, "-4"),
    ]);
    body += &groups(&[(0, "VERTEX"), (70, "128"), (71, "1"), (72, "3")]);
    body += &groups(&[(0, "SEQEND")]);
    let d = dxf::read(entities(&body).as_bytes()).unwrap();
    let m = d.mesh.unwrap();
    assert_eq!(m.triangles(), &[[0, 1, 2], [0, 2, 3]]);
    assert_eq!(
        m.face_attributes().column("layer").unwrap().data,
        ColumnData::Utf8(vec!["PF".into(), "PF".into()])
    );
    assert_eq!(
        d.warnings,
        vec!["skipped 1 faces with fewer than 3 corners"]
    );
}

#[test]
fn reads_cp1252_layer_names_and_units() {
    let mut text = groups(&[
        (0, "SECTION"),
        (2, "HEADER"),
        (9, "$ACADVER"),
        (1, "AC1015"),
        (9, "$INSUNITS"),
        (70, "2"),
        (0, "ENDSEC"),
    ])
    .into_bytes();
    text.extend(
        entities(&groups(&[
            (0, "POINT"),
            (8, "MIN\u{1}"),
            (10, "1"),
            (20, "2"),
            (30, "3"),
        ]))
        .into_bytes(),
    );
    let at = text.iter().position(|&b| b == 1).unwrap();
    text[at] = 0xC9;
    let d = dxf::read(&text).unwrap();
    assert_eq!(d.length_unit.as_deref(), Some("ft"));
    let ps = d.points.unwrap();
    assert_eq!(ps.coords(), &[[1.0, 2.0, 3.0]]);
    assert_eq!(
        ps.attributes().column("layer").unwrap().data,
        ColumnData::Utf8(vec!["MIN\u{C9}".into()])
    );
}

#[test]
fn reads_polygon_mesh_mesh_entity_lines_and_skips() {
    let mut body = groups(&[(0, "POLYLINE"), (70, "17"), (71, "2"), (72, "2")]);
    for p in [[0, 0], [0, 1], [1, 0], [1, 1]] {
        body += &groups(&[
            (0, "VERTEX"),
            (10, &p[0].to_string()),
            (20, &p[1].to_string()),
            (30, "5"),
            (70, "64"),
        ]);
    }
    body += &groups(&[(0, "SEQEND")]);
    body += &groups(&[
        (0, "MESH"),
        (8, "M"),
        (71, "2"),
        (72, "0"),
        (91, "0"),
        (92, "3"),
        (10, "10"),
        (20, "0"),
        (30, "0"),
        (10, "11"),
        (20, "0"),
        (30, "0"),
        (10, "10"),
        (20, "1"),
        (30, "0"),
        (93, "4"),
        (90, "3"),
        (90, "0"),
        (90, "1"),
        (90, "2"),
        (94, "0"),
    ]);
    body += &groups(&[
        (0, "LINE"),
        (10, "0"),
        (20, "0"),
        (30, "0"),
        (11, "3"),
        (21, "4"),
        (31, "0"),
    ]);
    body += &groups(&[
        (0, "LWPOLYLINE"),
        (90, "3"),
        (70, "1"),
        (38, "2"),
        (10, "0"),
        (20, "0"),
        (10, "1"),
        (20, "0"),
        (42, "0.5"),
        (10, "1"),
        (20, "1"),
        (210, "0"),
        (220, "0"),
        (230, "-1"),
    ]);
    body += &groups(&[(0, "INSERT"), (2, "BLK")]);
    body += &groups(&[(0, "INSERT"), (2, "BLK")]);
    let d = dxf::read(entities(&body).as_bytes()).unwrap();
    let m = d.mesh.unwrap();
    assert_eq!(m.n_faces(), 2 * 2 + 1);
    assert_eq!(m.vertices()[0], [0.0, 0.0, 5.0]);
    let poly = m.face_attributes().column("polygon").unwrap();
    assert_eq!(poly.data, ColumnData::Int64(vec![0, 0, 1, 1, 2]));
    let pl = d.polylines.unwrap();
    assert_eq!(pl.len(), 2);
    assert_eq!(pl.part(0), &[[0.0, 0.0, 0.0], [3.0, 4.0, 0.0]]);
    assert_eq!(
        pl.part(1),
        &[[-0.0, 0.0, -2.0], [-1.0, 0.0, -2.0], [-1.0, 1.0, -2.0]]
    );
    assert_eq!(pl.closed(), &[false, true]);
    assert_eq!(
        d.warnings,
        vec![
            "skipped 2 INSERT entities",
            "ignored nonzero bulge on 1 polyline vertices; arcs read as straight segments"
        ]
    );
}

#[test]
fn reads_binary_dxf() {
    let mut b = b"AutoCAD Binary DXF\r\n\x1a\0".to_vec();
    let s = |b: &mut Vec<u8>, code: i16, v: &str| {
        b.extend(code.to_le_bytes());
        b.extend(v.as_bytes());
        b.push(0);
    };
    let f = |b: &mut Vec<u8>, code: i16, v: f64| {
        b.extend(code.to_le_bytes());
        b.extend(v.to_le_bytes());
    };
    s(&mut b, 0, "SECTION");
    s(&mut b, 2, "ENTITIES");
    s(&mut b, 0, "3DFACE");
    s(&mut b, 8, "B");
    b.extend(62i16.to_le_bytes());
    b.extend(3i16.to_le_bytes());
    for (k, p) in [
        [0.1, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
    ]
    .iter()
    .enumerate()
    {
        for (a, &x) in p.iter().enumerate() {
            f(&mut b, (10 + k + 10 * a) as i16, x);
        }
    }
    s(&mut b, 0, "ENDSEC");
    s(&mut b, 0, "EOF");
    let d = dxf::read(&b).unwrap();
    let m = d.mesh.unwrap();
    assert_eq!(m.vertices()[0], [0.1, 0.0, 0.0]);
    assert_eq!(m.triangles(), &[[0, 1, 2], [0, 2, 3]]);
    assert_eq!(
        m.face_attributes().column("color").unwrap().data,
        ColumnData::Int64(vec![3, 3])
    );
    never_panics(&b, |x| {
        let _ = dxf::read(x);
    });
}

#[test]
fn malformed_files_error() {
    let cases: Vec<String> = vec![
        "  0\nSECTION\n  2\n".into(),
        "abc\nSECTION\n".into(),
        entities(&groups(&[(0, "3DFACE"), (10, "x")])),
        entities(&groups(&[
            (0, "POLYLINE"),
            (70, "64"),
            (0, "VERTEX"),
            (70, "128"),
            (71, "1"),
            (72, "2"),
            (73, "3"),
        ])),
        entities(&groups(&[
            (0, "POLYLINE"),
            (70, "16"),
            (71, "3"),
            (72, "3"),
            (0, "VERTEX"),
        ])),
        entities(&groups(&[(0, "MESH"), (92, "5"), (10, "0")])),
    ];
    for c in &cases {
        assert!(dxf::read(c.as_bytes()).is_err(), "{c:?}");
    }
    assert!(matches!(
        dxf::write(
            &Data {
                points: Some(
                    PointSet::new(vec![[0.0; 3]])
                        .unwrap()
                        .with_column("layer", Column::utf8(vec!["a\nb".into()]))
                        .unwrap()
                ),
                ..Data::default()
            },
            &Default::default()
        ),
        Err(Error::InvalidInput(_))
    ));
}

#[test]
fn damaged_files_never_panic() {
    let m = Mesh::new(
        vec![[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 1.0]],
        vec![[0, 1, 2], [1, 3, 2]],
    )
    .unwrap();
    let pl = Polylines::from_parts(
        vec![vec![[0.0; 3], [2.0, 2.0, 2.0], [2.0, 3.0, 2.0]]],
        vec![true],
        None,
    )
    .unwrap();
    let data = Data {
        mesh: Some(m),
        polylines: Some(pl),
        points: Some(PointSet::new(vec![[9.0; 3]]).unwrap()),
        ..Data::default()
    };
    for opts in [WriteOptions::default(), polyface()] {
        never_panics(&dxf::write(&data, &opts).unwrap(), |b| {
            let _ = dxf::read(b);
        });
    }
}

fn layered_polylines(pl: Polylines) -> Polylines {
    let n = pl.len();
    let pl = pl
        .with_column("layer", Column::utf8(vec!["S".into(); n]))
        .unwrap();
    let data = Data {
        polylines: Some(pl),
        ..Data::default()
    };
    dxf::read(&dxf::write(&data, &WriteOptions::default()).unwrap())
        .unwrap()
        .polylines
        .unwrap()
}

#[test]
fn multipart_feature_round_trips() {
    let pl = Polylines::from_parts(
        vec![
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            vec![[5.0, 0.0, 0.0], [6.0, 0.0, 0.0]],
            vec![[7.0, 0.0, 0.0], [8.0, 0.0, 0.0]],
            vec![[9.0, 0.0, 0.0], [9.5, 0.0, 0.0]],
        ],
        vec![false; 4],
        Some(vec![0, 0, 1, 1]),
    )
    .unwrap();
    let n = pl.len();
    let expect = pl
        .clone()
        .with_column("layer", Column::utf8(vec!["S".into(); n]))
        .unwrap();
    assert_eq!(layered_polylines(pl), expect);
}

#[test]
fn closed_part_with_repeated_end_round_trips() {
    let ring = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 0.0, 0.0],
    ];
    let pl = Polylines::from_parts(vec![ring], vec![true], None).unwrap();
    let back = layered_polylines(pl.clone());
    assert_eq!(back.n_vertices(), pl.n_vertices());
}

#[test]
fn empty_feature_is_refused() {
    let pl = Polylines::new(
        vec![[0.0; 3], [1.0; 3]],
        vec![0, 2],
        vec![0, 0, 1],
        vec![false],
        Table::new(2),
    )
    .unwrap();
    let data = Data {
        polylines: Some(pl),
        ..Data::default()
    };
    assert!(matches!(
        dxf::write(&data, &WriteOptions::default()),
        Err(Error::Unsupported(_))
    ));
}

fn point_layer(layer: &str) -> String {
    let ps = PointSet::new(vec![[1.0, 2.0, 3.0]])
        .unwrap()
        .with_column("layer", Column::utf8(vec![layer.to_string()]))
        .unwrap();
    let data = Data {
        points: Some(ps),
        ..Data::default()
    };
    let bytes = dxf::write(&data, &WriteOptions::default()).unwrap();
    assert!(bytes.is_ascii(), "an R2000 file is written as ASCII");
    let back = dxf::read(&bytes).unwrap();
    match &back
        .points
        .unwrap()
        .attributes()
        .column("layer")
        .unwrap()
        .data
    {
        ColumnData::Utf8(v) => v[0].clone(),
        _ => panic!("layer is not utf8"),
    }
}

#[test]
fn layer_names_round_trip_exactly() {
    for l in [
        " A ",
        "X\\U+0041",
        "caf\u{e9} \u{20AC}",
        "a\\b",
        "tab\there",
    ] {
        assert_eq!(point_layer(l), l);
    }
}

#[test]
fn mixed_case_unicode_escapes_are_both_decoded() {
    let body = "  0\nSECTION\n  2\nENTITIES\n  0\nPOINT\n  8\n\\u+0041\\U+0042\n 10\n1\n 20\n2\n 30\n3\n  0\nENDSEC\n  0\nEOF\n";
    let d = dxf::read(body.as_bytes()).unwrap();
    match &d.points.unwrap().attributes().column("layer").unwrap().data {
        ColumnData::Utf8(v) => assert_eq!(v[0], "AB"),
        _ => panic!("layer is not utf8"),
    }
}

#[test]
fn utf8_is_required_from_ac1021() {
    let file = |ver: &str| {
        let mut b = groups(&[
            (0, "SECTION"),
            (2, "HEADER"),
            (9, "$ACADVER"),
            (1, ver),
            (0, "ENDSEC"),
            (0, "SECTION"),
            (2, "ENTITIES"),
            (0, "POINT"),
        ])
        .into_bytes();
        b.extend_from_slice(b"  8\nL\xe9\n");
        b.extend(groups(&[(0, "ENDSEC"), (0, "EOF")]).into_bytes());
        b
    };
    let d = dxf::read(&file("AC1015")).unwrap();
    assert_eq!(
        d.points.unwrap().attributes().column("layer").unwrap().data,
        ColumnData::Utf8(vec!["L\u{e9}".into()])
    );
    assert!(matches!(
        dxf::read(&file("AC1021")),
        Err(Error::Format { .. })
    ));
}

#[test]
fn binary_mesh_is_skipped_with_a_warning() {
    let mut b = b"AutoCAD Binary DXF\r\n\x1a\0".to_vec();
    let s = |b: &mut Vec<u8>, code: i16, v: &str| {
        b.extend(code.to_le_bytes());
        b.extend(v.as_bytes());
        b.push(0);
    };
    s(&mut b, 0, "SECTION");
    s(&mut b, 2, "ENTITIES");
    s(&mut b, 0, "MESH");
    s(&mut b, 8, "M");
    s(&mut b, 0, "ENDSEC");
    s(&mut b, 0, "EOF");
    let d = dxf::read(&b).unwrap();
    assert!(d.mesh.is_none());
    assert_eq!(d.warnings, vec!["skipped 1 MESH entities".to_string()]);
}
