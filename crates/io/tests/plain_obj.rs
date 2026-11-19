mod plain_support;

use plain_support::*;
use proptest::prelude::*;
use wireframes_core::{Column, Error, PointSet, Polylines};
use wireframes_io::{Data, obj};

#[test]
fn fixture_with_negative_indices_and_forward_references() {
    let d = obj::read(&fixture("negative_forward.obj")).unwrap();
    let m = d.mesh.as_ref().unwrap();
    assert_eq!(m.n_vertices(), 5);
    assert_eq!(m.vertices()[4], [0.0, 0.0, 1.0]);
    assert_eq!(m.triangles(), &[[0, 1, 2], [0, 2, 3], [0, 1, 4]]);
    let fa = m.face_attributes();
    assert_eq!(fa.names(), &["object", "group", "polygon"]);
    assert_eq!(
        fa.column("object").unwrap(),
        &Column::utf8(vec!["slab".into(); 3])
    );
    assert_eq!(
        fa.column("group").unwrap(),
        &Column::utf8(vec!["top cap".into(), "top cap".into(), "side".into()])
    );
    assert_eq!(fa.column("polygon").unwrap(), &Column::i64(vec![0, 0, 1]));
    let pl = d.polylines.as_ref().unwrap();
    assert_eq!(pl.n_parts(), 1);
    assert_eq!(pl.closed(), &[true]);
    assert_eq!(
        pl.part(0),
        &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]]
    );
    assert_eq!(d.points.as_ref().unwrap().coords(), &[[1.0, 1.0, 0.0]]);
    assert!(d.warnings.is_empty(), "{:?}", d.warnings);
}

#[test]
fn fixture_reads_writes_and_reads_identically() {
    let d = obj::read(&fixture("negative_forward.obj")).unwrap();
    let again = obj::read(&obj::write(&d, &ascii()).unwrap()).unwrap();
    assert_eq!(again, d);
}

#[test]
fn polygons_write_back_as_polygons() {
    let d =
        obj::read(b"v 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 1 0\nv 2 2 0\nf 1 2 3 4 5\nf 1 2 3\n").unwrap();
    let text = String::from_utf8(obj::write(&d, &ascii()).unwrap()).unwrap();
    assert!(text.contains("f 1 2 3 4 5\nf 1 2 3\n"), "{text}");
}

#[test]
fn groups_default_only_when_some_group_exists() {
    let d = obj::read(b"v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\ng a  b\nf 1 3 2\n").unwrap();
    let g = d
        .mesh
        .as_ref()
        .unwrap()
        .face_attributes()
        .column("group")
        .unwrap();
    assert_eq!(g, &Column::utf8(vec!["default".into(), "a b".into()]));
    let d = obj::read(b"v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n").unwrap();
    assert_eq!(d.mesh.unwrap().face_attributes().num_columns(), 0);
}

#[test]
fn malformed_inputs_are_errors() {
    let cases: &[&[u8]] = &[
        b"v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 9\n",
        b"v 0 0 0\nv 1 0 0\nv 0 1 0\nf 0 1 2\n",
        b"v 0 0 0\nf -2 1 1\n",
        b"v 0 0\n",
        b"v nan 0 0\n",
        b"v inf 0 0\n",
        b"v 0 0 x\n",
        b"v 0 0 0\nf 1 a 1\n",
        b"v 0 0 0\nl 1 5\n",
        b"v 0 0 0\np 2\n",
        b"f 99999999999 1 1\n",
    ];
    for c in cases {
        let r = obj::read(c);
        assert!(
            matches!(r, Err(Error::Format { .. })),
            "{:?}: {r:?}",
            String::from_utf8_lossy(c)
        );
    }
}

#[test]
fn short_faces_and_unknown_statements_warn() {
    let d = obj::read(b"v 0 0 0\nv 1 0 0\nf 1 2\nl\ncurv 0 1\n").unwrap();
    assert_eq!(d.warnings.len(), 3, "{:?}", d.warnings);
    assert_eq!(d.mesh.unwrap().n_faces(), 0);
}

#[test]
fn unsupported_parts_and_columns() {
    let m = cube()
        .with_face_column("rock", Column::i64(vec![1; 12]))
        .unwrap();
    assert!(matches!(
        obj::write(&mesh_data(m.clone()), &ascii()),
        Err(Error::Unsupported(_))
    ));
    let out = obj::read(&obj::write(&mesh_data(m), &skipping(true)).unwrap()).unwrap();
    assert_eq!(out, mesh_data(cube()));
    let m = cube()
        .with_face_column("group", Column::i64(vec![1; 12]))
        .unwrap();
    assert!(obj::write(&mesh_data(m), &skipping(true)).is_err());
    let pl = Polylines::from_parts(
        vec![vec![[0.0; 3], [1.0; 3]], vec![[2.0; 3], [3.0; 3]]],
        vec![false, false],
        Some(vec![0, 0]),
    )
    .unwrap();
    let d = Data {
        polylines: Some(pl),
        ..Data::default()
    };
    assert!(matches!(
        obj::write(&d, &ascii()),
        Err(Error::Unsupported(_))
    ));
}

fn names() -> impl Strategy<Value = Option<String>> {
    let s = prop::sample::select(vec![
        "a",
        "rock type",
        "ore_2",
        "x y z",
        "",
        "zone#2",
        " pit ",
        "a  b",
        "%null",
        "%empty",
        "50%",
        "%41",
        "C:\\x\\",
        "default",
        "tab\there",
        "line\nbreak",
        "ção",
    ]);
    prop::option::weighted(0.9, s.prop_map(String::from))
}

fn name_column(v: Vec<Option<String>>) -> Column {
    let valid = v.iter().map(Option::is_some).collect();
    Column::utf8(v.into_iter().map(Option::unwrap_or_default).collect())
        .with_valid(valid)
        .unwrap()
}

proptest! {
    #[test]
    fn write_then_read_is_identity(
        (m, sizes) in polygon_mesh(12, 8),
        objects in prop::option::of(prop::collection::vec(names(), 8)),
        groups in prop::option::of(prop::collection::vec(names(), 8)),
        pl in prop::option::of(polylines()),
        ps in prop::option::of(points()),
    ) {
        let poly = m.face_attributes().column("polygon").cloned();
        let (vt, ft) = (m.vertex_attributes().clone(), wireframes_core::Table::new(m.n_faces()));
        let mut m = m.with_attributes(vt, ft).unwrap();
        if m.n_faces() > 0 {
            if let Some(o) = objects {
                m = m.with_face_column("object", name_column(per_polygon(&sizes, &o))).unwrap();
            }
            if let Some(g) = groups {
                m = m.with_face_column("group", name_column(per_polygon(&sizes, &g))).unwrap();
            }
        }
        if let Some(c) = poly {
            m = m.with_face_column("polygon", c).unwrap();
        }
        let d = Data { mesh: Some(m), polylines: pl, points: ps, ..Data::default() };
        let back = obj::read(&obj::write(&d, &ascii()).unwrap()).unwrap();
        prop_assert_eq!(back, d);
    }

    #[test]
    fn arbitrary_bytes_never_panic(b in prop::collection::vec(any::<u8>(), 0..400)) {
        let _ = obj::read(&b);
    }

    #[test]
    fn token_soup_never_panics(
        lines in prop::collection::vec(
            (prop::sample::select(vec!["v", "f", "l", "p", "o", "g", "vt", "#", "\\"]),
             prop::collection::vec(prop::sample::select(vec!["1", "-1", "0", "2/3", "/", "1e400", "x", "-7", "3//"]), 0..5)),
            0..20)
    ) {
        let text: String = lines.iter().map(|(k, a)| format!("{k} {}\n", a.join(" "))).collect();
        let _ = obj::read(text.as_bytes());
    }
}

#[test]
fn points_only_file_has_no_mesh() {
    let d = obj::read(b"v 1 2 3\nv 4 5 6\np 1 2\n").unwrap();
    assert!(d.mesh.is_none());
    assert_eq!(
        d.points,
        Some(PointSet::new(vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]).unwrap())
    );
    let empty = obj::read(b"# nothing\n").unwrap();
    assert_eq!(empty, Data::default());
}

fn tri() -> wireframes_core::Mesh {
    wireframes_core::Mesh::new(
        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        vec![[0, 1, 2]],
    )
    .unwrap()
}

fn group_round_trip(col: Column) {
    let d = mesh_data(tri().with_face_column("group", col).unwrap());
    let back = obj::read(&obj::write(&d, &ascii()).unwrap()).unwrap();
    assert_eq!(back, d);
}

#[test]
fn group_name_with_hash_round_trips() {
    group_round_trip(Column::utf8(vec!["zone#2".into()]));
}

#[test]
fn empty_group_name_round_trips() {
    group_round_trip(Column::utf8(vec![String::new()]));
}

#[test]
fn null_group_name_round_trips() {
    group_round_trip(
        Column::utf8(vec!["a".into()])
            .with_valid(vec![false])
            .unwrap(),
    );
}

#[test]
fn comment_ending_in_backslash_does_not_swallow_next_line() {
    let src = b"# exported from C:\\models\\\nv 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n";
    let d = obj::read(src).unwrap();
    assert_eq!(d.mesh.unwrap().n_vertices(), 3);
    let d = obj::read(b"v 0 0 \\\n0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n").unwrap();
    assert_eq!(d.mesh.unwrap().n_vertices(), 3);
}

#[test]
fn point_parts_round_trip() {
    let pl = Polylines::from_parts(
        vec![
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]],
            vec![[5.0, 5.0, 1.0]],
            vec![[6.0, 5.0, 1.0], [7.0, 5.0, 1.0]],
        ],
        vec![true, false, false],
        None,
    )
    .unwrap();
    let data = Data {
        polylines: Some(pl.clone()),
        ..Data::default()
    };
    let back = obj::read(&obj::write(&data, &ascii()).unwrap())
        .unwrap()
        .polylines
        .unwrap();
    assert_eq!(back.n_parts(), 3);
    for p in 0..3 {
        assert_eq!(
            (back.part(p), back.closed()[p]),
            (pl.part(p), pl.closed()[p])
        );
    }
}
