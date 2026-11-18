mod plain_support;

use drusa_core::{Column, Error, Mesh, Polylines};
use drusa_io::{Data, stl};
use plain_support::*;
use proptest::prelude::*;

#[test]
fn binary_whose_header_starts_with_solid() {
    let d = stl::read(&fixture("solid_header_binary.stl")).unwrap();
    let m = d.mesh.as_ref().unwrap();
    assert_eq!(m.n_vertices(), 4);
    assert_eq!(m.triangles(), &[[0, 1, 2], [0, 2, 3]]);
    assert_eq!(
        m.face_attributes().column("stl_attribute").unwrap(),
        &Column::i64(vec![7, 0])
    );
    assert!(d.warnings.is_empty());
}

#[test]
fn ascii_with_two_solids_mixed_case_and_crlf() {
    let d = stl::read(&fixture("two_solids.stl")).unwrap();
    let m = d.mesh.as_ref().unwrap();
    assert_eq!(m.n_vertices(), 4, "corners are welded across solids");
    assert_eq!(m.triangles(), &[[0, 1, 2], [0, 2, 3]]);
    assert_eq!(
        m.face_attributes().column("solid").unwrap(),
        &Column::utf8(vec!["first part".into(), "second".into()])
    );
}

#[test]
fn fixtures_read_write_read_identically() {
    let d = stl::read(&fixture("solid_header_binary.stl")).unwrap();
    assert_eq!(stl::read(&stl::write(&d, &binary()).unwrap()).unwrap(), d);
    let d = stl::read(&fixture("two_solids.stl")).unwrap();
    assert_eq!(stl::read(&stl::write(&d, &ascii()).unwrap()).unwrap(), d);
}

#[test]
fn welding_is_exact_and_maps_negative_zero() {
    let text = b"solid\nfacet normal 0 0 0\nouter loop\nvertex 0 0 0\nvertex 1 0 0\nvertex 0 1 0\nendloop\nendfacet\n\
facet normal 0 0 0\nouter loop\nvertex -0.0 0 0\nvertex 0 1 0\nvertex 0 1.0000000000000002 0\nendloop\nendfacet\nendsolid\n";
    let m = stl::read(text).unwrap().mesh.unwrap();
    assert_eq!(m.n_vertices(), 4);
    assert_eq!(m.triangles(), &[[0, 1, 2], [0, 2, 3]]);
    assert!(m.face_attributes().column("solid").is_none());
}

#[test]
fn trailing_bytes_fall_back_with_a_warning() {
    let mut b = stl::write(&mesh_data(cube()), &binary()).unwrap();
    b.extend_from_slice(b"junk");
    let d = stl::read(&b).unwrap();
    assert_eq!(d.warnings.len(), 1);
    assert_eq!(d.mesh.unwrap().n_faces(), 12);
}

#[test]
fn malformed_inputs_are_errors() {
    let mut truncated = stl::write(&mesh_data(cube()), &binary()).unwrap();
    truncated.truncate(300);
    let mut nan = stl::write(&mesh_data(cube()), &binary()).unwrap();
    nan[84 + 12..84 + 16].copy_from_slice(&f32::NAN.to_le_bytes());
    let cases: Vec<Vec<u8>> = vec![
        b"".to_vec(),
        b"garbage".to_vec(),
        truncated,
        nan,
        b"solid x\nfacet normal 0 0 1\nouter loop\nvertex 0 0 0\n".to_vec(),
        b"solid x\nfacet normal 0 0 1\nouter loop\nvertex 0 0 0\nvertex 1 0\nendloop\n".to_vec(),
        b"solid x\nfacet normal 0 0 1\nouter loop\nvertex 0 0 inf\nvertex 1 0 0\nvertex 0 1 0\nendloop\nendfacet\n".to_vec(),
        b"solid x\nfacet normal 0 0 1\nouter loop\nvertex 0 0 0\nvertex 1 0 0\nvertex 0 1 0\nendloop\n".to_vec(),
        b"solid x\nbogus\n".to_vec(),
    ];
    for c in &cases {
        let r = stl::read(c);
        assert!(
            matches!(r, Err(Error::Format { .. })),
            "{:?}: {r:?}",
            String::from_utf8_lossy(c)
        );
    }
}

#[test]
fn writer_rules() {
    let pl = Polylines::from_parts(vec![vec![[0.0; 3], [1.0; 3]]], vec![false], None).unwrap();
    let d = Data {
        polylines: Some(pl),
        ..Data::default()
    };
    assert!(matches!(
        stl::write(&d, &binary()),
        Err(Error::Unsupported(_))
    ));
    let far = Mesh::new(vec![[1e300, 0.0, 0.0], [0.0; 3], [1.0; 3]], vec![[0, 1, 2]]).unwrap();
    assert!(matches!(
        stl::write(&mesh_data(far.clone()), &binary()),
        Err(Error::InvalidInput(_))
    ));
    assert_eq!(
        stl::read(&stl::write(&mesh_data(far.clone()), &ascii()).unwrap()).unwrap(),
        mesh_data(far)
    );
    let m = cube()
        .with_face_column("stl_attribute", Column::i64(vec![70000; 12]))
        .unwrap();
    assert!(matches!(
        stl::write(&mesh_data(m.clone()), &binary()),
        Err(Error::InvalidInput(_))
    ));
    assert!(matches!(
        stl::write(&mesh_data(m.clone()), &ascii()),
        Err(Error::Unsupported(_))
    ));
    let back = stl::read(&stl::write(&mesh_data(m), &skipping(true)).unwrap()).unwrap();
    assert!(back.mesh.unwrap().face_attributes().num_columns() == 0);
    let empty = stl::read(&stl::write(&Data::default(), &binary()).unwrap()).unwrap();
    assert_eq!(empty, mesh_data(Mesh::empty()));
    let empty = stl::read(&stl::write(&Data::default(), &ascii()).unwrap()).unwrap();
    assert_eq!(empty, mesh_data(Mesh::empty()));
}

fn solid_names() -> impl Strategy<Value = Option<String>> {
    let s = prop::sample::select(vec![
        "a",
        "rock 2",
        "drusa ore",
        "x",
        "",
        "drusa",
        " pit ",
        "a  b",
        "%null",
        "%empty",
        "%64rusa",
        "50%",
        "tab\t",
        "line\nbreak",
        "ção",
        "endsolid",
    ]);
    prop::option::weighted(0.9, s.prop_map(String::from))
}

proptest! {
    #[test]
    fn binary_write_then_read_is_identity(m in stl_mesh(16), attr in prop::collection::vec(any::<u16>(), 16), tag in any::<bool>()) {
        let mut m = m;
        if tag {
            let mut a: Vec<i64> = attr[..m.n_faces()].iter().map(|&x| i64::from(x)).collect();
            a[0] = a[0].max(1);
            m = m.with_face_column("stl_attribute", Column::i64(a)).unwrap();
        }
        let d = mesh_data(m);
        prop_assert_eq!(stl::read(&stl::write(&d, &binary()).unwrap()).unwrap(), d);
    }

    #[test]
    fn ascii_write_then_read_is_identity(m in stl_mesh(16), names in prop::collection::vec(solid_names(), 16), tag in any::<bool>()) {
        let mut m = m;
        if tag {
            let v = &names[..m.n_faces()];
            let valid = v.iter().map(Option::is_some).collect();
            let col = Column::utf8(v.iter().map(|n| n.clone().unwrap_or_default()).collect())
                .with_valid(valid)
                .unwrap();
            m = m.with_face_column("solid", col).unwrap();
        }
        let d = mesh_data(m);
        prop_assert_eq!(stl::read(&stl::write(&d, &ascii()).unwrap()).unwrap(), d);
    }

    #[test]
    fn arbitrary_bytes_never_panic(b in prop::collection::vec(any::<u8>(), 0..400)) {
        let _ = stl::read(&b);
    }

    #[test]
    fn truncated_files_never_panic(m in stl_mesh(8), asc in any::<bool>(), cut in 0.0..1.0f64) {
        let opts = if asc { ascii() } else { binary() };
        let bytes = stl::write(&mesh_data(m), &opts).unwrap();
        let _ = stl::read(&bytes[..(bytes.len() as f64 * cut) as usize]);
    }
}

fn tri() -> Mesh {
    Mesh::new(
        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        vec![[0, 1, 2]],
    )
    .unwrap()
}

fn ascii_solid_round_trip(col: Column) {
    let d = mesh_data(tri().with_face_column("solid", col).unwrap());
    assert_eq!(stl::read(&stl::write(&d, &ascii()).unwrap()).unwrap(), d);
}

#[test]
fn ascii_solid_named_drusa_round_trips() {
    ascii_solid_round_trip(Column::utf8(vec!["drusa".into()]));
    ascii_solid_round_trip(Column::utf8(vec![String::new()]));
}

#[test]
fn ascii_solid_name_with_outer_spaces_round_trips() {
    ascii_solid_round_trip(Column::utf8(vec![" pit ".into()]));
}

#[test]
fn ascii_null_solid_name_round_trips() {
    ascii_solid_round_trip(
        Column::utf8(vec!["a".into()])
            .with_valid(vec![false])
            .unwrap(),
    );
}
