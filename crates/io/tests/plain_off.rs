mod plain_support;

use drusa_core::{Column, Error, Mesh};
use drusa_io::{Data, off};
use plain_support::*;
use proptest::prelude::*;

#[test]
fn coloured_off_with_comments_and_counts_on_next_line() {
    let d = off::read(&fixture("colored.off")).unwrap();
    let m = d.mesh.as_ref().unwrap();
    assert_eq!(m.n_vertices(), 5);
    assert_eq!(m.triangles(), &[[0, 1, 2], [0, 2, 3], [0, 4, 1]]);
    assert_eq!(m.face_attributes().names(), &["polygon"]);
    assert_eq!(
        m.face_attributes().column("polygon").unwrap(),
        &Column::i64(vec![0, 0, 1])
    );
    assert_eq!(m.vertex_attributes().num_columns(), 0);
}

#[test]
fn homogeneous_divides_by_w() {
    let d = off::read(&fixture("homogeneous.off")).unwrap();
    assert_eq!(
        d.mesh.unwrap().vertices(),
        &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0]]
    );
}

#[test]
fn keyword_variants() {
    let d = off::read(b"nOFF 3\n3 1 0\n0 0 0\n1 0 0\n0 1 0\n3 0 1 2\n").unwrap();
    assert_eq!(d.mesh.unwrap().n_faces(), 1);
    let d = off::read(b"STCNOFF 3 1 0\n0 0 0 0 0 1 1 1 1 0.5 0.5\n1 0 0 0 0 1 1 1 1 0.5 0.5\n0 1 0 0 0 1 1 1 1 0.5 0.5\n3 0 1 2\n").unwrap();
    assert_eq!(d.mesh.unwrap().vertices()[1], [1.0, 0.0, 0.0]);
}

#[test]
fn fixtures_read_write_read_identically() {
    for name in ["colored.off", "homogeneous.off"] {
        let d = off::read(&fixture(name)).unwrap();
        let again = off::read(&off::write(&d, &ascii()).unwrap()).unwrap();
        assert_eq!(again, d, "{name}");
    }
}

#[test]
fn polygons_write_back_as_polygons() {
    let d = off::read(&fixture("colored.off")).unwrap();
    let text = String::from_utf8(off::write(&d, &ascii()).unwrap()).unwrap();
    assert!(text.ends_with("4 0 1 2 3\n3 0 4 1\n"), "{text}");
}

#[test]
fn malformed_inputs() {
    let format: &[&[u8]] = &[
        b"",
        b"# only a comment\n",
        b"PFF\n3 1 0\n",
        b"XOFF\n",
        b"OFF\n3\n",
        b"OFF\n3 1 0\n0 0 0\n1 0 0\n",
        b"OFF\n3 1 0\n0 0 0\n1 0 0\n0 1\n3 0 1 2\n",
        b"OFF\n3 1 0\n0 0 0\n1 0 0\n0 1 0\n3 0 1 3\n",
        b"OFF\n3 1 0\n0 0 0\n1 0 0\n0 1 0\n3 0 1\n",
        b"OFF\n3 1 0\n0 0 0\n1 0 0\n0 1 0\n",
        b"OFF\n-3 1 0\n",
        b"OFF\n1 0 0\nnan 0 0\n",
        b"4OFF\n1 0 0\n1 1 1 0\n",
    ];
    for c in format {
        let r = off::read(c);
        assert!(
            matches!(r, Err(Error::Format { .. })),
            "{:?}: {r:?}",
            String::from_utf8_lossy(c)
        );
    }
    let unsupported: &[&[u8]] = &[
        b"OFF BINARY\n",
        b"OFF\nBINARY\n",
        b"nOFF 4\n1 0 0\n0 0 0 0\n",
    ];
    for c in unsupported {
        assert!(matches!(off::read(c), Err(Error::Unsupported(_))));
    }
}

#[test]
fn writer_refuses_what_off_cannot_hold() {
    let m = cube()
        .with_vertex_column("g", Column::f64(vec![0.0; 8]))
        .unwrap();
    assert!(matches!(
        off::write(&mesh_data(m.clone()), &ascii()),
        Err(Error::Unsupported(_))
    ));
    let back = off::read(&off::write(&mesh_data(m), &skipping(true)).unwrap()).unwrap();
    assert_eq!(back, mesh_data(cube()));
    let d = Data {
        points: Some(drusa_core::PointSet::new(vec![[0.0; 3]]).unwrap()),
        ..Data::default()
    };
    assert!(matches!(
        off::write(&d, &ascii()),
        Err(Error::Unsupported(_))
    ));
    assert_eq!(
        off::read(&off::write(&Data::default(), &ascii()).unwrap()).unwrap(),
        mesh_data(Mesh::empty())
    );
}

proptest! {
    #[test]
    fn write_then_read_is_identity((m, _) in polygon_mesh(12, 8)) {
        let d = mesh_data(m);
        let back = off::read(&off::write(&d, &ascii()).unwrap()).unwrap();
        prop_assert_eq!(back, d);
    }

    #[test]
    fn arbitrary_bytes_never_panic(b in prop::collection::vec(any::<u8>(), 0..400)) {
        let _ = off::read(&b);
    }

    #[test]
    fn truncated_files_never_panic((m, _) in polygon_mesh(8, 6), cut in 0.0..1.0f64) {
        let bytes = off::write(&mesh_data(m), &ascii()).unwrap();
        let _ = off::read(&bytes[..(bytes.len() as f64 * cut) as usize]);
    }
}
