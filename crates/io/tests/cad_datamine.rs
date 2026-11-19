mod cad_common;

use cad_common::*;
use proptest::prelude::*;
use wireframes_core::{Column, ColumnData, Error, Mesh, Table};
use wireframes_io::{Data, WriteOptions, datamine_csv};

fn utf8(v: &[&str]) -> Column {
    Column::utf8(v.iter().map(|s| s.to_string()).collect())
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    #[test]
    fn meshes_with_attributes_round_trip(
        verts in prop::collection::vec([coord(), coord(), coord()], 3..12),
        faces in prop::collection::vec([0u32..3, 0u32..3, 0u32..3], 1..10),
        grade in prop::collection::vec(prop::option::of(coord()), 12),
    ) {
        prop_assume!(grade.iter().all(|g| *g != Some(-1.0e30)));
        let n = verts.len();
        let mut m = Mesh::new(verts, faces.clone()).unwrap();
        let g: Vec<Option<f64>> = grade[..n].to_vec();
        let valid: Vec<bool> = g.iter().map(Option::is_some).collect();
        let mut col = Column::f64(g.iter().map(|x| x.unwrap_or(f64::NAN)).collect());
        if valid.iter().any(|v| !v) {
            col = col.with_valid(valid).unwrap();
        }
        let vt = Table::new(n).with_column("AU", col).unwrap();
        let rocks: Vec<String> = (0..faces.len()).map(|i| format!("r,{i} \"q\"")).collect();
        let ft = Table::new(faces.len())
            .with_column("ROCK", Column::utf8(rocks))
            .unwrap()
            .with_column("SID", Column::f64((0..faces.len()).map(|i| i as f64).collect()))
            .unwrap();
        m = m.with_attributes(vt, ft).unwrap();
        let (pt, tr) = datamine_csv::write(&mesh_data(m.clone()), &WriteOptions::default()).unwrap();
        let back = datamine_csv::read(&pt, &tr).unwrap();
        prop_assert_eq!(back.mesh, Some(m));
    }
}

#[test]
fn reads_semicolon_pairs_with_gapped_pids_and_type_rows() {
    let pt = "\u{feff}XP;YP;ZP;PID;AU\nN;N;N;N;N\n1,5;2;3;10;0,25\n4;5;6;20;-1.0E30\n7;8;9;5;\n";
    let tr = "TRISOL;PID1;PID2;PID3;LITH\n1;10;20;5;ox\n";
    let m = datamine_csv::read(pt.as_bytes(), tr.as_bytes())
        .unwrap()
        .mesh
        .unwrap();
    assert_eq!(
        m.vertices(),
        &[[1.5, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
    );
    assert_eq!(m.triangles(), &[[0, 1, 2]]);
    let au = m.vertex_attributes().column("AU").unwrap();
    assert_eq!(au.valid, Some(vec![true, false, false]));
    assert_eq!(
        m.face_attributes().column("LITH").unwrap().data,
        ColumnData::Utf8(vec!["ox".into()])
    );
    assert_eq!(m.face_attributes().names(), &["LITH"]);
}

#[test]
fn header_case_and_quotes_are_ignored() {
    let pt = " \"xp\" , yp ,Zp,pid\n0,0,0,1\n1,0,0,2\n0,1,0,3\n";
    let tr = "triangle,PID1,PID2,PID3\n1,1,2,3\n";
    let m = datamine_csv::read(pt.as_bytes(), tr.as_bytes())
        .unwrap()
        .mesh
        .unwrap();
    assert_eq!(m.n_faces(), 1);
}

#[test]
fn bad_pairs_error() {
    let pt = "XP,YP,ZP,PID\n0,0,0,1\n1,0,0,2\n0,1,0,3\n";
    let cases: &[(&str, &str)] = &[
        (pt, "TRIANGLE,PID1,PID2,PID3\n1,1,2,4\n"),
        ("XP,YP,ZP,PID\n0,0,0,1\n1,0,0,1\n", "PID1,PID2,PID3\n"),
        ("XP,YP,ZP,PID\n0,0,0,1.5\n", "PID1,PID2,PID3\n"),
        ("XP,YP,PID\n0,0,1\n", "PID1,PID2,PID3\n"),
        ("XP,YP,ZP,PID\n0,0,0,1\n0,x,0,2\n", "PID1,PID2,PID3\n"),
        ("XP,YP,ZP,PID\n0,0,0\n", "PID1,PID2,PID3\n"),
        (pt, "TRIANGLE,PID1,PID2\n1,1,2\n"),
        (pt, ""),
        (pt, "PID1,PID2,PID3\n1,2,-1.0E30\n"),
    ];
    for (p, t) in cases {
        assert!(
            datamine_csv::read(p.as_bytes(), t.as_bytes()).is_err(),
            "{p:?} {t:?}"
        );
    }
}

#[test]
fn writes_one_based_pids_and_refuses_clashes() {
    let m = Mesh::new(
        vec![[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        vec![[0, 1, 2]],
    )
    .unwrap();
    let (pt, tr) = datamine_csv::write(&mesh_data(m.clone()), &Default::default()).unwrap();
    assert_eq!(
        String::from_utf8(pt).unwrap(),
        "XP,YP,ZP,PID\n0.0,0.0,0.0,1\n1.0,0.0,0.0,2\n0.0,1.0,0.0,3\n"
    );
    assert_eq!(
        String::from_utf8(tr).unwrap(),
        "TRIANGLE,PID1,PID2,PID3\n1,1,2,3\n"
    );
    let bad = m
        .clone()
        .with_vertex_column("pid", Column::f64(vec![0.0; 3]))
        .unwrap();
    assert!(matches!(
        datamine_csv::write(&mesh_data(bad), &Default::default()),
        Err(Error::InvalidInput(_))
    ));
    let bad = m.with_face_column("TRI", utf8(&["a"])).unwrap();
    assert!(datamine_csv::write(&mesh_data(bad), &Default::default()).is_err());
    assert!(matches!(
        datamine_csv::write(&Data::default(), &Default::default()),
        Err(Error::Unsupported(_))
    ));
}

#[test]
fn damaged_pairs_never_panic() {
    let m = Mesh::new(
        vec![[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 1.0]],
        vec![[0, 1, 2], [1, 3, 2]],
    )
    .unwrap()
    .with_face_column("ROCK", utf8(&["a", "b"]))
    .unwrap();
    let (pt, tr) = datamine_csv::write(&mesh_data(m), &Default::default()).unwrap();
    never_panics(&pt, |b| {
        let _ = datamine_csv::read(b, &tr);
    });
    never_panics(&tr, |b| {
        let _ = datamine_csv::read(&pt, b);
    });
}

fn tri() -> Mesh {
    Mesh::new(
        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        vec![[0, 1, 2]],
    )
    .unwrap()
}

fn round(m: &Mesh) -> Mesh {
    let data = Data {
        mesh: Some(m.clone()),
        ..Data::default()
    };
    let (pt, tr) = datamine_csv::write(&data, &WriteOptions::default()).unwrap();
    datamine_csv::read(&pt, &tr).unwrap().mesh.unwrap()
}

#[test]
fn text_with_line_breaks_and_quotes_round_trips() {
    let m = tri()
        .with_face_column("note", utf8(&["a\nb\r\n\"c\", d;"]))
        .unwrap();
    assert_eq!(round(&m), m);
}

#[test]
fn column_name_with_semicolon_round_trips() {
    let m = tri()
        .with_vertex_column("a;b", Column::f64(vec![1.5, 2.5, 3.5]))
        .unwrap();
    assert_eq!(round(&m), m);
}

#[test]
fn unclosed_quote_errors() {
    let pt = b"XP,YP,ZP,PID,N\n0,0,0,1,\"x\n1,0,0,2,y\n0,1,0,3,z\n";
    let tr = b"TRIANGLE,PID1,PID2,PID3\n1,1,2,3\n";
    assert!(matches!(
        datamine_csv::read(pt, tr),
        Err(Error::Format { .. })
    ));
}

#[test]
fn every_triangle_id_alias_reads_and_none_becomes_a_column() {
    let pt = "XP,YP,ZP,PID\n0,0,0,1\n1,0,0,2\n0,1,0,3\n";
    for id in ["TRIANGLE", "TRIANG", "TRI", "TRISOL", "trisol"] {
        let tr = format!("{id},PID1,PID2,PID3\n7,1,2,3\n");
        let m = datamine_csv::read(pt.as_bytes(), tr.as_bytes())
            .unwrap()
            .mesh
            .unwrap();
        assert_eq!(m.triangles(), &[[0, 1, 2]], "{id}");
        assert!(m.face_attributes().names().is_empty(), "{id}");
        let (_, back) = datamine_csv::write(&mesh_data(m), &Default::default()).unwrap();
        assert!(back.starts_with(b"TRIANGLE,PID1,PID2,PID3\n"), "{id}");
    }
}
