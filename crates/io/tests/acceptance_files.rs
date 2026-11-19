use std::path::PathBuf;

use wireframes_core::Mesh;
use wireframes_io::{Data, WriteOptions, datamine_csv, dxf, off, read, vtu};

fn scratch(test: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "wireframes-io-accept-{}-{test}",
        std::process::id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn tetra() -> Mesh {
    Mesh::new(
        vec![
            [700_000.0, 7_500_000.0, 200.0],
            [700_010.5, 7_500_000.0, 200.0],
            [700_000.0, 7_500_010.5, 200.0],
            [700_000.0, 7_500_000.0, 210.5],
        ],
        vec![[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]],
    )
    .unwrap()
}

fn mesh_data(m: Mesh) -> Data {
    Data {
        mesh: Some(m),
        ..Data::default()
    }
}

#[test]
#[ignore = "acceptance: a dxf cut inside ENTITIES reads a truncated coordinate as a vertex"]
fn dxf_cut_inside_a_coordinate_is_an_error() {
    let bytes = dxf::write(&mesh_data(tetra()), &WriteOptions::default()).unwrap();
    let text = String::from_utf8(bytes).unwrap();
    let at = text.rfind("700010.5").unwrap() + 4;
    let cut = &text.as_bytes()[..at];
    let read = dxf::read(cut);
    assert!(
        read.is_err(),
        "a dxf without ENDSEC and EOF read as {:?}",
        read.map(|d| d.mesh.map(|m| m.vertices().to_vec()))
    );
}

#[test]
#[ignore = "acceptance: TRISOL is taken as the triangle id over TRIANGLE when it comes first, and its values are lost"]
fn datamine_trisol_beside_triangle_is_kept() {
    let pt = "XP,YP,ZP,PID\n0,0,0,1\n1,0,0,2\n0,1,0,3\n";
    let tr = "TRISOL,TRIANGLE,PID1,PID2,PID3\n7,1,1,2,3\n";
    let m = datamine_csv::read(pt.as_bytes(), tr.as_bytes())
        .unwrap()
        .mesh
        .unwrap();
    let names = m.face_attributes().names();
    assert!(
        names.iter().any(|n| n == "TRISOL"),
        "face columns are {names:?}; the TRISOL surface ids are gone"
    );
}

#[test]
fn datamine_triangle_before_trisol_keeps_trisol() {
    let pt = "XP,YP,ZP,PID\n0,0,0,1\n1,0,0,2\n0,1,0,3\n";
    let tr = "TRIANGLE,TRISOL,PID1,PID2,PID3\n1,7,1,2,3\n";
    let m = datamine_csv::read(pt.as_bytes(), tr.as_bytes())
        .unwrap()
        .mesh
        .unwrap();
    assert_eq!(m.face_attributes().names(), &["TRISOL"]);
}

#[test]
#[ignore = "acceptance: an error in _pt.csv is labelled with the _tr.csv path"]
fn datamine_error_names_the_file_that_holds_it() {
    let dir = scratch("datamine_error_names_the_file");
    std::fs::write(dir.join("pit_pt.csv"), "XP,YP,ZP,PID\n0,0,0,1\n1,0,0\n").unwrap();
    std::fs::write(dir.join("pit_tr.csv"), "TRIANGLE,PID1,PID2,PID3\n1,1,2,3\n").unwrap();
    let msg = read(dir.join("pit_tr.csv")).unwrap_err().to_string();
    assert!(
        msg.contains("pit_pt.csv:3"),
        "the bad row is line 3 of pit_pt.csv, but the error says {msg:?}"
    );
}

#[test]
#[ignore = "acceptance: a truncated off names neither the line nor the byte offset"]
fn truncated_off_names_the_line() {
    let off = "OFF\n4 4 0\n0 0 0\n1 0 0\n";
    let msg = off::read(off.as_bytes()).unwrap_err().to_string();
    assert!(msg.contains("line"), "error without a location: {msg:?}");
}

#[test]
#[ignore = "acceptance: a truncated appended vtu names neither the array offset nor the byte"]
fn truncated_vtu_names_the_byte_offset() {
    let bytes = vtu::write(&mesh_data(tetra()), &WriteOptions::default()).unwrap();
    let cut = &bytes[..bytes.len() - 80];
    let msg = vtu::read(cut).unwrap_err().to_string();
    assert!(msg.contains("byte"), "error without a location: {msg:?}");
}
