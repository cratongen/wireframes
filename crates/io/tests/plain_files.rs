mod plain_support;

use std::path::{Path, PathBuf};

use plain_support::*;
use wireframes_core::{Error, PointSet};
use wireframes_io::{
    Data, Format, ReadOptions, WriteOptions, canonicalize, detect, read, read_mesh, read_points,
    read_with, write, write_mesh, write_with,
};

fn scratch(test: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("wireframes-io-{}-{test}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

#[test]
fn detect_by_extension_ignoring_case() {
    let cases = [
        ("a.obj", Format::Obj),
        ("A.PLY", Format::Ply),
        ("x/y.Off", Format::Off),
        ("m.stl", Format::Stl),
        ("m.dxf", Format::Dxf),
        ("m.vtk", Format::Vtk),
        ("m.vtu", Format::Vtu),
        ("m.VTP", Format::Vtp),
        ("pit_tr.csv", Format::DatamineCsv),
        ("PIT_PT.CSV", Format::DatamineCsv),
    ];
    for (p, f) in cases {
        assert_eq!(detect(Path::new(p)).unwrap(), f, "{p}");
    }
    for p in ["a.csv", "a.txt", "noext", "a.obj.bak"] {
        assert!(
            matches!(detect(Path::new(p)), Err(Error::Unsupported(_))),
            "{p}"
        );
    }
}

#[test]
fn plain_formats_round_trip_through_files() {
    let dir = scratch("round");
    for name in ["c.obj", "c.ply", "c.off", "c.stl", "C.STL"] {
        let p = dir.join(name);
        let want = canonicalize(&cube(), detect(&p).unwrap());
        write(&p, &mesh_data(cube())).unwrap();
        assert_eq!(read(&p).unwrap(), mesh_data(want.clone()), "{name}");
        write_mesh(&p, &cube(), &ascii()).unwrap();
        assert_eq!(read_mesh(&p).unwrap(), want, "{name}");
    }
    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn read_options_label_every_object() {
    let dir = scratch("labels");
    let p = dir.join("m.obj");
    std::fs::write(&p, "v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\nv 5 5 5\np 4\n").unwrap();
    let opts = ReadOptions {
        crs: Some("EPSG:31983".into()),
        length_unit: Some("m".into()),
    };
    let d = read_with(&p, &opts).unwrap();
    let m = d.mesh.unwrap();
    assert_eq!(m.crs.as_deref(), Some("EPSG:31983"));
    assert_eq!(m.length_unit.as_deref(), Some("m"));
    let ps = d.points.unwrap();
    assert_eq!(
        ps,
        PointSet::new(vec![[5.0, 5.0, 5.0]])
            .unwrap()
            .with_crs(Some("EPSG:31983".into()))
            .with_length_unit(Some("m".into()))
    );
    assert_eq!(read(&p).unwrap().mesh.unwrap().crs, None);
    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn missing_parts_and_files_are_errors() {
    let dir = scratch("missing");
    let p = dir.join("pts.obj");
    std::fs::write(&p, "v 1 2 3\np 1\n").unwrap();
    assert!(
        matches!(read_mesh(&p), Err(Error::InvalidInput(m)) if m.contains("file holds no mesh"))
    );
    assert_eq!(read_points(&p).unwrap().len(), 1);
    assert!(matches!(read(dir.join("absent.ply")), Err(Error::Io(_))));
    let bad = dir.join("bad.off");
    std::fs::write(&bad, "OFF\n1 0 0\n0 0\n").unwrap();
    match read(&bad) {
        Err(Error::Format { path, .. }) => assert_eq!(path.as_deref(), Some(bad.as_path())),
        r => panic!("{r:?}"),
    }
    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn datamine_pair_is_found_from_either_name() {
    let dir = scratch("pair");
    let d = mesh_data(cube());
    write_with(dir.join("Pit_TR.csv"), &d, &WriteOptions::default()).unwrap();
    assert!(dir.join("Pit_PT.csv").exists());
    let a = read(dir.join("Pit_TR.csv")).unwrap();
    let b = read(dir.join("Pit_PT.csv")).unwrap();
    assert_eq!(a, b);
    assert_eq!(a.mesh.unwrap().n_faces(), 12);
    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn unsupported_parts_are_refused_before_writing() {
    let dir = scratch("refuse");
    let p = dir.join("x.stl");
    let d = Data {
        points: Some(PointSet::new(vec![[0.0; 3]]).unwrap()),
        ..Data::default()
    };
    assert!(matches!(write(&p, &d), Err(Error::Unsupported(_))));
    assert!(!p.exists());
    std::fs::remove_dir_all(dir).unwrap();
}
