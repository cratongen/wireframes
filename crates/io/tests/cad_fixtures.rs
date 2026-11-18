use drusa_core::ColumnData;
use drusa_io::{Data, DxfEntity, WriteOptions, datamine_csv, dxf, vtk, vtu};

fn fixture(name: &str) -> Vec<u8> {
    let path = format!("{}/tests/fixtures/{name}", env!("CARGO_MANIFEST_DIR"));
    std::fs::read(&path).unwrap_or_else(|e| panic!("{path}: {e}"))
}

fn content(d: &Data) -> Data {
    Data {
        warnings: Vec::new(),
        ..d.clone()
    }
}

fn opts() -> [WriteOptions; 2] {
    [
        WriteOptions::default(),
        WriteOptions {
            ascii: true,
            zlib: true,
            dxf_entity: DxfEntity::Polyface,
            ..WriteOptions::default()
        },
    ]
}

fn read_write_read(
    name: &str,
    read: impl Fn(&[u8]) -> drusa_core::Result<Data>,
    write: impl Fn(&Data, &WriteOptions) -> drusa_core::Result<Vec<u8>>,
) -> Data {
    let first = read(&fixture(name)).unwrap();
    for o in opts() {
        let again = read(&write(&first, &o).unwrap()).unwrap();
        assert_eq!(content(&again), content(&first), "{name} with {o:?}");
    }
    first
}

fn through_legacy_vtk(d: &Data) {
    for o in opts() {
        let again = vtk::read(&vtk::write(d, &o).unwrap()).unwrap();
        assert_eq!(&again, d, "legacy vtk with {o:?}");
    }
}

#[test]
fn dxf_polyface_with_negative_indices() {
    let d = read_write_read("polyface_negative.dxf", dxf::read, dxf::write);
    assert_eq!(d.mesh.unwrap().triangles(), &[[0, 1, 2], [0, 2, 3]]);
}

#[test]
fn dxf_cp1252_layer_names() {
    let d = read_write_read("cp1252_layers.dxf", dxf::read, dxf::write);
    assert_eq!(d.length_unit.as_deref(), Some("ft"));
    assert_eq!(
        d.points.unwrap().attributes().column("layer").unwrap().data,
        ColumnData::Utf8(vec!["MIN\u{C9}".into()])
    );
    assert_eq!(
        d.polylines
            .unwrap()
            .attributes()
            .column("layer")
            .unwrap()
            .data,
        ColumnData::Utf8(vec!["Caf\u{e9} \u{20AC}".into()])
    );
}

#[test]
fn datamine_gapped_semicolon_pair() {
    let (pt, tr) = (
        fixture("gapped_semicolon_pt.csv"),
        fixture("gapped_semicolon_tr.csv"),
    );
    let first = datamine_csv::read(&pt, &tr).unwrap();
    let (pt2, tr2) = datamine_csv::write(&first, &WriteOptions::default()).unwrap();
    assert_eq!(datamine_csv::read(&pt2, &tr2).unwrap(), first);
    assert_eq!(first.mesh.unwrap().n_vertices(), 3);
}

#[test]
fn vtp_with_polys_strips_and_cell_data() {
    let d = read_write_read("polys_strips.vtp", vtu::read, vtu::write_vtp);
    through_legacy_vtk(&d);
    assert_eq!(d.mesh.unwrap().n_faces(), 3);
    assert_eq!(d.polylines.unwrap().len(), 1);
}

#[test]
fn multi_piece_vtu() {
    let d = read_write_read("two_pieces.vtu", vtu::read, vtu::write);
    through_legacy_vtk(&d);
    assert_eq!(d.grid.unwrap().offsets, vec![0, 4, 8]);
}

#[test]
fn polyhedron_vtu() {
    let first = vtu::read(&fixture("polyhedron.vtu")).unwrap();
    for o in opts() {
        let again = vtu::read(&vtu::write(&first, &o).unwrap()).unwrap();
        assert_eq!(again, first);
    }
    assert!(first.grid.unwrap().faces.is_some());
}
