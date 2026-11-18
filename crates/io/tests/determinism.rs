use drusa_core::fixtures;
use drusa_io::{Data, DxfEntity, WriteOptions, datamine_csv, dxf, obj, off, ply, stl, vtk, vtu};

fn round_trips() -> Vec<Vec<u8>> {
    let data = Data {
        mesh: Some(fixtures::pit_surface(40, 7)),
        ..Data::default()
    };
    let bin = WriteOptions::default();
    let ascii = WriteOptions {
        ascii: true,
        ..WriteOptions::default()
    };
    let polyface = WriteOptions {
        dxf_entity: DxfEntity::Polyface,
        ..WriteOptions::default()
    };
    let zlib = WriteOptions {
        zlib: true,
        ..WriteOptions::default()
    };
    let (pt, tr) = datamine_csv::write(&data, &bin).unwrap();
    let mut out = vec![
        stl::write(&data, &bin).unwrap(),
        stl::write(&data, &ascii).unwrap(),
        ply::write(&data, &bin).unwrap(),
        ply::write(&data, &ascii).unwrap(),
        obj::write(&data, &bin).unwrap(),
        off::write(&data, &bin).unwrap(),
        dxf::write(&data, &bin).unwrap(),
        dxf::write(&data, &polyface).unwrap(),
        vtk::write(&data, &bin).unwrap(),
        vtu::write(&data, &zlib).unwrap(),
        pt.clone(),
        tr.clone(),
    ];
    let reread = [
        stl::read(&out[0]).unwrap(),
        stl::read(&out[1]).unwrap(),
        ply::read(&out[2]).unwrap(),
        ply::read(&out[3]).unwrap(),
        obj::read(&out[4]).unwrap(),
        off::read(&out[5]).unwrap(),
        dxf::read(&out[6]).unwrap(),
        dxf::read(&out[7]).unwrap(),
        vtk::read(&out[8]).unwrap(),
        vtu::read(&out[9]).unwrap(),
        datamine_csv::read(&pt, &tr).unwrap(),
    ];
    for d in &reread {
        let m = d.mesh.as_ref().unwrap();
        out.push(
            m.vertices()
                .iter()
                .flatten()
                .flat_map(|v| v.to_bits().to_le_bytes())
                .collect(),
        );
        out.push(
            m.triangles()
                .iter()
                .flatten()
                .flat_map(|i| i.to_le_bytes())
                .collect(),
        );
    }
    out
}

fn in_pool(threads: usize) -> Vec<Vec<u8>> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap()
        .install(round_trips)
}

#[test]
fn read_and_write_do_not_depend_on_thread_count() {
    assert_eq!(in_pool(1), in_pool(4));
}
