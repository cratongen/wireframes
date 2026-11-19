use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::time::Duration;
use wireframes_core::{Result, fixtures};
use wireframes_io::{Data, WriteOptions, datamine_csv, dxf, obj, off, ply, stl, vtk, vtu};

type Writer = fn(&Data, &WriteOptions) -> Result<Vec<u8>>;
type Reader = fn(&[u8]) -> Result<Data>;

fn read_write(c: &mut Criterion) {
    let data = Data {
        mesh: Some(fixtures::torus(500, 100)),
        ..Data::default()
    };
    assert_eq!(data.mesh.as_ref().unwrap().n_faces(), 100_000);
    let opts = WriteOptions::default();
    let formats: [(&str, Writer, Reader); 7] = [
        ("obj", obj::write, obj::read),
        ("ply", ply::write, ply::read),
        ("off", off::write, off::read),
        ("stl", stl::write, stl::read),
        ("dxf", dxf::write, dxf::read),
        ("vtk", vtk::write, vtk::read),
        ("vtu", vtu::write, vtu::read),
    ];
    let mut g = c.benchmark_group("torus_100k");
    g.sample_size(10).measurement_time(Duration::from_secs(5));
    for (name, write, read) in formats {
        let bytes = write(&data, &opts).unwrap();
        assert_eq!(read(&bytes).unwrap().mesh.unwrap().n_faces(), 100_000);
        g.bench_function(format!("{name}_write"), |b| {
            b.iter(|| write(black_box(&data), &opts).unwrap())
        });
        g.bench_function(format!("{name}_read"), |b| {
            b.iter(|| read(black_box(&bytes)).unwrap())
        });
    }
    let (pt, tr) = datamine_csv::write(&data, &opts).unwrap();
    assert_eq!(
        datamine_csv::read(&pt, &tr)
            .unwrap()
            .mesh
            .unwrap()
            .n_faces(),
        100_000
    );
    g.bench_function("datamine_csv_write", |b| {
        b.iter(|| datamine_csv::write(black_box(&data), &opts).unwrap())
    });
    g.bench_function("datamine_csv_read", |b| {
        b.iter(|| datamine_csv::read(black_box(&pt), black_box(&tr)).unwrap())
    });
    g.finish();
}

criterion_group!(benches, read_write);
criterion_main!(benches);
