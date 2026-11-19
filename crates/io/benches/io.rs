use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use wireframes_core::fixtures;
use wireframes_io::{Data, WriteOptions, ply, stl};

fn binary_read_write(c: &mut Criterion) {
    let data = Data {
        mesh: Some(fixtures::pit_surface(501, 1)),
        ..Data::default()
    };
    assert_eq!(data.mesh.as_ref().unwrap().triangles().len(), 500_000);
    let opts = WriteOptions::default();
    let stl_bytes = stl::write(&data, &opts).unwrap();
    let ply_bytes = ply::write(&data, &opts).unwrap();
    let mut g = c.benchmark_group("pit_500k");
    g.sample_size(20);
    g.bench_function("stl_write", |b| {
        b.iter(|| stl::write(black_box(&data), &opts).unwrap())
    });
    g.bench_function("stl_read", |b| {
        b.iter(|| stl::read(black_box(&stl_bytes)).unwrap())
    });
    g.bench_function("ply_write", |b| {
        b.iter(|| ply::write(black_box(&data), &opts).unwrap())
    });
    g.bench_function("ply_read", |b| {
        b.iter(|| ply::read(black_box(&ply_bytes)).unwrap())
    });
    g.finish();
}

criterion_group!(benches, binary_read_write);
criterion_main!(benches);
