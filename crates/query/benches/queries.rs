use criterion::{Criterion, criterion_group, criterion_main};
use drusa_core::fixtures;
use drusa_query::{Blocks, SolidIndex, Trace};
use std::hint::black_box;

fn queries(c: &mut Criterion) {
    let s = SolidIndex::new(&fixtures::torus(64, 32)).unwrap();
    let points = fixtures::random_points(1, 10_000, [-1.5; 3], [1.5; 3]);
    c.bench_function("locate_10k", |b| {
        b.iter(|| s.locate_all(black_box(&points)))
    });
    c.bench_function("signed_distance_10k", |b| {
        b.iter(|| s.signed_distance(black_box(&points)))
    });
    let blocks = Blocks {
        centroids: fixtures::random_points(2, 2_000, [-1.5; 3], [1.5; 3]),
        size: [0.1; 3],
        rotation: None,
    };
    c.bench_function("block_fraction_2k", |b| {
        b.iter(|| s.block_fraction(black_box(&blocks)).unwrap())
    });
    let xyz = fixtures::random_points(3, 200, [-1.5; 3], [1.5; 3]);
    let depth = (0..xyz.len()).map(|i| i as f64).collect();
    let trace = Trace { xyz, depth };
    c.bench_function("trace_intervals_200", |b| {
        b.iter(|| s.trace_intervals(black_box(&trace)).unwrap())
    });
}

criterion_group!(benches, queries);
criterion_main!(benches);
