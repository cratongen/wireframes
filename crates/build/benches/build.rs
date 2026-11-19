use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::time::Duration;
use wireframes_build::{DtmOptions, dtm};
use wireframes_core::fixtures;

fn build(c: &mut Criterion) {
    let points: Vec<[f64; 3]> =
        fixtures::random_points(4, 100_000, [0.0; 3], [1000.0, 1000.0, 1.0])
            .into_iter()
            .map(|[x, y, z]| [x, y, 100.0 + 0.05 * x - 0.02 * y + z])
            .collect();
    let opts = DtmOptions::default();
    assert_eq!(dtm(&points, &opts).unwrap().n_vertices(), 100_000);
    let mut g = c.benchmark_group("points_100k");
    g.sample_size(10).measurement_time(Duration::from_secs(10));
    g.bench_function("dtm", |b| {
        b.iter(|| dtm(black_box(&points), &opts).unwrap())
    });
    g.finish();
}

criterion_group!(benches, build);
criterion_main!(benches);
