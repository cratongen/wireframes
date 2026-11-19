use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::time::Duration;
use wireframes_core::fixtures;
use wireframes_query::{Blocks, SolidIndex, Trace};

fn block_model(n: [usize; 3], min: [f64; 3], max: [f64; 3]) -> Blocks {
    let size = [0, 1, 2].map(|k| (max[k] - min[k]) / n[k] as f64);
    let mut centroids = Vec::with_capacity(n[0] * n[1] * n[2]);
    for k in 0..n[2] {
        for j in 0..n[1] {
            for i in 0..n[0] {
                let ijk = [i, j, k];
                centroids.push([0, 1, 2].map(|a| min[a] + (ijk[a] as f64 + 0.5) * size[a]));
            }
        }
    }
    Blocks {
        centroids,
        size,
        rotation: None,
    }
}

fn queries(c: &mut Criterion) {
    let torus = fixtures::torus(500, 100);
    assert_eq!(torus.n_faces(), 100_000);
    let s = SolidIndex::new(&torus).unwrap();
    let points = fixtures::random_points(1, 100_000, [-2.6, -2.6, -0.6], [2.6, 2.6, 0.6]);
    let mut g = c.benchmark_group("torus_100k");
    g.sample_size(10).measurement_time(Duration::from_secs(10));
    g.bench_function("solid_index", |b| {
        b.iter(|| SolidIndex::new(black_box(&torus)).unwrap())
    });
    g.bench_function("locate_100k_points", |b| {
        b.iter(|| s.locate_all(black_box(&points)))
    });
    g.bench_function("signed_distance_100k_points", |b| {
        b.iter(|| s.signed_distance(black_box(&points)))
    });
    let blocks = block_model([100, 100, 50], [-2.5, -2.5, -0.5], [2.5, 2.5, 0.5]);
    let f = s.block_fraction(&blocks).unwrap();
    let block_volume = blocks.size.iter().product::<f64>();
    let total: f64 = f.iter().sum::<f64>() * block_volume;
    let exact = wireframes_ops::volume(&torus).unwrap();
    assert!((total - exact).abs() <= 1e-9 * exact, "{total} vs {exact}");
    g.bench_function("block_fraction_100x100x50", |b| {
        b.iter(|| s.block_fraction(black_box(&blocks)).unwrap())
    });
    let xyz = fixtures::random_points(3, 200, [-2.6, -2.6, -0.6], [2.6, 2.6, 0.6]);
    let depth = (0..xyz.len()).map(|i| i as f64).collect();
    let trace = Trace { xyz, depth };
    g.bench_function("trace_intervals_200", |b| {
        b.iter(|| s.trace_intervals(black_box(&trace)).unwrap())
    });
    g.finish();
}

criterion_group!(benches, queries);
criterion_main!(benches);
