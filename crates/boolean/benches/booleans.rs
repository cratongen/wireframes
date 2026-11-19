use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::time::Duration;
use wireframes_boolean::{BooleanOptions, Op, arrangement_size, boolean, slice_by_plane};
use wireframes_core::{Mesh, fixtures};

const SIZES: [(&str, u32, u32); 3] = [("1k", 25, 20), ("10k", 100, 50), ("100k", 500, 100)];
const SHIFT: [f64; 3] = [0.37, 0.29, 0.11];

fn pair(major: u32, minor: u32) -> (Mesh, Mesh) {
    let a = fixtures::torus(major, minor);
    let b = fixtures::translated(&a, SHIFT);
    (a, b)
}

fn booleans(c: &mut Criterion) {
    let opts = BooleanOptions::default();
    let mut g = c.benchmark_group("torus_pair");
    g.sample_size(10).measurement_time(Duration::from_secs(10));
    for (name, major, minor) in SIZES {
        let (a, b) = pair(major, minor);
        assert_eq!(a.n_faces(), 2 * (major * minor) as usize);
        assert!(arrangement_size(&a, &b).unwrap() > a.n_faces() + b.n_faces());
        g.bench_with_input(
            BenchmarkId::new("arrangement", name),
            &(&a, &b),
            |bn, (a, b)| bn.iter(|| arrangement_size(black_box(a), black_box(b)).unwrap()),
        );
        for (op_name, op) in [("union", Op::Union), ("difference", Op::Difference)] {
            assert!(boolean(&a, &b, op, &opts).unwrap().mesh.n_faces() > 0);
            g.bench_with_input(BenchmarkId::new(op_name, name), &(&a, &b), |bn, (a, b)| {
                bn.iter(|| boolean(black_box(a), black_box(b), op, &opts).unwrap())
            });
        }
        g.bench_with_input(BenchmarkId::new("validate", name), &a, |bn, a| {
            bn.iter(|| black_box(a).validate(0.0).unwrap())
        });
    }
    let (s, _) = pair(500, 100);
    let offsets: Vec<f64> = (0..50).map(|k| -0.49 + k as f64 * 0.02).collect();
    assert_eq!(
        slice_by_plane(&s, [0.0; 3], [0.0, 0.0, 1.0], &offsets)
            .unwrap()
            .n_parts(),
        100
    );
    g.bench_function("slice_50_planes/100k", |bn| {
        bn.iter(|| slice_by_plane(black_box(&s), [0.0; 3], [0.0, 0.0, 1.0], &offsets).unwrap())
    });
    g.finish();
}

criterion_group!(benches, booleans);
criterion_main!(benches);
