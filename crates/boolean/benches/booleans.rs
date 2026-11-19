use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use wireframes_boolean::{BooleanOptions, Op, boolean, slice_by_plane};
use wireframes_core::fixtures;

fn booleans(c: &mut Criterion) {
    let a = fixtures::sphere(4);
    let b = fixtures::translated(&a, [0.4, 0.3, 0.2]);
    let opts = BooleanOptions::default();
    for (name, op) in [
        ("union_sphere_5k", Op::Union),
        ("difference_sphere_5k", Op::Difference),
    ] {
        c.bench_function(name, |bn| {
            bn.iter(|| boolean(black_box(&a), black_box(&b), op, &opts).unwrap())
        });
    }
    let s = fixtures::sphere(5);
    let offsets: Vec<f64> = (0..50).map(|k| -0.98 + k as f64 * 0.04).collect();
    c.bench_function("slice_by_plane_sphere_20k_50", |bn| {
        bn.iter(|| slice_by_plane(black_box(&s), [0.0; 3], [0.0, 0.0, 1.0], &offsets).unwrap())
    });
    c.bench_function("validate_sphere_20k", |bn| {
        bn.iter(|| black_box(&s).validate(0.0).unwrap())
    });
}

criterion_group!(benches, booleans);
criterion_main!(benches);
