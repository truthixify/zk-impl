use criterion::{Criterion, criterion_group, criterion_main};
mod circuit_benchmarks;

use circuit_benchmarks::circuit_benchmarks;

criterion_group!(
    name = polynomials;
    config = Criterion::default().sample_size(10).configure_from_args();
    targets = circuit_benchmarks
);
criterion_main!(polynomials);
