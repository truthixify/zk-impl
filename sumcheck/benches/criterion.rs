use criterion::{Criterion, criterion_group, criterion_main};

mod sumcheck_benchmarks;

use sumcheck_benchmarks::sumcheck_benchmarks;

criterion_group!(
    name = sumcheck;
    config = Criterion::default().sample_size(10).configure_from_args();
    targets = sumcheck_benchmarks
);
criterion_main!(sumcheck);
