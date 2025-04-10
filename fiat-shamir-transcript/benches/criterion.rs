use criterion::{Criterion, criterion_group, criterion_main};
mod fiat_shamir_transcript_benchmarks;

use fiat_shamir_transcript_benchmarks::fiat_shamir_transcript_benchmarks;

criterion_group!(
    name = fiat_shamir_transcript;
    config = Criterion::default().sample_size(100).configure_from_args();
    targets = fiat_shamir_transcript_benchmarks
);
criterion_main!(fiat_shamir_transcript);
