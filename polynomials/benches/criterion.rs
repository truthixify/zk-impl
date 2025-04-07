use criterion::{Criterion, criterion_group, criterion_main};
mod univariate_benchmarks;
mod multilinear_benchmarks;

use univariate_benchmarks::{
    dense::dense_univariate_polynomial_benchmarks, sparse::sparse_univariate_polynomial_benchmarks,
};
use multilinear_benchmarks::{
    sparse::sparse_multilinear_polynomial_benchmarks,
};

criterion_group!(
    name = polynomials;
    config = Criterion::default().sample_size(100).configure_from_args();
    targets = sparse_multilinear_polynomial_benchmarks //dense_univariate_polynomial_benchmarks, sparse_univariate_polynomial_benchmarks
);
criterion_main!(polynomials);
