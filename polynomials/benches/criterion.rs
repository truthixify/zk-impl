use criterion::{Criterion, criterion_group, criterion_main};
mod multilinear_benchmarks;
mod univariate_benchmarks;

use multilinear_benchmarks::{
    dense::dense_multilinear_polynomial_benchmarks,
    sparse::sparse_multilinear_polynomial_benchmarks,
};
use univariate_benchmarks::{
    dense::dense_univariate_polynomial_benchmarks, sparse::sparse_univariate_polynomial_benchmarks,
};

criterion_group!(
    name = polynomials;
    config = Criterion::default().sample_size(10).configure_from_args();
    targets = dense_multilinear_polynomial_benchmarks, sparse_multilinear_polynomial_benchmarks, dense_univariate_polynomial_benchmarks, sparse_univariate_polynomial_benchmarks
);
criterion_main!(polynomials);
