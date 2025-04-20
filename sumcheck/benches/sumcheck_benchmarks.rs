use ark_bls12_381::Fq;
use ark_ff::UniformRand;
use criterion::{Criterion, black_box};
use polynomials::multilinear::MultilinearPolynomial;
use sumcheck::sumcheck_over_multilinear::{prove, verify};

pub fn sumcheck_benchmarks(c: &mut Criterion) {
    let mut rng = rand::thread_rng();
    let mut group = c.benchmark_group("sumcheck");

    let poly_size = 1 << 12; // 2^12 = 4096 values
    let polynomial =
        MultilinearPolynomial::new((0..poly_size).map(|_| Fq::rand(&mut rng)).collect());
    let sum: Fq = polynomial.evals_slice().iter().copied().sum();

    group.bench_function("sumcheck prove", |b| {
        b.iter(|| black_box(prove(&polynomial, sum)))
    });

    let proof = prove(&polynomial, sum);

    group.bench_function("sumcheck verify", |b| {
        b.iter(|| black_box(verify(&polynomial, &proof)))
    });
}
