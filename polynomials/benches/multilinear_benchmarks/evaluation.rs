use ark_bls12_381::Fq;
use ark_ff::UniformRand;
use criterion::{Criterion, black_box};
use polynomials::multilinear::evaluation::MultilinearPolynomial;
use rand::thread_rng;

fn sample_poly(num_vars: usize) -> MultilinearPolynomial<Fq> {
    let mut rng = thread_rng();
    let num_evals = 1 << num_vars;
    let evals = (0..num_evals)
        .map(|_| Fq::rand(&mut rng))
        .collect::<Vec<_>>();

    MultilinearPolynomial::new(evals)
}

pub fn evaluation_form_multilinear_polynomial_benchmarks(c: &mut Criterion) {
    let mut group = c.benchmark_group("multilinear polynomials in evaluation form");
    let num_vars = 10; // 2^10 = 1024 evaluations
    let poly = sample_poly(num_vars);
    let poly2 = sample_poly(num_vars);
    let mut rng = thread_rng();

    group.bench_function("scalar multiplication", |b| {
        let scalar = Fq::rand(&mut rng);
        b.iter(|| {
            black_box(poly.scalar_mul(scalar));
        });
    });

    group.bench_function("full evaluation", |b| {
        let point: Vec<Fq> = (0..num_vars).map(|_| Fq::rand(&mut rng)).collect();
        b.iter(|| {
            black_box(poly.evaluate(&point));
        });
    });

    group.bench_function("partial evaluation (fix 5 vars)", |b| {
        let fixed: Vec<(Fq, usize)> = (0..5).map(|i| (Fq::rand(&mut rng), i)).collect();
        b.iter(|| {
            black_box(poly.partial_evaluate(&fixed));
        });
    });

    group.bench_function("partial evaluation (fix all vars)", |b| {
        let fixed: Vec<(Fq, usize)> = (0..num_vars).map(|i| (Fq::rand(&mut rng), i)).collect();
        b.iter(|| {
            black_box(poly.partial_evaluate(&fixed));
        });
    });

    group.bench_function("polynomial addition", |b| {
        b.iter(|| {
            let _sum = MultilinearPolynomial {
                evals: poly
                    .evals
                    .iter()
                    .zip(&poly2.evals)
                    .map(|(a, b)| *a + *b)
                    .collect(),
            };
            black_box(_sum);
        });
    });

    group.finish();
}
