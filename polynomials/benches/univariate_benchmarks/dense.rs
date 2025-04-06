use ark_bls12_381::Fq;
use ark_ff::UniformRand;
use criterion::{Criterion, black_box};
use polynomials::univariate::dense::DenseUnivariatePolynomial;

fn sample_poly() -> DenseUnivariatePolynomial<Fq> {
    let mut rng = rand::thread_rng();
    let mut coeffs = Vec::with_capacity(100);
    for _ in 0..100 {
        coeffs.push(Fq::rand(&mut rng));
    }

    DenseUnivariatePolynomial::new(coeffs)
}
fn sample_poly_2() -> DenseUnivariatePolynomial<Fq> {
    let mut rng = rand::thread_rng();
    let mut coeffs = Vec::with_capacity(100);
    for _ in 0..100 {
        coeffs.push(Fq::rand(&mut rng));
    }

    DenseUnivariatePolynomial::new(coeffs)
}

pub fn dense_univariate_polynomial_benchmarks(c: &mut Criterion) {
    let mut rng = rand::thread_rng();
    let mut group = c.benchmark_group("univariate dense polynomials");
    let poly = sample_poly();
    let poly_2 = sample_poly_2();

    group.bench_function("polynomial degree", |b| {
        b.iter(|| {
            black_box(poly.degree());
        })
    });

    group.bench_function("polynomial scalar multiplication", |b| {
        let scalar = Fq::rand(&mut rng);
        b.iter(|| black_box(poly.scalar_mul(scalar)))
    });

    group.bench_function("polynomial evaluation", |b| {
        b.iter(|| {
            let x = Fq::rand(&mut rng);
            black_box(poly.evaluate(x));
        })
    });

    group.bench_function("polynomial addition", |b| {
        b.iter(|| black_box(&poly + &poly_2))
    });

    group.bench_function("polynomial multiplication", |b| {
        b.iter(|| black_box(&poly * &poly_2))
    });

    group.bench_function("polynomial interpolation", |b| {
        b.iter(|| {
            let mut xs = Vec::with_capacity(100);
            let mut ys = Vec::with_capacity(100);
            for _ in 0..100 {
                xs.push(Fq::rand(&mut rng));
                ys.push(Fq::rand(&mut rng));
            }
            black_box(DenseUnivariatePolynomial::interpolate(xs, ys))
        });
    });
}
