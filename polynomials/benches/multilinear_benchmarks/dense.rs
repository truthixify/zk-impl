use ark_bls12_381::Fq;
use ark_ff::UniformRand;
use criterion::{Criterion, black_box};
use polynomials::multilinear::dense::DenseMultilinearPolynomial;

fn sample_poly() -> DenseMultilinearPolynomial<Fq> {
    let mut rng = rand::thread_rng();
    let mut terms = Vec::with_capacity(100);

    for _ in 0..100 {
        terms.push(Fq::rand(&mut rng));
    }

    DenseMultilinearPolynomial::new_with_coefficients(terms, 100)
}
fn sample_poly_2() -> DenseMultilinearPolynomial<Fq> {
    let mut rng = rand::thread_rng();
    let mut terms = Vec::with_capacity(100);

    for _ in 0..100 {
        terms.push(Fq::rand(&mut rng));
    }

    DenseMultilinearPolynomial::new_with_coefficients(terms, 100)
}

pub fn dense_multilinear_polynomial_benchmarks(c: &mut Criterion) {
    let mut rng = rand::thread_rng();
    let mut group = c.benchmark_group("multilinear dense polynomials");
    let poly = sample_poly();
    let poly_2 = sample_poly_2();

    group.bench_function("polynomial scalar multiplication", |b| {
        let scalar = Fq::rand(&mut rng);
        b.iter(|| black_box(poly.scalar_mul(scalar)))
    });

    group.bench_function("polynomial evaluation", |b| {
        b.iter(|| {
            let mut points = Vec::with_capacity(100);
            for i in 0..100 {
                points.push((Fq::rand(&mut rng), i));
            }
            black_box(poly.evaluate(&points));
        })
    });

    group.bench_function("polynomial partial evaluation", |b| {
        b.iter(|| {
            let mut partial_terms = Vec::with_capacity(10);
            for i in 0..10 {
                partial_terms.push((Fq::rand(&mut rng), i));
            }
            black_box(poly.partial_evaluate(&partial_terms));
        })
    });

    group.bench_function("polynomial addition", |b| {
        b.iter(|| black_box(&poly + &poly_2))
    });

    // group.bench_function("polynomial multiplication", |b| {
    //     b.iter(|| black_box(&poly * &poly_2))
    // });

    // group.bench_function("polynomial interpolation", |b| {
    //     b.iter(|| {
    //         let mut xs = Vec::with_capacity(100);
    //         let mut ys = Vec::with_capacity(100);
    //         for _ in 0..100 {
    //             xs.push(Fq::rand(&mut rng));
    //             ys.push(Fq::rand(&mut rng));
    //         }
    //         black_box(DenseMultilinearPolynomial::interpolate(xs, ys))
    //     });
    // });
}
