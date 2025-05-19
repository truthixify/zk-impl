use ark_bls12_381::Fq;
use ark_ff::UniformRand;
use criterion::{Criterion, black_box};
use polynomials::{
    composed::{ProductPolynomial, SumPolynomial},
    multilinear::MultilinearPolynomial,
};
use sumcheck::{prove, verify};

// Generate a synthetic SumPolynomial for testing
fn setup_polynomial(num_evals: usize) -> SumPolynomial<Fq> {
    let mut rng = rand::thread_rng();
    let mut products = Vec::new();

    for _ in 0..2 {
        let mut evals1 = Vec::new();
        let mut evals2 = Vec::new();

        for _ in 0..num_evals {
            let eval1 = Fq::rand(&mut rng);
            let eval2 = Fq::rand(&mut rng);

            evals1.push(eval1);
            evals2.push(eval2);
        }

        let poly1 = MultilinearPolynomial::new(evals1);
        let poly2 = MultilinearPolynomial::new(evals2);

        products.push(ProductPolynomial::new(vec![poly1, poly2]));
    }

    SumPolynomial::new(products)
}

pub fn sumcheck_benchmarks(c: &mut Criterion) {
    let mut group = c.benchmark_group("sumcheck");

    let sum_polynomial = setup_polynomial(16);
    let (claimed_sum, round_polys) = prove(sum_polynomial.clone());

    group.bench_function("sumcheck prove", |b| {
        b.iter(|| black_box(prove(sum_polynomial.clone())))
    });

    group.bench_function("sumcheck verify", |b| {
        b.iter(|| {
            black_box(verify(
                sum_polynomial.clone(),
                claimed_sum,
                round_polys.clone(),
            ))
        })
    });
}
