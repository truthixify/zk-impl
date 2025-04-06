use ark_bls12_381::Fq;
use ark_ff::UniformRand;
use criterion::{Criterion, black_box};
use shamir_secret_sharing::sss_with_password::{recover_secret, shares};

pub fn sss_with_password_benchmarks(c: &mut Criterion) {
    let mut rng = rand::thread_rng();
    let mut group = c.benchmark_group("shamir secret sharing with password");

    group.bench_function("shares", |b| {
        b.iter(|| {
            let secret = Fq::rand(&mut rng);
            let password = Fq::from(299792458);
            let num_shares = 100;
            let threshold = 50;
            black_box(shares(secret, password, num_shares, threshold));
        })
    });

    group.bench_function("recover secret", |b| {
        b.iter(|| {
            let secret = Fq::rand(&mut rng);
            let password = Fq::from(299792458);
            let num_shares = 100;
            let threshold = 50;
            let shares = shares(secret, password, num_shares, threshold);
            black_box(recover_secret(shares, password));
        })
    });
}
