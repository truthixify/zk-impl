use ark_bls12_381::Fq;
use ark_ff::UniformRand;
use circuit::{Circuit, Gate, Layer, Op};
use criterion::{Criterion, black_box};

fn build_sample_circuit() -> Circuit<Fq> {
    let layer2 = Layer::new(vec![
        Gate::new(Op::Add, 0, 0, 1),
        Gate::new(Op::Mul, 1, 2, 3),
        Gate::new(Op::Add, 2, 4, 5),
        Gate::new(Op::Mul, 3, 6, 7),
    ]);

    let layer1 = Layer::new(vec![
        Gate::new(Op::Mul, 0, 0, 1),
        Gate::new(Op::Add, 1, 2, 3),
    ]);

    let layer0 = Layer::new(vec![Gate::new(Op::Add, 0, 0, 1)]);

    Circuit::new(vec![layer0, layer1, layer2])
}

pub fn circuit_benchmarks(c: &mut Criterion) {
    let mut rng = rand::thread_rng();
    let mut group = c.benchmark_group("circuit");
    let mut circuit = build_sample_circuit();
    let mut input = Vec::with_capacity(20);

    for _ in 0..20 {
        input.push(Fq::rand(&mut rng));
    }

    group.bench_function("circuit_run_8_inputs", |b| {
        b.iter(|| {
            black_box(circuit.evaluate(input.clone()));
        });
    });

    group.bench_function("mle generation layer 2", |b| {
        b.iter(|| {
            circuit.add_i_and_mul_i_polynomials(2);
        });
    });

    group.bench_function("mle generation layer 1", |b| {
        b.iter(|| {
            circuit.add_i_and_mul_i_polynomials(1);
        });
    });

    group.bench_function("mle generation layer 0", |b| {
        b.iter(|| {
            circuit.add_i_and_mul_i_polynomials(0);
        });
    });
}
