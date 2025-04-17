use ark_bls12_381::Fq;
use ark_ff::UniformRand;
use circuit::{Circuit, Gate, Layer, Op};
use criterion::{Criterion, black_box};
use rand::Rng;

fn build_sample_circuit(num_of_layers: usize) -> Circuit<Fq> {
    let input_size = 1 << num_of_layers;
    let mut i = input_size;
    let mut layers: Vec<Layer<Fq>> = vec![];

    while i > 1 {
        let mut layer = vec![];

        for j in (0..i).step_by(2) {
            let rng = rand::thread_rng().gen_range(0..1);
            let gate: Gate;

            if rng == 1 {
                gate = Gate::new(Op::Add, j / 2, j, j + 1);
            } else {
                gate = Gate::new(Op::Mul, j / 2, j, j + 1);
            }

            layer.push(gate);
        }

        layers.push(Layer::new(layer));

        i = i / 2;
    }

    layers.reverse();
    Circuit::new(layers)
}

pub fn circuit_benchmarks(c: &mut Criterion) {
    let mut rng = rand::thread_rng();
    let num_of_layers = 10;
    let input_size = 1 << (num_of_layers + 1);
    let mut group = c.benchmark_group("circuit");
    let mut circuit = build_sample_circuit(num_of_layers);
    let mut input = Vec::with_capacity(input_size);

    for _ in 0..input_size {
        input.push(Fq::rand(&mut rng));
    }

    group.bench_function(format!("circuit evaluate {} inputs", input_size), |b| {
        b.iter(|| {
            black_box(circuit.evaluate(input.clone()));
        });
    });

    for i in (0..num_of_layers).rev() {
        group.bench_function(format!("mle generation layer {}", i + 1), |b| {
            b.iter(|| {
                black_box(circuit.add_i_and_mul_i_polynomials(i as usize));
            });
        });
    }
}
