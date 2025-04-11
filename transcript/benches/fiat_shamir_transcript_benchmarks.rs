use ark_bls12_381::Fq;
use ark_ff::UniformRand;
use criterion::{Criterion, black_box};
use sha3::Keccak256;
use transcript::Transcript;

pub fn fiat_shamir_transcript_benchmarks(c: &mut Criterion) {
    let mut rng = rand::thread_rng();
    let mut group = c.benchmark_group("fiat-shamir transcript");

    let mut transcript = Transcript::<Fq, Keccak256>::new();
    let element = Fq::rand(&mut rng);

    group.bench_function("append_field_element", |b| {
        b.iter(|| {
            black_box(transcript.append_field_element(&element));
        });
    });

    let mut transcript = Transcript::<Fq, Keccak256>::new();
    let element = Fq::rand(&mut rng);

    transcript.append_field_element(&element);

    group.bench_function("sample_field_element", |b| {
        b.iter(|| {
            black_box(transcript.sample_field_element());
        });
    });

    let mut transcript = Transcript::<Fq, Keccak256>::new();
    transcript.append_field_element(&Fq::rand(&mut rng));

    group.bench_function("sample_n_field_elements_100", |b| {
        b.iter(|| {
            black_box(transcript.sample_n_field_elements(100));
        });
    });
}
