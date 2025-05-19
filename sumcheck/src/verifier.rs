use ark_ff::PrimeField;
use polynomials::{composed::SumPolynomial, univariate::DenseUnivariatePolynomial};
use sha3::Keccak256;
use transcript::Transcript;

pub fn partial_verify<F: PrimeField>(
    transcript: &mut Transcript<F, Keccak256>,
    claimed_sum: F,
    round_polynomials: Vec<DenseUnivariatePolynomial<F>>,
) -> (bool, F, Vec<F>) {
    if round_polynomials.is_empty() {
        return (false, claimed_sum, vec![]);
    }

    transcript.append_field_element(&claimed_sum);
    let mut current_sum: F = claimed_sum;
    let mut challenges: Vec<F> = Vec::new();

    for round_polynomial in round_polynomials {
        let p_0 = round_polynomial.evaluate(F::ZERO);
        let p_1 = round_polynomial.evaluate(F::ONE);

        if current_sum != p_0 + p_1 {
            println!(
                "cs: {}, p_0: {}, p_1: {}, fal: {}",
                current_sum,
                p_0,
                p_1,
                claimed_sum == p_0 + p_1
            );

            return (false, current_sum, challenges);
        }

        transcript.append(&round_polynomial.to_bytes());

        let challenge = transcript.sample_field_element();

        current_sum = round_polynomial.evaluate(challenge);
        challenges.push(challenge);
    }

    (true, current_sum, challenges)
}

pub fn verify<F: PrimeField>(
    sum_polynomial: SumPolynomial<F>,
    claimed_sum: F,
    round_polynomials: Vec<DenseUnivariatePolynomial<F>>,
) -> bool {
    let mut transcript: Transcript<F, Keccak256> = Transcript::new();

    transcript.append(&sum_polynomial.to_bytes());

    let (is_partially_verified, claimed_sum, challenges) =
        partial_verify(&mut transcript, claimed_sum, round_polynomials);

    if !is_partially_verified {
        return false;
    }

    let derived_sum = sum_polynomial.evaluate(&challenges);

    // Perform oracle check
    claimed_sum == derived_sum
}
