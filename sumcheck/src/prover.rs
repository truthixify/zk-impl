use ark_ff::PrimeField;
use polynomials::{composed::SumPolynomial, univariate::DenseUnivariatePolynomial};
use sha3::Keccak256;
use transcript::Transcript;

pub fn partial_prove<F: PrimeField>(
    mut sum_polynomial: SumPolynomial<F>,
    transcript: &mut Transcript<F, Keccak256>,
) -> (F, Vec<DenseUnivariatePolynomial<F>>, Vec<F>) {
    let claimed_sum = sum_polynomial
        .element_wise_add()
        .evals_slice()
        .into_iter()
        .sum();
    let n_vars = sum_polynomial.n_vars();
    let mut round_polynomials = Vec::with_capacity(n_vars);
    let mut challenges = Vec::with_capacity(n_vars);

    transcript.append_field_element(&claimed_sum);

    for _ in 0..n_vars {
        let num_evals = sum_polynomial.degree() + 1;
        let mut evals = Vec::with_capacity(num_evals);

        for i in 0..num_evals {
            let point = F::from(i as u64);
            let partial_polynomial = sum_polynomial.partial_evaluate(point, 0);

            let eval: F = partial_polynomial
                .element_wise_add()
                .evals_slice()
                .iter()
                .sum();

            evals.push(eval);
        }

        let round_polynomial = DenseUnivariatePolynomial::interpolate_y(evals);

        transcript.append(&round_polynomial.to_bytes());
        round_polynomials.push(round_polynomial);

        let challenge = transcript.sample_field_element();
        challenges.push(challenge);

        sum_polynomial = sum_polynomial.partial_evaluate(challenge, 0);
    }

    (claimed_sum, round_polynomials, challenges)
}

pub fn prove<F: PrimeField>(
    sum_polynomial: SumPolynomial<F>,
) -> (F, Vec<DenseUnivariatePolynomial<F>>, Vec<F>) {
    let mut transcript: Transcript<F, Keccak256> = Transcript::new();

    transcript.append(&sum_polynomial.to_bytes());

    partial_prove(sum_polynomial, &mut transcript)
}
