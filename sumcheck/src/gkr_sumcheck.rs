use ark_ff::PrimeField;
use sha3::Keccak256;
use polynomials::composed::SumPolynomial;
use polynomials::univariate::DenseUnivariatePolynomial;
use transcript::Transcript;

#[derive(Debug)]
pub struct SumcheckProof<F: PrimeField> {
    claimed_sum: F,
    round_polynomials: Vec<DenseUnivariatePolynomial<F>>,
    challenges: Vec<F>,
}

pub fn prove<F: PrimeField>(claimed_sum: F, sum_polynomial: SumPolynomial<F>) -> SumcheckProof<F> {
    let mut transcript: Transcript<F, Keccak256> = Transcript::new();

    transcript.append_field_element(&claimed_sum);
    transcript.append(&sum_polynomial.to_bytes());

    prover_partial(claimed_sum, sum_polynomial, &mut transcript)
}

pub fn prover_partial<F: PrimeField>(claimed_sum: F, sum_polynomial: SumPolynomial<F>, transcript: &mut Transcript<F, Keccak256>) -> SumcheckProof<F> {
    transcript.append_field_element(&claimed_sum);

    let mut round_polynomials = vec![];
    let mut challenges = vec![];
    let mut current_polynomial = sum_polynomial.clone();

    for i in 0..sum_polynomial.n_vars() {
        let round_polynomial_evals = get_round_polynomial(sum_polynomial.clone());
        let xs = (0..sum_polynomial.degree() + 1).map(|x| F::from(x as u64)).collect::<Vec<F>>();
        let univariate_polynomial = DenseUnivariatePolynomial::interpolate(&xs, &round_polynomial_evals);

        transcript.append(&univariate_polynomial.to_bytes());
        round_polynomials.push(univariate_polynomial);

        let challenge = transcript.sample_field_element();

        current_polynomial = current_polynomial.partial_evaluate(challenge, 0);
        challenges.push(challenge);
    }

    SumcheckProof {
        claimed_sum,
        round_polynomials,
        challenges,
    }
}

fn get_round_polynomial<F: PrimeField>(polynomial: SumPolynomial<F>) -> Vec<F> {
    let num_evals = polynomial.degree() + 1;
    let mut evals = Vec::with_capacity(num_evals);

    for i in 0..num_evals {
        let point = F::from(i as u64);
        let partial_polynomial = polynomial.partial_evaluate(point, 0);

        let eval = partial_polynomial.element_wise_add().evals_slice().iter().sum();

        evals.push(eval);
    }

    evals
}
