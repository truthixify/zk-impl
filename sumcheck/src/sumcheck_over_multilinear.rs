use ark_ff::{BigInteger, PrimeField};
use polynomials::multilinear::MultilinearPolynomial;
use polynomials::univariate::DenseUnivariatePolynomial as UnivariatePolynomial;
use sha3::Keccak256;
use transcript::Transcript;

#[derive(Debug)]
pub struct SumcheckProof<F: PrimeField> {
    claimed_sum: F,
    round_polynomials: Vec<UnivariatePolynomial<F>>,
}

impl<F: PrimeField> SumcheckProof<F> {
    pub fn new(claimed_sum: F, round_polynomials: Vec<UnivariatePolynomial<F>>) -> Self {
        SumcheckProof {
            claimed_sum,
            round_polynomials,
        }
    }
}

pub fn prove<F: PrimeField>(
    polynomial: &MultilinearPolynomial<F>,
    claimed_sum: F,
) -> SumcheckProof<F> {
    let mut round_polynomials = vec![];

    let mut transcript: Transcript<F, Keccak256> = Transcript::new();
    transcript.append_field_element(&claimed_sum);
    transcript.append(&polynomial.to_bytes());

    let mut polynomial = polynomial.clone();

    for _ in 0..polynomial.n_vars() {
        let round_polynomial = skip_one_and_sum_over_boolean_hypercube(&polynomial);

        transcript.append(
            &round_polynomial
                .coefficients_slice()
                .iter()
                .flat_map(|coeff| coeff.into_bigint().to_bytes_be())
                .collect::<Vec<_>>(),
        );

        round_polynomials.push(round_polynomial);

        let challenge = transcript.sample_field_element();

        polynomial = polynomial.partial_evaluate(challenge, 0);
    }

    SumcheckProof::new(claimed_sum, round_polynomials)
}

pub fn verify<F: PrimeField>(
    polynomial: &MultilinearPolynomial<F>,
    proof: &SumcheckProof<F>,
) -> bool {
    if proof.round_polynomials.len() != polynomial.n_vars() {
        return false;
    }

    let mut transcript: Transcript<F, Keccak256> = Transcript::new();
    transcript.append_field_element(&proof.claimed_sum);
    transcript.append(&polynomial.to_bytes());

    let mut claimed_sum = proof.claimed_sum;
    let mut challenges = vec![];

    for round_polynomial in &proof.round_polynomials {
        let p_0 = round_polynomial.evaluate(F::ZERO);
        let p_1 = round_polynomial.evaluate(F::ONE);

        if claimed_sum != p_0 + p_1 {
            return false;
        }

        transcript.append(&round_polynomial.to_bytes());

        let challenge = transcript.sample_field_element();

        claimed_sum = round_polynomial.evaluate(challenge);
        challenges.push(challenge);
    }

    // perform oracle check
    if claimed_sum != polynomial.evaluate(&challenges) {
        return false;
    }

    true
}

pub fn skip_one_and_sum_over_boolean_hypercube<F: PrimeField>(
    polynomial: &MultilinearPolynomial<F>,
) -> UnivariatePolynomial<F> {
    let evals = polynomial.evals_slice();
    let (left_half, right_half) = evals.split_at(evals.len() / 2);
    let f_0 = left_half.iter().sum();
    let f_1 = right_half.iter().sum();

    UnivariatePolynomial::interpolate_y(vec![f_0, f_1])
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::Fq as ArkField;
    use field_tracker::{Ft, end_tscope, print_summary, start_tscope};

    type Fq = Ft!(ArkField);

    fn fq(x: u64) -> Fq {
        Fq::from(x)
    }

    #[test]
    fn test_sumcheck() {
        start_tscope!("sumcheck");
        let polynomial: MultilinearPolynomial<Fq> = MultilinearPolynomial::new(vec![
            fq(0),
            fq(0),
            fq(0),
            fq(3),
            fq(0),
            fq(0),
            fq(2),
            fq(5),
        ]);

        let proof = prove(&polynomial, fq(10));
        assert!(verify(&polynomial, &proof));
        end_tscope!();
        print_summary!();
    }

    #[test]
    fn test_sumcheck_valid_proof() {
        let polynomial: MultilinearPolynomial<Fq> = MultilinearPolynomial::new(vec![
            fq(0),
            fq(0),
            fq(0),
            fq(3),
            fq(0),
            fq(0),
            fq(2),
            fq(5),
        ]);
        let proof = prove(&polynomial, fq(10));
        assert!(verify(&polynomial, &proof));
    }

    #[test]
    fn test_sumcheck_invalid_sum() {
        let polynomial: MultilinearPolynomial<Fq> = MultilinearPolynomial::new(vec![
            fq(0),
            fq(0),
            fq(0),
            fq(3),
            fq(0),
            fq(0),
            fq(2),
            fq(5),
        ]);
        let proof = prove(&polynomial, fq(9)); // incorrect claimed sum
        assert!(!verify(&polynomial, &proof));
    }

    #[test]
    fn test_sumcheck_invalid_polynomial() {
        let poly_correct = MultilinearPolynomial::new(vec![
            fq(0),
            fq(0),
            fq(0),
            fq(3),
            fq(0),
            fq(0),
            fq(2),
            fq(5),
        ]);
        let poly_wrong = MultilinearPolynomial::new(vec![
            fq(0),
            fq(0),
            fq(0),
            fq(3),
            fq(0),
            fq(0),
            fq(2),
            fq(4),
        ]);
        let proof = prove(&poly_correct, fq(10));
        assert!(!verify(&poly_wrong, &proof));
    }
}
