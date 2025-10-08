use ark_ff::PrimeField;
use circuit::{Circuit, Gate, Layer, Op};
use polynomials::{
    composed::SumPolynomial, multilinear::MultilinearPolynomial,
    univariate::DenseUnivariatePolynomial,
};
use sha3::Keccak256;
use sumcheck::{partial_prove, partial_verify, prove as sumcheck_prove, verify as sumcheck_verify};
use transcript::Transcript;

pub struct GKRProofResult<F: PrimeField> {
    pub claimed_sum: F,
    pub output_layer: Vec<F>,
    pub proofs: Vec<(F, Vec<DenseUnivariatePolynomial<F>>, Vec<F>)>,
    pub wb_evals: Vec<F>,
    pub wc_evals: Vec<F>,
}

// pub fn prove<F: PrimeField>(circuit: &mut Circuit<F>, input: Vec<F>) -> GKRProofResult<F> {
// let circuit_eval = circuit.evaluate(input);
// let mut transcript: Transcript<F, Keccak256> = Transcript::new();
// let mut wb_evals = Vec::new();
// let mut wc_evals = Vec::new();
// let mut alpha = F::ZERO;
// let mut beta = F::ONE;
// let mut rbs = Vec::new();
// let mut rcs = Vec::new();
// let mut w_0 = circuit.w_i_polynomial(0);

// if w_0.evals_slice().len() == 1 {
//     let mut w_0_evals = w_0.evals_slice().to_vec();

//     w_0_evals.push(F::ZERO);

//     w_0 = MultilinearPolynomial::new(w_0_evals);
// }

// transcript.append(&w_0.to_bytes());
// let challenge = transcript.sample_field_element(); // r_0
// let mut claimed_sum = w_0.evaluate(&vec![challenge]); // D(r_0) = m_0

// for layer_index in 0..circuit.layers.len() {
//     let (add_i_rbc, mul_i_rbc) = circuit.add_i_and_mul_i_polynomials(layer_index);
//     let add_i_bc = add_i_rbc.partial_evaluate(challenge, 0);
//     let mul_i_bc = mul_i_rbc.partial_evaluate(challenge, 0);

//     // let (add_i_bc, mul_i_bc) = if layer_index == 0 {
//     //     (
//     //         add_i_rbc.partial_evaluate(challenge, 0),
//     //         mul_i_rbc.partial_evaluate(challenge, 0),
//     //     )
//     // } else {
//     // };
// }

// todo!()
// }
