use ark_ff::PrimeField;
use polynomials::multilinear::MultilinearPolynomial;
use std::marker::PhantomData;

#[derive(Debug)]
pub enum Op {
    Add,
    Mul,
}

#[derive(Debug)]
pub struct Gate {
    op: Op,
    output: usize,
    left_index: usize,
    right_index: usize,
}

impl Gate {
    pub fn new(op: Op, output: usize, left_index: usize, right_index: usize) -> Self {
        Gate {
            op,
            output,
            left_index,
            right_index,
        }
    }

    pub fn eval_gate<F: PrimeField>(&self, layer_eval: &[F]) -> F {
        let left_val = layer_eval[self.left_index];
        let right_val = layer_eval[self.right_index];

        match self.op {
            Op::Add => left_val + right_val,
            Op::Mul => left_val * right_val,
        }
    }
}

#[derive(Debug)]
pub struct Layer<F: PrimeField> {
    gates: Vec<Gate>,
    _phantom: PhantomData<F>,
}

impl<F: PrimeField> Layer<F> {
    pub fn new(gates: Vec<Gate>) -> Self {
        Self {
            gates,
            _phantom: PhantomData,
        }
    }

    pub fn num_layer_vars(&self) -> usize {
        let layer_index = self.layer_index();

        if layer_index == 0 {
            3
        } else {
            3 * layer_index + 2
        }
    }

    pub fn layer_index(&self) -> usize {
        self.gates.len().ilog2() as usize
    }

    pub fn add_i_and_mul_i_polynomials(
        &self,
    ) -> (MultilinearPolynomial<F>, MultilinearPolynomial<F>) {
        let num_boolean_hypercube_evals = 1 << self.num_layer_vars();
        let layer_index = self.layer_index();

        let mut add_i_evals = vec![F::ZERO; num_boolean_hypercube_evals];
        let mut mul_i_evals = vec![F::ZERO; num_boolean_hypercube_evals];

        for gate in &self.gates {
            let postional_index =
                get_positional_index(layer_index, gate.output, gate.left_index, gate.right_index);

            match gate.op {
                Op::Add => add_i_evals[postional_index] = F::ONE,
                Op::Mul => mul_i_evals[postional_index] = F::ONE,
            }
        }

        (
            MultilinearPolynomial::new(add_i_evals),
            MultilinearPolynomial::new(mul_i_evals),
        )
    }
}

#[derive(Debug)]
pub struct Circuit<F: PrimeField> {
    layers: Vec<Layer<F>>,
    layer_evals: Vec<Vec<F>>,
}

impl<F: PrimeField> Circuit<F> {
    pub fn new(layers: Vec<Layer<F>>) -> Self {
        assert!(
            !layers.is_empty(),
            "Circuit must contain at least one layer"
        );

        let layer_evals = vec![vec![]; 1 << (layers.len() - 1)];

        Circuit {
            layers,
            layer_evals,
        }
    }

    pub fn evaluate(&mut self, initial_layer_eval: Vec<F>) -> Vec<F> {
        let mut current_layer_eval = initial_layer_eval;
        let mut resultant_evals = vec![];

        resultant_evals.push(current_layer_eval.clone());

        for layer in self.layers.iter().rev() {
            let max_layer_index = layer
                .gates
                .iter()
                .map(|gate| gate.output)
                .max()
                .unwrap_or(0);
            let mut evals = vec![F::ZERO; max_layer_index + 1];

            for gate in layer.gates.iter() {
                let current_gate_eval = gate.eval_gate(&current_layer_eval);

                evals[gate.output] += current_gate_eval;
            }

            current_layer_eval = evals;
            resultant_evals.push(current_layer_eval.clone());
        }

        resultant_evals.reverse();
        self.layer_evals = resultant_evals.clone();

        (&resultant_evals[0]).to_vec()
    }

    pub fn add_i_and_mul_i_polynomials(
        &self,
        layer_index: usize,
    ) -> (MultilinearPolynomial<F>, MultilinearPolynomial<F>) {
        self.layers[layer_index].add_i_and_mul_i_polynomials()
    }

    pub fn w_i_polynomial(&self, layer_index: usize) -> MultilinearPolynomial<F> {
        assert!(
            layer_index < self.layer_evals.len(),
            "Layer index cannot be greater than total number of layers"
        );

        MultilinearPolynomial::new(self.layer_evals[layer_index].clone())
    }
}

pub fn get_positional_index(
    layer_index: usize,
    output_index: usize,
    left_index: usize,
    right_index: usize,
) -> usize {
    let output_padded_bin = format!("{:0>width$b}", output_index, width = layer_index);
    let left_padded_bin = format!("{:0>width$b}", left_index, width = layer_index + 1);
    let right_padded_bin = format!("{:0>width$b}", right_index, width = layer_index + 1);

    let sum = output_padded_bin + &left_padded_bin + &right_padded_bin;

    usize::from_str_radix(&sum, 2).unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::Fq;
    use ark_ff::Field;

    fn fq(val: u64) -> Fq {
        Fq::from(val)
    }

    #[test]
    fn test_circuit_evaluation() {
        // Initial input values for the circuit (bottom layer)
        // These will feed into the gates of the next layer up
        let input = vec![fq(1), fq(2), fq(3), fq(4)];

        // --- Layer 2 (Bottom Layer): 4 inputs -> 2 outputs ---

        // Gate 0: Add input[0] + input[1] = 1 + 2 = 3 → stored at index 0 of layer 2 output
        let layer2_gate1 = Gate::new(Op::Add, 0, 0, 1);

        // Gate 1: Mul input[2] * input[3] = 3 * 4 = 12 → stored at index 1 of layer 2 output
        let layer_2gate2 = Gate::new(Op::Mul, 1, 2, 3);

        // Layer 2 has two gates, producing [3, 12]
        let layer2 = Layer::new(vec![layer_2gate2, layer2_gate1]);

        // --- Layer 1 (Middle Layer): 2 inputs → 1 output ---

        // Gate: Add layer2_output[0] + layer2_output[1] = 3 + 12 = 15 → output of entire circuit
        let layer1_gate1 = Gate::new(Op::Add, 0, 0, 1);

        // Layer 1 has one gate, producing [15]
        let layer1 = Layer::new(vec![layer1_gate1]);

        // --- Construct the circuit ---
        // The layers are in reverse order because the circuit runs from top to bottom (layer0 is final output)
        let mut circuit = Circuit::<Fq>::new(vec![layer1, layer2]);

        // --- Run the circuit on the given input ---
        // It evaluates layer2 first (from inputs), then layer1, and returns the final result
        let result = circuit.evaluate(input);

        // --- Expected evaluations for each layer ---
        let expected_layers_evaluation = vec![
            vec![fq(15)],                     // Final output from layer1
            vec![fq(3), fq(12)],              // Output from layer2 (from input)
            vec![fq(1), fq(2), fq(3), fq(4)], // Original input
        ];

        // --- Check that intermediate values match expected layer evaluations ---
        assert_eq!(circuit.layer_evals, expected_layers_evaluation);

        // --- Check that the final result is as expected ---
        assert_eq!(result[0], fq(15));
    }

    #[test]
    fn test_gate_eval_add_and_mul() {
        let layer_eval = vec![fq(2), fq(3)];

        let add_gate = Gate::new(Op::Add, 0, 0, 1);
        let mul_gate = Gate::new(Op::Mul, 0, 0, 1);

        assert_eq!(add_gate.eval_gate(&layer_eval), fq(5));
        assert_eq!(mul_gate.eval_gate(&layer_eval), fq(6));
    }

    #[test]
    fn test_get_positional_index() {
        // layer_index = 2, binary lengths: output = 2, left = 3, right = 3
        let idx = get_positional_index(2, 2, 3, 4);

        // binary: output=10, left=011, right=100 => combined = 10011100
        // binary 10011100 = decimal 156
        assert_eq!(idx, 156);
    }

    #[test]
    fn test_add_i_and_mul_i_polynomials() {
        let add_gate = Gate::new(Op::Add, 0, 0, 1);
        let mul_gate = Gate::new(Op::Mul, 1, 1, 2);
        let layer = Layer::<Fq>::new(vec![add_gate, mul_gate]);

        let (add_poly, mul_poly) = layer.add_i_and_mul_i_polynomials();

        let add_count = add_poly
            .evals_slice()
            .iter()
            .filter(|&&x| x == Fq::ONE)
            .count();
        let mul_count = mul_poly
            .evals_slice()
            .iter()
            .filter(|&&x| x == Fq::ONE)
            .count();

        assert_eq!(add_count, 1);
        assert_eq!(mul_count, 1);
    }

    #[test]
    fn test_w_i_polynomial_returns_correct_layer_eval() {
        let input = vec![fq(1), fq(1), fq(1), fq(1)];
        let gate1 = Gate::new(Op::Add, 0, 0, 1);
        let gate2 = Gate::new(Op::Mul, 1, 2, 3);
        let layer = Layer::new(vec![gate1, gate2]);

        let mut circuit = Circuit::<Fq>::new(vec![layer]);
        circuit.evaluate(input.clone());

        let poly = circuit.w_i_polynomial(1);

        assert_eq!(poly.evals_slice(), input);
    }

    #[test]
    fn test_circuit_evaluation_add_mul_combo() {
        let input = vec![fq(1), fq(2), fq(3), fq(4)];

        let gate_add = Gate::new(Op::Add, 0, 0, 1); // 1 + 2 = 3
        let gate_mul = Gate::new(Op::Mul, 1, 2, 3); // 3 * 4 = 12
        let layer2 = Layer::new(vec![gate_add, gate_mul]);

        let gate_final = Gate::new(Op::Add, 0, 0, 1); // 3 + 12 = 15
        let layer1 = Layer::new(vec![gate_final]);

        let mut circuit = Circuit::<Fq>::new(vec![layer1, layer2]);

        let result = circuit.evaluate(input);

        let expected_layers = vec![
            vec![fq(15)],                     // output layer
            vec![fq(3), fq(12)],              // mid layer
            vec![fq(1), fq(2), fq(3), fq(4)], // input
        ];

        assert_eq!(circuit.layer_evals, expected_layers);
        assert_eq!(result, vec![fq(15)]);
    }

    #[test]
    fn test_circuit_with_single_layer_add_only() {
        let input = vec![fq(5), fq(7)];

        let gate = Gate::new(Op::Add, 0, 0, 1); // 5 + 7 = 12
        let layer = Layer::new(vec![gate]);

        let mut circuit = Circuit::<Fq>::new(vec![layer]);

        let result = circuit.evaluate(input);

        assert_eq!(result, vec![fq(12)]);
    }

    #[test]
    fn test_circuit_with_single_layer_mul_only() {
        let input = vec![fq(6), fq(2)];

        let gate = Gate::new(Op::Mul, 0, 0, 1); // 6 * 2 = 12
        let layer = Layer::new(vec![gate]);

        let mut circuit = Circuit::<Fq>::new(vec![layer]);

        let result = circuit.evaluate(input);

        assert_eq!(result, vec![fq(12)]);
    }

    #[test]
    #[should_panic(expected = "Circuit must contain at least one layer")]
    fn test_empty_circuit_no_layers_fails() {
        let _ = Circuit::<Fq>::new(vec![]);
    }

    #[test]
    fn test_invalid_layer_index_panics() {
        let input = vec![fq(1), fq(1)];

        let gate = Gate::new(Op::Add, 0, 0, 1);
        let layer = Layer::new(vec![gate]);

        let mut circuit = Circuit::<Fq>::new(vec![layer]);
        circuit.evaluate(input);

        let result = std::panic::catch_unwind(|| {
            circuit.w_i_polynomial(100); // way out of bounds
        });

        assert!(result.is_err());
    }

    #[test]
    fn test_multiple_layer_circuit_with_mixed_ops() {
        // Layer 3 input
        let input = vec![fq(2), fq(3), fq(4), fq(5)];

        // Layer 2:
        // Gate 0: Add 0 + 1 = 2 + 3 = 5
        // Gate 1: Mul 2 * 3 = 4 * 5 = 20
        let gate1 = Gate::new(Op::Add, 0, 0, 1);
        let gate2 = Gate::new(Op::Mul, 1, 2, 3);
        let layer2 = Layer::new(vec![gate1, gate2]);

        // Layer 1:
        // Gate 0: Mul of 5 * 20 = 100
        let gate3 = Gate::new(Op::Mul, 0, 0, 1);
        let layer1 = Layer::new(vec![gate3]);

        let mut circuit = Circuit::<Fq>::new(vec![layer1, layer2]);

        let output = circuit.evaluate(input);

        assert_eq!(output, vec![fq(100)]);
    }

    #[test]
    fn test_circuit_with_8_inputs_and_mixed_operations_per_layer() {
        let input = vec![fq(1), fq(2), fq(3), fq(4), fq(5), fq(6), fq(7), fq(8)];

        // --- Layer 2: 8 → 4, mixed ops ---
        let layer2_gate0 = Gate::new(Op::Add, 0, 0, 1); // 1 + 2 = 3
        let layer2_gate1 = Gate::new(Op::Mul, 1, 2, 3); // 3 * 4 = 12
        let layer2_gate2 = Gate::new(Op::Add, 2, 4, 5); // 5 + 6 = 11
        let layer2_gate3 = Gate::new(Op::Mul, 3, 6, 7); // 7 * 8 = 56
        let layer2 = Layer::new(vec![layer2_gate0, layer2_gate1, layer2_gate2, layer2_gate3]);

        // --- Layer 1: 4 → 2, mixed ops ---
        let layer1_gate0 = Gate::new(Op::Mul, 0, 0, 1); // 3 * 12 = 36
        let layer1_gate1 = Gate::new(Op::Add, 1, 2, 3); // 11 + 56 = 67
        let layer1 = Layer::new(vec![layer1_gate0, layer1_gate1]);

        // --- Layer 0: 2 → 1, mixed ops ---
        let layer0_gate0 = Gate::new(Op::Add, 0, 0, 1); // 36 + 67 = 103
        let layer0 = Layer::new(vec![layer0_gate0]);

        let mut circuit = Circuit::<Fq>::new(vec![layer0, layer1, layer2]);
        let result = circuit.evaluate(input.clone());

        // Final output
        assert_eq!(result, vec![fq(103)]);

        // Intermediate layer checks
        assert_eq!(circuit.layer_evals[2], vec![fq(3), fq(12), fq(11), fq(56)]);
        assert_eq!(circuit.layer_evals[1], vec![fq(36), fq(67)]);
        assert_eq!(circuit.layer_evals[0], vec![fq(103)]);

        // MLEs for layer 2 (2 Add, 2 Mul)
        let (add_poly2, mul_poly2) = circuit.add_i_and_mul_i_polynomials(2);
        assert_eq!(
            add_poly2
                .evals_slice()
                .iter()
                .filter(|&&x| x == Fq::ONE)
                .count(),
            2
        );
        assert_eq!(
            mul_poly2
                .evals_slice()
                .iter()
                .filter(|&&x| x == Fq::ONE)
                .count(),
            2
        );

        // MLEs for layer 1 (1 Add, 1 Mul)
        let (add_poly1, mul_poly1) = circuit.add_i_and_mul_i_polynomials(1);
        assert_eq!(
            add_poly1
                .evals_slice()
                .iter()
                .filter(|&&x| x == Fq::ONE)
                .count(),
            1
        );
        assert_eq!(
            mul_poly1
                .evals_slice()
                .iter()
                .filter(|&&x| x == Fq::ONE)
                .count(),
            1
        );

        // MLE for layer 0 (1 Add)
        let (add_poly0, mul_poly0) = circuit.add_i_and_mul_i_polynomials(0);
        assert_eq!(
            add_poly0
                .evals_slice()
                .iter()
                .filter(|&&x| x == Fq::ONE)
                .count(),
            1
        );
        assert_eq!(
            mul_poly0
                .evals_slice()
                .iter()
                .filter(|&&x| x == Fq::ONE)
                .count(),
            0
        );
    }
}
