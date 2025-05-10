use crate::gate::{Gate, Op};
use ark_ff::PrimeField;
use polynomials::multilinear::MultilinearPolynomial;
use std::marker::PhantomData;

#[derive(Debug)]
pub struct Layer<F: PrimeField> {
    pub gates: Vec<Gate>,
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

    #[test]
    fn test_get_positional_index() {
        // layer_index = 2, binary lengths: output = 2, left = 3, right = 3
        let idx = get_positional_index(2, 2, 3, 4);

        // binary: output=10, left=011, right=100 => combined = 10011100
        // binary 10011100 = decimal 156
        assert_eq!(idx, 156);
    }
}
