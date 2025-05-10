use ark_ff::PrimeField;

#[derive(Debug)]
pub enum Op {
    Add,
    Mul,
}

#[derive(Debug)]
pub struct Gate {
    pub op: Op,
    pub output: usize,
    pub left_index: usize,
    pub right_index: usize,
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
