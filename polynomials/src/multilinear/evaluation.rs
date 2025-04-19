use ark_ff::{BigInteger, PrimeField};

#[derive(Clone, Debug, PartialEq)]
pub struct MultilinearPolynomial<F: PrimeField> {
    evals: Vec<F>,
}

impl<F: PrimeField> MultilinearPolynomial<F> {
    pub fn new(evals: Vec<F>) -> Self {
        assert!(
            evals.len().is_power_of_two(),
            "Number of evaluations must be a power of two"
        );

        Self { evals }
    }

    pub fn n_vars(&self) -> usize {
        self.evals.len().ilog2() as usize
    }

    pub fn evals_slice(&self) -> &[F] {
        &self.evals
    }

    pub fn scalar_mul(&self, scalar: F) -> Self {
        Self {
            evals: self.evals.iter().map(|&x| x * scalar).collect(),
        }
    }

    pub fn evaluate(&self, points: &[F]) -> F {
        assert_eq!(
            points.len(),
            self.n_vars(),
            "Number of points must match number of variables"
        );

        self.partial_evaluate_many_vars(&points.iter().map(|&x| (x, 0)).collect::<Vec<_>>())
            .evals[0]
    }

    pub fn partial_evaluate(&self, point: F, var_index: usize) -> Self {
        self.partial_evaluate_many_vars(&[(point, var_index)])
    }

    pub fn partial_evaluate_many_vars(&self, points: &[(F, usize)]) -> Self {
        assert!(
            points.len() <= self.n_vars(),
            "Number of points must be less than or equal to number of variables"
        );

        let mut evals = self.evals.clone();
        let mut current_n_vars = self.n_vars();

        let mut points_sorted = points.to_vec();
        points_sorted.sort_by_key(|&(_, idx)| std::cmp::Reverse(idx));

        for &(value, var_index) in &points_sorted {
            assert!(
                var_index < current_n_vars,
                "Variable index {} out of bounds (max {})",
                var_index,
                current_n_vars
            );

            // For fixing variable at `var_index`, we collapse the dimension corresponding to that variable.
            // The evaluations are ordered lexicographically, so fixing a variable means interpolating
            // between pairs of points that differ only in that variable.
            //
            // `stride` is the distance between those two points in the evaluation vector.
            // `chunk_size` defines the region where those pairs appear repeatedly.
            //
            // For example, if we fix variable `i` in an `n`-var polynomial,
            //   then stride = 2^(n - i - 1)
            // This ensures that for each fixed variable, we walk over the correct (a, b) pairs
            // corresponding to that dimension in the hypercube, and interpolate across them.
            let stride = 1 << (current_n_vars - var_index - 1);
            let chunk_size = stride << 1; // 2 chunks of size stride (stride << 1 = stride * 2)
            let mut new_evals = Vec::with_capacity(evals.len() / 2);

            // impl 1
            // this was faster in benchmarks even though it had one more loop
            for chunk in evals.chunks(chunk_size) {
                for i in 0..stride {
                    let y1 = chunk[i];
                    let y2 = chunk[i + stride];
                    // linear interpolation: (1 - x) * a + x * b = a + (b - a) * x
                    let term = if value.is_zero() {
                        y1
                    } else if value.is_one() {
                        y2
                    } else {
                        y1 + (y2 - y1) * value
                    };

                    new_evals.push(term);
                }
            }

            // impl 2
            // let mut i = 0;

            // while i < evals.len() {
            //     let y1 = evals[i];
            //     let y2 = evals[i + stride];

            //     // linear interpolation: (1 - x) * a + x * b = a + (b - a) * x
            //     new_evals.push(y1 + (y2 - y1) * value);

            //     i += 1;

            //     if i % chunk_size == stride {
            //         i += stride;
            //     }
            // }

            evals = new_evals;
            current_n_vars -= 1;
        }

        MultilinearPolynomial::new(evals)
    }

    pub fn tensor_add(&self, other: &Self) -> Self {
        assert_eq!(
            self.evals.len(),
            other.evals.len(),
            "Polynomials must have the same number of evaluations"
        );

        let evals = self
            .evals
            .iter()
            .zip(other.evals.iter())
            .map(|(x, y)| *x + *y)
            .collect();

        Self { evals }
    }

    pub fn tensor_mul(&self, other: &Self) -> Self {
        assert_eq!(
            self.evals.len(),
            other.evals.len(),
            "Polynomials must have the same number of evaluations"
        );

        let evals = self
            .evals
            .iter()
            .zip(other.evals.iter())
            .map(|(x, y)| *x * *y)
            .collect();

        Self { evals }
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        self.evals
            .iter()
            .flat_map(|el| el.into_bigint().to_bytes_be())
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::Fq;
    use ark_ff::UniformRand;

    fn fq(x: u64) -> Fq {
        Fq::from(x)
    }

    #[test]
    fn test_evaluate() {
        let evaluated_values = vec![Fq::from(0), Fq::from(0), Fq::from(3), Fq::from(8)];
        let polynomial = MultilinearPolynomial::new(evaluated_values);
        let values = vec![Fq::from(6), Fq::from(2)];

        assert_eq!(polynomial.evaluate(&values), Fq::from(78));
    }

    #[test]
    fn test_partial_evaluate() {
        // 3-variable polynomial over (a, b, c)
        // Evaluation order: (a, b, c) in lex: 000, 001, 010, 011, 100, 101, 110, 111
        let poly_abc = MultilinearPolynomial::new(vec![
            fq(1),
            fq(3),
            fq(5),
            fq(7),
            fq(2),
            fq(4),
            fq(6),
            fq(8),
        ]);

        // Fix c = 0 → reduces to polynomial over (a, b)
        // Retains values at c=0: indices 0, 2, 4, 6
        assert_eq!(
            poly_abc.partial_evaluate(fq(0), 2),
            MultilinearPolynomial::new(vec![fq(1), fq(5), fq(2), fq(6)])
        );

        // Fix b = 1 → reduces to polynomial over (a, c)
        // Indices with b=1: (0,1,0)=2, (0,1,1)=3, (1,1,0)=6, (1,1,1)=7 → [5, 7, 6, 8]
        assert_eq!(
            poly_abc.partial_evaluate(fq(1), 1),
            MultilinearPolynomial::new(vec![fq(5), fq(7), fq(6), fq(8)])
        );

        // Fix a = 1, c = 1 → indices: (1,0,1)=5, (1,1,1)=7 → [4, 8]
        assert_eq!(
            poly_abc.partial_evaluate_many_vars(&[(fq(1), 0), (fq(1), 2)]),
            MultilinearPolynomial::new(vec![fq(4), fq(8)])
        );

        // Fix a = 0, b = 0 → (0,0,0)=1, (0,0,1)=3
        assert_eq!(
            poly_abc.partial_evaluate_many_vars(&[(fq(0), 0), (fq(0), 1)]),
            MultilinearPolynomial::new(vec![fq(1), fq(3)])
        );

        // Fix a = 0, b = 1, c = 0 → single point: (0,1,0)=2 → value 5
        assert_eq!(
            poly_abc.partial_evaluate_many_vars(&[(fq(0), 0), (fq(1), 1), (fq(0), 2)]),
            MultilinearPolynomial::new(vec![fq(5)])
        );

        // -------------------------------
        // New test: 4-variable polynomial over (a, b, c, d)
        // Evaluation order: (a, b, c, d) lex order → 16 entries
        let poly_abcd = MultilinearPolynomial::new(vec![
            fq(1),
            fq(2),
            fq(3),
            fq(4),
            fq(5),
            fq(6),
            fq(7),
            fq(8),
            fq(9),
            fq(10),
            fq(11),
            fq(12),
            fq(13),
            fq(14),
            fq(15),
            fq(16),
        ]);

        // Fix d = 0 → keep every even index
        assert_eq!(
            poly_abcd.partial_evaluate(fq(0), 3),
            MultilinearPolynomial::new(vec![
                fq(1),
                fq(3),
                fq(5),
                fq(7),
                fq(9),
                fq(11),
                fq(13),
                fq(15)
            ])
        );

        // Fix b = 1, c = 0 → pick (a,1,0,d): indices 4,5,12,13 → values [5,6,13,14]
        assert_eq!(
            poly_abcd.partial_evaluate_many_vars(&[(fq(1), 1), (fq(0), 2)]),
            MultilinearPolynomial::new(vec![fq(5), fq(6), fq(13), fq(14)])
        );

        // Fix a = 0, b = 1, c = 1 → indices: (0,1,1,0)=6, (0,1,1,1)=7 → [7,8]
        assert_eq!(
            poly_abcd.partial_evaluate_many_vars(&[(fq(0), 0), (fq(1), 1), (fq(1), 2)]),
            MultilinearPolynomial::new(vec![fq(7), fq(8)])
        );

        // Fix all: a=1, b=0, c=1, d=0 → (1,0,1,0) = index 10 → value 11
        assert_eq!(
            poly_abcd.partial_evaluate_many_vars(&[(fq(1), 0), (fq(0), 1), (fq(1), 2), (fq(0), 3)]),
            MultilinearPolynomial::new(vec![fq(11)])
        );
    }

    #[test]
    fn test_new_and_n_vars() {
        let evals = vec![fq(0), fq(1), fq(2), fq(3)];
        let poly = MultilinearPolynomial::new(evals.clone());

        assert_eq!(poly.evals, evals);
        assert_eq!(poly.n_vars(), 2);
    }

    #[test]
    #[should_panic(expected = "Number of evaluations must be a power of two")]
    fn test_new_invalid_length() {
        let evals = vec![fq(0), fq(1), fq(2)];
        let _ = MultilinearPolynomial::new(evals); // Should panic
    }

    #[test]
    fn test_scalar_mul() {
        let evals = vec![fq(1), fq(2), fq(3), fq(4)];
        let poly = MultilinearPolynomial::new(evals);
        let result = poly.scalar_mul(fq(2));
        assert_eq!(result.evals, vec![fq(2), fq(4), fq(6), fq(8)]);
    }

    #[test]
    fn test_evaluate_constant_polynomial() {
        let evals = vec![fq(7), fq(7)];
        let poly = MultilinearPolynomial::new(evals);
        assert_eq!(poly.evaluate(&[fq(0)]), fq(7));
        assert_eq!(poly.evaluate(&[fq(1)]), fq(7));
    }

    #[test]
    fn test_partial_evaluation_stepwise() {
        let evals = vec![fq(1), fq(2), fq(3), fq(4), fq(5), fq(6), fq(7), fq(8)];

        let poly = MultilinearPolynomial::new(evals);

        let a = fq(1);
        let b = fq(0);
        let c = fq(1);

        // Full evaluation directly
        let full_eval = poly.evaluate(&[a, b, c]);

        // Step-by-step partial evaluation
        let partial1 = poly.partial_evaluate(a, 0);
        let partial2 = partial1.partial_evaluate(b, 0);
        let final_eval = partial2.evaluate(&[c]);

        assert_eq!(full_eval, final_eval);
    }

    #[test]
    fn test_partial_evaluate_randomized() {
        let mut rng = rand::thread_rng();
        let num_vars = 10;
        let num_evals = 1 << num_vars;

        let evals: Vec<Fq> = (0..num_evals).map(|_| Fq::rand(&mut rng)).collect();
        let poly = MultilinearPolynomial::new(evals);

        let assignment: Vec<Fq> = (0..num_vars).map(|_| Fq::rand(&mut rng)).collect();
        let full_eval = poly.evaluate(&assignment);

        for fixed_vars in 1..=num_vars {
            let fixed: Vec<(Fq, usize)> = (0..fixed_vars).map(|i| (assignment[i], i)).collect();

            let partial = poly.partial_evaluate_many_vars(&fixed);
            let remaining_assignment: Vec<Fq> = (0..num_vars)
                .filter(|i| i >= &fixed_vars)
                .map(|i| assignment[i])
                .collect();
            let final_eval = partial.evaluate(&remaining_assignment);

            assert_eq!(full_eval, final_eval);
        }
    }

    #[test]
    #[should_panic(expected = "Number of points must match number of variables")]
    fn test_evaluate_invalid_input_length() {
        let evals = vec![fq(1), fq(2), fq(3), fq(4)];
        let poly = MultilinearPolynomial::new(evals);
        let _ = poly.evaluate(&[fq(1)]); // Not enough points
    }

    #[test]
    fn test_tensor_add() {
        let poly1 = MultilinearPolynomial::new(vec![fq(1), fq(2), fq(3), fq(4)]);
        let poly2 = MultilinearPolynomial::new(vec![fq(5), fq(6), fq(7), fq(8)]);
        let result = poly1.tensor_add(&poly2);

        assert_eq!(result.evals, vec![fq(6), fq(8), fq(10), fq(12)]);
    }

    #[test]
    #[should_panic(expected = "Polynomials must have the same number of evaluations")]
    fn test_tensor_add_invalid_length() {
        let poly1 = MultilinearPolynomial::new(vec![fq(1), fq(2), fq(3), fq(4)]);
        let poly2 = MultilinearPolynomial::new(vec![fq(5), fq(6)]);
        let _ = poly1.tensor_add(&poly2);
    }

    #[test]
    fn test_tensor_product() {
        let poly1 = MultilinearPolynomial::new(vec![fq(1), fq(2), fq(3), fq(4)]);
        let poly2 = MultilinearPolynomial::new(vec![fq(5), fq(6), fq(7), fq(8)]);
        let result = poly1.tensor_mul(&poly2);

        assert_eq!(result.evals, vec![fq(5), fq(12), fq(21), fq(32)]);
    }

    #[test]
    #[should_panic(expected = "Polynomials must have the same number of evaluations")]
    fn test_tensor_mul_invalid_length() {
        let poly1 = MultilinearPolynomial::new(vec![fq(1), fq(2), fq(3), fq(4)]);
        let poly2 = MultilinearPolynomial::new(vec![fq(5), fq(6)]);
        let _ = poly1.tensor_mul(&poly2);
    }
}
