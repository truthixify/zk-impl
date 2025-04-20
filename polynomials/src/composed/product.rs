use crate::multilinear::MultilinearPolynomial;
use ark_ff::PrimeField;

#[derive(Debug, Clone, PartialEq)]
pub struct ProductPolynomial<F: PrimeField> {
    pub polynomials: Vec<MultilinearPolynomial<F>>,
}

impl<F: PrimeField> ProductPolynomial<F> {
    pub fn new(polynomials: Vec<MultilinearPolynomial<F>>) -> Self {
        let n_vars = polynomials[0].n_vars();

        assert!(
            polynomials.iter().all(|poly| poly.n_vars() == n_vars),
            "All polynomials in product polynomial must have the same number of variable"
        );

        Self { polynomials }
    }

    pub fn degree(&self) -> usize {
        self.polynomials.len()
    }

    pub fn evaluate(&self, points: &[F]) -> F {
        self.polynomials
            .iter()
            .map(|poly| poly.evaluate(points))
            .product()
    }

    pub fn partial_evaluate_many_vars(&self, points: &[(F, usize)]) -> Self {
        let polynomials = self
            .polynomials
            .iter()
            .map(|poly| poly.partial_evaluate_many_vars(points))
            .collect();

        Self::new(polynomials)
    }

    pub fn partial_evaluate(&self, point: F, var_index: usize) -> Self {
        self.partial_evaluate_many_vars(&[(point, var_index)])
    }

    pub fn element_wise_mul(&self) -> MultilinearPolynomial<F> {
        assert!(
            self.polynomials.len() > 1,
            "At least two polynomials are needed for multiplication"
        );

        let init = self.polynomials[0].clone();

        self.polynomials
            .iter()
            .skip(1)
            .fold(init, |acc, curr| acc.tensor_mul(curr))
    }

    pub fn reduce(&self) -> Vec<F> {
        self.element_wise_mul().evals_slice().to_vec()
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        self.polynomials
            .iter()
            .flat_map(|poly| poly.to_bytes())
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::Fq;

    fn fq(val: u64) -> Fq {
        Fq::from(val)
    }

    fn create_multilinear_poly(evals: Vec<u64>) -> MultilinearPolynomial<Fq> {
        let coeffs: Vec<Fq> = evals.iter().copied().map(fq).collect();
        MultilinearPolynomial::new(coeffs)
    }

    #[test]
    fn test_new_valid_polys() {
        let polys = vec![
            create_multilinear_poly(vec![1, 2, 3, 4]),
            create_multilinear_poly(vec![5, 6, 7, 8]),
        ];
        let pp = ProductPolynomial::new(polys);

        assert_eq!(pp.degree(), 2);
    }

    #[test]
    #[should_panic(
        expected = "All polynomials in product polynomial must have the same number of variable"
    )]
    fn test_new_inconsistent_vars() {
        let p1 = create_multilinear_poly(vec![1, 2, 3, 4]);
        let p2 = create_multilinear_poly(vec![1, 2, 3, 4, 5, 6, 7, 8]);
        ProductPolynomial::new(vec![p1, p2]);
    }

    #[test]
    fn test_evaluate() {
        let p1 = create_multilinear_poly(vec![1, 2, 3, 4]);
        let p2 = create_multilinear_poly(vec![5, 6, 7, 8]);
        let pp = ProductPolynomial::new(vec![p1.clone(), p2.clone()]);
        let point = vec![fq(1), fq(0)]; // Evaluate both polys at (1, 0)
        let expected = p1.evaluate(&point) * p2.evaluate(&point);
        let result = pp.evaluate(&point);

        assert_eq!(expected, result);
    }

    #[test]
    fn test_partial_evaluate() {
        let p1 = create_multilinear_poly(vec![1, 2, 3, 4]);
        let p2 = create_multilinear_poly(vec![5, 6, 7, 8]);
        let pp = ProductPolynomial::new(vec![p1.clone(), p2.clone()]);

        let partial_points = vec![(fq(1), 0)]; // Fix variable x₀ = 1
        let evaluated_pp = pp.partial_evaluate_many_vars(&partial_points);

        for (original, evaluated) in [p1, p2].iter().zip(evaluated_pp.polynomials.iter()) {
            let expected = original.partial_evaluate_many_vars(&partial_points);

            assert_eq!(expected, *evaluated);
        }
    }

    #[test]
    fn test_element_wise_mul() {
        let poly1 = create_multilinear_poly(vec![1, 2, 3, 4]);
        let poly2 = create_multilinear_poly(vec![2, 3, 4, 5]);
        let poly3 = create_multilinear_poly(vec![1, 1, 1, 1]);

        let product = ProductPolynomial::new(vec![poly1.clone(), poly2.clone(), poly3.clone()]);
        let result = product.element_wise_mul();

        // Multiply all three polynomials element-wise
        let expected = create_multilinear_poly(vec![1 * 2 * 1, 2 * 3 * 1, 3 * 4 * 1, 4 * 5 * 1]); // [2, 6, 12, 20]

        assert_eq!(result, expected);
    }

    #[test]
    #[should_panic(expected = "At least two polynomials are needed for multiplication")]
    fn test_element_wise_mul_panics_on_single_poly() {
        let poly = create_multilinear_poly(vec![1, 2, 3, 4]);
        let product = ProductPolynomial::new(vec![poly]);

        product.element_wise_mul(); // Should panic
    }

    #[test]
    fn test_to_bytes() {
        let p1 = create_multilinear_poly(vec![1, 2, 3, 4]);
        let p2 = create_multilinear_poly(vec![5, 6, 7, 8]);
        let pp = ProductPolynomial::new(vec![p1.clone(), p2.clone()]);
        let expected_bytes: Vec<u8> = p1.to_bytes().into_iter().chain(p2.to_bytes()).collect();

        assert_eq!(pp.to_bytes(), expected_bytes);
    }
}
