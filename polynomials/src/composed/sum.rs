use super::product::ProductPolynomial;
use ark_ff::PrimeField;

pub struct SumPolynomial<F: PrimeField> {
    product_polynomials: Vec<ProductPolynomial<F>>,
}

impl<F: PrimeField> SumPolynomial<F> {
    pub fn new(product_polynomials: Vec<ProductPolynomial<F>>) -> Self {
        let n_vars = product_polynomials[0].polynomials[0].n_vars();

        assert!(
            product_polynomials.iter().all(|prod_poly| prod_poly
                .polynomials
                .iter()
                .all(|poly| poly.n_vars() == n_vars)),
            "All polynomials in sum polynomial must have the same number of variable"
        );

        Self {
            product_polynomials,
        }
    }

    pub fn degree(&self) -> usize {
        self.product_polynomials[0].degree()
    }

    pub fn evaluate(&self, points: &[F]) -> F {
        self.product_polynomials
            .iter()
            .map(|prod_poly| prod_poly.evaluate(points))
            .sum()
    }

    pub fn partial_evaluate(&self, points: &[(F, usize)]) -> Self {
        let product_polynomials = self
            .product_polynomials
            .iter()
            .map(|prod_poly| prod_poly.partial_evaluate(points))
            .collect();

        Self::new(product_polynomials)
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        self.product_polynomials
            .iter()
            .flat_map(|prod_poly| prod_poly.to_bytes())
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::multilinear::MultilinearPolynomial;
    use ark_bls12_381::Fq;

    fn fq(val: u64) -> Fq {
        Fq::from(val)
    }

    fn create_multilinear_poly(values: &[u64]) -> MultilinearPolynomial<Fq> {
        let evals = values.iter().copied().map(fq).collect();
        MultilinearPolynomial::new(evals)
    }

    fn create_product_poly(polys: &[&[u64]]) -> ProductPolynomial<Fq> {
        let multilinears = polys
            .iter()
            .map(|vals| create_multilinear_poly(vals))
            .collect();
        ProductPolynomial::new(multilinears)
    }

    #[test]
    fn test_new_valid_sum_poly() {
        let prod1 = create_product_poly(&[&[1, 2, 3, 4], &[5, 6, 7, 8]]);
        let prod2 = create_product_poly(&[&[9, 10, 11, 12], &[13, 14, 15, 16]]);
        let sum_poly = SumPolynomial::new(vec![prod1, prod2]);

        assert_eq!(sum_poly.degree(), 2); // 2 polynomials in each product
    }

    #[test]
    #[should_panic(
        expected = "All polynomials in sum polynomial must have the same number of variable"
    )]
    fn test_new_mismatched_vars() {
        let prod1 = create_product_poly(&[&[1, 2, 3, 4]]);
        let prod2 = create_product_poly(&[&[1, 2, 3, 4, 5, 6, 7, 8]]); // 3 vars (2^3 = 8)

        SumPolynomial::new(vec![prod1, prod2]);
    }

    #[test]
    fn test_evaluate_sum_poly() {
        let prod1 = create_product_poly(&[&[1, 2, 3, 4]]); // degree 1
        let prod2 = create_product_poly(&[&[5, 6, 7, 8]]); // degree 1
        let sum_poly = SumPolynomial::new(vec![prod1.clone(), prod2.clone()]);
        let point = &[fq(1), fq(0)];
        let expected = prod1.evaluate(point) + prod2.evaluate(point);

        assert_eq!(sum_poly.evaluate(point), expected);
    }

    #[test]
    fn test_partial_evaluate_sum_poly() {
        let prod1 = create_product_poly(&[&[1, 2, 3, 4]]);
        let prod2 = create_product_poly(&[&[5, 6, 7, 8]]);
        let sum_poly = SumPolynomial::new(vec![prod1.clone(), prod2.clone()]);

        let partial = sum_poly.partial_evaluate(&[(fq(1), 0)]);

        for (original, reduced) in sum_poly
            .product_polynomials
            .iter()
            .zip(partial.product_polynomials.iter())
        {
            let expected = original.partial_evaluate(&[(fq(1), 0)]);

            assert_eq!(*reduced, expected);
        }
    }

    #[test]
    fn test_to_bytes_sum_poly() {
        let prod1 = create_product_poly(&[&[1, 2, 3, 4]]);
        let prod2 = create_product_poly(&[&[5, 6, 7, 8]]);
        let sum_poly = SumPolynomial::new(vec![prod1.clone(), prod2.clone()]);

        let expected: Vec<u8> = prod1
            .to_bytes()
            .into_iter()
            .chain(prod2.to_bytes())
            .collect();

        assert_eq!(sum_poly.to_bytes(), expected);
    }
}
