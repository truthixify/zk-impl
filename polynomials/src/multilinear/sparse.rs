use ark_ff::PrimeField;
use std::iter::{Product, Sum};
use std::ops::{Add, Mul};

#[derive(Debug, Clone, PartialEq)]
pub struct SparseMultilinearPolynomial<F: PrimeField> {
    pub terms: Vec<(F, usize)>,
    pub n_vars: usize,
}

impl<F: PrimeField> SparseMultilinearPolynomial<F> {
    pub fn new(terms: Vec<(F, usize)>, n_vars: usize) -> Self {
        Self { terms, n_vars }
    }

    pub fn n_vars(&self) -> usize {
        self.n_vars
    }

    pub fn scalar_mul(&self, scalar: F) -> Self {
        let new_terms = self
            .terms
            .iter()
            .map(|(coeff, monomial_index)| (coeff.mul(scalar), *monomial_index))
            .collect();

        Self {
            terms: new_terms,
            n_vars: self.n_vars,
        }
    }

    pub fn evaluate(&self, point: &[F]) -> F {
        assert_eq!(
            point.len(),
            self.n_vars,
            "point length must be equal to n_vars"
        );

        self.terms
            .iter()
            .map(|(coeff, monomial_index)| {
                let mut result = F::ONE;
                for i in 0..self.n_vars {
                    if monomial_index & (1 << i) != 0 {
                        result = result.mul(point[i]);
                    }
                }
                coeff.mul(result)
            })
            .sum()
    }

    pub fn partial_evaluate(&self, partial_terms: &[(F, usize)]) -> Self {
        let mut new_terms = Vec::new();

        for (coeff, monomial_index) in &self.terms {
            let mut new_coeff = *coeff;
            let mut new_monomial_index = *monomial_index;
            for &(partial_coeff, partial_monomial_index) in partial_terms {
                if monomial_index & (1 << partial_monomial_index) != 0 {
                    new_coeff = new_coeff.mul(partial_coeff);
                    new_monomial_index &= !(1 << partial_monomial_index);
                }
            }
            new_terms.push((new_coeff, new_monomial_index));
        }

        SparseMultilinearPolynomial::new(new_terms, self.n_vars)
    }

    fn basis(point: &[u8], val: F) -> Self {
        let n_vars = point.len();
        let mut poly = SparseMultilinearPolynomial::new(vec![(val, 0)], n_vars);

        for (j, &x) in point.iter().enumerate() {
            let basis_term = if x == 0 {
                // (1 - x_j) = 1 - x_j = represented by (1, 0) + (-1, 1 << j)
                SparseMultilinearPolynomial::new(vec![(F::ONE, 0), (-F::ONE, 1 << j)], n_vars)
            } else {
                // x_j = represented by (1, 1 << j)
                SparseMultilinearPolynomial::new(vec![(F::ONE, 1 << j)], n_vars)
            };

            poly = &poly * &basis_term;
        }

        poly
    }

    pub fn interpolate(points: &[Vec<u8>], values: &[F]) -> Self {
        assert_eq!(points.len(), values.len());

        let n_vars = points[0].len();

        let mut interpolated_polynomial = SparseMultilinearPolynomial::new(vec![], n_vars);

        for (i, point) in points.iter().enumerate() {
            interpolated_polynomial = &interpolated_polynomial + &Self::basis(point, values[i]);
        }

        interpolated_polynomial
    }
}

impl<F: PrimeField> Add for &SparseMultilinearPolynomial<F> {
    type Output = SparseMultilinearPolynomial<F>;

    fn add(self, rhs: Self) -> Self::Output {
        assert_eq!(self.n_vars, rhs.n_vars, "n_vars must be equal");

        let mut summed_terms = Vec::new();
        let mut i = 0;
        let mut j = 0;

        while i < self.terms.len() && j < rhs.terms.len() {
            let (coeff1, monomial_index1) = self.terms[i];
            let (coeff2, monomial_index2) = rhs.terms[j];

            if monomial_index1 == monomial_index2 {
                summed_terms.push((coeff1.add(coeff2), monomial_index1));
                i += 1;
                j += 1;
            } else if monomial_index1 < monomial_index2 {
                summed_terms.push((coeff1, monomial_index1));
                i += 1;
            } else {
                summed_terms.push((coeff2, monomial_index2));
                j += 1;
            }
        }

        while i < self.terms.len() {
            summed_terms.push(self.terms[i]);
            i += 1;
        }

        while j < rhs.terms.len() {
            summed_terms.push(rhs.terms[j]);
            j += 1;
        }

        summed_terms.retain(|&(coeff, _)| coeff != F::ZERO);

        SparseMultilinearPolynomial::new(summed_terms, self.n_vars)
    }
}

impl<F: PrimeField> Sum for SparseMultilinearPolynomial<F> {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(SparseMultilinearPolynomial::new(vec![], 0), |acc, x| {
            &acc + &x
        })
    }
}

impl<F: PrimeField> Mul for &SparseMultilinearPolynomial<F> {
    type Output = SparseMultilinearPolynomial<F>;
    fn mul(self, rhs: Self) -> Self::Output {
        let mut product_polynomial =
            SparseMultilinearPolynomial::new(vec![(F::ZERO, 0)], self.n_vars);

        for (coeff1, monomial_index1) in &self.terms {
            for (coeff2, monomial_index2) in &rhs.terms {
                assert!(
                    monomial_index1 & monomial_index2 == 0,
                    "monomial indices must not overlap"
                );
                let new_coeff = coeff1.mul(*coeff2);
                let new_monomial_index = monomial_index1 | monomial_index2;

                let poly = SparseMultilinearPolynomial::new(
                    vec![(new_coeff, new_monomial_index)],
                    self.n_vars,
                );

                product_polynomial = &product_polynomial + &poly;
            }
        }

        product_polynomial
            .terms
            .sort_by_key(|&(_, monomial_index)| monomial_index);

        product_polynomial
    }
}

impl<F: PrimeField> Product for SparseMultilinearPolynomial<F> {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(SparseMultilinearPolynomial::new(vec![], 0), |acc, x| {
            &acc * &x
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::Fq;

    fn fq(val: u64) -> Fq {
        Fq::from(val)
    }

    #[test]
    fn test_scalar_mul() {
        let n_vars = 2;

        // f(x, y) = 3xy + 2x + 1
        let poly = SparseMultilinearPolynomial::new(
            vec![(fq(3), 0b11), (fq(2), 0b01), (fq(1), 0b00)],
            n_vars,
        );

        let scalar = fq(4);
        let scaled_poly = poly.scalar_mul(scalar);

        // Expected terms after multiplying by 4:
        let expected_terms = vec![(fq(12), 0b11), (fq(8), 0b01), (fq(4), 0b00)];
        let expected_poly = SparseMultilinearPolynomial::new(expected_terms, 2);

        assert_eq!(scaled_poly, expected_poly);
        assert_eq!(scaled_poly.n_vars, 2);
    }

    #[test]
    fn test_addition_with_overlap() {
        let n_vars = 2;

        let poly1 = SparseMultilinearPolynomial::new(vec![(fq(3), 0b01)], n_vars);
        let poly2 = SparseMultilinearPolynomial::new(vec![(fq(2), 0b01)], n_vars);
        let expected = SparseMultilinearPolynomial::new(vec![(fq(5), 0b01)], n_vars);

        assert_eq!(&poly1 + &poly2, expected);
    }

    #[test]
    #[should_panic(expected = "n_vars must be equal")]
    fn test_addition_mismatched_n_vars_panics() {
        let poly1 = SparseMultilinearPolynomial::new(vec![(fq(1), 0b00)], 2);
        let poly2 = SparseMultilinearPolynomial::new(vec![(fq(1), 0b00)], 3);
        let _ = &poly1 + &poly2;
    }

    #[test]
    fn test_multiplication() {
        let n_vars = 2;

        // 2x + y
        let poly1 = SparseMultilinearPolynomial::new(vec![(fq(2), 0b001), (fq(1), 0b010)], n_vars);
        // 4z
        let poly2 = SparseMultilinearPolynomial::new(vec![(fq(4), 0b100)], 2);
        // (2x + y)*4z = 8xz + 4yz
        let expected =
            SparseMultilinearPolynomial::new(vec![(fq(8), 0b101), (fq(4), 0b110)], n_vars);

        assert_eq!(&poly1 * &poly2, expected);
    }

    #[test]
    fn test_evaluate() {
        // f(x, y) = 3xy + 2x + 4
        let poly = SparseMultilinearPolynomial::new(
            vec![
                (fq(3), 0b11), // 3xy
                (fq(2), 0b01), // 2x
                (fq(4), 0b00), // 4
            ],
            2,
        );
        let point = vec![fq(2), fq(5)];
        // 3*2*5 + 2*2 + 4 = 30 + 4 + 4 = 38
        assert_eq!(poly.evaluate(&point), fq(38));
    }

    #[test]
    fn test_partial_evaluate() {
        let n_vars = 3;

        // f(x, y, z) = 2xyz + 5xz + 1
        let poly = SparseMultilinearPolynomial::new(
            vec![
                (fq(2), 0b111), // xyz
                (fq(5), 0b101), // xz
                (fq(1), 0b000), // constant
            ],
            3,
        );

        // Set y = 3 (index 1), z = 2 (index 2)
        let partial_terms = vec![(fq(3), 1), (fq(2), 2)];
        let partially_evaluated = poly.partial_evaluate(&partial_terms);

        // Expected:
        // 2xyz → 2*3*2 * x = 12x => (12, 0b001)
        // 5xz → 5*2 * x = 10x => (10, 0b001)
        // 1 stays => (1, 0b000)
        let expected_terms = vec![(fq(12), 0b001), (fq(10), 0b001), (fq(1), 0b000)];
        let expected_poly = SparseMultilinearPolynomial::new(expected_terms, n_vars);

        assert_eq!(partially_evaluated, expected_poly);
    }

    #[test]
    fn test_partial_evaluate_all_vars() {
        let n_vars = 2;

        // f(x, y) = 3xy + 2x + 4
        let poly = SparseMultilinearPolynomial::new(
            vec![(fq(3), 0b11), (fq(2), 0b01), (fq(4), 0b00)],
            n_vars,
        );
        let partial_terms = vec![(fq(2), 0), (fq(3), 1)];
        let partially_evaluated = poly.partial_evaluate(&partial_terms);

        // 3xy => 3*2*3 = 18 → constant
        // 2x => 2*2 = 4 → constant
        // 4 → constant
        let expected_terms = vec![(fq(18), 0b00), (fq(4), 0b00), (fq(4), 0b00)];
        let expected_poly = SparseMultilinearPolynomial::new(expected_terms, n_vars);

        assert_eq!(partially_evaluated, expected_poly);
    }

    #[test]
    fn test_interpolate() {
        // Boolean hypercube points for 3 variables
        let points = vec![
            vec![0, 0, 0],
            vec![0, 0, 1],
            vec![0, 1, 0],
            vec![0, 1, 1],
            vec![1, 0, 0],
            vec![1, 0, 1],
            vec![1, 1, 0],
            vec![1, 1, 1],
        ];

        // Corresponding values of the function f = 2ab + 3bc
        let values = vec![fq(0), fq(0), fq(0), fq(3), fq(0), fq(0), fq(2), fq(5)];

        let poly = SparseMultilinearPolynomial::interpolate(&points, &values);
        let expected_terms = vec![(fq(2), 0b011), (fq(3), 0b110)];
        let expected_n_vars = 3;
        let expected_poly = SparseMultilinearPolynomial::new(expected_terms, expected_n_vars);

        assert_eq!(poly, expected_poly);
    }
}
