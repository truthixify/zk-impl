use ark_ff::PrimeField;
use std::ops::{Add, Mul};

#[derive(Debug, Clone, PartialEq)]
pub struct DenseMultilinearPolynomial<F: PrimeField> {
    coefficients: Vec<F>,
    n_vars: usize,
}

impl<F: PrimeField> DenseMultilinearPolynomial<F> {
    pub fn new(n_vars: usize) -> Self {
        let coefficients = vec![F::ZERO; 1 << n_vars];

        Self::new_with_coefficients(coefficients, n_vars)
    }

    pub fn new_with_coefficients(coefficients: Vec<F>, n_vars: usize) -> Self {
        Self {
            coefficients,
            n_vars,
        }
    }

    fn unit_poly(n_vars: usize) -> Self {
        let mut coeffs = vec![F::ZERO; 1 << n_vars];
        coeffs[0] = F::ONE;

        DenseMultilinearPolynomial::new_with_coefficients(coeffs, n_vars)
    }

    pub fn n_vars(&self) -> usize {
        self.n_vars
    }

    pub fn coefficients_slice(&self) -> &[F] {
        &self.coefficients
    }

    pub fn scalar_mul(&self, scalar: F) -> Self {
        let mut new_coeffs = Vec::with_capacity(self.coefficients.len());

        for coeff in &self.coefficients {
            new_coeffs.push(coeff.mul(scalar));
        }

        DenseMultilinearPolynomial::new_with_coefficients(new_coeffs, self.n_vars)
    }

    pub fn evaluate(&self, point: &[(F, u8)]) -> F {
        assert_eq!(
            point.len(),
            self.n_vars,
            "point length must be equal to number of variables"
        );

        self.coefficients
            .iter()
            .enumerate()
            .map(|(i, coeff)| {
                let mut result = F::ONE;

                for (val, pos) in point.iter() {
                    if i & (1 << pos) != 0 {
                        result = result.mul(val);
                    }
                }

                coeff.mul(result)
            })
            .sum()
    }

    pub fn partial_evaluate(&self, partial_terms: &[(F, usize)]) -> Self {
        let mut new_coeffs = vec![F::ZERO; 1 << (self.n_vars - partial_terms.len())];

        for (i, &coeff) in self.coefficients.iter().enumerate() {
            let mut new_coeff = coeff;
            let mut new_index = i;

            for &(val, pos) in partial_terms {
                if i & (1 << pos) != 0 {
                    new_coeff *= val;
                    new_index ^= 1 << pos; // Remove variable from index
                }
            }

            // Collapse monomial into lower-dimensional space
            let mut collapsed_index = 0;
            let mut shift = 0;

            for j in 0..self.n_vars {
                if !partial_terms.iter().any(|&(_, pos)| pos == j) {
                    if (new_index >> j) & 1 != 0 {
                        collapsed_index |= 1 << shift;
                    }

                    shift += 1;
                }
            }

            new_coeffs[collapsed_index] += new_coeff;
        }

        DenseMultilinearPolynomial::new_with_coefficients(
            new_coeffs,
            self.n_vars - partial_terms.len(),
        )
    }

    fn basis(point: &[u8]) -> Self {
        let n_vars = point.len();
        let mut poly = Self::unit_poly(n_vars);

        for (i, &x) in point.iter().enumerate() {
            let mut coeffs = vec![F::ZERO; 1 << n_vars];
            let basis_term = if x == 0 {
                // (1 - x_j) = 1 - x_j
                coeffs[0] = F::ONE;
                coeffs[1 << i] = -F::ONE;
                DenseMultilinearPolynomial::new_with_coefficients(coeffs, n_vars)
            } else {
                // x_j
                coeffs[1 << i] = F::ONE;
                DenseMultilinearPolynomial::new_with_coefficients(coeffs, n_vars)
            };

            poly = poly * basis_term;
        }

        poly
    }

    pub fn interpolate(points: &[Vec<u8>], values: &[F]) -> Self {
        assert_eq!(points.len(), values.len());

        let n_vars = points[0].len();

        let mut interpolated_polynomial = DenseMultilinearPolynomial::new(n_vars);

        for (i, point) in points.iter().enumerate() {
            interpolated_polynomial =
                &interpolated_polynomial + &Self::basis(point).scalar_mul(values[i]);
        }

        interpolated_polynomial
    }
}

impl<F: PrimeField> Add for DenseMultilinearPolynomial<F> {
    type Output = DenseMultilinearPolynomial<F>;

    fn add(self, rhs: Self) -> Self::Output {
        assert_eq!(
            self.n_vars, rhs.n_vars,
            "polynomials must have the same number of variables"
        );

        let new_coeffs = self
            .coefficients
            .iter()
            .zip(rhs.coefficients.iter())
            .map(|(a, b)| *a + *b)
            .collect::<Vec<_>>();

        DenseMultilinearPolynomial::new_with_coefficients(new_coeffs, self.n_vars)
    }
}

impl<F: PrimeField> Add for &DenseMultilinearPolynomial<F> {
    type Output = DenseMultilinearPolynomial<F>;

    fn add(self, rhs: Self) -> Self::Output {
        assert_eq!(
            self.n_vars, rhs.n_vars,
            "polynomials must have the same number of variables"
        );

        let new_coeffs = self
            .coefficients
            .iter()
            .zip(rhs.coefficients.iter())
            .map(|(a, b)| *a + *b)
            .collect::<Vec<_>>();

        DenseMultilinearPolynomial::new_with_coefficients(new_coeffs, self.n_vars)
    }
}

impl<F: PrimeField> Mul for DenseMultilinearPolynomial<F> {
    type Output = DenseMultilinearPolynomial<F>;

    fn mul(self, rhs: Self) -> Self::Output {
        assert_eq!(
            self.n_vars, rhs.n_vars,
            "polynomials must have the same number of variables"
        );

        let mut product_polynomial = DenseMultilinearPolynomial::new(self.n_vars);

        for (index1, coeff1) in self.coefficients.iter().enumerate() {
            let mut coeffs = vec![F::ZERO; 1 << self.n_vars];

            for (index2, coeff2) in rhs.coefficients.iter().enumerate() {
                if !coeff1.is_zero() && !coeff2.is_zero() {
                    assert!(index1 & index2 == 0, "monomial indices must not overlap");
                }

                let new_index = index1 | index2;
                let new_coeff = *coeff1 * coeff2;

                coeffs[new_index] += new_coeff;
            }

            let poly = DenseMultilinearPolynomial::new_with_coefficients(coeffs, self.n_vars);

            product_polynomial = product_polynomial + poly;
        }

        product_polynomial
    }
}

impl<F: PrimeField> Mul for &DenseMultilinearPolynomial<F> {
    type Output = DenseMultilinearPolynomial<F>;

    fn mul(self, rhs: Self) -> Self::Output {
        assert_eq!(
            self.n_vars, rhs.n_vars,
            "polynomials must have the same number of variables"
        );

        let mut product_polynomial = DenseMultilinearPolynomial::new(self.n_vars);

        for (index1, coeff1) in self.coefficients.iter().enumerate() {
            let mut coeffs = vec![F::ZERO; 1 << self.n_vars];

            for (index2, coeff2) in rhs.coefficients.iter().enumerate() {
                if !coeff1.is_zero() && !coeff2.is_zero() {
                    assert!(index1 & index2 == 0, "monomial indices must not overlap");
                }

                let new_index = index1 | index2;
                let new_coeff = *coeff1 * coeff2;

                coeffs[new_index] += new_coeff;
            }

            let poly = DenseMultilinearPolynomial::new_with_coefficients(coeffs, self.n_vars);

            product_polynomial = product_polynomial + poly;
        }

        product_polynomial
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
        // f(x, y) = 3xy + 2x + 1
        let poly =
            DenseMultilinearPolynomial::new_with_coefficients(vec![fq(1), fq(2), fq(0), fq(3)], 2);

        let scalar = fq(4);
        let scaled_poly = poly.scalar_mul(scalar);

        // Expected terms after multiplying by 4:
        let expected_terms = vec![fq(4), fq(8), fq(0), fq(12)];
        let expected_poly = DenseMultilinearPolynomial::new_with_coefficients(expected_terms, 2);

        assert_eq!(scaled_poly, expected_poly);
        assert_eq!(scaled_poly.n_vars, 2);
    }

    #[test]
    fn test_addition() {
        let n_vars = 2;

        // f(x, y) = 2x + 3y + 6xy
        let poly1 = DenseMultilinearPolynomial::new_with_coefficients(
            vec![fq(0), fq(2), fq(3), fq(6)],
            n_vars,
        );
        // g(x, y) = 1 + 17x
        let poly2 = DenseMultilinearPolynomial::new_with_coefficients(
            vec![fq(1), fq(17), fq(0), fq(0)],
            n_vars,
        );
        // f(x, y) = 2x + 3y + 6xy + 1 + 17x = 1 + 19x + 3y + 6xy
        let expected = DenseMultilinearPolynomial::new_with_coefficients(
            vec![fq(1), fq(19), fq(3), fq(6)],
            n_vars,
        );

        assert_eq!(poly1 + poly2, expected);
    }

    #[test]
    fn test_multiplication() {
        let n_vars = 3;

        // 2x + y
        let poly1 = DenseMultilinearPolynomial::new_with_coefficients(
            vec![fq(0), fq(2), fq(1), fq(0), fq(0), fq(0), fq(0), fq(0)],
            n_vars,
        );
        // 4z
        let poly2 = DenseMultilinearPolynomial::new_with_coefficients(
            vec![fq(0), fq(0), fq(0), fq(0), fq(4), fq(0), fq(0), fq(0)],
            n_vars,
        );
        // (2x + y)*4z = 8xz + 4yz
        let expected = DenseMultilinearPolynomial::new_with_coefficients(
            vec![fq(0), fq(0), fq(0), fq(0), fq(0), fq(8), fq(4), fq(0)],
            n_vars,
        );

        assert_eq!(poly1 * poly2, expected);
    }

    #[test]
    fn test_evaluate() {
        // f(x, y) = 3xy + 2x + 4
        let poly = DenseMultilinearPolynomial::new_with_coefficients(
            vec![
                fq(4), // 4
                fq(2), // 2x
                fq(0), // 0y
                fq(3), // 3xy
            ],
            2,
        );
        // Evaluate at x = 2, y = 5
        let point = vec![(fq(2), 0), (fq(5), 1)];
        // 3*2*5 + 2*2 + 4 = 30 + 4 + 4 = 38
        assert_eq!(poly.evaluate(&point), fq(38));
    }

    #[test]
    fn test_partial_evaluate() {
        // f(x, y, z) = 2xyz + 5xz + 1
        let poly = DenseMultilinearPolynomial::new_with_coefficients(
            vec![
                fq(1), // 1
                fq(0), // 0x
                fq(0), // 0y
                fq(0), // 2xy
                fq(0), // 0z
                fq(5), // 5xz
                fq(0), // 0yz
                fq(2), // 0xyz
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
        let expected_terms = vec![fq(1), fq(22)];
        let expected_poly = DenseMultilinearPolynomial::new_with_coefficients(expected_terms, 1);

        assert_eq!(partially_evaluated, expected_poly);
    }

    #[test]
    fn test_partial_evaluate_all_vars() {
        // f(x, y) = 3xy + 2x + 4
        let poly =
            DenseMultilinearPolynomial::new_with_coefficients(vec![fq(4), fq(2), fq(0), fq(3)], 2);
        let partial_terms = vec![(fq(2), 0), (fq(3), 1)];
        let partially_evaluated = poly.partial_evaluate(&partial_terms);

        // 3xy => 3*2*3 = 18 → constant
        // 2x => 2*2 = 4 → constant
        // 4 → constant
        // 18 + 4 + 4 = 26
        let expected_terms = vec![fq(26)];
        let expected_poly = DenseMultilinearPolynomial::new_with_coefficients(expected_terms, 0);

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

        let poly = DenseMultilinearPolynomial::interpolate(&points, &values);
        let expected_terms = vec![fq(0), fq(0), fq(0), fq(2), fq(0), fq(0), fq(3), fq(0)];
        let expected_n_vars = 3;
        let expected_poly =
            DenseMultilinearPolynomial::new_with_coefficients(expected_terms, expected_n_vars);

        assert_eq!(poly, expected_poly);
    }
}
