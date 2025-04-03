use std::iter::{Product, Sum};
use std::ops::{Add, Mul};

#[derive(Debug, Clone, PartialEq)]
pub struct SparseUnivariatePolynomial {
    terms: Vec<(f64, usize)>,
}

impl SparseUnivariatePolynomial {
    fn new(terms: Vec<(f64, usize)>) -> Self {
        Self { terms }
    }

    fn degree(&self) -> usize {
        match self.terms.iter().max_by_key(|&(_, exp)| exp) {
            Some((_, degree)) => *degree,
            None => 0,
        }
    }

    fn scalar_mul(&self, scalar: f64) -> Self {
        let new_terms = self
            .terms
            .iter()
            .map(|(coeff, exp)| (coeff * scalar, *exp))
            .collect::<Vec<(f64, usize)>>();

        Self { terms: new_terms }
    }

    fn basis(x: f64, interpolating_set: &Vec<f64>) -> Self {
        //  numerator
        let numerators = interpolating_set
            .iter()
            .filter(|&val| *val != x)
            .map(|x_prime| Self::new(vec![(-x_prime, 0), (1.0, 1)]))
            .product::<SparseUnivariatePolynomial>();

        // denominator
        let denominator = 1.0 / numerators.evaluate(x);

        numerators.scalar_mul(denominator)
    }

    fn evaluate(&self, x: f64) -> f64 {
        self.terms
            .iter()
            .map(|(coeff, exp)| coeff * x.powf(*exp as f64))
            .sum()
    }

    fn interpolate(xs: Vec<f64>, ys: Vec<f64>) -> Self {
        assert_eq!(xs.len(), ys.len());

        xs.iter()
            .zip(ys.iter())
            .map(|(x, y)| Self::basis(*x, &xs).scalar_mul(*y))
            .sum()
    }
}

impl Add for SparseUnivariatePolynomial {
    type Output = SparseUnivariatePolynomial;

    fn add(self, rhs: Self) -> Self::Output {
        let (bigger_poly, smaller_poly) = if self.degree() < rhs.degree() {
            (rhs.clone(), self)
        } else {
            (self.clone(), rhs)
        };
        let min_len = smaller_poly.terms.len();

        let min_terms = bigger_poly
            .terms
            .iter()
            .take(min_len)
            .zip(smaller_poly.terms.iter())
            .flat_map(|(&(coeff1, exp1), &(coeff2, exp2))| match exp1.cmp(&exp2) {
                std::cmp::Ordering::Less => {
                    Some((coeff1, exp1)).into_iter().chain(Some((coeff2, exp2)))
                }
                std::cmp::Ordering::Greater => {
                    Some((coeff2, exp2)).into_iter().chain(Some((coeff1, exp1)))
                }
                std::cmp::Ordering::Equal => (coeff1 + coeff2 != 0.0)
                    .then_some((coeff1 + coeff2, exp1))
                    .into_iter()
                    .chain(None),
            });

        let rest = bigger_poly
            .terms
            .iter()
            .skip(min_len)
            .map(|&(coeff, exp)| (coeff, exp));

        let summed_terms = min_terms.chain(rest).collect();

        SparseUnivariatePolynomial::new(summed_terms)
    }
}

impl Mul for SparseUnivariatePolynomial {
    type Output = SparseUnivariatePolynomial;

    fn mul(self, rhs: Self) -> Self::Output {
        // mul for sparse
        let mut result = vec![];
        for (coeff1, exp1) in self.terms {
            for (coeff2, exp2) in &rhs.terms {
                result.push((coeff1 * coeff2, exp1 + exp2));
            }
        }

        SparseUnivariatePolynomial { terms: result }
    }
}

impl Sum for SparseUnivariatePolynomial {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut result = SparseUnivariatePolynomial::new(vec![(0.0, 0)]);

        for poly in iter {
            result = result + poly;
        }

        result
    }
}

impl Product for SparseUnivariatePolynomial {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut result = SparseUnivariatePolynomial::new(vec![(1.0, 0)]);

        for poly in iter {
            result = result * poly
        }

        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_degree() {
        // f(x) = 1 + 2x + 3x^2
        let poly = SparseUnivariatePolynomial::new(vec![(1.0, 0), (2.0, 1), (3.0, 2)]);

        assert_eq!(poly.degree(), 2);
    }

    #[test]
    fn test_evaluation() {
        let poly = SparseUnivariatePolynomial::new(vec![(1.0, 0), (2.0, 1), (3.0, 2)]);

        assert_eq!(poly.evaluate(2.0), 17.0);
    }

    #[test]
    fn test_scalar_mul() {
        let poly = SparseUnivariatePolynomial::new(vec![(1.0, 0), (2.0, 1), (3.0, 2)]);
        let expected_result = SparseUnivariatePolynomial::new(vec![(2.0, 0), (4.0, 1), (6.0, 2)]);

        assert_eq!(poly.scalar_mul(2.0), expected_result);
    }

    #[test]
    fn test_addition() {
        let expected_result = SparseUnivariatePolynomial::new(vec![
            (4.0, 0),
            (6.0, 1),
            (3.0, 2),
            (5.0, 11),
            (-1.0, 120),
        ]);
        // f(x) = 1 + 2x + 3x^2
        let poly_1 = SparseUnivariatePolynomial::new(vec![(1.0, 0), (2.0, 1), (3.0, 2)]);
        // f(x) = 3 + 4x + 5x^11
        let poly_2 =
            SparseUnivariatePolynomial::new(vec![(3.0, 0), (4.0, 1), (5.0, 11), (-1.0, 120)]);

        assert_eq!(poly_1 + poly_2, expected_result);
    }

    #[test]
    fn test_multiplication() {
        // f(x) = 5 + 2x^2
        let poly_1 = SparseUnivariatePolynomial::new(vec![(5.0, 0), (2.0, 2)]);
        // f(x) = 6 + 2x
        let poly_2 = SparseUnivariatePolynomial::new(vec![(6.0, 0), (2.0, 1)]);
        let expected_result =
            SparseUnivariatePolynomial::new(vec![(30.0, 0), (10.0, 1), (12.0, 2), (4.0, 3)]);

        assert_eq!(poly_1 * poly_2, expected_result);
    }

    #[test]
    fn test_interpolation() {
        // f(x) = 2x
        // [(2, 4), (4, 8)]

        let interpolated_poly =
            SparseUnivariatePolynomial::interpolate(vec![2.0, 4.0], vec![4.0, 8.0]);
        let expected_result = SparseUnivariatePolynomial::new(vec![(2.0, 1)]);

        assert_eq!(interpolated_poly, expected_result);
    }
}
