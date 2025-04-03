use std::iter::{Product, Sum};
use std::ops::{Add, Mul};

// dense polynomial
#[derive(Debug, Clone, PartialEq)]
pub struct DenseUnivariatePolynomial {
    // 1 coefficient for each power of x
    coefficients: Vec<f64>,
}

impl DenseUnivariatePolynomial {
    fn new(coefficients: Vec<f64>) -> Self {
        Self { coefficients }
    }

    fn degree(&self) -> usize {
        self.coefficients.len() - 1
    }

    fn scalar_mul(&self, scalar: f64) -> Self {
        DenseUnivariatePolynomial {
            coefficients: self
                .coefficients
                .iter()
                .map(|coeff| coeff * scalar)
                .collect(),
        }
    }

    fn basis(x: f64, interpolating_set: &Vec<f64>) -> Self {
        //  numerator
        let numerators = interpolating_set
            .iter()
            .filter(|&val| *val != x)
            .map(|x_prime| Self::new(vec![-x_prime, 1.0]))
            .product::<DenseUnivariatePolynomial>();

        // denominator
        let denominator = 1.0 / numerators.evaluate(x);

        numerators.scalar_mul(denominator)
    }

    fn evaluate(&self, x: f64) -> f64 {
        // for loop method
        // let mut eval = 0.0;
        // let mut current_x = 1.0;

        // for i in 0..self.coefficients.len() {
        //     eval += self.coefficients[i] * current_x;
        //     current_x *= x;
        // }

        // eval
        // c1 + c2*x + c3*x^2 = c1 + x*(c2 + c3*x)
        self.coefficients
            .iter()
            .rev()
            .cloned()
            .reduce(|acc, curr| acc * x + curr)
            .expect("Something went wrong when evaluationg polynomial")

        // iterator method
        // self.coefficients
        //     .iter()
        //     .enumerate()
        //     .map(|(power, &coeff)| coeff * x.powf(power as f64))
        //     .sum()
    }

    fn interpolate(xs: Vec<f64>, ys: Vec<f64>) -> Self {
        assert_eq!(xs.len(), ys.len());
        // dot product between the ys and the lagrange basis

        xs.iter()
            .zip(ys.iter())
            .map(|(x, y)| Self::basis(*x, &xs).scalar_mul(*y))
            .sum()
    }
}

impl Mul for DenseUnivariatePolynomial {
    type Output = DenseUnivariatePolynomial;

    fn mul(self, rhs: Self) -> Self::Output {
        // mul for dense
        let new_degree = self.degree() + rhs.degree();
        let mut result = vec![0.0; new_degree + 1];
        for i in 0..self.coefficients.len() {
            for j in 0..rhs.coefficients.len() {
                result[i + j] += self.coefficients[i] * rhs.coefficients[j]
            }
        }

        DenseUnivariatePolynomial {
            coefficients: result,
        }
    }
}

impl Add for DenseUnivariatePolynomial {
    type Output = DenseUnivariatePolynomial;

    fn add(self, rhs: Self) -> Self::Output {
        let (mut bigger_poly, smaller_poly) = if self.degree() < rhs.degree() {
            (rhs.clone(), self)
        } else {
            (self.clone(), rhs)
        };

        let _ = bigger_poly
            .coefficients
            .iter_mut()
            .zip(smaller_poly.coefficients.iter())
            .map(|(b_coeff, s_coeff)| *b_coeff += s_coeff)
            .collect::<()>();

        bigger_poly
    }
}

impl Sum for DenseUnivariatePolynomial {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut result = DenseUnivariatePolynomial::new(vec![0.0]);

        for poly in iter {
            result = result + poly;
        }

        result
    }
}

impl Product for DenseUnivariatePolynomial {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut result = DenseUnivariatePolynomial::new(vec![1.0]);

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
        let poly = DenseUnivariatePolynomial::new(vec![1.0, 2.0, 3.0]);

        assert_eq!(poly.degree(), 2);
    }

    #[test]
    fn test_evaluation() {
        let poly = DenseUnivariatePolynomial::new(vec![1.0, 2.0, 3.0]);

        assert_eq!(poly.evaluate(2.0), 17.0);
    }

    #[test]
    fn test_scalar_mul() {
        let poly = DenseUnivariatePolynomial::new(vec![1.0, 2.0, 3.0]);
        let expected_result = DenseUnivariatePolynomial::new(vec![2.0, 4.0, 6.0]);

        assert_eq!(poly.scalar_mul(2.0), expected_result);
    }

    #[test]
    fn test_addition() {
        let expected_result = DenseUnivariatePolynomial::new(vec![
            4.0, 6.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0,
        ]);
        let poly_1 = DenseUnivariatePolynomial::new(vec![1.0, 2.0, 3.0]);
        let poly_2 = DenseUnivariatePolynomial::new(vec![
            3.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0,
        ]);

        assert_eq!(poly_1 + poly_2, expected_result);
    }

    #[test]
    fn test_multiplication() {
        // f(x) = 5 + 2x^2
        let poly_1 = DenseUnivariatePolynomial::new(vec![5.0, 0.0, 2.0]);
        // f(x) = 6 + 2x
        let poly_2 = DenseUnivariatePolynomial::new(vec![6.0, 2.0]);
        let expected_result = DenseUnivariatePolynomial::new(vec![30.0, 10.0, 12.0, 4.0]);

        assert_eq!(poly_1 * poly_2, expected_result);
    }

    #[test]
    fn test_interpolation() {
        // f(x) = 2x
        // [(2, 4), (4, 8)]

        let interpolated_poly =
            DenseUnivariatePolynomial::interpolate(vec![2.0, 4.0], vec![4.0, 8.0]);
        let expected_result = DenseUnivariatePolynomial::new(vec![0.0, 2.0]);

        assert_eq!(interpolated_poly, expected_result);
    }
}
