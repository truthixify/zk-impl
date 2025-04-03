use core::num;
use std::{
    iter::{Product, Sum},
    ops::{Add, Mul},
};

// dense polynomial
#[derive(Debug, Clone, PartialEq)]
struct UnivariatePolynomial {
    // 1 coefficient for each power of x
    coefficients: Vec<f64>,
}

impl UnivariatePolynomial {
    fn new(coefficients: Vec<f64>) -> Self {
        Self { coefficients }
    }

    fn degree(&self) -> usize {
        self.coefficients.len() - 1
    }

    fn scalar_mul(&self, scalar: f64) -> Self {
        UnivariatePolynomial {
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
            .product::<UnivariatePolynomial>();

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
            .unwrap()

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

impl Mul for UnivariatePolynomial {
    type Output = UnivariatePolynomial;

    fn mul(self, rhs: Self) -> Self::Output {
        // mul for dense
        let new_degree = self.degree() + rhs.degree();
        let mut result = vec![0.0; new_degree + 1];
        for i in 0..self.coefficients.len() {
            for j in 0..rhs.coefficients.len() {
                result[i + j] += self.coefficients[i] * rhs.coefficients[j]
            }
        }

        UnivariatePolynomial {
            coefficients: result,
        }
    }
}

impl Add for UnivariatePolynomial {
    type Output = UnivariatePolynomial;

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

impl Sum for UnivariatePolynomial {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut result = UnivariatePolynomial::new(vec![0.0]);

        for poly in iter {
            result = result + poly;
        }

        result
    }
}

impl Product for UnivariatePolynomial {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut result = UnivariatePolynomial::new(vec![1.0]);

        for poly in iter {
            result = result * poly
        }

        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn poly_with_coefficient(coefficients: Vec<f64>) -> UnivariatePolynomial {
        UnivariatePolynomial { coefficients }
    }

    fn poly_1() -> UnivariatePolynomial {
        // f(x) = 1 + 2x + 3x^2
        poly_with_coefficient(vec![1.0, 2.0, 3.0])
    }

    fn poly_2() -> UnivariatePolynomial {
        // f(x) = 3 + 4x + 5x^11
        poly_with_coefficient(vec![
            3.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0,
        ])
    }

    #[test]
    fn test_degree() {
        let poly = poly_1();

        assert_eq!(poly.degree(), 2);
    }

    #[test]
    fn test_evaluation() {
        let poly = poly_1();

        assert_eq!(poly.evaluate(2.0), 17.0);
    }

    #[test]
    fn test_scalar_mul() {
        let poly = poly_1();
        let expected_result = UnivariatePolynomial::new(vec![2.0, 4.0, 6.0]);

        assert_eq!(poly.scalar_mul(2.0), expected_result);
    }

    #[test]
    fn test_addition() {
        let expected_result = poly_with_coefficient(vec![
            4.0, 6.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0,
        ]);

        assert_eq!(poly_1() + poly_2(), expected_result);
    }

    #[test]
    fn test_multiplication() {
        // f(x) = 5 + 2x^2
        let poly_1 = poly_with_coefficient(vec![5.0, 0.0, 2.0]);
        // f(x) = 6 + 2x
        let poly_2 = poly_with_coefficient(vec![6.0, 2.0]);
        let expected_result = poly_with_coefficient(vec![30.0, 10.0, 12.0, 4.0]);

        assert_eq!(poly_1 * poly_2, expected_result);
    }

    #[test]
    fn test_interpolation() {
        // f(x) = 2x
        // [(2, 4), (4, 8)]

        let interpolated_poly = UnivariatePolynomial::interpolate(vec![2.0, 4.0], vec![4.0, 8.0]);
        let expected_result = UnivariatePolynomial::new(vec![0.0, 2.0]);

        assert_eq!(interpolated_poly, expected_result);
    }
}
