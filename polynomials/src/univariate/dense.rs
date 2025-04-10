use ark_ff::PrimeField;
use std::iter::{Product, Sum};
use std::ops::{Add, Mul};

// dense polynomial
#[derive(Debug, Clone, PartialEq)]
pub struct DenseUnivariatePolynomial<F: PrimeField> {
    // 1 coefficient for each power of x
    pub coefficients: Vec<F>,
}

impl<F: PrimeField> DenseUnivariatePolynomial<F> {
    pub fn new(coefficients: Vec<F>) -> Self {
        Self { coefficients }
    }

    pub fn degree(&self) -> usize {
        self.coefficients.len() - 1
    }

    pub fn scalar_mul(&self, scalar: F) -> Self {
        DenseUnivariatePolynomial {
            coefficients: self
                .coefficients
                .iter()
                .map(|coeff| coeff.mul(scalar))
                .collect(),
        }
    }

    pub fn basis(x: F, interpolating_set: &[F]) -> Self {
        //  numerator
        let numerators = interpolating_set
            .iter()
            .filter(|&val| *val != x)
            .map(|x_prime| Self::new(vec![x_prime.neg(), F::ONE]))
            .product::<DenseUnivariatePolynomial<F>>();

        // denominator
        let denominator = F::ONE.div(numerators.evaluate(x));

        numerators.scalar_mul(denominator)
    }

    pub fn evaluate(&self, x: F) -> F {
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

    pub fn interpolate(xs: &[F], ys: &[F]) -> Self {
        assert_eq!(xs.len(), ys.len());
        // dot product between the ys and the lagrange basis

        xs.iter()
            .zip(ys.iter())
            .map(|(x, y)| Self::basis(*x, xs).scalar_mul(*y))
            .sum()
    }
}

impl<F: PrimeField> Mul for &DenseUnivariatePolynomial<F> {
    type Output = DenseUnivariatePolynomial<F>;

    fn mul(self, rhs: Self) -> Self::Output {
        // mul for dense
        let new_degree = self.degree() + rhs.degree();
        let mut result = vec![F::ZERO; new_degree + 1];
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

impl<F: PrimeField> Product for DenseUnivariatePolynomial<F> {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut result = DenseUnivariatePolynomial::new(vec![F::ONE]);

        for poly in iter {
            result = &result * &poly
        }

        result
    }
}

impl<F: PrimeField> Add for &DenseUnivariatePolynomial<F> {
    type Output = DenseUnivariatePolynomial<F>;

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

impl<F: PrimeField> Sum for DenseUnivariatePolynomial<F> {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut result = DenseUnivariatePolynomial::new(vec![F::ZERO]);

        for poly in iter {
            result = &result + &poly;
        }

        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::Fq;

    fn test_poly() -> DenseUnivariatePolynomial<Fq> {
        let coeffs = vec![Fq::from(1), Fq::from(2), Fq::from(3)];

        DenseUnivariatePolynomial::new(coeffs)
    }

    #[test]
    fn test_degree() {
        let poly = test_poly();

        assert_eq!(poly.degree(), 2);
    }

    #[test]
    fn test_evaluation() {
        let poly = test_poly();

        assert_eq!(poly.evaluate(Fq::from(2)), Fq::from(17));
    }

    #[test]
    fn test_scalar_mul() {
        let poly = test_poly();
        let expected_result =
            DenseUnivariatePolynomial::new(vec![Fq::from(2), Fq::from(4), Fq::from(6)]);

        assert_eq!(poly.scalar_mul(Fq::from(2)), expected_result);
    }

    #[test]
    fn test_addition() {
        let poly_1 = test_poly();
        let poly_2 = DenseUnivariatePolynomial::new(vec![
            Fq::from(3),
            Fq::from(4),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(5),
        ]);
        let expected_result = DenseUnivariatePolynomial::new(vec![
            Fq::from(4),
            Fq::from(6),
            Fq::from(3),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(5),
        ]);

        assert_eq!(&poly_1 + &poly_2, expected_result);
    }

    #[test]
    fn test_multiplication() {
        // f(x) = 5 + 2x^2
        let poly_1 = DenseUnivariatePolynomial::new(vec![Fq::from(5), Fq::from(0), Fq::from(2)]);
        // f(x) = 6 + 2x
        let poly_2 = DenseUnivariatePolynomial::new(vec![Fq::from(6), Fq::from(2)]);
        let expected_result = DenseUnivariatePolynomial::new(vec![
            Fq::from(30),
            Fq::from(10),
            Fq::from(12),
            Fq::from(4),
        ]);

        assert_eq!(&poly_1 * &poly_2, expected_result);
    }

    #[test]
    fn test_interpolation() {
        // f(x) = 2x
        // [(2, 4), (4, 8)]

        let interpolated_poly = DenseUnivariatePolynomial::interpolate(
            &[Fq::from(2), Fq::from(4)],
            &[Fq::from(4), Fq::from(8)],
        );
        let expected_result = DenseUnivariatePolynomial::new(vec![Fq::from(0), Fq::from(2)]);

        assert_eq!(interpolated_poly, expected_result);
    }
}
