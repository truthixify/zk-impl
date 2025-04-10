use ark_ff::PrimeField;
use std::cmp::Ordering;
use std::iter::{Product, Sum};
use std::ops::{Add, Mul};

#[derive(Debug, Clone, PartialEq)]
pub struct SparseUnivariatePolynomial<F: PrimeField> {
    terms: Vec<(F, usize)>,
}

impl<F: PrimeField> SparseUnivariatePolynomial<F> {
    pub fn new(terms: Vec<(F, usize)>) -> Self {
        Self { terms }
    }

    pub fn degree(&self) -> usize {
        match self.terms.iter().max_by_key(|&(_, exp)| exp) {
            Some((_, degree)) => *degree,
            None => 0,
        }
    }

    pub fn scalar_mul(&self, scalar: F) -> Self {
        let new_terms = self
            .terms
            .iter()
            .map(|(coeff, exp)| (coeff.mul(scalar), *exp))
            .collect::<Vec<(F, usize)>>();

        Self { terms: new_terms }
    }

    pub fn basis(x: F, interpolating_set: &[F]) -> Self {
        //  numerator
        let numerators = interpolating_set
            .iter()
            .filter(|&val| *val != x)
            .map(|x_prime| Self::new(vec![(x_prime.neg(), 0), (F::ONE, 1)]))
            .product::<SparseUnivariatePolynomial<F>>();

        // denominator
        let denominator = F::ONE.div(numerators.evaluate(x));

        numerators.scalar_mul(denominator)
    }

    pub fn evaluate(&self, x: F) -> F {
        self.terms
            .iter()
            .map(|(coeff, exp)| coeff.mul(x.pow([*exp as u64])))
            .sum()
    }

    pub fn interpolate(xs: &[F], ys: &[F]) -> Self {
        assert_eq!(xs.len(), ys.len());

        xs.iter()
            .zip(ys.iter())
            .map(|(x, y)| Self::basis(*x, xs).scalar_mul(*y))
            .sum()
    }
}

impl<F: PrimeField> Add for &SparseUnivariatePolynomial<F> {
    type Output = SparseUnivariatePolynomial<F>;

    fn add(self, rhs: Self) -> Self::Output {
        // let (bigger_poly, smaller_poly) = if self.degree() < rhs.degree() {
        //     (rhs.clone(), self)
        // } else {
        //     (self.clone(), rhs)
        // };
        // let min_len = smaller_poly.terms.len();

        // let min_terms = bigger_poly
        //     .terms
        //     .iter()
        //     .take(min_len)
        //     .zip(smaller_poly.terms.iter())
        //     .flat_map(|(&(coeff1, exp1), &(coeff2, exp2))| match exp1.cmp(&exp2) {
        //         std::cmp::Ordering::Less => {
        //             Some((coeff1, exp1)).into_iter().chain(Some((coeff2, exp2)))
        //         }
        //         std::cmp::Ordering::Greater => {
        //             Some((coeff2, exp2)).into_iter().chain(Some((coeff1, exp1)))
        //         }
        //         std::cmp::Ordering::Equal => (coeff1 + coeff2 != F::ZERO)
        //             .then_some((coeff1 + coeff2, exp1))
        //             .into_iter()
        //             .chain(None),
        //     });

        // let rest = bigger_poly
        //     .terms
        //     .iter()
        //     .skip(min_len)
        //     .map(|&(coeff, exp)| (coeff, exp));

        // let summed_terms = min_terms.chain(rest).collect();
        let mut summed_terms = vec![];
        let mut i = 0;
        let mut j = 0;

        while i < self.terms.len() && j < rhs.terms.len() {
            let (coeff1, exp1) = self.terms[i];
            let (coeff2, exp2) = rhs.terms[j];

            match exp1.cmp(&exp2) {
                Ordering::Equal => {
                    summed_terms.push((coeff1.add(coeff2), exp1));
                    i += 1;
                    j += 1;
                }
                Ordering::Less => {
                    summed_terms.push((coeff1, exp1));
                    i += 1;
                }
                Ordering::Greater => {
                    summed_terms.push((coeff2, exp2));
                    j += 1;
                }
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

        SparseUnivariatePolynomial::new(summed_terms)
    }
}

impl<F: PrimeField> Mul for &SparseUnivariatePolynomial<F> {
    type Output = SparseUnivariatePolynomial<F>;

    fn mul(self, rhs: Self) -> Self::Output {
        // mul for sparse
        let mut result = vec![];
        for (coeff1, exp1) in &self.terms {
            for (coeff2, exp2) in &rhs.terms {
                result.push((*coeff1 * coeff2, exp1 + exp2));
            }
        }

        SparseUnivariatePolynomial { terms: result }
    }
}

impl<F: PrimeField> Sum for SparseUnivariatePolynomial<F> {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut result = SparseUnivariatePolynomial::new(vec![(F::ZERO, 0)]);

        for poly in iter {
            result = &result + &poly;
        }

        result
    }
}

impl<F: PrimeField> Product for SparseUnivariatePolynomial<F> {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut result = SparseUnivariatePolynomial::new(vec![(F::ONE, 0)]);

        for poly in iter {
            result = &result * &poly
        }

        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::Fq;

    fn test_poly() -> SparseUnivariatePolynomial<Fq> {
        let coeffs = vec![(Fq::from(1), 0), (Fq::from(2), 1), (Fq::from(3), 2)];

        SparseUnivariatePolynomial::new(coeffs)
    }

    #[test]
    fn test_degree() {
        // f(x) = 1 + 2x + 3x^2
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
        let expected_result = SparseUnivariatePolynomial::new(vec![
            (Fq::from(2), 0),
            (Fq::from(4), 1),
            (Fq::from(6), 2),
        ]);

        assert_eq!(poly.scalar_mul(Fq::from(2)), expected_result);
    }

    #[test]
    fn test_addition() {
        let expected_result = SparseUnivariatePolynomial::new(vec![
            (Fq::from(4), 0),
            (Fq::from(6), 1),
            (Fq::from(3), 2),
            (Fq::from(5), 11),
            (Fq::from(-1), 120),
        ]);
        // f(x) = 1 + 2x + 3x^2
        let poly_1 = test_poly();
        // f(x) = 3 + 4x + 5x^11
        let poly_2 = SparseUnivariatePolynomial::new(vec![
            (Fq::from(3), 0),
            (Fq::from(4), 1),
            (Fq::from(5), 11),
            (Fq::from(-1), 120),
        ]);

        assert_eq!(&poly_1 + &poly_2, expected_result);
    }

    #[test]
    fn test_multiplication() {
        // f(x) = 5 + 2x^2
        let poly_1 = SparseUnivariatePolynomial::new(vec![(Fq::from(5), 0), (Fq::from(2), 2)]);
        // f(x) = 6 + 2x
        let poly_2 = SparseUnivariatePolynomial::new(vec![(Fq::from(6), 0), (Fq::from(2), 1)]);
        let expected_result = SparseUnivariatePolynomial::new(vec![
            (Fq::from(30), 0),
            (Fq::from(10), 1),
            (Fq::from(12), 2),
            (Fq::from(4), 3),
        ]);

        assert_eq!(&poly_1 * &poly_2, expected_result);
    }

    #[test]
    fn test_interpolation() {
        // f(x) = 2x
        // [(2, 4), (4, 8)]

        let interpolated_poly = SparseUnivariatePolynomial::interpolate(
            &[Fq::from(2), Fq::from(4)],
            &[Fq::from(4), Fq::from(8)],
        );
        let expected_result = SparseUnivariatePolynomial::new(vec![(Fq::from(2), 1)]);

        assert_eq!(interpolated_poly, expected_result);
    }
}
