use ark_ff::PrimeField;
use std::ops::{Add, Mul};

#[derive(Debug, Clone, PartialEq)]
pub struct Matrix<F: PrimeField> {
    rep: Vec<Vec<F>>,
}

impl<F: PrimeField> Matrix<F> {
    pub fn new(rep: Vec<Vec<F>>) -> Self {
        Matrix { rep }
    }

    pub fn nrows(&self) -> usize {
        self.rep.len()
    }

    pub fn ncols(&self) -> usize {
        self.rep[0].len()
    }

    pub fn scalar_mul(&self, scalar: F) -> Self {
        let new_rep = self
            .rep
            .iter()
            .map(|row| row.iter().map(|x| x.mul(scalar)).collect())
            .collect();

        Matrix::new(new_rep)
    }

    pub fn transpose(&self) -> Self {
        let mut transposed = vec![vec![]; self.ncols()];
        for row in &self.rep {
            for (j, &val) in row.iter().enumerate() {
                transposed[j].push(val);
            }
        }

        Matrix::new(transposed)
    }

    pub fn add_matrices(&self, other: &Self) -> Self {
        assert_eq!(
            self.nrows(),
            other.nrows(),
            "The two matrices must have the same number of rows"
        );
        assert_eq!(
            self.ncols(),
            other.ncols(),
            "The two matrices must have the same number of columns"
        );

        let new_rep = self
            .rep
            .iter()
            .zip(&other.rep)
            .map(|(row_a, row_b)| row_a.iter().zip(row_b).map(|(&a, &b)| a + b).collect())
            .collect();

        Matrix::new(new_rep)
    }

    pub fn mul_matrices(&self, other: &Self) -> Self {
        assert_eq!(
            self.ncols(),
            other.nrows(),
            "Inner dimensions must match for multiplication"
        );

        let mut new_rep = vec![vec![F::ZERO; other.ncols()]; self.nrows()];
        for i in 0..self.nrows() {
            for j in 0..other.ncols() {
                for k in 0..self.ncols() {
                    new_rep[i][j] = new_rep[i][j] + self.rep[i][k] * other.rep[k][j];
                }
            }
        }

        Matrix::new(new_rep)
    }
}

impl<F: PrimeField> Add for Matrix<F> {
    type Output = Matrix<F>;

    fn add(self, rhs: Self) -> Self::Output {
        self.add_matrices(&rhs)
    }
}

impl<F: PrimeField> Add for &Matrix<F> {
    type Output = Matrix<F>;

    fn add(self, rhs: Self) -> Self::Output {
        self.add_matrices(&rhs)
    }
}

impl<F: PrimeField> Mul for Matrix<F> {
    type Output = Matrix<F>;

    fn mul(self, rhs: Self) -> Self::Output {
        self.mul_matrices(&rhs)
    }
}

impl<F: PrimeField> Mul for &Matrix<F> {
    type Output = Matrix<F>;

    fn mul(self, rhs: Self) -> Self::Output {
        self.mul_matrices(&rhs)
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
    fn test_nrows_and_ncols() {
        let m = Matrix::new(vec![vec![fq(1), fq(2), fq(3)], vec![fq(4), fq(5), fq(6)]]);

        assert_eq!(m.nrows(), 2);
        assert_eq!(m.ncols(), 3);
    }

    #[test]
    fn test_addition() {
        let a = Matrix::new(vec![vec![fq(1), fq(2)], vec![fq(3), fq(4)]]);
        let b = Matrix::new(vec![vec![fq(5), fq(6)], vec![fq(7), fq(8)]]);
        let expected = Matrix::new(vec![vec![fq(6), fq(8)], vec![fq(10), fq(12)]]);

        assert_eq!(a + b, expected);
    }

    #[test]
    fn test_transpose() {
        let m = Matrix::new(vec![vec![fq(1), fq(2), fq(3)], vec![fq(4), fq(5), fq(6)]]);
        let expected = Matrix::new(vec![
            vec![fq(1), fq(4)],
            vec![fq(2), fq(5)],
            vec![fq(3), fq(6)],
        ]);

        assert_eq!(m.transpose(), expected);
    }

    #[test]
    #[should_panic(expected = "The two matrices must have the same number of rows")]
    fn test_addition_row_mismatch_panics() {
        let a = Matrix::new(vec![vec![fq(1), fq(2)]]);
        let b = Matrix::new(vec![vec![fq(1), fq(2)], vec![fq(3), fq(4)]]);

        let _ = a + b;
    }

    #[test]
    #[should_panic(expected = "The two matrices must have the same number of columns")]
    fn test_addition_column_mismatch_panics() {
        let a = Matrix::new(vec![vec![fq(1), fq(2)]]);
        let b = Matrix::new(vec![vec![fq(1)]]);

        let _ = a + b;
    }

    #[test]
    fn test_scalar_multiplication() {
        let m = Matrix::new(vec![vec![fq(1), fq(2)], vec![fq(3), fq(4)]]);
        let expected = Matrix::new(vec![vec![fq(2), fq(4)], vec![fq(6), fq(8)]]);

        assert_eq!(m.scalar_mul(fq(2)), expected);
    }

    #[test]
    fn test_matrix_multiplication() {
        let a = Matrix::new(vec![vec![fq(1), fq(2)], vec![fq(3), fq(4)]]);
        let b = Matrix::new(vec![vec![fq(5), fq(6)], vec![fq(7), fq(8)]]);
        let expected = Matrix::new(vec![vec![fq(19), fq(22)], vec![fq(43), fq(50)]]);

        assert_eq!(a * b, expected);
    }

    #[test]
    #[should_panic(expected = "Inner dimensions must match for multiplication")]
    fn test_matrix_mul_invalid_dimensions() {
        let a = Matrix::new(vec![vec![fq(1), fq(2)]]);
        let b = Matrix::new(vec![vec![fq(1), fq(2)]]);

        let _ = a * b;
    }
}
