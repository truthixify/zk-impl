use ark_ff::PrimeField;
use matrix::Matrix;

mod matrix;

pub struct Freivald<F: PrimeField> {
    x: Vec<F>,
}

impl<F: PrimeField> Freivald<F> {
    fn new(array_size: usize) -> Self {
        // Generate random number
        // Populate vector with values r^i for i=0..matrix_size
        // Return freivald value with this vector as its x value
        let mut rng = rand::thread_rng();
        let r = F::rand(&mut rng);
        let x = (0..array_size).map(|i| r.pow([i as u64])).collect();

        Self { x }
    }

    pub fn verify(&self, matrix_a: Matrix<F>, matrix_b: Matrix<F>, supposed_ab: Matrix<F>) -> bool {
        assert!(
            check_matrix_dimensions(&matrix_a, &matrix_b, &supposed_ab),
            "Inner dimensions must match for multiplication"
        );

        // Check if a * b * x == c * x
        let x = Matrix::new(vec![self.x.clone()]).transpose();

        matrix_a * (matrix_b * x.clone()) == &supposed_ab * &x
    }

    // utility function to not have to instantiate Freivalds if you just want to make one
    // verification.
    pub fn verify_once(matrix_a: Matrix<F>, matrix_b: Matrix<F>, supposed_ab: Matrix<F>) -> bool {
        let freivald = Freivald::new(supposed_ab.nrows());
        freivald.verify(matrix_a, matrix_b, supposed_ab)
    }
}

pub fn check_matrix_dimensions<F: PrimeField>(
    matrix_a: &Matrix<F>,
    matrix_b: &Matrix<F>,
    supposed_ab: &Matrix<F>,
) -> bool {
    // Check if dimensions of making matrix_a * matrix_b matches values in supposed_ab.
    // If it doesn't we know its not the correct result independently of matrix contents
    let a_m = matrix_a.nrows();
    let a_n = matrix_a.ncols();
    let b_n = matrix_b.nrows();
    let b_p = matrix_b.ncols();
    let c_m = supposed_ab.nrows();
    let c_p = supposed_ab.ncols();

    a_m == c_m && a_n == b_n && b_p == c_p
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::Fq;

    fn fq(n: u64) -> Fq {
        Fq::from(n)
    }

    #[test]
    fn test_freivald_verify_success() {
        let a = Matrix::new(vec![vec![fq(1), fq(2)], vec![fq(3), fq(4)]]);

        let b = Matrix::new(vec![vec![fq(5), fq(6)], vec![fq(7), fq(8)]]);

        // A * B = [[19, 22], [43, 50]]
        let ab = Matrix::new(vec![vec![fq(19), fq(22)], vec![fq(43), fq(50)]]);

        let freivald = Freivald::new(2);
        assert!(freivald.verify(a.clone(), b.clone(), ab.clone()));
        assert!(Freivald::verify_once(a, b, ab));
    }

    #[test]
    #[should_panic(expected = "Inner dimensions must match for multiplication")]
    fn test_freivald_verify_fail() {
        let a = Matrix::new(vec![vec![fq(1), fq(0)], vec![fq(0), fq(1)]]);

        let b = Matrix::new(vec![vec![fq(2), fq(3)]]);

        // Incorrect result
        let wrong_ab = Matrix::new(vec![vec![fq(0), fq(0)], vec![fq(0), fq(0)]]);

        let freivald = Freivald::new(2);
        assert!(!freivald.verify(a.clone(), b.clone(), wrong_ab.clone()));
        assert!(!Freivald::verify_once(a, b, wrong_ab));
    }
}
