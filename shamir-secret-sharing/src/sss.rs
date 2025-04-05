use ark_ff::PrimeField;
use polynomials::univariate::dense::DenseUnivariatePolynomial;

pub fn shares<F: PrimeField>(secret: F, num_shares: u64, threshold: u64) -> Vec<(F, F)> {
    let mut shares: Vec<(F, F)> = Vec::new();
    let mut rng = rand::thread_rng();
    let mut coeffs = (1..threshold)
        .map(|_| F::rand(&mut rng))
        .collect::<Vec<F>>();

    coeffs.splice(0..0, [secret]);

    let poly = DenseUnivariatePolynomial::new(coeffs);

    for i in 1..num_shares {
        shares.push((F::from(i), poly.evaluate(F::from(i))));
    }

    shares
}

pub fn recover_secret<F: PrimeField>(shares: Vec<(F, F)>) -> F {
    let mut xs: Vec<F> = Vec::new();
    let mut ys: Vec<F> = Vec::new();

    for share in shares {
        xs.push(share.0);
        ys.push(share.1);
    }

    let poly = DenseUnivariatePolynomial::interpolate(xs, ys);
    let secret = poly.evaluate(F::from(0));

    secret
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::Fq;

    #[test]
    fn test_recover_secret() {
        let secret = Fq::from(1729);
        let threshold = 4;
        let num_of_shares = 10;

        let shares = shares(secret, num_of_shares, threshold);

        let recovered_secret = recover_secret(shares);

        assert_eq!(recovered_secret, secret);
    }

    #[test]
    fn test_recover_wrong_secret_fails() {
        let secret = Fq::from(220284);
        let threshold = 4;
        let num_of_shares = 10;

        let shares = shares(secret, num_of_shares, threshold);

        let recovered_secret = recover_secret(shares);

        assert_ne!(recovered_secret, Fq::from(10));
    }
}
