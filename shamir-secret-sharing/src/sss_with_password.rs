use ark_ff::PrimeField;
use polynomials::univariate::dense::DenseUnivariatePolynomial;

pub fn shares<F: PrimeField>(
    secret: F,
    password: F,
    num_shares: u64,
    threshold: u64,
) -> Vec<(F, F)> {
    let mut shares: Vec<(F, F)> = Vec::new();
    let mut rng = rand::thread_rng();

    loop {
        let (mut xs, mut ys) = (1..threshold)
            .map(|i| (F::from(i), F::rand(&mut rng)))
            .unzip::<_, _, Vec<F>, Vec<F>>();

        xs.splice(0..0, [password]);
        ys.splice(0..0, [secret]);

        let poly = DenseUnivariatePolynomial::interpolate(&xs, &ys);

        if poly.degree() == (threshold - 1) as usize {
            for i in 1..num_shares {
                shares.push((F::from(i), poly.evaluate(F::from(i))));
            }
            break;
        }
    }

    shares
}

pub fn recover_secret<F: PrimeField>(shares: Vec<(F, F)>, password: F) -> F {
    let mut xs: Vec<F> = Vec::new();
    let mut ys: Vec<F> = Vec::new();

    for share in shares {
        xs.push(share.0);
        ys.push(share.1);
    }

    let poly = DenseUnivariatePolynomial::interpolate(&xs, &ys);
    let secret = poly.evaluate(F::from(password));

    secret
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::Fq;
    use ark_ff::AdditiveGroup;

    #[test]
    fn test_recover_secret_with_password() {
        let secret = Fq::from(1729);
        let password = Fq::from(123);
        let threshold = 4;
        let num_of_shares = 10;

        let shares = shares(secret, password, num_of_shares, threshold);

        let recovered_secret = recover_secret(shares, password);

        assert_eq!(recovered_secret, secret);
    }

    #[test]
    fn test_recover_with_password_wrong_secret_fails() {
        let secret = Fq::from(220284);
        let password = Fq::from(123);
        let threshold = 4;
        let num_of_shares = 10;

        let shares = shares(secret, password, num_of_shares, threshold);

        let recovered_secret = recover_secret(shares, password);

        assert_ne!(recovered_secret, Fq::from(10));
    }

    #[test]
    fn test_recover_with_wrong_password_fails() {
        let secret = Fq::from(220284);
        let password = Fq::from(123);
        let threshold = 4;
        let num_of_shares = 10;

        let shares = shares(secret, password, num_of_shares, threshold);

        let recovered_secret = recover_secret(shares, Fq::ZERO);

        assert_ne!(recovered_secret, secret);
    }
}
