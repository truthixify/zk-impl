use ark_ff::PrimeField;

pub struct ReedSolomonFingerprint<F: PrimeField> {
    r: F,
    v: F,
}

pub fn hash<F: PrimeField>(data_a: &[F]) -> ReedSolomonFingerprint<F> {
    assert!(
        F::MODULUS.gt(&F::BigInt::from(data_a.len() as u64)),
        "Length of input data is greater than modulus of the prime field."
    );

    let mut rng = rand::thread_rng();
    let r = F::rand(&mut rng);
    let v = data_a
        .iter()
        .enumerate()
        .map(|(index, x)| *x * r.pow([index as u64]))
        .sum();

    ReedSolomonFingerprint { r, v }
}

pub fn verify<F: PrimeField>(data_b: &[F], rsf: ReedSolomonFingerprint<F>) -> bool {
    let eval_b = data_b
        .iter()
        .enumerate()
        .map(|(index, x)| *x * rsf.r.pow([index as u64]))
        .sum();

    rsf.v == eval_b
}

#[cfg(test)]
mod tests {
    use crate::{hash, verify};
    use ark_bls12_381::Fq;
    use rand::Rng;

    fn fq(val: u64) -> Fq {
        Fq::from(val)
    }

    #[test]
    fn test_fingerprint_with_small_data() {
        let data = vec![fq(1), fq(2), fq(3), fq(4), fq(5)];
        let fingerprint = hash(&data);
        assert!(verify(&data, fingerprint));
    }

    #[test]
    fn test_fingerprint_with_large_data() {
        let data: Vec<Fq> = (0..1000).map(fq).collect();
        let fingerprint = hash(&data);
        assert!(verify(&data, fingerprint));
    }

    #[test]
    fn test_fingerprint_invalid_data() {
        let data_a = vec![fq(1), fq(2), fq(3), fq(4), fq(5)];
        let mut data_b = data_a.clone();
        data_b[2] = fq(100); // mutate one value

        let fingerprint = hash(&data_a);
        assert!(!verify(&data_b, fingerprint));
    }

    #[test]
    fn test_fingerprint_with_zeros() {
        let data = vec![fq(0); 500];
        let fingerprint = hash(&data);
        assert!(verify(&data, fingerprint));
    }

    #[test]
    fn test_fingerprint_with_random_data() {
        let mut rng = rand::thread_rng();
        let data: Vec<Fq> = (0..256).map(|_| fq(rng.gen_range(0..10000))).collect();
        let fingerprint = hash(&data);
        assert!(verify(&data, fingerprint));
    }
}
