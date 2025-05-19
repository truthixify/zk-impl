pub mod prover;
pub mod verifier;

pub use prover::*;
pub use verifier::*;

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::Fq;
    use polynomials::{
        composed::{ProductPolynomial, SumPolynomial},
        multilinear::MultilinearPolynomial,
    };

    fn fq(x: i64) -> Fq {
        Fq::from(x)
    }

    fn poly1a() -> MultilinearPolynomial<Fq> {
        MultilinearPolynomial::new(vec![
            fq(6),
            fq(9),
            fq(7),
            fq(6),
            fq(9),
            fq(12),
            fq(10),
            fq(9),
            fq(7),
            fq(10),
            fq(8),
            fq(7),
            fq(6),
            fq(9),
            fq(7),
            fq(6),
        ])
    }

    fn poly1b() -> MultilinearPolynomial<Fq> {
        MultilinearPolynomial::new(vec![
            fq(0),
            fq(0),
            fq(0),
            fq(0),
            fq(0),
            fq(0),
            fq(0),
            fq(0),
            fq(0),
            fq(0),
            fq(0),
            fq(2),
            fq(0),
            fq(0),
            fq(0),
            fq(0),
        ])
    }

    fn poly2a() -> MultilinearPolynomial<Fq> {
        MultilinearPolynomial::new(vec![
            fq(9),
            fq(18),
            fq(12),
            fq(9),
            fq(18),
            fq(36),
            fq(24),
            fq(18),
            fq(12),
            fq(24),
            fq(16),
            fq(12),
            fq(9),
            fq(18),
            fq(12),
            fq(9),
        ])
    }

    fn poly2b() -> MultilinearPolynomial<Fq> {
        MultilinearPolynomial::new(vec![
            fq(0),
            fq(-1),
            fq(0),
            fq(0),
            fq(0),
            fq(0),
            fq(0),
            fq(0),
            fq(0),
            fq(0),
            fq(0),
            fq(0),
            fq(0),
            fq(0),
            fq(0),
            fq(0),
        ])
    }

    fn prod_poly1() -> ProductPolynomial<Fq> {
        ProductPolynomial::new(vec![poly1a(), poly1b()])
    }

    fn prod_poly2() -> ProductPolynomial<Fq> {
        ProductPolynomial::new(vec![poly2a(), poly2b()])
    }

    fn sum_poly() -> SumPolynomial<Fq> {
        SumPolynomial::new(vec![prod_poly1(), prod_poly2()])
    }

    // This test is from lambdaclass blog: https://blog.lambdaclass.com/gkr-protocol-a-step-by-step-example/
    #[test]
    fn test_full_sumcheck() {
        let (claimed_sum, round_polys) = prove(sum_poly());

        assert!(verify(sum_poly(), claimed_sum, round_polys))
    }

    // This test is from Sir Casweeney: https://github.com/casweeney/zk-cryptography-research-implementations/blob/main/sumcheck_protocol/src/gkr_sumcheck/sumcheck_gkr_protocol.rs
    #[test]
    fn test_prover_and_verifier() {
        let poly1a = MultilinearPolynomial::new(vec![fq(0), fq(0), fq(0), fq(2)]);
        let poly2a = MultilinearPolynomial::new(vec![fq(0), fq(0), fq(0), fq(3)]);
        let product_poly1 = ProductPolynomial::new(vec![poly1a, poly2a]);

        let poly1b = MultilinearPolynomial::new(vec![fq(0), fq(0), fq(0), fq(2)]);
        let poly2b = MultilinearPolynomial::new(vec![fq(0), fq(0), fq(0), fq(3)]);
        let product_poly2 = ProductPolynomial::new(vec![poly1b, poly2b]);

        let sum_polynomial = SumPolynomial::new(vec![product_poly1, product_poly2]);

        let (claimed_sum, round_polys) = prove(sum_polynomial.clone());

        let verified = verify(sum_polynomial.clone(), claimed_sum, round_polys);

        assert_eq!(verified, true);
    }
}
