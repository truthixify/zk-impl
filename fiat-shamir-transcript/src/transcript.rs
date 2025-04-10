use ark_ff::{BigInteger, PrimeField};
use sha3::Digest;
use sha3::digest::FixedOutputReset;
use std::marker::PhantomData;

#[derive(Debug)]
pub struct Transcript<F, H> {
    hasher: H,
    field_elements: PhantomData<F>,
}

impl<F: PrimeField, H: Clone + Digest + FixedOutputReset> Transcript<F, H> {
    pub fn new() -> Self {
        Transcript {
            hasher: H::new(),
            field_elements: PhantomData,
        }
    }

    pub fn append(&mut self, data: &[u8]) {
        Digest::update(&mut self.hasher, data);
    }

    pub fn append_field_element(&mut self, element: &F) {
        self.append(&element.into_bigint().to_bytes_be());
    }

    pub fn sample_field_element(&mut self) -> F {
        let hash = &self.hasher.finalize_reset();

        Digest::update(&mut self.hasher, hash);

        F::from_be_bytes_mod_order(hash)
    }

    pub fn sample_n_elements(&mut self, n: usize) -> Vec<F> {
        (0..n).map(|_| self.sample_field_element()).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::Fq; use ark_ff::ToConstraintField;
    // Example field type
    use sha3::Keccak256;

    #[test]
    fn test_basic_challenge_generation() {
        let mut transcript = Transcript::<Fq, Keccak256>::new();

        let data = b"test_data";

        transcript.append(data);

        let challenge = transcript.sample_field_element();

        // Generate expected challenge with Keccak256 manually
        let expected_challenge_bytes = Keccak256::digest(data);
        let expected_challenge = Fq::from_be_bytes_mod_order(&expected_challenge_bytes);

        assert_eq!(
            challenge, expected_challenge,
            "Challenge does not match expected value"
        );
    }

    #[test]
    fn test_append_multiple_field_elements() {
        let mut transcript = Transcript::<Fq, Keccak256>::new();

        let field_element1: Fq = Fq::from(12345u64);
        let field_element2: Fq = Fq::from(67890u64);

        transcript.append_field_element(&field_element1);
        transcript.append_field_element(&field_element2);

        let challenge = transcript.sample_field_element();

        // Generate expected challenge with Keccak256 manually
        let data1 = field_element1.into_bigint().to_bytes_be();
        let data2 = field_element2.into_bigint().to_bytes_be();
        let mut keccak = Keccak256::new();
        keccak.update(data1);
        keccak.update(data2);
        let expected_challenge_bytes = keccak.finalize();
        let expected_challenge = Fq::from_be_bytes_mod_order(&expected_challenge_bytes);

        assert_eq!(
            challenge, expected_challenge,
            "Challenge does not match expected value after appending multiple elements"
        );
    }

    #[test]
    fn test_reuse_transcript() {
        let mut transcript1 = Transcript::<Fq, Keccak256>::new();
        let mut transcript2 = Transcript::<Fq, Keccak256>::new();

        let data = b"test_data";

        transcript1.append(data);
        transcript2.append(data);

        let challenge1 = transcript1.sample_field_element();
        let challenge2 = transcript2.sample_field_element();

        assert_eq!(
            challenge1, challenge2,
            "Challenges should match for identical inputs"
        );
    }

    #[test]
    fn test_empty_input() {
        let mut transcript = Transcript::<Fq, Keccak256>::new();

        let empty_data: &[u8] = b"";

        transcript.append(empty_data);

        let challenge = transcript.sample_field_element();

        // Generate expected challenge for empty input
        let expected_challenge_bytes = Keccak256::digest(empty_data);
        let expected_challenge = Fq::from_be_bytes_mod_order(&expected_challenge_bytes);

        assert_eq!(
            challenge, expected_challenge,
            "Challenge does not match expected value for empty input"
        );
    }

    #[test]
    fn test_field_element_encoding() {
        let mut transcript = Transcript::<Fq, Keccak256>::new();

        let field_element: Fq = Fq::from(12345u64);

        transcript.append_field_element(&field_element);

        let challenge = transcript.sample_field_element();

        // Generate expected challenge with Keccak256 manually
        let field_element_bytes = field_element.into_bigint().to_bytes_be();
        let expected_challenge_bytes = Keccak256::digest(&field_element_bytes);
        let expected_challenge = Fq::from_be_bytes_mod_order(&expected_challenge_bytes);

        assert_eq!(
            challenge, expected_challenge,
            "Challenge does not match expected value for field element encoding"
        );
    }

    #[test]
    fn test_sample_n_elements() {
        let mut transcript: Transcript<Fq, Keccak256> = Transcript::new();
        
        // Test sampling 5 elements
        let elements = transcript.sample_n_elements(5);

        // Ensure the returned vector has the correct length
        assert_eq!(elements.len(), 5, "The number of sampled elements should be 5");

        // Ensure all elements are valid field elements
        for element in &elements {
            assert!(element.to_field_elements().is_some(), "The element should be a valid field element");
        }
    }

    // Optional: Test that elements are unique (if that's important)
    #[test]
    fn test_sample_n_elements_unique() {
        let mut transcript: Transcript<Fq, Keccak256> = Transcript::new();

        // Sample 100 elements
        let elements = transcript.sample_n_elements(100);

        // Check that the sampled elements are unique
        let mut seen = std::collections::HashSet::new();
        for element in elements {
            assert!(seen.insert(element), "Found duplicate element");
        }
    }
}
