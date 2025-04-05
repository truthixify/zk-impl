# Polynomial Library in Rust

This repository contains implementations of univariate polynomials over finite fields in Rust, with two primary data structures:

- **Dense Polynomial Representation** (`dense.rs`)
- **Sparse Polynomial Representation** (`sparse.rs`)

These implementations leverage the `ark_ff` crate from the [arkworks ecosystem](https://arkworks.rs/) for finite field arithmetic.

## Files and Data Structures

### Dense Representation (`dense.rs`)
- Represents polynomials in dense form, storing a coefficient for each power of the polynomial explicitly.
- Efficient for operations involving mostly non-zero terms.

#### Main Features:
- Polynomial evaluation using Horner's method.
- Addition, multiplication, scalar multiplication, and polynomial interpolation.
- Supports standard arithmetic operations via Rust traits (`Add`, `Mul`, `Sum`, `Product`).

### Sparse Representation (`sparse.rs`)
- Represents polynomials by explicitly storing only non-zero terms as `(coefficient, exponent)` pairs.
- Efficient for polynomials with many zero coefficients or high-degree sparse polynomials.

#### Main Features:
- Efficient polynomial evaluation optimized for sparse polynomials.
- Arithmetic operations similar to dense representation but optimized for sparse data.
- Implements interpolation and arithmetic via Rust traits (`Add`, `Mul`, `Sum`, `Product`).

## Usage

These modules can be imported into your Rust project for algebraic computations. Here's a quick example:

```rust
use ark_bls12_381::Fq;
use crate::dense::DenseUnivariatePolynomial;

let poly = DenseUnivariatePolynomial::new(vec![Fq::from(1), Fq::from(2), Fq::from(3)]);
println!("Polynomial evaluation at x=2: {:?}", poly.evaluate(Fq::from(2)));
```

Replace `DenseUnivariatePolynomial` with `SparseUnivariatePolynomial` to utilize sparse representation.

## Testing

Tests are provided within each file. Run the tests using Cargo:

```bash
cargo test
```

## Dependencies

- [ark_ff](https://docs.rs/ark-ff/latest/ark_ff/): For finite field arithmetic.
- [ark_bls12_381](https://docs.rs/ark-bls12-381/latest/ark_bls12_381/): Used in tests to demonstrate polynomial arithmetic with BLS12-381 curve field elements.

## Contributing

Contributions are welcome! Feel free to submit issues and pull requests to enhance the polynomial operations or to fix bugs.

## License

This project is licensed under the MIT License.