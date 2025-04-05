# Secret Sharing Schemes in Rust

This repository provides implementations of Secret Sharing Schemes (SSS) and Password-Protected Secret Sharing Schemes using Rust. These implementations leverage polynomial interpolation over finite fields, utilizing the Arkworks ecosystem for finite field operations.

## Files and Implementations

### Basic Secret Sharing (`sss.rs`)

This implementation provides a classic Shamir's Secret Sharing Scheme without password protection:

- **`shares()`**: Generates shares from a secret.
- **`recover_secret()`**: Recovers the secret from the given shares.

**Features:**
- Based on polynomial interpolation using random polynomials.
- Customizable threshold and share number.

### Password-Protected Secret Sharing (`sss_with_password.rs`)

This enhanced implementation adds a password-based protection layer:

- **`shares()`**: Generates shares based on a secret and a password.
- **`recover_secret()`**: Recovers the secret using provided shares and the password.

**Features:**
- Password-based polynomial interpolation for additional security.
- Validates polynomial degree consistency during share generation.

## Usage Example

Here's a simple example demonstrating basic secret sharing:

```rust
use ark_bls12_381::Fq;
use sss::shares;

let secret = Fq::from(1234);
let num_shares = 5;
let threshold = 3;

let shares_generated = shares(secret, num_shares, threshold);
println!("Generated shares: {:?}", shares_generated);
```

For password-protected sharing, include the password parameter:

```rust
use ark_bls12_381::Fq;
use sss_with_password::shares;

let secret = Fq::from(1234);
let password = Fq::from(4321);
let num_shares = 5;
let threshold = 3;

let shares_generated = shares(secret, password, num_shares, threshold);
println!("Generated shares with password: {:?}", shares_generated);
```

## Testing

Run tests to verify implementation correctness:

```bash
cargo test
```

## Dependencies

- [Arkworks (ark_ff, ark_bls12_381)](https://arkworks.rs/): Finite field arithmetic.
- Polynomial implementations from the local `polynomials` crate.

## Contributing

Feel free to open issues or submit pull requests to improve this project.

## License

This project is licensed under the MIT License.

