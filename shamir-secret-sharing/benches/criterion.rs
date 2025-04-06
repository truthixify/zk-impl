use criterion::{Criterion, criterion_group, criterion_main};
mod sss;
mod sss_with_password;

use sss::sss_benchmarks;
use sss_with_password::sss_with_password_benchmarks;

criterion_group!(
    name = sss;
    config = Criterion::default().sample_size(100).configure_from_args();
    targets = sss_benchmarks, sss_with_password_benchmarks
);
criterion_main!(sss);
