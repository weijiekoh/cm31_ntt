use std::hint::black_box;
use criterion::{criterion_group, criterion_main, Criterion};
use cm31_ntt::ntt::ntt_radix_8;
use cm31_ntt::cm31::CF;
use num_traits::Zero;
use rand::Rng;
use rand_chacha::ChaCha8Rng;
use rand_chacha::rand_core::SeedableRng;

fn benchmark(c: &mut Criterion) {
    let mut rng = ChaCha8Rng::seed_from_u64(0);
    let mut inputs = [CF::zero(); 8];
    for i in 0..8 {
        inputs[i] = rng.r#gen();
    }

    let w = CF::root_of_unity_8(0).unwrap();
    let w_neg_1 = w.mul_neg_1();

    c.bench_function("ntt_radix_8", |b| b.iter(|| ntt_radix_8(black_box(inputs), black_box(w), black_box(w_neg_1))));
}

criterion_group!(
    name = benches;
    config = Criterion::default()
        .sample_size(10000)
        .warm_up_time(std::time::Duration::new(5,0))
        .measurement_time(std::time::Duration::new(20,0));
    targets = benchmark
);
criterion_main!(benches);
