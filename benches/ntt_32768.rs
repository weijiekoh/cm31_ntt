use std::hint::black_box;
use criterion::{criterion_group, criterion_main, Criterion};
use cm31_ntt::ntt::{ntt_radix_8, get_root_of_unity};
use cm31_ntt::cm31::CF;
use num_traits::Zero;
use rand::Rng;
use rand_chacha::ChaCha8Rng;
use rand_chacha::rand_core::SeedableRng;

fn benchmark(c: &mut Criterion) {
    let n = 32768;
    let mut rng = ChaCha8Rng::seed_from_u64(0);
    let mut inputs = vec![CF::zero(); 32768];
    for i in 0..n {
        inputs[i] = rng.r#gen();
    }

    let w = get_root_of_unity(n);
    let w8 = CF::root_of_unity_8(0).unwrap();
    let w8_neg_1 = w8.mul_neg_1();

    c.bench_function("ntt_radix_8", |b| b.iter(|| ntt_radix_8(black_box(inputs.clone()), w, w8, w8_neg_1)));
}

criterion_group!(
    name = benches;
    config = Criterion::default()
        .sample_size(10)
        .warm_up_time(std::time::Duration::new(5,0))
        .measurement_time(std::time::Duration::new(20,0));
    targets = benchmark
);
criterion_main!(benches);
