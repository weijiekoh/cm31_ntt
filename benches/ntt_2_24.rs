use std::hint::black_box;
use criterion::{criterion_group, criterion_main, Criterion};
use cm31_ntt::ntt::{ntt_radix_8, ntt_radix_8_precomputed, get_root_of_unity, precompute_twiddles};
use cm31_ntt::cm31::CF;
use num_traits::Zero;
use rand::Rng;
use rand_chacha::ChaCha8Rng;
use rand_chacha::rand_core::SeedableRng;

fn bench_2_24(c: &mut Criterion) {
    let radix = 8;

    let n = (8usize).pow(8);
    let w = get_root_of_unity(n as usize);
    let twiddles = precompute_twiddles(n as usize, w, radix);

    let mut group = c.benchmark_group("NTT (2^24)");

    let mut rng = ChaCha8Rng::seed_from_u64(0);

    let mut f = vec![CF::zero(); n];
    for i in 0..n {
        f[i] = rng.r#gen();
    }

    let w8 = CF::root_of_unity_8(0).unwrap();
    let w8_neg_1 = w8.mul_neg_1();


    group.bench_function(format!("size {n} without precomputation"), |b| {
        b.iter(|| {
            let f_clone = f.clone();
            ntt_radix_8(black_box(f_clone), w, w8, w8_neg_1);
        })
    });
    group.bench_function(format!("size {n} with precomputation"), |b| {
        b.iter(|| {
            let f_clone = f.clone();
            ntt_radix_8_precomputed(black_box(f_clone), &twiddles, w8, w8_neg_1);
        })
    });
    group.finish();
}

criterion_group! {
    name = benches;
    config = Criterion::default().sample_size(10);
    targets = bench_2_24
}
criterion_main!(benches);

