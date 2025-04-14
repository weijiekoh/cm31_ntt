use std::hint::black_box;
use criterion::{criterion_group, criterion_main, Criterion};
use cm31_ntt::ntt::ntt;
use cm31_ntt::ntt_utils::get_root_of_unity;
use cm31_ntt::ntt_unoptimised::{ntt_radix_8_precomputed, precompute_twiddles};
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

    let mut arr = Box::new([CF::zero(); 16777216]);
    let mut g = vec![CF::zero(); 16777216];
    for i in 0..n {
        arr[i] = rng.r#gen();
        g[i] = arr[i];
    }

    group.bench_function(format!("size {n} with precomputation (in-place)"), |b| {
        b.iter(|| {
            ntt::<16777216>(black_box(&mut arr), &twiddles);
        })
    });
    group.bench_function(format!("size {n} with precomputation (Vec)"), |b| {
        b.iter(|| {
            let g_clone = g.clone();
            ntt_radix_8_precomputed(black_box(g_clone), &twiddles);
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
