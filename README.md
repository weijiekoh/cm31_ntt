# cm31_ntt

## Implementations of:

- [x] M31 field arithmetic
- [x] M31 field arithmetic using redundant representation [x] Complex M31 field arithmetic (using the redundant representation of M31s)
- [x] NTT (radix-8)
- [x] Benchmarks
- [ ] Optimsations


## Benchmarks

To run benchmarks:

```bash
rustup default nightly
cargo bench
```

At the time of writing, the results on a Raspberry Pi 5 are:

```
     Running benches/ntt_32768.rs (target/release/deps/ntt_32768-380c493e74327e27)
Gnuplot not found, using plotters backend
ntt_radix_8             time:   [19.014 ms 19.032 ms 19.059 ms]
                        change: [+102.58% +102.92% +103.27%] (p = 0.00 < 0.05)
                        Performance has regressed.
Found 3 outliers among 10 measurements (30.00%)
  2 (20.00%) low mild
  1 (10.00%) high mild

     Running benches/ntt_block_8.rs (target/release/deps/ntt_block_8-3f25dbef695374ee)
Gnuplot not found, using plotters backend
ntt_block_8             time:   [194.70 ns 194.70 ns 194.70 ns]
                        change: [+157.39% +157.56% +157.68%] (p = 0.00 < 0.05)
                        Performance has regressed.
Found 719 outliers among 10000 measurements (7.19%)
  2 (0.02%) low severe
  14 (0.14%) low mild
  395 (3.95%) high mild
  308 (3.08%) high severe
```

### Precomputation of twiddle factors

At the time of writing, the results on a Raspberry Pi 5 are:

| NTT size | Without precomputation | With precomputation |
|-|-|-|
| 512      |         |         |
| 4096     |         |         |
| 32768    |         |         |
