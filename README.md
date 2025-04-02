# cm31_ntt

## Implementations of:

- [x] M31 field arithmetic
- [x] M31 field arithmetic using redundant representation
- [x] Complex M31 field arithmetic (using the redundant representation of M31s)
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
     Running benches/ntt_32768.rs (target/release/deps/ntt_32768-28e6d3fe29b7c7a3)
Gnuplot not found, using plotters backend
Benchmarking ntt_radix_8: Collecting 10 samples in estimated 20.145 s (2145 iteration
ntt_radix_8             time:   [9.3906 ms 9.4047 ms 9.4104 ms]
                        change: [+4796428% +4803849% +4810947%] (p = 0.00 < 0.05)
                        Performance has regressed.

     Running benches/ntt_block_8.rs (target/release/deps/ntt_block_8-eb678ddc29dce925)
Gnuplot not found, using plotters backend
Benchmarking ntt_block_8: Collecting 10000 samples in estimated 22.649 s (300M iterat
ntt_block_8             time:   [75.632 ns 75.653 ns 75.675 ns]
Found 450 outliers among 10000 measurements (4.50%)
  338 (3.38%) high mild
  112 (1.12%) high severe
```
