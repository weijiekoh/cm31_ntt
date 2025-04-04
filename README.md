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


### Radix-8 NTT block 
At the time of writing, the results on a Raspberry Pi 5 are:

```
     Running benches/ntt_32768.rs (target/release/deps/ntt_32768-380c493e74327e27)
Gnuplot not found, using plotters backend
ntt_radix_8             time:   [19.014 ms 19.032 ms 19.059 ms]
                        change: [+102.58% +102.92% +103.27%] (p = 0.00 < 0.05)
Found 3 outliers among 10 measurements (30.00%)
  2 (20.00%) low mild
  1 (10.00%) high mild

     Running benches/ntt_block_8.rs (target/release/deps/ntt_block_8-3f25dbef695374ee)
Gnuplot not found, using plotters backend
ntt_block_8             time:   [194.70 ns 194.70 ns 194.70 ns]
                        change: [+157.39% +157.56% +157.68%] (p = 0.00 < 0.05)
Found 719 outliers among 10000 measurements (7.19%)
  2 (0.02%) low severe
  14 (0.14%) low mild
  395 (3.95%) high mild
  308 (3.08%) high severe
```

### Radix-8 NTT with and without twiddle factor precomputation

At the time of writing, the results on a Raspberry Pi 5 are:

```
NTT/size 512 without precomputation
                        time:   [74.598 µs 74.604 µs 74.609 µs]
                        change: [+1.5666% +1.5796% +1.5968%] (p = 0.00 < 0.05)
Found 4 outliers among 100 measurements (4.00%)
  2 (2.00%) high mild
  2 (2.00%) high severe
NTT/size 512 with precomputation
                        time:   [71.535 µs 71.541 µs 71.547 µs]
                        change: [-0.7891% -0.7731% -0.7566%] (p = 0.00 < 0.05)
                        Change within noise threshold.
Found 6 outliers among 100 measurements (6.00%)
  6 (6.00%) high mild
NTT/size 4096 without precomputation
                        time:   [812.02 µs 812.13 µs 812.25 µs]
                        change: [-0.5310% -0.5021% -0.4625%] (p = 0.00 < 0.05)
                        Change within noise threshold.
Found 5 outliers among 100 measurements (5.00%)
  2 (2.00%) high mild
  3 (3.00%) high severe
NTT/size 4096 with precomputation
                        time:   [672.43 µs 672.51 µs 672.60 µs]
                        change: [+0.7131% +0.8691% +0.9617%] (p = 0.00 < 0.05)
                        Change within noise threshold.
Found 10 outliers among 100 measurements (10.00%)
  6 (6.00%) high mild
  4 (4.00%) high severe
NTT/size 32768 without precomputation
                        time:   [7.9226 ms 7.9231 ms 7.9236 ms]
                        change: [-2.3829% -2.3645% -2.3478%] (p = 0.00 < 0.05)
Found 12 outliers among 100 measurements (12.00%)
  6 (6.00%) high mild
  6 (6.00%) high severe
NTT/size 32768 with precomputation
                        time:   [6.2222 ms 6.2230 ms 6.2240 ms]
                        change: [+0.3462% +0.3679% +0.3897%] (p = 0.00 < 0.05)
                        Change within noise threshold.
Found 7 outliers among 100 measurements (7.00%)
  4 (4.00%) high mild
  3 (3.00%) high severe
Benchmarking NTT/size 262144 without precomputation: Warming up for 3.0000 s
NTT/size 262144 without precomputation
                        time:   [91.309 ms 91.396 ms 91.534 ms]
                        change: [+0.9257% +1.0216% +1.1868%] (p = 0.00 < 0.05)
                        Change within noise threshold.
Found 8 outliers among 100 measurements (8.00%)
  7 (7.00%) low mild
  1 (1.00%) high severe
Benchmarking NTT/size 262144 with precomputation: Warming up for 3.0000 s
NTT/size 262144 with precomputation
                        time:   [69.033 ms 69.060 ms 69.089 ms]
                        change: [+0.8209% +0.8727% +0.9229%] (p = 0.00 < 0.05)
                        Change within noise threshold.
```
