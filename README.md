# cm31_ntt

## Implementations of:

- [x] M31 field arithmetic
- [x] M31 field arithmetic using redundant representation
- [x] Complex M31 field arithmetic (using the redundant representation of M31s)
- [x] NTT
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
Benchmarking ntt_radix_8: Collecting 10000 samples in estimated 29.314 s (150M iterat
ntt_radix_8             time:   [195.39 ns 195.39 ns 195.40 ns]
                        change: [+0.0323% +0.0495% +0.0811%] (p = 0.03 < 0.05)       
                        Change within noise threshold.                               
Found 538 outliers among 10000 measurements (5.38%)                                  
  1 (0.01%) low severe                                                               
  31 (0.31%) low mild                 
  277 (2.77%) high mild                
  229 (2.29%) high severe
```
