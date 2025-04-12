use crate::cm31::{CF, W_8};
use num_traits::Zero;
use crate::ntt_utils::*;
//use crate::ntt_utils::ntt_block_8;

/// A radix-8 NTT butterfly.
#[inline]
pub fn ntt_block_8_in_place<const N: usize>(
    read_from: &mut Box<[CF; N]>,
    write_to: &mut Box<[CF; N]>,
    read_indices: &[usize; 8],
    write_indices: &[usize; 8],
    wt: CF,
    wt2: CF,
    wt3: CF,
    wt4: CF,
    wt5: CF,
    wt6: CF,
    wt7: CF,
) {
    // Refer to Yuval's Radix 8 DIT diagram.
    // 1st columm of black dots: a0-a8
    // 2nd columm of black dots: b0-b8
    // 3nd columm of black dots: res[0]-res[8]
    //
    let t0 = read_from[read_indices[0]];
    let t1 = read_from[read_indices[1]] * wt;
    let t2 = read_from[read_indices[2]] * wt2;
    let t3 = read_from[read_indices[3]] * wt3;
    let t4 = read_from[read_indices[4]] * wt4;
    let t5 = read_from[read_indices[5]] * wt5;
    let t6 = read_from[read_indices[6]] * wt6;
    let t7 = read_from[read_indices[7]] * wt7;

    // Column 1
    let a0 = t0 + t4;
    let a1 = t0 - t4;
    let a2 = t2 + t6;
    let a3 = t2 - t6;
    let a4 = t1 + t5;
    let a5 = t1 - t5;
    let a6 = t3 + t7;
    let a7 = t3 - t7;

    // Column 2
    let a3_j = a3.mul_j();
    let a7_j = a7.mul_j();

    let b0 = a0 + a2;
    let b1 = a0 - a2;
    let b2 = a1 + a3_j;
    let b3 = a1 - a3_j;
    let b4 = a4 + a6;
    let b5 = a4 - a6;
    let b6 = a5 + a7_j;
    let b7 = a5 - a7_j;

    // Column 3
    let b5_j = b5.mul_j();
    let b7_j = b7.mul_j();
    let b6_w8 = b6 * W_8;
    let b7_j_w8 = b7_j * W_8;

    write_to[write_indices[0]] = b0 + b4;
    write_to[write_indices[4]] = b0 - b4;
    write_to[write_indices[2]] = b1 + b5_j;
    write_to[write_indices[6]] = b1 - b5_j;
    write_to[write_indices[1]] = b2 + b6_w8;
    write_to[write_indices[5]] = b2 - b6_w8;
    write_to[write_indices[3]] = b3 + b7_j_w8;
    write_to[write_indices[7]] = b3 - b7_j_w8;
}

/// An in-place radix-8 NTT with precomputed twiddles.
pub fn ntt<const N: usize>(
    f: &mut Box<[CF; N]>,
    twiddles: &Vec<CF>,
) {
    debug_assert!(N >= 8, "N must be at least 8");
    debug_assert!(is_power_of_8(N as u32), "N must be a power of 8");

    let mut scratch = Box::new([CF::zero(); N]);

    fn do_ntt<const N: usize>(
        f: &mut Box<[CF; N]>,
        scratch: &mut Box<[CF; N]>,
        twiddles: &Vec<CF>,
        offset: usize,
        stride: usize,
        n: usize,
        //w: CF,
        depth: usize,
        overall_transform_size: usize,
    ) {
        if n == 1 {
            return;
        }

        let m = n / 8;

        for r in 0..8 {
            do_ntt(f, scratch, twiddles, offset + r * stride, stride * 8, m, depth + 1, overall_transform_size);
        }

        for i in 0..n {
            scratch[i] = f[offset + i * stride];
        }

        let level_size = 1 + 7 * m;
        let lvl_offset = level_offset(overall_transform_size, depth);
        let level_twiddles = &twiddles[lvl_offset..lvl_offset + level_size];

        for k in 0..m {
            let base_idx = 1 + 7 * k;
            let wt = level_twiddles[base_idx];
            let wt2 = level_twiddles[base_idx + 1];
            let wt3 = level_twiddles[base_idx + 2];
            let wt4 = level_twiddles[base_idx + 3];
            let wt5 = level_twiddles[base_idx + 4];
            let wt6 = level_twiddles[base_idx + 5];
            let wt7 = level_twiddles[base_idx + 6];

            ntt_block_8_in_place(
                scratch,
                f,
                &[
                    0 + 8 * k,
                    1 + 8 * k,
                    2 + 8 * k,
                    3 + 8 * k,
                    4 + 8 * k,
                    5 + 8 * k,
                    6 + 8 * k,
                    7 + 8 * k
                ],
                &[
                    offset + (k + 0 * m) * stride,
                    offset + (k + 1 * m) * stride,
                    offset + (k + 2 * m) * stride,
                    offset + (k + 3 * m) * stride,
                    offset + (k + 4 * m) * stride,
                    offset + (k + 5 * m) * stride,
                    offset + (k + 6 * m) * stride,
                    offset + (k + 7 * m) * stride
                ],
                wt,
                wt2,
                wt3,
                wt4,
                wt5,
                wt6,
                wt7,
            );
        }
    }

    do_ntt(&mut *f, &mut scratch, twiddles, 0, 1, N, /*w,*/ 0, N);
}


#[cfg(test)]
pub mod tests {
    use crate::ntt::*;
    use crate::ntt_unoptimised::*;
    use crate::cm31::CF;
    use num_traits::Zero;
    use rand::Rng;
    use rand_chacha::ChaCha8Rng;
    use rand_chacha::rand_core::SeedableRng;

    #[test]
    fn test_ntt() {
        let mut rng = ChaCha8Rng::seed_from_u64(0);

        let n = 512;
        let w = get_root_of_unity(n);
        let radix = 8;
        let twiddles = precompute_twiddles(n, w, radix);

        for _ in 0..1 {
            let mut f = Box::new([CF::zero(); 512]);
            for i in 0..n {
                f[i] = rng.r#gen();
            }

            let expected = ntt_radix_8(f.clone().to_vec(), w);

            ntt::<512>(&mut f, &twiddles);

            let is_correct = f.to_vec() == expected;
            assert!(is_correct);
        }
    }
}
