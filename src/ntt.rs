use crate::cm31::{CF, gen_roots_of_unity};
use crate::rm31::RF;
use num_traits::{Zero, One};
use num_traits::pow::Pow;

/// A radix-4 NTT butterfly.
pub fn ntt_block_4(f: [CF; 4], w: CF, w_neg_1: CF) -> [CF; 4] {
    debug_assert_eq!(f.len(), 4);
    debug_assert_eq!(-w, w_neg_1);
    let mut res = [CF::zero(); 4];

    let a0 = f[0] + f[2];
    let a1 = f[0] + f[2].mul_neg_1();
    let a2 = f[1] + f[3];
    let a3 = f[1] + f[3].mul_neg_1();

    res[0] = a0 + a2;
    res[2] = a0 + a2.mul_neg_1();
    res[1] = a1 + a3.mul_j();
    res[3] = a1 + a3.mul_j().mul_neg_1();

    res
}

/// A radix-8 NTT butterfly.
pub fn ntt_block_8(f: [CF; 8], w: CF, w_neg_1: CF) -> [CF; 8] {
    debug_assert_eq!(f.len(), 8);
    debug_assert_eq!(-w, w_neg_1);
    let mut res = [CF::zero(); 8];

    // Refer to Yuval's Radix 8 DIT diagram.
    // 1st columm of black dots: a0-a8
    // 2nd columm of black dots: b0-b8
    // 3nd columm of black dots: res[0]-res[8]

    // Column 1
    let a0 = f[0] + f[4];
    let a1 = f[0] + f[4].mul_neg_1();
    let a2 = f[2] + f[6];
    let a3 = f[2] + f[6].mul_neg_1();
    let a4 = f[1] + f[5];
    let a5 = f[1] + f[5].mul_neg_1();
    let a6 = f[3] + f[7];
    let a7 = f[3] + f[7].mul_neg_1();

    // Column 2
    let b0 = a0 + a2;
    let b1 = a0 + a2.mul_neg_1();
    let b2 = a1 + a3.mul_j();
    let b3 = a1 + a3.mul_j().mul_neg_1();
    let b4 = a4 + a6;
    let b5 = a4 + a6.mul_neg_1();
    let b6 = a5 + a7.mul_j();
    let b7 = a5 + a7.mul_j().mul_neg_1();

    let b5_j = b5.mul_j();

    // Column 3
    res[0] = b0 + b4;
    res[4] = b0 + b4.mul_neg_1();
    res[2] = b1 + b5_j;
    res[6] = b1 + b5_j.mul_neg_1();
    res[1] = b2 + b6 * w;
    res[5] = b2 + b6 * w_neg_1;
    res[3] = b3 + b7.mul_j() * w;
    res[7] = b3 + b7.mul_j() * w_neg_1;

    res
}

pub fn get_root_of_unity(n: usize) -> CF {
    assert!(n.is_power_of_two(), "n must be a power of 2");
    let roots_of_unity = gen_roots_of_unity((n as f64).log2() as usize);
    roots_of_unity[roots_of_unity.len() - 1]
}

pub fn naive_ntt(f: Vec<CF>) -> Vec<CF> {
    let n = f.len();
    let wn = get_root_of_unity(n);
    let mut res = vec![CF::zero(); f.len()];
    for i in 0..n {
        for j in 0..n {
            res[i] += f[j] * wn.pow(i * j);
        }
    }

    res
}

pub fn naive_intt(f: Vec<CF>) -> Vec<CF> {
    let n = f.len();
    let wn = get_root_of_unity(n);

    let mut res = vec![CF::zero(); n];
    for i in 0..n {
        for j in 0..n {
            res[i] += f[j] * wn.pow(i * j).try_inverse().unwrap();
        }
    }
    for i in 0..n {
        res[i] = res[i].mul_by_f(RF::new(n as u32).try_inverse().unwrap());
    }

    res
}

pub fn naive_ntt_radix_4(f: [CF; 4]) -> [CF; 4] {
    naive_ntt(f.to_vec()).try_into().unwrap()
}

pub fn naive_intt_radix_4(f: [CF; 4]) -> [CF; 4] {
    naive_intt(f.to_vec()).try_into().unwrap()
}

pub fn naive_ntt_radix_8(f: [CF; 8]) -> [CF; 8] {
    naive_ntt(f.to_vec()).try_into().unwrap()
}

pub fn naive_intt_radix_8(f: [CF; 8]) -> [CF; 8] {
    naive_intt(f.to_vec()).try_into().unwrap()
}

pub fn naive_ntt_16(f: [CF; 16]) -> [CF; 16] {
    naive_ntt(f.to_vec()).try_into().unwrap()
}

pub fn naive_intt_16(f: [CF; 16]) -> [CF; 16] {
    naive_intt(f.to_vec()).try_into().unwrap()
}


// TODO: deprecate
pub fn ntt_16(f: [CF; 16]) -> [CF; 16] {
    let mut res = [CF::zero(); 16];
    let w8 = CF::root_of_unity_8(0).unwrap();
    let w16 = w8.try_sqrt().unwrap();

    // Split input into even and odd indices
    let mut f_even = [CF::zero(); 8];
    let mut f_odd = [CF::zero(); 8];
    for i in 0..8 {
        f_even[i] = f[2*i];
        f_odd[i] = f[2*i + 1];
    }

    // Compute 8-point NTTs on even and odd parts
    let ntt_even = ntt_block_8(f_even, w8, w8.mul_neg_1());
    let ntt_odd = ntt_block_8(f_odd, w8, w8.mul_neg_1());

    // Combine the results with twiddle factors
    for i in 0..8 {
        let twiddle = w16.pow(i);
        res[i] = ntt_even[i] + twiddle * ntt_odd[i];
        res[i + 8] = ntt_even[i] + twiddle.mul_neg_1() * ntt_odd[i];
    }

    res
}

// TODO: deprecate
pub fn ntt_32(f: [CF; 32]) -> [CF; 32] {
    let mut res = [CF::zero(); 32];
    let w8 = CF::root_of_unity_8(0).unwrap();
    let w16 = w8.try_sqrt().unwrap();
    let w32 = w16.try_sqrt().unwrap();
    
    // Split input into 4 parts for radix-4 approach using 8-point NTTs
    let mut f_0 = [CF::zero(); 8];
    let mut f_1 = [CF::zero(); 8];
    let mut f_2 = [CF::zero(); 8];
    let mut f_3 = [CF::zero(); 8];
    
    // Distribute elements to the 4 parts (stride of 4)
    for i in 0..8 {
        f_0[i] = f[4*i];
        f_1[i] = f[4*i + 1];
        f_2[i] = f[4*i + 2];
        f_3[i] = f[4*i + 3];
    }
    
    // Compute 8-point NTTs on each part
    let ntt_0 = ntt_block_8(f_0, w8, w8.mul_neg_1());
    let ntt_1 = ntt_block_8(f_1, w8, w8.mul_neg_1());
    let ntt_2 = ntt_block_8(f_2, w8, w8.mul_neg_1());
    let ntt_3 = ntt_block_8(f_3, w8, w8.mul_neg_1());
    
    // Combine the results with twiddle factors
    for i in 0..8 {
        // w^0 = 1 (no need to explicitly calculate)
        let w_i_1 = w32.pow(1 * i);      // w^i
        let w_i_2 = w32.pow(2 * i);      // w^(2i)
        let w_i_3 = w32.pow(3 * i);      // w^(3i)
        
        // First 8 points (k=0)
        res[i] = ntt_0[i] + w_i_1 * ntt_1[i] + w_i_2 * ntt_2[i] + w_i_3 * ntt_3[i];
        
        // Second 8 points (k=1)
        let k1_idx = i + 8;
        res[k1_idx] = ntt_0[i] + w_i_1 * w32.pow(8) * ntt_1[i] + 
                       w_i_2 * w32.pow(16) * ntt_2[i] + w_i_3 * w32.pow(24) * ntt_3[i];
        
        // Third 8 points (k=2)
        let k2_idx = i + 16;
        res[k2_idx] = ntt_0[i] + w_i_1 * w32.pow(16) * ntt_1[i] + 
                       w_i_2 * w32.pow(32) * ntt_2[i] + w_i_3 * w32.pow(48) * ntt_3[i];
        
        // Fourth 8 points (k=3)
        let k3_idx = i + 24;
        res[k3_idx] = ntt_0[i] + w_i_1 * w32.pow(24) * ntt_1[i] + 
                       w_i_2 * w32.pow(48) * ntt_2[i] + w_i_3 * w32.pow(72) * ntt_3[i];
    }
    
    res
}

pub fn ntt_radix_2(f: Vec<CF>, w: CF) -> Vec<CF> {
    let n = f.len();
    assert!(n >= 2, "n must be at least 2");
    assert!(n.is_power_of_two(), "n must be a power of 2");

    fn do_ntt(f: Vec<CF>, w: CF) -> Vec<CF> {
        let n = f.len();
        if n == 1 {
            return f;
        }

        // Divide
        let mut f_even = vec![CF::zero(); n/2];
        let mut f_odd = vec![CF::zero(); n/2];
        for i in 0..n/2 {
            f_even[i] = f[2*i];
            f_odd[i] = f[2*i + 1];
        }

        // Recurse
        let ntt_even = do_ntt(f_even, w.pow(2));
        let ntt_odd = do_ntt(f_odd, w.pow(2));

        // Combine
        let mut res = vec![CF::zero(); n];

        let mut wk = CF::one();
        for i in 0..n/2 {
            res[i] = ntt_even[i] + wk * ntt_odd[i];
            res[i + n/2] = ntt_even[i] + wk.mul_neg_1() * ntt_odd[i];
            wk = wk * w;
        }

        res
    }
    do_ntt(f, w)
}

fn is_power_of_8(n: u32) -> bool {
    if n == 0 {
        return false;
    }

    let mut num = n;
    while num % 8 == 0 {
        num /= 8;
    }

    num == 1
}

pub fn ntt_radix_8(f: Vec<CF>, w: CF, w8: CF, w8_neg_1: CF) -> Vec<CF> {
    let n = f.len();
    assert!(n >= 8, "n must be at least 8");
    assert!(is_power_of_8(n as u32), "n must be a power of 8");

    fn do_ntt(f: Vec<CF>, w: CF, w8: CF, w8_neg_1: CF) -> Vec<CF> {
        let n = f.len();
        if n == 1 {
            return f;
        }

        // n is divisible by 8.
        let m = n / 8;

        // Partition f into eight subarrays.
        let mut a0 = vec![CF::zero(); m];
        let mut a1 = vec![CF::zero(); m];
        let mut a2 = vec![CF::zero(); m];
        let mut a3 = vec![CF::zero(); m];
        let mut a4 = vec![CF::zero(); m];
        let mut a5 = vec![CF::zero(); m];
        let mut a6 = vec![CF::zero(); m];
        let mut a7 = vec![CF::zero(); m];
        for i in 0..m {
            a0[i] = f[8 * i];
            a1[i] = f[8 * i + 1];
            a2[i] = f[8 * i + 2];
            a3[i] = f[8 * i + 3];
            a4[i] = f[8 * i + 4];
            a5[i] = f[8 * i + 5];
            a6[i] = f[8 * i + 6];
            a7[i] = f[8 * i + 7];
        }

        // Recurse
        let ntt_a0 = do_ntt(a0, w.pow(8), w8, w8_neg_1);
        let ntt_a1 = do_ntt(a1, w.pow(8), w8, w8_neg_1);
        let ntt_a2 = do_ntt(a2, w.pow(8), w8, w8_neg_1);
        let ntt_a3 = do_ntt(a3, w.pow(8), w8, w8_neg_1);
        let ntt_a4 = do_ntt(a4, w.pow(8), w8, w8_neg_1);
        let ntt_a5 = do_ntt(a5, w.pow(8), w8, w8_neg_1);
        let ntt_a6 = do_ntt(a6, w.pow(8), w8, w8_neg_1);
        let ntt_a7 = do_ntt(a7, w.pow(8), w8, w8_neg_1);

        let mut res = vec![CF::zero(); n];
        for k in 0..m {
            // Twiddle factors
            let wt   = w.pow(k);
            let wt2  = wt * wt;
            let wt3  = wt2 * wt;
            let wt4  = wt3 * wt;
            let wt5  = wt4 * wt;
            let wt6  = wt5 * wt;
            let wt7  = wt6 * wt;

            // Apply twiddle factors to each of the eight sub-transform outputs.
            let t0 = ntt_a0[k];
            let t1 = wt * ntt_a1[k];
            let t2 = wt2 * ntt_a2[k];
            let t3 = wt3 * ntt_a3[k];
            let t4 = wt4 * ntt_a4[k];
            let t5 = wt5 * ntt_a5[k];
            let t6 = wt6 * ntt_a6[k];
            let t7 = wt7 * ntt_a7[k];

            // Collect the eight values.
            let ts = [t0, t1, t2, t3, t4, t5, t6, t7];

            // Use the provided 8-point butterfly.
            let butterfly = ntt_block_8(ts, w8, w8_neg_1);

            // Write the butterfly result into the correct positions.
            for r in 0..8 {
                res[k + r * m] = butterfly[r];
            }
        }
        res
    }

    do_ntt(f, w, w8, w8_neg_1)
}


#[cfg(test)]
pub mod tests {
    use crate::ntt::*;
    use crate::cm31::CF;
    use num_traits::{Zero, Pow};
    use rand::Rng;
    use rand_chacha::ChaCha8Rng;
    use rand_chacha::rand_core::SeedableRng;

    // Schoolbook multiplication
    fn naive_poly_mul_2(f1: [CF; 2], f2: [CF; 2]) -> [CF; 3] {
        let mut res = [CF::zero(); 3];
        for i in 0..2 {
            for j in 0..2 {
                res[i + j] += f1[i] * f2[j];
            }
        }
        res
    }

    // Schoolbook multiplication
    fn naive_poly_mul_4(f1: [CF; 4], f2: [CF; 4]) -> [CF; 7] {
        let mut res = [CF::zero(); 7];
        for i in 0..4 {
            for j in 0..4 {
                res[i + j] += f1[i] * f2[j];
            }
        }
        res
    }

    fn naive_poly_mul_8(f1: [CF; 8], f2: [CF; 8]) -> [CF; 15] {
        let mut res = [CF::zero(); 15];
        for i in 0..8 {
            for j in 0..8 {
                res[i + j] += f1[i] * f2[j];
            }
        }
        res
    }

    #[test]
    fn test_naive_poly_mul_4() {
        let f1 = [CF::new(1, 0), CF::new(2, 0), CF::new(3, 0), CF::new(4, 0)];
        let f3 = [CF::new(1, 0), CF::new(3, 0), CF::new(5, 0), CF::new(7, 0)];
        let res = naive_poly_mul_4(f1, f3);
        let expected = [CF::new(1, 0), CF::new(5, 0), CF::new(14, 0), CF::new(30, 0), CF::new(41, 0), CF::new(41, 0), CF::new(28, 0)];
        assert_eq!(res, expected);
    }

    #[test]
    fn test_naive_ntt_4() {
        // Test that the naive NTT and inverse NTT functions are the opposite of each other.
        let mut rng = ChaCha8Rng::seed_from_u64(0);
        for _ in 0..128 {
            let mut f = [CF::zero(); 4];
            for i in 0..4 {
                f[i] = rng.r#gen();
            }

            let res = naive_ntt_radix_4(f);
            let ires = naive_intt_radix_4(res);

            assert_eq!(ires.to_vec(), f);
        }
    }

    #[test]
    fn test_naive_ntt_8() {
        // Test that the naive NTT and inverse NTT functions are the opposite of each other.
        let mut rng = ChaCha8Rng::seed_from_u64(0);
        for _ in 0..128 {
            let mut f = [CF::zero(); 8];
            for i in 0..8 {
                f[i] = rng.r#gen();
            }

            let res = naive_ntt_radix_8(f);
            let ires = naive_intt_radix_8(res);

            assert_eq!(ires, f);
        }
    }

    #[test]
    fn test_ntt_4_by_property() {
        // Test the correctness of the native NTT and inverse NTT functions.
        let mut rng = ChaCha8Rng::seed_from_u64(0);

        for _ in 0..128 {
            let mut poly1 = [CF::zero(); 2];
            let mut poly2 = [CF::zero(); 2];
            for i in 0..2 {
                poly1[i] = rng.r#gen();
                poly2[i] = rng.r#gen();
            }

            let mut poly1_padded = [CF::zero(); 4];
            let mut poly2_padded = [CF::zero(); 4];

            for i in 0..2 {
                poly1_padded[i] = poly1[i];
                poly2_padded[i] = poly2[i];
            }

            let poly1_ntt = naive_ntt_radix_4(poly1_padded);
            let poly2_ntt = naive_ntt_radix_4(poly2_padded);
            let mut product_ntt = [CF::zero(); 4];

            for i in 0..4 {
                product_ntt[i] = poly1_ntt[i] * poly2_ntt[i];
            }
            let product_poly = naive_intt_radix_4(product_ntt);
            let expected_product = naive_poly_mul_2(poly1, poly2);

            for i in 0..expected_product.len() {
                assert_eq!(product_poly[i], expected_product[i]);
            }
        }
    }

    #[test]
    fn test_naive_ntt_8_by_property() {
        // Test the correctness of the native NTT and inverse NTT functions.
        let mut rng = ChaCha8Rng::seed_from_u64(0);

        for _ in 0..128 {
            let mut poly1 = [CF::zero(); 4];
            let mut poly2 = [CF::zero(); 4];
            for i in 0..4 {
                poly1[i] = rng.r#gen();
                poly2[i] = rng.r#gen();
            }

            let mut poly1_padded = [CF::zero(); 8];
            let mut poly2_padded = [CF::zero(); 8];

            for i in 0..4 {
                poly1_padded[i] = poly1[i];
                poly2_padded[i] = poly2[i];
            }

            let poly1_ntt = naive_ntt_radix_8(poly1_padded);
            let poly2_ntt = naive_ntt_radix_8(poly2_padded);
            let mut product_ntt = [CF::zero(); 8];

            for i in 0..8 {
                product_ntt[i] = poly1_ntt[i] * poly2_ntt[i];
            }
            let product_poly = naive_intt_radix_8(product_ntt);
            let expected_product = naive_poly_mul_4(poly1, poly2);

            for i in 0..expected_product.len() {
                assert_eq!(product_poly[i], expected_product[i]);
            }
        }
    }

    #[test]
    fn test_ntt_16() {
        // Test that the naive NTT and inverse NTT functions are the opposite of each other.
        let mut rng = ChaCha8Rng::seed_from_u64(0);
        for _ in 0..128 {
            let mut f = [CF::zero(); 16];
            for i in 0..16 {
                f[i] = rng.r#gen();
            }

            let res = naive_ntt_16(f);
            let res2 = ntt_16(f);
            let ires = naive_intt_16(res);

            assert_eq!(ires, f);
            assert_eq!(res2, res);
        }
    }
    
    #[test]
    fn test_ntt_32() {
        // Test the 32-point NTT implementation against a naive implementation
        let mut rng = ChaCha8Rng::seed_from_u64(0);
        for _ in 0..8 {  // Reduced iterations due to size of computation
            let mut f = [CF::zero(); 32];
            for i in 0..32 {
                f[i] = rng.r#gen();
            }
            
            // Create a naive 32-point NTT for comparison
            let w8 = CF::root_of_unity_8(0).unwrap();
            let w16 = w8.try_sqrt().unwrap();
            let w32 = w16.try_sqrt().unwrap();
            
            let mut naive_result = [CF::zero(); 32];
            for i in 0..32 {
                for j in 0..32 {
                    naive_result[i] += f[j] * w32.pow(i * j);
                }
            }
            
            // Compare with our optimized implementation
            let result = ntt_32(f);
            assert_eq!(result, naive_result);
        }
    }

    #[test]
    fn test_ntt_16_by_property() {
        // Test the correctness of the native NTT and inverse NTT functions.
        let mut rng = ChaCha8Rng::seed_from_u64(0);

        for _ in 0..128 {
            let mut poly1 = [CF::zero(); 8];
            let mut poly2 = [CF::zero(); 8];
            for i in 0..8 {
                poly1[i] = rng.r#gen();
                poly2[i] = rng.r#gen();
            }

            let mut poly1_padded = [CF::zero(); 16];
            let mut poly2_padded = [CF::zero(); 16];

            for i in 0..8 {
                poly1_padded[i] = poly1[i];
                poly2_padded[i] = poly2[i];
            }

            let poly1_ntt = ntt_16(poly1_padded);
            let poly2_ntt = ntt_16(poly2_padded);
            let mut product_ntt = [CF::zero(); 16];

            for i in 0..16 {
                product_ntt[i] = poly1_ntt[i] * poly2_ntt[i];
            }
            let product_poly = naive_intt_16(product_ntt);
            let expected_product = naive_poly_mul_8(poly1, poly2);

            for i in 0..expected_product.len() {
                assert_eq!(product_poly[i], expected_product[i]);
            }
        }
    }

    #[test]
    fn test_ntt_block_4() {
        // Test the radix-4 NTT butterfly.
        let w = CF::root_of_unity_4(0).unwrap();
        let w_neg_1 = w.mul_neg_1();
        let mut rng = ChaCha8Rng::seed_from_u64(0);
        for _ in 0..1024 {
            let mut poly = [CF::zero(); 4];
            for j in 0..4 {
                poly[j] = rng.r#gen();
            }

            let naive = naive_ntt_radix_4(poly);
            let res = ntt_block_4(poly, w, w_neg_1);
            assert_eq!(naive, res);

            let ires = naive_intt_radix_4(res);
            assert_eq!(ires, poly);
        }
    }

    #[test]
    fn test_ntt_block_8() {
        // Test the radix-8 NTT butterfly.
        let w = CF::root_of_unity_8(0).unwrap();
        let w_neg_1 = w.mul_neg_1();
        let mut rng = ChaCha8Rng::seed_from_u64(0);
        for _ in 0..1024 {
            let mut poly = [CF::zero(); 8];
            for j in 0..8 {
                poly[j] = rng.r#gen();
            }

            let naive = naive_ntt_radix_8(poly);
            let res = ntt_block_8(poly, w, w_neg_1);
            assert_eq!(naive, res);

            let ires = naive_intt_radix_8(res);
            assert_eq!(ires, poly);
        }
    }

    #[test]
    fn test_ntt_radix_2() {
        for i in 2..10 {
            let n = 1 << i;
            let wn = get_root_of_unity(n);
            let mut rng = ChaCha8Rng::seed_from_u64(0);
            let mut f = vec![CF::zero(); n];
            for i in 0..n {
                f[i] = rng.r#gen();
            }
            let res = ntt_radix_2(f.clone(), wn);
            let naive_res = naive_ntt(f.clone());
            assert_eq!(res, naive_res);
        }
    }

    #[test]
    fn test_ntt_radix_8() {
        let n = 64;
        let mut rng = ChaCha8Rng::seed_from_u64(0);
        let mut f = vec![CF::zero(); n];
        for i in 0..n {
            f[i] = rng.r#gen();
        }
        let w = get_root_of_unity(n);
        let w8 = CF::root_of_unity_8(0).unwrap();
        let w8_neg_1 = w8.mul_neg_1();
        let res = ntt_radix_8(f.clone(), w, w8, w8_neg_1);
        let naive_res = naive_ntt(f.clone());
        assert_eq!(res, naive_res);
    }
}
