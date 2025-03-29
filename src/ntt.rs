use crate::cm31::CF;
use crate::rm31::RF;
use num_traits::Zero;
use num_traits::pow::Pow;

/// A radix-8 NTT butterfly.
// TODO: optimise by deferring reduction.
pub fn ntt_radix_8(f: [CF; 8], w: CF, w_neg_1: CF) -> [CF; 8] {
    debug_assert_eq!(f.len(), 8);
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

pub fn naive_ntt_radix_8(f: [CF; 8]) -> [CF; 8] {
    assert_eq!(f.len(), 8);

    let mut res = [CF::zero(); 8];

    let w8 = CF::root_of_unity_8(0).unwrap();
    
    for i in 0..8 {
        for j in 0..8 {
            res[i] += f[j] * w8.pow(i * j);
        }
    }

    res
}

pub fn naive_intt_radix_8(f: [CF; 8]) -> [CF; 8] {
    assert_eq!(f.len(), 8);

    let mut res = [CF::zero(); 8];

    let w8 = CF::root_of_unity_8(0).unwrap();
    
    for i in 0..8 {
        for j in 0..8 {
            res[i] += f[j] * w8.pow(i * j).try_inverse().unwrap();
        }
    }

    for i in 0..8 {
        res[i] = res[i].mul_by_f(RF::new(8).try_inverse().unwrap());
    }

    res
}

#[cfg(test)]
pub mod tests {
    use crate::ntt::{ntt_radix_8, naive_ntt_radix_8, naive_intt_radix_8};
    use crate::cm31::CF;
    use num_traits::{Zero, Pow};
    use rand::Rng;
    use rand_chacha::ChaCha8Rng;
    use rand_chacha::rand_core::SeedableRng;

    // Schoolbook multiplication
    fn naive_poly_mul(f1: [CF; 4], f2: [CF; 4]) -> [CF; 7] {
        let mut res = [CF::zero(); 7];
        for i in 0..4 {
            for j in 0..4 {
                res[i + j] += f1[i] * f2[j];
            }
        }
        res
    }

    #[test]
    fn test_naive_poly_mul() {
        let f1 = [CF::new(1, 0), CF::new(2, 0), CF::new(3, 0), CF::new(4, 0)];
        let f3 = [CF::new(1, 0), CF::new(3, 0), CF::new(5, 0), CF::new(7, 0)];
        let res = naive_poly_mul(f1, f3);
        let expected = [CF::new(1, 0), CF::new(5, 0), CF::new(14, 0), CF::new(30, 0), CF::new(41, 0), CF::new(41, 0), CF::new(28, 0)];
        assert_eq!(res, expected);
    }

    #[test]
    fn test_naive_ntt() {
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
    fn test_naive_ntt_by_property() {
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
            let expected_product = naive_poly_mul(poly1, poly2);

            for i in 0..expected_product.len() {
                assert_eq!(product_poly[i], expected_product[i]);
            }
        }
    }

    #[test]
    fn test_ntt_radix_8() {
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
            let res = ntt_radix_8(poly, w, w_neg_1);
            assert_eq!(naive, res);
        }
    }

    #[test]
    fn test_opts() {
        // Test the optimized functions.
        let w = CF::root_of_unity_8(0).unwrap();
        let j = w.pow(2);
        let neg_1 = w.pow(4);

        let v: CF = CF::new(0x12345678, 0x87654321);
        let v_neg_1 = v.mul_neg_1();
        assert_eq!(v * neg_1, v_neg_1);

        let v_j = v.mul_j();
        assert_eq!(v * j, v_j);
    }
}
