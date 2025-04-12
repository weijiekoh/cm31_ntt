use crate::cm31::{
    CF,
    gen_roots_of_unity,
    W_8,
};
use crate::rm31::RF;
use num_traits::Zero;
use num_traits::pow::Pow;

//pub(crate) fn log_8(n: usize) -> usize {
    //let mut log = 0;
    //let mut m = n;
    //while m > 1 {
        //m /= 8;
        //log += 1;
    //}
    //log
//}

pub(crate) fn is_power_of_8(n: u32) -> bool {
    if n == 0 {
        return false;
    }

    let mut num = n;
    while num % 8 == 0 {
        num /= 8;
    }

    num == 1
}


/// Helper function to get the n-th root of unity in CF.
pub fn get_root_of_unity(n: usize) -> CF {
    assert!(n.is_power_of_two(), "n must be a power of 2");
    let roots_of_unity = gen_roots_of_unity((n as f64).log2() as usize);
    roots_of_unity[roots_of_unity.len() - 1]
}

pub(crate) fn level_offset(overall_transform_size: usize, d: usize) -> usize {
    let mut offset = 0;
    let mut current = overall_transform_size;
    for _ in 0..d {
        offset += 1 + 7 * (current / 8);
        current /= 8;
    }
    offset
}

/// A radix-4 NTT butterfly.
pub fn ntt_block_4(f: [CF; 4]) -> [CF; 4] {
    debug_assert_eq!(f.len(), 4);
    let mut res = [CF::zero(); 4];

    let a0 = f[0] + f[2];
    let a1 = f[0] - f[2];
    let a2 = f[1] + f[3];
    let a3 = f[1] - f[3];

    let a3_j = a3.mul_j();

    res[0] = a0 + a2;
    res[2] = a0 - a2;
    res[1] = a1 + a3_j;
    res[3] = a1 - a3_j;

    res
}

/// A radix-8 NTT butterfly.
#[inline]
pub fn ntt_block_8(f: [CF; 8]) -> [CF; 8] {
    debug_assert_eq!(f.len(), 8);
    let mut res = [CF::zero(); 8];

    // Refer to Yuval's Radix 8 DIT diagram.
    // 1st columm of black dots: a0-a8
    // 2nd columm of black dots: b0-b8
    // 3nd columm of black dots: res[0]-res[8]

    // Column 1
    let a0 = f[0] + f[4];
    let a1 = f[0] - f[4];
    let a2 = f[2] + f[6];
    let a3 = f[2] - f[6];
    let a4 = f[1] + f[5];
    let a5 = f[1] - f[5];
    let a6 = f[3] + f[7];
    let a7 = f[3] - f[7];

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

    res[0] = b0 + b4;
    res[4] = b0 - b4;
    res[2] = b1 + b5_j;
    res[6] = b1 - b5_j;
    res[1] = b2 + b6_w8;
    res[5] = b2 - b6_w8;
    res[3] = b3 + b7_j_w8;
    res[7] = b3 - b7_j_w8;

    res
}

/// A helper function to perform the NTT in a very simple but unoptimised O(n^2) way to test for
/// correctness of the optimised NTT functions.
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

/// A helper function to perform the inverse NTT in a very simple but unoptimised O(n^2) way to
/// help test for correctness of the NTT.
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
