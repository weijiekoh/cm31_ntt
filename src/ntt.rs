use crate::cm31::{
    CF,
    gen_roots_of_unity,
    W_8,
};
use crate::rm31::RF;
use num_traits::{Zero, One};
use num_traits::pow::Pow;

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

/// A radix-2 NTT.
/// @param f The coefficients of the polynomial to be transformed.
/// @param w The n-th root of unity, where n is the length of f.
/// @return The transformed polynomial in evaluation form.
pub fn ntt_radix_2(f: Vec<CF>, w: CF) -> Vec<CF> {
    let n = f.len();
    assert!(n >= 2, "n must be at least 2");
    assert!(n.is_power_of_two(), "n must be a power of 2");

    fn do_ntt(f: Vec<CF>, w: CF) -> Vec<CF> {
        let n = f.len();

        // Base case
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
            // Perform a radix-2 butterfly
            res[i] = ntt_even[i] + wk * ntt_odd[i];
            res[i + n/2] = ntt_even[i] - wk * ntt_odd[i];
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

/// Performs a radix-8 NTT on the polynomial f.
/// @param f The coefficients of the polynomial to be transformed.
/// @param w The n-th root of unity, where n is the length of f.
/// @param w8 The 8th root of unity.
/// @param w8_neg_1 The 8th root of unity multiplied by -1.
/// @return The transformed polynomial in evaluation form.
pub fn ntt_radix_8(f: Vec<CF>, w: CF) -> Vec<CF> {
    let n = f.len();
    debug_assert!(n >= 8, "n must be at least 8");
    debug_assert!(is_power_of_8(n as u32), "n must be a power of 8");

    fn do_ntt(f: Vec<CF>, w: CF) -> Vec<CF> {
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
            let i_8 = i * 8;
            a0[i] = f[i_8];
            a1[i] = f[i_8 + 1];
            a2[i] = f[i_8 + 2];
            a3[i] = f[i_8 + 3];
            a4[i] = f[i_8 + 4];
            a5[i] = f[i_8 + 5];
            a6[i] = f[i_8 + 6];
            a7[i] = f[i_8 + 7];
        }

        let w_pow_8 = w.pow(8);
        // Recurse
        let ntt_a0 = do_ntt(a0, w_pow_8);
        let ntt_a1 = do_ntt(a1, w_pow_8);
        let ntt_a2 = do_ntt(a2, w_pow_8);
        let ntt_a3 = do_ntt(a3, w_pow_8);
        let ntt_a4 = do_ntt(a4, w_pow_8);
        let ntt_a5 = do_ntt(a5, w_pow_8);
        let ntt_a6 = do_ntt(a6, w_pow_8);
        let ntt_a7 = do_ntt(a7, w_pow_8);

        let mut res = vec![CF::zero(); n];
        for k in 0..m {
            // TODO:
            // - Support NTTs that are not a power of 8. Do the mixed-radix NTT. Start with 128
            // - which is 2 x 64. Make sure that our mixed-radix NTTs are done in radix 8 stages,
            // - followed by radix 4, then radix 2
            // - Check if we do any unnecessary reductions

            let wt   = w.pow(k);
            let wt2  = wt  * wt;
            let wt3  = wt2 * wt;
            let wt4  = wt3 * wt;
            let wt5  = wt4 * wt;
            let wt6  = wt5 * wt;
            let wt7  = wt6 * wt;

            // Apply twiddle factors
            let t0 = ntt_a0[k];
            let t1 = wt  * ntt_a1[k];
            let t2 = wt2 * ntt_a2[k];
            let t3 = wt3 * ntt_a3[k];
            let t4 = wt4 * ntt_a4[k];
            let t5 = wt5 * ntt_a5[k];
            let t6 = wt6 * ntt_a6[k];
            let t7 = wt7 * ntt_a7[k];

            let ts = [t0, t1, t2, t3, t4, t5, t6, t7];

            // Use the provided 8-point butterfly.
            let butterfly = ntt_block_8(ts);

            // Write the butterfly result into the correct positions.
            for r in 0..8 {
                res[k + r * m] = butterfly[r];
            }
        }
        res
    }

    do_ntt(f, w)
}

/// Precomputes twiddle factors for a given size `n`.
/// @param n The size of the NTT (must be a power of 2).
/// @param w The nth root of unity
/// @param radix The butterfly size
/// @return A vector of precomputed twiddle factors.
pub fn precompute_twiddles(n: usize, w: CF, radix: usize) -> Vec<CF> {
    assert!(n.is_power_of_two(), "n must be a power of 2");
    assert!(n >= radix, "n must be at least as large as radix");
    assert!(radix.is_power_of_two(), "radix must be a power of 2");

    let mut twiddles = Vec::new();
    let mut current_n = n;
    let mut current_w = w;
    
    while current_n > 1 {
        let m = current_n / radix;
        let next_w = current_w.pow(radix);
        twiddles.push(next_w.reduce());
        
        for k in 0..m {
            let base = current_w.pow(k);
            let mut factor = CF::one();
            for _r in 1..radix {
                factor = factor * base;
                twiddles.push(factor.reduce());
            }
        }
        
        current_n /= radix;
        current_w = next_w;
    }
    
    twiddles
}

/// Precomputes twiddle factors needed for a stride-2 combination stage of an NTT.
/// @param n The size of the full NTT
/// @return Vector of w^i factors for i in 0..n/2
pub fn precompute_twiddles_stride2(n: usize) -> Vec<CF> {
    assert!(n.is_power_of_two(), "n must be a power of 2");
    assert!(n >= 2, "n must be at least 2");
    
    let w = get_root_of_unity(n);
    
    // Precompute w^i for i in 0..n/2
    let mut w_powers = Vec::with_capacity(n/2);
    let mut w_i = CF::one();
    
    for _ in 0..n/2 {
        w_powers.push(w_i);
        
        // Update for next iteration
        w_i = w_i * w;
    }
    
    w_powers
}

/// Precomputes twiddle factors needed for a stride-4 combination stage of an NTT.
/// @param n The size of the full NTT
/// @return Vector of [w^i, w^(2i), w^(3i)] arrays for i in 0..n/4
pub fn precompute_twiddles_stride4(n: usize) -> Vec<[CF; 3]> {
    assert!(n.is_power_of_two(), "n must be a power of 2");
    assert!(n >= 4, "n must be at least 4");
    assert!(n % 4 == 0, "n must be divisible by 4");
    
    let w = get_root_of_unity(n);
    
    // Precompute w^i, w^(2i), w^(3i) for i in 0..n/4
    let subn = n / 4;
    let mut w_powers = Vec::with_capacity(subn);
    let mut w_i = CF::one();
    
    for _ in 0..subn {
        let w_2i = w_i * w_i;       // w^(2i)
        let w_3i = w_2i * w_i;      // w^(3i)
        w_powers.push([w_i, w_2i, w_3i]);
        
        // Update for next iteration
        w_i = w_i * w;
    }
    
    w_powers
}

fn level_offset(overall_transform_size: usize, d: usize) -> usize {
    let mut offset = 0;
    let mut current = overall_transform_size;
    for _ in 0..d {
        offset += 1 + 7 * (current / 8);
        current /= 8;
    }
    offset
}

/// Performs a radix‑8 NTT using precomputed twiddle factors in `twiddles`.
/// The twiddle table must have been generated for the overall transform size.
/// This implementation uses a helper that, based on the recursion depth,
/// computes the correct offsets into the flat table.
pub fn ntt_radix_8_precomputed(
    f: Vec<CF>,
    twiddles: &Vec<CF>,
) -> Vec<CF> {
    let n = f.len();

    /// Recursive helper for the precomputed radix‑8 NTT.
    fn do_ntt_precomputed(
        f: Vec<CF>,
        twiddles: &Vec<CF>,
        depth: usize,
        overall_transform_size: usize,
    ) -> Vec<CF> {
        let n = f.len();
        if n == 1 { return f; }

        let m = n / 8;

        // Compute the starting offset for the current recursion level.
        let offset = level_offset(overall_transform_size, depth);

        // Block size for this level.
        let block_size = 1 + 7 * m;
        let level_twiddles = &twiddles[offset .. offset + block_size];

        // Partition the input into eight subarrays of length m.
        let mut a = vec![vec![CF::zero(); m]; 8];
        for j in 0..8 {
            for i in 0..m {
                a[j][i] = f[i * 8 + j];
            }
        }

        // Recurse
        let ntt_a: Vec<Vec<CF>> = a
            .into_iter()
            .map(|sub_f| do_ntt_precomputed(sub_f, twiddles, depth + 1, overall_transform_size))
            .collect();

        let mut res = vec![CF::zero(); n];
        for k in 0..m {
            let base_index = 1 + k * 7;
            let wt   = level_twiddles[base_index];
            let wt2  = level_twiddles[base_index + 1];
            let wt3  = level_twiddles[base_index + 2];
            let wt4  = level_twiddles[base_index + 3];
            let wt5  = level_twiddles[base_index + 4];
            let wt6  = level_twiddles[base_index + 5];
            let wt7  = level_twiddles[base_index + 6];

            let t0 = ntt_a[0][k];
            let t1 = wt  * ntt_a[1][k];
            let t2 = wt2 * ntt_a[2][k];
            let t3 = wt3 * ntt_a[3][k];
            let t4 = wt4 * ntt_a[4][k];
            let t5 = wt5 * ntt_a[5][k];
            let t6 = wt6 * ntt_a[6][k];
            let t7 = wt7 * ntt_a[7][k];
            let ts = [t0, t1, t2, t3, t4, t5, t6, t7];

            let butterfly = ntt_block_8(ts);
            for r in 0..8 {
                res[k + r * m] = butterfly[r];
            }
        }
        res
    }

    do_ntt_precomputed(f, twiddles, 0, n)
}

/// Performs an NTT on inputs whose length is of the form (8^k) * 2
/// Uses ntt_block_8 for the initial transforms, followed by a radix-2 combination stage
/// @param f The coefficients of the polynomial to be transformed.
/// @return The transformed polynomial in evaluation form.
pub fn ntt_8_stride_2(f: Vec<CF>) -> Vec<CF> {
    let n = f.len();
    // Input length must be at least 16 and of the form (8^k) * 2
    assert!(n >= 16, "n must be at least 16");
    assert!(n % 2 == 0, "n must be divisible by 2");
    
    // Check if n/2 is a power of 8
    let half_n = n / 2;
    let mut temp = half_n;
    while temp % 8 == 0 {
        temp /= 8;
    }
    assert_eq!(temp, 1, "n must be of the form (8^k) * 2");
    
    // Split input into even and odd parts
    let mut f_even = vec![CF::zero(); half_n];
    let mut f_odd = vec![CF::zero(); half_n];
    
    for i in 0..half_n {
        f_even[i] = f[2*i];
        f_odd[i] = f[2*i + 1];
    }
    
    // Get the primitive nth root of unity and its square
    let w = get_root_of_unity(n);
    let w_squared = w.pow(2); // This is the (n/2)th root of unity
    
    // Perform radix-8 NTT on each half
    let ntt_even = ntt_radix_8(f_even, w_squared);
    let ntt_odd = ntt_radix_8(f_odd, w_squared);
    
    // Combine using radix-2 butterfly operations
    let mut res = vec![CF::zero(); n];
    
    let mut w_k = CF::one();
    for i in 0..half_n {
        // Perform radix-2 butterfly operations
        res[i] = ntt_even[i] + w_k * ntt_odd[i];
        res[i + half_n] = ntt_even[i] - w_k * ntt_odd[i];
        
        // Update twiddle factor
        w_k = w_k * w;
    }
    
    res
}

/// Performs an NTT on inputs whose length is of the form (8^k) * 2 using precomputed twiddle factors
/// Uses ntt_radix_8_precomputed for the initial transforms, followed by a radix-2 combination stage
/// @param f The coefficients of the polynomial to be transformed.
/// @param precomputed_r8_twiddles The precomputed twiddle factors for radix-8 NTTs on half the input size.
/// @param precomputed_stride2_twiddles Precomputed twiddle factors for the stride-2 stage.
/// @return The transformed polynomial in evaluation form.
pub fn ntt_8_stride_2_precomputed(
    f: Vec<CF>,
    precomputed_r8_twiddles: &Vec<CF>,
    precomputed_stride2_twiddles: &Vec<CF>
) -> Vec<CF> {
    let n = f.len();
    // Input length must be at least 16 and of the form (8^k) * 2
    assert!(n >= 16, "n must be at least 16");
    assert!(n % 2 == 0, "n must be divisible by 2");
    
    // Check if n/2 is a power of 8
    let half_n = n / 2;
    let mut temp = half_n;
    while temp % 8 == 0 {
        temp /= 8;
    }
    assert_eq!(temp, 1, "n must be of the form (8^k) * 2");
    
    // Split input into even and odd parts
    let mut f_even = vec![CF::zero(); half_n];
    let mut f_odd = vec![CF::zero(); half_n];
    
    for i in 0..half_n {
        f_even[i] = f[2*i];
        f_odd[i] = f[2*i + 1];
    }
    
    // Perform radix-8 NTT on each half using precomputed twiddles
    let ntt_even = ntt_radix_8_precomputed(f_even, precomputed_r8_twiddles);
    let ntt_odd = ntt_radix_8_precomputed(f_odd, precomputed_r8_twiddles);
    
    // Combine using radix-2 butterfly operations with precomputed twiddles
    let mut res = vec![CF::zero(); n];
    
    for i in 0..half_n {
        // Get the precomputed twiddle factor
        let w_i = precomputed_stride2_twiddles[i];
        
        // Perform radix-2 butterfly operations
        res[i] = ntt_even[i] + w_i * ntt_odd[i];
        res[i + half_n] = ntt_even[i] - w_i * ntt_odd[i];
    }
    
    res
}

/// Performs an NTT on inputs whose length is of the form (8^k) * 4
/// Uses ntt_block_8 for the initial transforms, followed by a radix-4 combination stage
/// @param f The coefficients of the polynomial to be transformed.
/// @return The transformed polynomial in evaluation form.
pub fn ntt_8_stride_4(f: Vec<CF>) -> Vec<CF> {
    let n = f.len();
    // Input length must be at least 32 and of the form (8^k) * 4
    assert!(n >= 32, "n must be at least 32");
    assert!(n % 4 == 0, "n must be divisible by 4");
    
    // Check if n/4 is a power of 8
    let quarter_n = n / 4;
    let mut temp = quarter_n;
    while temp % 8 == 0 {
        temp /= 8;
    }
    assert_eq!(temp, 1, "n must be of the form (8^k) * 4");
    
    // Split input into 4 parts with stride 4
    let mut f0 = vec![CF::zero(); quarter_n];
    let mut f1 = vec![CF::zero(); quarter_n];
    let mut f2 = vec![CF::zero(); quarter_n];
    let mut f3 = vec![CF::zero(); quarter_n];
    
    for i in 0..quarter_n {
        f0[i] = f[4*i];
        f1[i] = f[4*i + 1];
        f2[i] = f[4*i + 2];
        f3[i] = f[4*i + 3];
    }
    
    // Get the primitive nth root of unity and its powers
    let w = get_root_of_unity(n);
    let w_4 = w.pow(4); // This is the (n/4)th root of unity
    
    // Perform radix-8 NTT on each quarter
    let ntt_f0 = ntt_radix_8(f0, w_4);
    let ntt_f1 = ntt_radix_8(f1, w_4);
    let ntt_f2 = ntt_radix_8(f2, w_4);
    let ntt_f3 = ntt_radix_8(f3, w_4);
    
    // Combine using radix-4 butterfly operations
    let mut res = vec![CF::zero(); n];
    
    for i in 0..quarter_n {
        // Calculate twiddle factors for this position
        let w_i = w.pow(i);                   // w^i
        let w_2i = w_i * w_i;                 // w^(2i)
        let w_3i = w_2i * w_i;                // w^(3i)
        
        // Apply twiddle factors to the NTT results
        let t0 = ntt_f0[i];
        let t1 = ntt_f1[i] * w_i;
        let t2 = ntt_f2[i] * w_2i;
        let t3 = ntt_f3[i] * w_3i;
        
        // Perform radix-4 butterfly operations
        let a0 = t0 + t2;
        let a1 = t0 - t2;
        let a2 = t1 + t3;
        let a3 = t1 - t3;
        
        // Using j = sqrt(-1)
        let a3_j = a3.mul_j();
        
        // Final combination using the radix-4 pattern
        res[i]                = a0 + a2;
        res[i + quarter_n]    = a1 + a3_j;
        res[i + 2 * quarter_n] = a0 - a2;
        res[i + 3 * quarter_n] = a1 - a3_j;
    }
    
    res
}

/// Performs an NTT on inputs whose length is of the form (8^k) * 4 using precomputed twiddle factors
/// Uses ntt_radix_8_precomputed for the initial transforms, followed by a radix-4 combination stage
/// @param f The coefficients of the polynomial to be transformed.
/// @param precomputed_r8_twiddles The precomputed twiddle factors for radix-8 NTTs on quarter the input size.
/// @param precomputed_stride4_twiddles Precomputed twiddle factors for the stride-4 stage.
/// @return The transformed polynomial in evaluation form.
pub fn ntt_8_stride_4_precomputed(
    f: Vec<CF>,
    precomputed_r8_twiddles: &Vec<CF>,
    precomputed_stride4_twiddles: &Vec<[CF; 3]>
) -> Vec<CF> {
    let n = f.len();
    // Input length must be at least 32 and of the form (8^k) * 4
    assert!(n >= 32, "n must be at least 32");
    assert!(n % 4 == 0, "n must be divisible by 4");
    
    // Check if n/4 is a power of 8
    let quarter_n = n / 4;
    let mut temp = quarter_n;
    while temp % 8 == 0 {
        temp /= 8;
    }
    assert_eq!(temp, 1, "n must be of the form (8^k) * 4");
    
    // Split input into 4 parts with stride 4
    let mut f0 = vec![CF::zero(); quarter_n];
    let mut f1 = vec![CF::zero(); quarter_n];
    let mut f2 = vec![CF::zero(); quarter_n];
    let mut f3 = vec![CF::zero(); quarter_n];
    
    for i in 0..quarter_n {
        f0[i] = f[4*i];
        f1[i] = f[4*i + 1];
        f2[i] = f[4*i + 2];
        f3[i] = f[4*i + 3];
    }
    
    // Perform radix-8 NTT on each quarter using precomputed twiddles
    let ntt_f0 = ntt_radix_8_precomputed(f0, precomputed_r8_twiddles);
    let ntt_f1 = ntt_radix_8_precomputed(f1, precomputed_r8_twiddles);
    let ntt_f2 = ntt_radix_8_precomputed(f2, precomputed_r8_twiddles);
    let ntt_f3 = ntt_radix_8_precomputed(f3, precomputed_r8_twiddles);
    
    // Combine using radix-4 butterfly operations with precomputed twiddles
    let mut res = vec![CF::zero(); n];
    
    for i in 0..quarter_n {
        // Get precomputed twiddle factors for this position
        let [w_i, w_2i, w_3i] = precomputed_stride4_twiddles[i];
        
        // Apply twiddle factors to the NTT results
        let t0 = ntt_f0[i];
        let t1 = ntt_f1[i] * w_i;
        let t2 = ntt_f2[i] * w_2i;
        let t3 = ntt_f3[i] * w_3i;
        
        // Perform radix-4 butterfly operations
        let a0 = t0 + t2;
        let a1 = t0 - t2;
        let a2 = t1 + t3;
        let a3 = t1 - t3;
        
        // Using j = sqrt(-1)
        let a3_j = a3.mul_j();
        
        // Final combination using the radix-4 pattern
        res[i]                = a0 + a2;
        res[i + quarter_n]    = a1 + a3_j;
        res[i + 2 * quarter_n] = a0 - a2;
        res[i + 3 * quarter_n] = a1 - a3_j;
    }
    
    res
}

#[cfg(test)]
pub mod tests {
    use crate::ntt::*;
    use crate::cm31::CF;
    use num_traits::Zero;
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

    #[test]
    fn test_naive_poly_mul_4() {
        let f1 = [CF::new(1, 0), CF::new(2, 0), CF::new(3, 0), CF::new(4, 0)];
        let f3 = [CF::new(1, 0), CF::new(3, 0), CF::new(5, 0), CF::new(7, 0)];
        let res = naive_poly_mul_4(f1, f3);
        let expected = [CF::new(1, 0), CF::new(5, 0), CF::new(14, 0), CF::new(30, 0), CF::new(41, 0), CF::new(41, 0), CF::new(28, 0)];
        assert_eq!(res, expected);
    }

    #[test]
    fn test_naive_ntt() {
        let mut rng = ChaCha8Rng::seed_from_u64(0);
        let sizes = [4, 8, 16, 32, 64, 128];
        for size in sizes.iter() {
            let mut f = vec![CF::zero(); *size];
            for i in 0..*size {
                f[i] = rng.r#gen();
            }

            let res = naive_ntt(f.clone());
            let ires = naive_intt(res.clone());

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

            let poly1_ntt = naive_ntt(poly1_padded.to_vec());
            let poly2_ntt = naive_ntt(poly2_padded.to_vec());
            let mut product_ntt = [CF::zero(); 4];

            for i in 0..4 {
                product_ntt[i] = poly1_ntt[i] * poly2_ntt[i];
            }
            let product_poly = naive_intt(product_ntt.to_vec());
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

            let poly1_ntt = naive_ntt(poly1_padded.to_vec());
            let poly2_ntt = naive_ntt(poly2_padded.to_vec());
            let mut product_ntt = [CF::zero(); 8];

            for i in 0..8 {
                product_ntt[i] = poly1_ntt[i] * poly2_ntt[i];
            }
            let product_poly = naive_intt(product_ntt.to_vec());
            let expected_product = naive_poly_mul_4(poly1, poly2);

            for i in 0..expected_product.len() {
                assert_eq!(product_poly[i], expected_product[i]);
            }
        }
    }

    #[test]
    fn test_ntt_block_4() {
        // Test the radix-4 NTT butterfly.
        let mut rng = ChaCha8Rng::seed_from_u64(0);
        for _ in 0..1024 {
            let mut poly = [CF::zero(); 4];
            for j in 0..4 {
                poly[j] = rng.r#gen();
            }

            let naive = naive_ntt(poly.to_vec());
            let res = ntt_block_4(poly);
            assert_eq!(naive, res);

            let ires = naive_intt(res.to_vec());
            assert_eq!(ires, poly);
        }
    }

    #[test]
    fn test_ntt_block_8() {
        // Test the radix-8 NTT butterfly.
        let mut rng = ChaCha8Rng::seed_from_u64(0);
        for _ in 0..1024 {
            let mut poly = [CF::zero(); 8];
            for j in 0..8 {
                poly[j] = rng.r#gen();
            }

            let naive = naive_ntt(poly.to_vec());
            let res = ntt_block_8(poly);
            assert_eq!(naive, res);

            let ires = naive_intt(res.to_vec());
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
        let mut rng = ChaCha8Rng::seed_from_u64(0);

        let n = 8 * 8 * 8;
        let w = get_root_of_unity(n);

        for _ in 0..4 {
            let mut f = vec![CF::zero(); n];
            for i in 0..n {
                f[i] = rng.r#gen();
            }
            let res = ntt_radix_8(f.clone(), w);
            let naive_res = naive_ntt(f.clone());
            assert_eq!(res, naive_res);
        }
    }

    #[test]
    fn test_ntt_radix_8_precomputed() {
        let n = 512;
        let radix = 8;
        let w = get_root_of_unity(n);

        let mut rng = ChaCha8Rng::seed_from_u64(0);

        let twiddles = precompute_twiddles(n, w, radix);

        for _ in 0..4 {
            let mut f = vec![CF::zero(); n];
            for i in 0..n {
                f[i] = rng.r#gen();
            }
            let res = ntt_radix_8_precomputed(f.clone(), &twiddles);

            let naive_res = naive_ntt(f);
            assert_eq!(res, naive_res);
        }
    }

    #[test]
    fn test_ntt_8_stride_2() {
        // Test with 1024 = 8^3 * 2 elements
        let n = 1024;
        let mut rng = ChaCha8Rng::seed_from_u64(0);
        
        // Generate random input
        let mut f = vec![CF::zero(); n];
        for i in 0..n {
            f[i] = rng.r#gen();
        }
        
        // Compute NTT using our implementation
        let res = ntt_8_stride_2(f.clone());
        
        // Compute NTT using naive approach for comparison
        let naive_res = naive_ntt(f.clone());
        
        // Verify correctness
        assert_eq!(res, naive_res);
        
        // Also test with 16 = 8^1 * 2 elements (smallest valid size)
        let n_small = 16;
        let mut f_small = vec![CF::zero(); n_small];
        for i in 0..n_small {
            f_small[i] = rng.r#gen();
        }
        
        let res_small = ntt_8_stride_2(f_small.clone());
        let naive_res_small = naive_ntt(f_small.clone());
        
        assert_eq!(res_small, naive_res_small);
    }
    
    #[test]
    fn test_ntt_8_stride_2_precomputed() {
        // Test for different sizes: 16, 128, 1024
        let sizes = [16, 128, 1024];
        let radix = 8;
        let mut rng = ChaCha8Rng::seed_from_u64(0);
        
        for &n in &sizes {
            let half_n = n / 2;
            
            // Precompute twiddle factors for half_n-point radix-8 NTTs
            let whalf = get_root_of_unity(half_n);
            let r8_twiddles = precompute_twiddles(half_n, whalf, radix);
            
            // Precompute twiddle factors for stride-2 combination stage
            let stride2_twiddles = precompute_twiddles_stride2(n);
            
            // Generate random input
            let mut f = vec![CF::zero(); n];
            for i in 0..n {
                f[i] = rng.r#gen();
            }
            
            // Compute NTT using the precomputed implementation
            let res_precomputed = ntt_8_stride_2_precomputed(f.clone(), &r8_twiddles, &stride2_twiddles);
            
            // Compare with regular non-precomputed implementation
            let res_regular = ntt_8_stride_2(f.clone());
            assert_eq!(res_precomputed, res_regular);
            
            // Compare with naive approach
            let naive_res = naive_ntt(f);
            assert_eq!(res_precomputed, naive_res);
        }
    }
    
    #[test]
    fn test_ntt_8_stride_4() {
        // Test with 32, 256, 2048 elements (32 = 8^1 * 4, 256 = 8^2 * 4, 2048 = 8^3 * 4)
        let sizes = [32, 256, 2048];
        let mut rng = ChaCha8Rng::seed_from_u64(0);
        
        for &n in &sizes {
            // Generate random input
            let mut f = vec![CF::zero(); n];
            for i in 0..n {
                f[i] = rng.r#gen();
            }
            
            // Compute NTT using our implementation
            let res = ntt_8_stride_4(f.clone());
            
            // Compute NTT using naive approach for comparison
            let naive_res = naive_ntt(f.clone());
            
            // Verify correctness
            assert_eq!(res, naive_res);
        }
    }
    
    #[test]
    fn test_ntt_8_stride_4_precomputed() {
        // Test with 32, 256, 2048 elements
        let sizes = [32, 256, 2048];
        let radix = 8;
        let mut rng = ChaCha8Rng::seed_from_u64(0);
        
        for &n in &sizes {
            let quarter_n = n / 4;
            
            // Precompute twiddle factors for quarter_n-point radix-8 NTTs
            let wquarter = get_root_of_unity(quarter_n);
            let r8_twiddles = precompute_twiddles(quarter_n, wquarter, radix);
            
            // Precompute twiddle factors for stride-4 combination stage
            let stride4_twiddles = precompute_twiddles_stride4(n);
            
            // Generate random input
            let mut f = vec![CF::zero(); n];
            for i in 0..n {
                f[i] = rng.r#gen();
            }
            
            // Compute NTT using the precomputed implementation
            let res_precomputed = ntt_8_stride_4_precomputed(f.clone(), &r8_twiddles, &stride4_twiddles);
            
            // Compare with regular non-precomputed implementation
            let res_regular = ntt_8_stride_4(f.clone());
            assert_eq!(res_precomputed, res_regular);
            
            // Compare with naive approach
            let naive_res = naive_ntt(f);
            assert_eq!(res_precomputed, naive_res);
        }
    }
}
