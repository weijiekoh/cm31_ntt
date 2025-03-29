/// Complex M31 field arithmetic.

use crate::rm31::{ RF, reduce, P, P_64 };
use std::ops::{ Add, AddAssign, Sub, SubAssign, Neg, Mul, MulAssign };
use core::fmt::Display;
use std::convert::{ From, Into };
use num_traits::Zero;
use num_traits::identities::One;
use num_traits::pow::Pow;
use rand::distributions::{Distribution, Standard};
use rand::Rng;

#[derive(Copy, Clone, Debug)]
pub struct CF {
    pub(crate) a: RF, // the real part
    pub(crate) b: RF, // the imaginary part
}

impl CF {
    pub const fn new(real: u32, imag: u32) -> CF {
        CF { a: RF::new(real), b: RF::new(imag) }
    }

    pub const fn real(self) -> RF {
        self.a
    }

    pub const fn imag(self) -> RF {
        self.b
    }

    pub fn mul_by_f(self, f: RF) -> CF {
        CF { 
            a: f * self.real(),
            b: f * self.imag(),
        }
    }

    /*
    pub fn add_base_f(self, f: RF) -> CF {
        // TODO: check if this is correct
        CF { 
            a: self.real() + f,
            b: self.imag(),
        }
    }
    */

    pub fn try_inverse(&self) -> Option<Self> {
        if self.a.val == 0 && self.b.val == 0 {
            return None;
        }

        // TODO: optimise by deferring reduction.
        let a2b2 = (self.a * self.a + self.b * self.b).reduce();
        if a2b2.is_zero() {
            return None;
        }

        let a2b2_inv = a2b2.try_inverse().unwrap().reduce();
        debug_assert!((a2b2 * a2b2_inv).reduce() == RF::new(1));

        let neg_b = self.b.neg();
        let a_neg_b = CF { a: self.a, b: neg_b };

        let result = a_neg_b.mul_by_f(a2b2_inv);
        Some(result)
    }

    /// Returns the 8th root of unity. Since there are 4 options for this value, select the one you
    /// want using the input `i`.
    /// The 8th root of unity should be (+-2^15, +-2^15)
    /// Let v = 2^15
    /// The options denoted by `i` are:
    /// 0. ( v,  v)
    /// 1. ( v, -v)
    /// 2. (-v,  v)
    /// 3. (-v, -v)
    pub fn root_of_unity_8(i: u32) -> Result<CF, String> {
        assert!(i < 4);
        let v = 2u32.pow(15);
        let neg_v = P - v;
        // i = 0: (v, v)
        // i = 1: (v, -v)
        // i = 2: (-v, v)
        // i = 3: (-v, -v)
        if i == 0 {
            return Ok(CF::new(v, v));
        }
        if i == 1 {
            return Ok(CF::new(v, neg_v));
        }
        if i == 2 {
            return Ok(CF::new(neg_v, v));
        }
        if i == 3 {
            return Ok(CF::new(neg_v, neg_v));
        }
        panic!("i must be 0, 1, 2 or 3");
    }

    #[inline]
    pub fn mul_neg_1(self) -> Self {
        // TODO: figure out if we can work with non-reduced values.
        debug_assert!(self.a.val < P_64);
        debug_assert!(self.b.val < P_64);
        let c = RF::new(P - (self.a.val as u32));
        let d = RF::new(P - (self.b.val as u32));
        CF { a: c, b: d }
    }

    #[inline]
    pub fn mul_j(self) -> Self {
        // TODO: figure out if we can work with non-reduced values.
        debug_assert!(self.a.val < P_64);
        debug_assert!(self.b.val < P_64);
        CF { a: RF::new(0) - self.b, b: self.a }
    }

    #[inline]
    pub fn mul_j_neg_1(self) -> Self {
        let a = self.a;
        let b = self.b;
        let c = RF::new(0);
        let d = RF::new(P - 1);

        let ac = (a * c).reduce();
        let bd = (b * d).reduce();
        let real = ac - bd;
        let imag = ((a + b).reduce() * (c + d).reduce() - ac).reduce() - bd;

        CF { a: real.reduce(), b: imag.reduce() }
    }
}

impl Zero for CF {
    #[inline]
    fn zero() -> CF {
        CF::new(0, 0)
    }
    #[inline]
    fn is_zero(&self) -> bool {
        self.a.is_zero() && self.b.is_zero()
    }
}

impl One for CF {
    #[inline]
    fn one() -> CF {
        CF::new(1, 0)
    }
}

impl Into<CF> for u32 {
    #[inline]
    /// Converts a u32 into a CF where the real part is the specified u32, and the imaginary part
    /// is 0
    fn into(self) -> CF {
        CF::new(self, 0)
    }
}

impl Into<CF> for (u32, u32) {
    #[inline]
    fn into(self) -> CF {
        CF { a: RF::new(self.0), b: RF::new(self.1) }
    }
}

impl From<CF> for (u32, u32) {
    #[inline]
    fn from(f: CF) -> (u32, u32) {
        (reduce(f.a.val) as u32, reduce(f.b.val) as u32)
    }
}

impl Add for CF {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        let a = self.a;
        let b = self.b;
        let c = rhs.a;
        let d = rhs.b;

        CF { a: (a + c).reduce(), b: (b + d).reduce() }
    }
}

impl AddAssign for CF {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl Sub for CF {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        let a = self.a;
        let b = self.b;
        let c = rhs.a;
        let d = rhs.b;

        CF { a: (a - c).reduce(), b: (b - d).reduce() }
    }
}

impl SubAssign for CF {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl Mul for CF {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: Self) -> Self::Output {
        // TODO: optimise by deferring reduction. (One must be very careful when doing so.)
        // To do so we need add_without_reduction and mul_without_reduction functions, and also
        // prove that the results remain within bounds, and do a final reduction on the final
        // result (real, imag).
        let a = self.a;
        let b = self.b;
        let c = rhs.a;
        let d = rhs.b;

        let ac = (a * c).reduce();
        let bd = (b * d).reduce();
        let real = ac - bd;
        let imag = ((a + b).reduce() * (c + d).reduce() - ac).reduce() - bd;

        CF { a: real.reduce(), b: imag.reduce() }
    }
}

impl MulAssign for CF {
    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl Neg for CF {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output {
        CF { a: -self.a, b: -self.b }
    }
}

impl PartialEq for CF {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.a == other.a && self.b == other.b
    }
}

impl Eq for CF {}

impl Display for CF {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} + {}i", self.a, self.b)
    }
}

impl Pow<usize> for CF {
    type Output = CF;

    #[inline]
    fn pow(self, exp: usize) -> Self::Output {
        let mut result = CF::one();
        let mut base = self;
        let mut exp = exp;
        while exp > 0 {
            if exp % 2 == 1 {
                result *= base;
            }
            base *= base;
            exp /= 2;
        }
        result
    }
}

impl Distribution<CF> for Standard {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> CF {
        CF {
            a: rng.r#gen(),
            b: rng.r#gen(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::m31::P;
    use rand_chacha::ChaCha8Rng;
    use rand_chacha::rand_core::SeedableRng;


    #[test]
    fn test_new() {
        assert_eq!(CF::new(0, 0).a.val, 0);
        assert_eq!(CF::new(1, 0).a.val, 1);
        assert_eq!(CF::new(0, 1).b.val, 1);
    }

    #[test]
    fn test_one() {
        assert_eq!(CF::one().a.val, 1);
        assert_eq!(CF::one().b.val, 0);
    }

    #[test]
    fn test_into() {
        let x: CF = 0u32.into();
        assert_eq!(x, CF::new(0, 0));
        
        let x: CF = 1u32.into();
        assert_eq!(x, CF::new(1, 0));

        let x: CF = (1u32, 2u32).into();
        assert_eq!(x, CF::new(1, 2));

        let y: (u32, u32) = x.into();
        assert_eq!(y, (1, 2));
    }

    #[test]
    fn test_add() {
        assert_eq!(CF::new(1, 0) + CF::new(1, 0), CF::new(2, 0));
        assert_eq!(CF::new(0, 1) + CF::new(0, 1), CF::new(0, 2));
        assert_eq!(CF::new(P - 1, P - 1) + CF::new(1, 2), CF::new(0, 1));
    }

    #[test]
    fn test_neg() {
        let x = CF::new(1, 2);
        assert_eq!(-x, CF::new(P - 1, P - 2));
    }

    #[test]
    fn test_mul() {
        assert_eq!(
            CF::new(2, 2) * CF::new(4, 5),
            CF::new(P - 2, 18)
        );
    }

    #[test]
    fn test_inverse() {
        let mut rng = ChaCha8Rng::seed_from_u64(0);
        for _ in 0..1024 {
            let x: CF = rng.r#gen();
            let x_inv = CF::try_inverse(&x).unwrap();
            assert_eq!(x * x_inv, CF::new(1, 0));
        }
    }

    #[test]
    fn test_pow() {
        let mut rng = ChaCha8Rng::seed_from_u64(0);
        for _ in 0..128 {
            let x: CF = rng.r#gen();
            let mut r = CF::one();
            for i in 0..1024 {
                assert_eq!(r, x.pow(i));
                r *= x;
            }
        }
    }

    #[test]
    fn test_w8() {
        fn do_test(w8: CF) {
            // w8 ^ 8 should equal 1, and not 1 for w8 ^ 1..7
            let r1 = w8;
            let r2 = w8 * w8;
            let r3 = r2 * w8;
            let r4 = r3 * w8;
            let r5 = r4 * w8;
            let r6 = r5 * w8;
            let r7 = r6 * w8;
            let r8 = r7 * w8;

            let one = CF::one();
            assert_ne!(r1, one);
            assert_ne!(r2, one);
            assert_ne!(r3, one);
            assert_ne!(r4, one);
            assert_ne!(r5, one);
            assert_ne!(r6, one);
            assert_ne!(r7, one);
            assert_eq!(r8, one);
        }
        for i in 0..4 {
            do_test(CF::root_of_unity_8(i).unwrap());
        }
    }
}
