/// Complex M31 field arithmetic.

use crate::rm31::{ RF, reduce };
use std::ops::{ Add, AddAssign, Sub, SubAssign, Neg, Mul, MulAssign };
use core::fmt::Display;
use std::convert::{ From, Into };
use num_traits::Zero;

#[derive(Copy, Clone, Debug)]
pub struct CF {
    pub(crate) a: RF, // the real part
    pub(crate) b: RF, // the imaginary part
}

impl CF {
    pub const fn new(a: u32, b: u32) -> CF {
        CF { a: RF::new(a), b: RF::new(b) }
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

    fn try_inverse(&self) -> Option<Self> {
        if self.a.val == 0 && self.b.val == 0 {
            return None;
        }

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

        CF { a: a + c, b: b + d }
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

        CF { a: a - c, b: b - d }
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
        let a = self.a;
        let b = self.b;
        let c = rhs.a;
        let d = rhs.b;

        let ac = a * c;
        let bd = b * d;
        let real = ac - bd;
        let imag = (a + b) * (c + d) - ac - bd;

        CF { a: real, b: imag }
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::m31::P;

    #[test]
    fn test_new() {
        assert_eq!(CF::new(0, 0).a.val, 0);
        assert_eq!(CF::new(1, 0).a.val, 1);
        assert_eq!(CF::new(0, 1).b.val, 1);
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
        let x = CF::new(3, 4);
        let x_inv = CF::try_inverse(&x).unwrap();

        assert_eq!(x * x_inv, CF::new(1, 0));
    }
}
