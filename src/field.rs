//! 有限体 GF(p) 上の演算を提供する。
//!
//! Groth16 実装の Layer 1。多項式・楕円曲線・QAPなどの上位レイヤーがこの上に構築される。
//!
//! ## 主要型
//! - [`FieldElement`]: 法 `p` の元。`Add`, `Sub`, `Mul`, `Div` を実装。
//!
//! ## 制約
//! - [`FieldElement::sqrt`] は `p ≡ 3 (mod 4)` の素数でのみ計算する。
//!   それ以外は `None` を返す（Tonelli-Shanks 法は未実装）。

use num_bigint::BigInt;
use std::fmt;
use std::ops::{Add, Div, Mul, Sub};

/// 有限体 GF(p) 上の元を表す。
/// 
/// 内部的に `value` は `0 <= value < p` の範囲に正規化される。
/// `p` は素数を想定しているが、構造体側ではチェックしない（呼び出し側の責務）。
/// 
/// # 例
/// 
/// ```text
/// let a = FieldElement::new(3, 7);
/// let b = FieldElement::new(5, 7);
/// let sum = &a + &b; // 8 mod 7 = 1
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FieldElement {
    pub value: BigInt, // 値
    pub p: BigInt,     // 法となる素数
}

impl FieldElement {
    /// 法 `p` のもとで `value` を正規化した `FieldElement` を生成する。
    ///
    /// 負の値や `p` を超える値も自動で `0 <= v < p` に丸める。
    /// `value` と `p` は `Into<BigInt>` を満たす任意の型を受け付ける
    /// （`i32`, `i64`, `BigInt` など）。
    ///
    /// # 例
    ///
    /// ```text
    /// let a = FieldElement::new(-1, 7);  // value == 6
    /// let b = FieldElement::new(10, 7);  // value == 3
    /// ```
    pub fn new(value: impl Into<BigInt>, p: impl Into<BigInt>) -> Self {
        let value = value.into();
        let p = p.into();
        // 値が 0 <= value < p の範囲に収まるように正規化
        // Rust の % は余りを求める演算子であるため、負数を割ると結果がマイナスになる。
        // そこで、「負の数」を「正の整数」に無理やり引き戻す。
        let normalized_value = ((value % &p) + &p) % &p;
        FieldElement {
            value: normalized_value,
            p,
        }
    }

    /// 逆元 a^-1 mod p を求める。0 の場合は None を返す。
    /// 
    /// 内部的には `BigInt::modinv` を使い、拡張ユークリッド法で計算する。
    /// `p` が素数なら、フェルマーの小定理 `a^[p-2] ≡ a^{-1} (mod p)` でも
    /// 同じ結果が得られるが、拡張ユークリッド法の方が速い。
    pub fn inverse(&self) -> Option<Self> {
        let inv_value = self.value.modinv(&self.p)?;
        Some(FieldElement::new(inv_value, self.p.clone()))
    }

    /// 割り算 `a / b = a * b^{-1}`。`b == 0` のとき panic する。
    pub fn div(&self, other: &Self) -> Self {
        self * &other.inverse().expect("division by zero")
    }

    /// `self^exponent mod p` を計算する。
    ///
    /// 繰り返し二乗法（square-and-multiply）で `O(log exponent)` 時間。
    pub fn pow(&self, exponent: impl Into<BigInt>) -> Self {
        let mut res = FieldElement::new(BigInt::from(1), self.p.clone());
        let mut base = self.clone();
        let mut exp = exponent.into();

        let zero = BigInt::from(0);
        let two = BigInt::from(2);

        while exp > zero {
            // 指数の最下位ビットが1（奇数）なら、現在の base を結果に掛ける
            if &exp % &two != zero {
                res = &res * &base;
            }
            // base を二乗する（base = base^2）
            base = &base * &base;

            // 指数を半分にする（右シフト）
            exp = &exp / &two;
        }
        res
    }

    /// モジュラ平方根 `√self mod p` を返す。
    ///
    /// - 平方剰余でない場合: `None`
    /// - `p % 4 != 3` の場合: `None`（Tonelli-Shanks 法は未実装）
    /// - それ以外: `Some(root)` を返す。`root` と `p - root` の 2 つの根のうち
    ///   どちらか一方が返る（どちらかは保証しない）。
    ///
    /// # アルゴリズム
    ///
    /// `p ≡ 3 (mod 4)` のとき、`a^((p+1)/4) mod p` が候補となる。
    /// 検算 (`root^2 == self`) で平方剰余かどうかを判定する。
    pub fn sqrt(&self) -> Option<Self> {
        // 1. 定数の準備
        let three = BigInt::from(3);
        let four = BigInt::from(4);
        let one = BigInt::from(1);

        // 2. 素数の型チェック（p % 4 == 3 か？）
        // Tonelli-Shanks 法は未実装のため、上記以外の素数では None を返す
        if &self.p % &four != three {
            return None;
        }

        // 3. 指数の計算: exponent = (p + 1) / 4
        let exponent = (&self.p + &one) / &four;

        // 4. 候補の計算: root = self^exponent
        let root = self.pow(exponent);

        // 5. 検算: root * root = self に戻ることを確認する
        // 戻らない場合は平方剰余でない（ルートが存在しない）ことを表す
        if &root * &root == *self {
            Some(root)
        } else {
            None
        }
    }
}

/// `&a + &b`: 加法。法 `p` が異なる場合は panic する。
impl<'a, 'b> Add<&'b FieldElement> for &'a FieldElement {
    type Output = FieldElement;

    fn add(self, other: &'b FieldElement) -> FieldElement {
        assert_eq!(self.p, other.p, "異なる標数の体では計算できません");
        FieldElement::new(&self.value + &other.value, self.p.clone())
    }
}

/// `&a - &b`: 減法。法 `p` が異なる場合は panic する。
impl<'a, 'b> Sub<&'b FieldElement> for &'a FieldElement {
    type Output = FieldElement;

    fn sub(self, other: &'b FieldElement) -> FieldElement {
        assert_eq!(self.p, other.p, "異なる標数の体では計算できません");
        FieldElement::new(&self.value - &other.value, self.p.clone())
    }
}

/// `&a * &b`: 乗法。法 `p` が異なる場合は panic する。
impl<'a, 'b> Mul<&'b FieldElement> for &'a FieldElement {
    type Output = FieldElement;

    fn mul(self, other: &'b FieldElement) -> FieldElement {
        assert_eq!(self.p, other.p, "異なる標数の体では計算できません");
        FieldElement::new(&self.value * &other.value, self.p.clone())
    }
}

/// `&a / &b`: 除法 = `&a * &b.inverse()`。
/// 法 `p` が異なる場合 / `b == 0` の場合は panic する。
impl<'a, 'b> Div<&'b FieldElement> for &'a FieldElement {
    type Output = FieldElement;

    fn div(self, other: &'b FieldElement) -> FieldElement {
        assert_eq!(self.p, other.p, "異なる標数の体では計算できません");
        // 有限体の割り算は a * (bの逆元)
        self * &other.inverse().expect("division by zero")
    }
}

/// `"value mod p"` 形式で表示する。
impl fmt::Display for FieldElement {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // "value mod p" という形式で表示するルールを定義
        write!(f, "{} mod {}", self.value, self.p)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fe(value: i64, p: i64) -> FieldElement {
        FieldElement::new(value, p)
    }

    #[test]
    fn new_normalizes_negative_value() {
        // -1 mod 7 == 6
        assert_eq!(fe(-1, 7).value, BigInt::from(6));
    }

    #[test]
    fn new_normalizes_value_larger_than_p() {
        // 10 mod 7 == 3
        assert_eq!(fe(10, 7).value, BigInt::from(3));
    }

    #[test]
    fn add_sub_mul_basic() {
        // F_7 上で
        assert_eq!((&fe(3, 7) + &fe(5, 7)).value, BigInt::from(1));
        assert_eq!((&fe(2, 7) - &fe(5, 7)).value, BigInt::from(4));
        assert_eq!((&fe(3, 7) * &fe(5, 7)).value, BigInt::from(1));
    }

    #[test]
    fn div_basic() {
        // F_7: 6 / 3 == 2
        assert_eq!((&fe(6, 7) / &fe(3, 7)).value, BigInt::from(2));
    }

    #[test]
    #[should_panic(expected = "異なる標数")]
    fn add_with_different_modulus_panics() {
        let _ = &fe(1, 7) + &fe(1, 11);
    }

    #[test]
    fn inverse_of_zero_returns_none() {
        assert!(fe(0, 7).inverse().is_none());
    }

    #[test]
    fn inverse_times_self_is_one() {
        // F_7 で 3 の逆元は 5（3*5 = 15 = 1 mod 7）
        let a = fe(3, 7);
        let inv = a.inverse().unwrap();
        assert_eq!((&a * &inv).value, BigInt::from(1));
    }

    #[test]
    fn pow_basic_cases() {
        let a = fe(3, 7);
        assert_eq!(a.pow(BigInt::from(0)).value, BigInt::from(1)); // a^0 == 1
        assert_eq!(a.pow(BigInt::from(1)).value, BigInt::from(3)); // a^1 == a
                                                                   // フェルマーの小定理: a^(p-1) == 1
        assert_eq!(a.pow(BigInt::from(6)).value, BigInt::from(1));
    }

    #[test]
    fn sqrt_quadratic_residue() {
        // F_7 で 4 の平方根は 2 または 5（5 = -2 mod 7）
        let root = fe(4, 7).sqrt().unwrap();
        assert!(root.value == BigInt::from(2) || root.value == BigInt::from(5));
    }

    #[test]
    fn sqrt_non_residue_returns_none() {
        // F_7 で 3 は平方剰余ではない
        assert!(fe(3, 7).sqrt().is_none());
    }

    #[test]
    fn sqrt_unsupported_prime_returns_none() {
        // p = 5 は 5 % 4 == 1 なので未対応
        assert!(fe(4, 5).sqrt().is_none());
    }

    #[test]
    fn display_format() {
        assert_eq!(format!("{}", fe(3, 7)), "3 mod 7");
    }
}
