//! 有限体係数の多項式と、その演算を提供する。
//!
//! Groth16 実装の Layer 1（数学的基盤）。R1CS から QAP への変換や、
//! 証明生成における h(x) = (A·B - C)/Z(x) の計算で使われる。
//!
//! ## 主要型
//! - [`Polynomial`]: [`FieldElement`] を係数とする dense 表現。
//!   `Add`, `Sub`, `Mul`, `Div` を実装。
//!
//! ## 主要メソッド
//! - [`Polynomial::evaluate`][]: ホーナー法で多項式を評価
//! - [`Polynomial::div_rem`][]: 多項式の長除法（商と余りを返す）
//! - [`Polynomial::lagrange_interpolation`][]: x = 0, 1, 2, ... の点列からラグランジュ補間

use crate::field::FieldElement;
use num_bigint::BigInt;
use std::ops::{Add, Div, Mul, Sub};

/// 有限体係数の多項式を dense 表現で保持する。
///
/// `coefficients[i]` が x^i の係数。例： `[1, 2, 3]` は `1 + 2x + 3x^2` を表す。
/// 末尾の 0 係数は [`Polynomial::new`] で自動的に取り除かれるため、
/// 意味的な次数と `coefficients.len() - 1` は常に一致する。
///
/// # 例
///
/// ```text
/// let p = Polynomial::new(vec![
///     FieldElement::new(1, 7),
///     FieldElement::new(2, 7),
/// ]); // 1 + 2x in F_7
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Polynomial {
    // coefficients[i] が x^i の係数
    // Dense（密）表現 を採用
    pub coefficients: Vec<FieldElement>,
}

impl Polynomial {
    /// 係数列から多項式を生成する。末尾の 0 係数は自動的に削除される。
    ///
    /// 全係数が 0 のときは `[0]` を残す（空 vec にはならない）。
    /// 空 vec を渡した場合は空のままとなる。
    ///
    /// # 例
    //
    /// ```text
    /// [1, 2, 0, 0] → [1, 2]  (1 + 2x)
    /// [0, 0, 0]    → [0]     (定数 0)
    ///```
    pub fn new(mut coefficients: Vec<FieldElement>) -> Self {
        while coefficients.len() > 1 && coefficients.last().unwrap().value == BigInt::from(0) {
            coefficients.pop();
        }
        Polynomial { coefficients }
    }

    /// 多項式の次数を返す。
    ///
    /// 定数 `c` の次数は 0、空多項式 (`coefficients.is_empty()`) の場合も 0 を返す。
    /// 両者を区別したい場合は `coefficients.is_empty()` で判定すること。
    pub fn degree(&self) -> usize {
        if self.coefficients.is_empty() {
            return 0;
        }
        self.coefficients.len() - 1
    }

    /// 多項式が 0 多項式かどうかを返す。
    ///
    /// [`Polynomial::new`] の正規化ルール（全 0 のとき `[0]` を残す）に依存。
    /// よって「係数 1 個 かつ それが 0」という単純判定を行う。
    pub fn is_zero(&self) -> bool {
        self.coefficients.len() == 1 && self.coefficients[0].value == BigInt::from(0)
    }

    /// 与えられた `x` で多項式を評価し、`P(x)` を返す。
    ///
    /// ホーナー法で実装しており、係数長 `n` に対して計算量は `O(n)`。
    ///
    /// # 例
    ///
    /// ```text
    /// // P(x) = 1 + 2x in F_7;  P(3) = 7 ≡ 0 (mod 7)
    /// let p = Polynomial::new(vec![
    ///     FieldElement::new(1, 7),
    ///     FieldElement::new(2, 7),
    /// ]);
    /// assert_eq!(p.evaluate(&FieldElement::new(3, 7)).value, BigInt::from(0));
    /// ```
    pub fn evaluate(&self, x: &FieldElement) -> FieldElement {
        let mut result = FieldElement::new(BigInt::from(0), x.p.clone());
        for coeff in self.coefficients.iter().rev() {
            result = &(&result * x) + coeff;
        }
        result
    }

    /// 多項式の長除法を行い、`(quotient, remainder)` を返す。
    ///
    /// 結果は不変式 `self == divisor * quotient + remainder` を満たし、
    /// `remainder` の次数は `divisor` より厳密に小さい。
    /// 被除数の次数が除数より小さいときは `(0, self)` を返す。
    ///
    /// # Panics
    ///
    /// `divisor` が 0 多項式の場合 panic する。
    pub fn div_rem(&self, divisor: &Polynomial) -> (Polynomial, Polynomial) {
        let p = self.coefficients[0].p.clone();

        // 0 で割ろうとした場合はパニック
        if divisor.is_zero() {
            panic!("0多項式で割ることはできません");
        }

        // 被除数の次数が除数より低い場合、商は 0、余りは被除数自身
        if self.degree() < divisor.degree() {
            return (
                Polynomial::new(vec![FieldElement::new(BigInt::from(0), p.clone())]),
                self.clone(),
            );
        }

        let mut quotient_coeffs = vec![
            FieldElement::new(BigInt::from(0), p.clone());
            self.degree() - divisor.degree() + 1
        ];
        let mut remainder = self.clone();

        // 長除法のメインループ
        while remainder.degree() >= divisor.degree() && !remainder.is_zero() {
            let deg_r = remainder.degree();
            let deg_d = divisor.degree();

            // a. 最高次の項同士の割り算
            let leading_r = remainder.coefficients.last().unwrap();
            let leading_d = divisor.coefficients.last().unwrap();
            let ratio = leading_r.div(leading_d); // FieldElement の割り算

            // 次数の差
            let deg_diff = deg_r - deg_d;
            quotient_coeffs[deg_diff] = ratio.clone();

            // b. 減算用の多項式（ratio * x^deg_diff * divisor）を作成
            let mut sub_coeffs = vec![
                FieldElement::new(BigInt::from(0), p.clone());
                deg_diff + divisor.coefficients.len()
            ];
            for (i, coeff) in divisor.coefficients.iter().enumerate() {
                sub_coeffs[i + deg_diff] = coeff * &ratio;
            }
            let sub_poly = Polynomial::new(sub_coeffs);

            // c. 余りから引く
            remainder = &remainder - &sub_poly;
        }

        (Polynomial::new(quotient_coeffs), remainder)
    }

    /// `y_values[i]` を `x = i` での値とする多項式を補間して返す。
    ///
    /// n 点から n-1 次以下の多項式が一意に定まる。計算量は `O(n^2)`。
    /// 補間点は `x = 0, 1, 2, ..., n-1` に固定（QAP 構築に最適化）。
    ///
    /// # 例
    ///
    /// ```text
    /// // y_i = (i+1)^2 mod 7  → [1, 4, 2]
    /// let y = vec![fe(1), fe(4), fe(2)];
    /// let p = Polynomial::lagrange_interpolation(&y);
    /// assert_eq!(p.evaluate(&fe(2)), fe(2));
    /// ```
    pub fn lagrange_interpolation(y_values: &[FieldElement]) -> Polynomial {
        if y_values.is_empty() {
            return Polynomial::new(vec![]);
        }

        // 素数 p を取得（計算に必要）
        let p = y_values[0].p.clone();

        // 合計用の多項式（最初は 0）
        let mut total_poly = Polynomial::new(vec![FieldElement::new(BigInt::from(0), p.clone())]);

        let num_points = y_values.len();

        // 各点 x_i = 0, 1, 2 ... についてループする
        for i in 0..num_points {
            let y_i = &y_values[i];

            // y_i が 0 なら計算しても結果は 0 なのでスキップ（高速化）
            // ただし厳密には基底計算が必要だが、結果に寄与しないのでOK
            if y_i.value == BigInt::from(0) {
                continue;
            }

            // 基底多項式 L_i(x) の作成
            // 分子（Numerator）： (x - x0)(xi - x1)...
            // 分母（Denominator）： (xi - x0)(xi - x1)...
            let mut numerator =
                Polynomial::new(vec![FieldElement::new(BigInt::from(1), p.clone())]);
            let mut denominator = FieldElement::new(BigInt::from(1), p.clone());

            let xi = FieldElement::new(BigInt::from(i), p.clone());

            for j in 0..num_points {
                // 自分自身はスキップ
                if i == j {
                    continue;
                }

                let xj = FieldElement::new(BigInt::from(j), p.clone());

                // 分子に (x - xj) をかける
                // (x - xj) という多項式は、係数が [-xj, 1]
                // つまり [xj * -1, 1]
                let zero = FieldElement::new(BigInt::from(0), p.clone());
                let neg_xj = &zero - &xj;
                let one = FieldElement::new(BigInt::from(1), p.clone());
                let term = Polynomial::new(vec![neg_xj, one]);
                numerator = &numerator * &term; // 多項式の掛け算

                // 分母に (xi - xj) をかける
                let diff = &xi - &xj;
                denominator = &denominator * &diff; // スカラーの掛け算
            }

            // 分母の逆数を計算して、分子にかける（割り算の代わり）
            let denom_inv = denominator
                .inverse()
                .expect("xi - xj is non-zero by construction (i != j");
            let basis_poly = numerator.scale(&denom_inv);

            // 高さをあわせて合計に足す： total += y_i * basis_poly
            let weighted_poly = basis_poly.scale(y_i);
            total_poly = &total_poly + &weighted_poly;
        }

        total_poly
    }

    /// 全係数に `factor` を掛けたスカラー倍多項式を返す。
    pub fn scale(&self, factor: &FieldElement) -> Polynomial {
        let new_coeffs = self.coefficients.iter().map(|c| c * factor).collect();
        Polynomial::new(new_coeffs)
    }
}

/// 多項式の加算: 同じ次数の係数同士を加算する。
impl<'a, 'b> Add<&'b Polynomial> for &'a Polynomial {
    type Output = Polynomial;

    fn add(self, other: &'b Polynomial) -> Polynomial {
        // 1. 両方とも空なら、空を返す
        if self.coefficients.is_empty() && other.coefficients.is_empty() {
            return Polynomial::new(vec![]);
        }

        // 2. p を安全に取得する
        // self が空なら other から取得する
        let p = if !self.coefficients.is_empty() {
            self.coefficients[0].p.clone()
        } else {
            other.coefficients[0].p.clone()
        };

        let max_len = std::cmp::max(self.coefficients.len(), other.coefficients.len());
        let mut res_coeffs = Vec::with_capacity(max_len);

        for i in 0..max_len {
            let zero = FieldElement::new(num_bigint::BigInt::from(0), p.clone());
            let a = self.coefficients.get(i).unwrap_or(&zero);
            let b = other.coefficients.get(i).unwrap_or(&zero);

            // ここで参照同士の足し算
            res_coeffs.push(a + b);
        }

        Polynomial::new(res_coeffs)
    }
}

/// 多項式の減算: 同じ次数の係数同士を減算する。
impl<'a, 'b> Sub<&'b Polynomial> for &'a Polynomial {
    type Output = Polynomial;

    fn sub(self, other: &'b Polynomial) -> Polynomial {
        let max_len = std::cmp::max(self.coefficients.len(), other.coefficients.len());
        let mut res_coeffs = Vec::with_capacity(max_len);
        let p = self.coefficients[0].p.clone();

        for i in 0..max_len {
            let zero = FieldElement::new(BigInt::from(0), p.clone());
            let a = self.coefficients.get(i).unwrap_or(&zero);
            let b = other.coefficients.get(i).unwrap_or(&zero);
            res_coeffs.push(a - b);
        }

        Polynomial::new(res_coeffs)
    }
}

/// 多項式の乗算: 各係数を畳み込んで `i + j` 次の項に集約する（計算量 `O(n*m)`）。
impl<'a, 'b> Mul<&'b Polynomial> for &'a Polynomial {
    type Output = Polynomial;

    fn mul(self, other: &'b Polynomial) -> Polynomial {
        let p = self.coefficients[0].p.clone();
        // どちらの多項式にも含まれている 0次のオフセットを、重複して数えないように調整
        let new_len = self.coefficients.len() + other.coefficients.len() - 1;
        let mut res_coeffs = vec![FieldElement::new(BigInt::from(0), p.clone()); new_len];

        for i in 0..self.coefficients.len() {
            for j in 0..other.coefficients.len() {
                let product = &self.coefficients[i] * &other.coefficients[j];
                res_coeffs[i + j] = &res_coeffs[i + j] + &product;
            }
        }

        Polynomial::new(res_coeffs)
    }
}

impl<'a, 'b> Div<&'b Polynomial> for &'a Polynomial {
    type Output = Polynomial;

    fn div(self, other: &'b Polynomial) -> Polynomial {
        let (q, _r) = self.div_rem(other);
        q
    }
}

impl std::fmt::Display for Polynomial {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.coefficients.is_empty() {
            return write!(f, "0");
        }

        let s = self
            .coefficients
            .iter()
            .enumerate()
            .rev()
            .filter(|(_, coeff)| coeff.value != BigInt::from(0) || self.degree() == 0)
            .map(|(i, coeff)| {
                if i == 0 {
                    format!("{}", coeff.value) // 定数項
                } else if i == 1 {
                    format!("{}x", coeff.value) // 1次の項
                } else {
                    format!("{}x^{}", coeff.value, i) // 2次以上の項
                }
            })
            .collect::<Vec<_>>()
            .join(" + ");

        write!(f, "{}", if s.is_empty() { "0".to_string() } else { s })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const P: i64 = 7;

    fn fe(v: i64) -> FieldElement {
        FieldElement::new(v, P)
    }

    fn poly(coeffs: &[i64]) -> Polynomial {
        Polynomial::new(coeffs.iter().map(|&c| fe(c)).collect())
    }

    #[test]
    fn new_strips_trailing_zeros() {
        // [1, 2, 0, 0] -> [1, 2]
        let p = poly(&[1, 2, 0, 0]);
        assert_eq!(p.coefficients, vec![fe(1), fe(2)]);
    }

    #[test]
    fn new_keeps_single_zero_for_zero_polynomial() {
        let p = poly(&[0, 0, 0]);
        assert_eq!(p.coefficients, vec![fe(0)]);
    }

    #[test]
    fn degree_basic_cases() {
        assert_eq!(poly(&[5]).degree(), 0);
        assert_eq!(poly(&[1, 2]).degree(), 1);
        assert_eq!(poly(&[1, 0, 3]).degree(), 2);
    }

    #[test]
    fn evaluate_uses_horner() {
        // P(x) = 1 + 2x; P(3) = 7 ≡ 0 (mod 7)
        let p = poly(&[1, 2]);
        assert_eq!(p.evaluate(&fe(3)), fe(0));
    }

    #[test]
    fn add_handles_different_lengths() {
        // (1 + 2x) + (3 + x^2) = 4 + 2x + x^2
        let a = poly(&[1, 2]);
        let b = poly(&[3, 0, 1]);
        assert_eq!((&a + &b).coefficients, vec![fe(4), fe(2), fe(1)]);
    }

    #[test]
    fn sub_normalizes_negative_results() {
        // (3 + x^2) - (1 + 2x) = 2 - 2x + x^2 ≡ 2 + 5x + x^2 (mod 7)
        let a = poly(&[3, 0, 1]);
        let b = poly(&[1, 2]);
        assert_eq!((&a - &b).coefficients, vec![fe(2), fe(5), fe(1)]);
    }

    #[test]
    fn mul_basic() {
        // (1 + x)(1 - x) = 1 - x^2 ≡ 1 + 6x^2 (mod 7)
        let a = poly(&[1, 1]);
        let b = poly(&[1, -1]);
        assert_eq!((&a * &b).coefficients, vec![fe(1), fe(0), fe(6)]);
    }

    #[test]
    fn div_rem_exact_division() {
        // (x^2 - 1) / (x - 1) = x + 1, remainder 0
        let dividend = poly(&[-1, 0, 1]);
        let divisor = poly(&[-1, 1]);
        let (q, r) = dividend.div_rem(&divisor);
        assert_eq!(q.coefficients, vec![fe(1), fe(1)]);
        assert_eq!(r.coefficients, vec![fe(0)]);
    }

    #[test]
    fn div_rem_with_remainder() {
        // (x^2 + 1) / x = x, remainder 1
        let dividend = poly(&[1, 0, 1]);
        let divisor = poly(&[0, 1]);
        let (q, r) = dividend.div_rem(&divisor);
        assert_eq!(q.coefficients, vec![fe(0), fe(1)]);
        assert_eq!(r.coefficients, vec![fe(1)]);
    }

    #[test]
    fn div_rem_dividend_smaller_than_divisor() {
        // (x + 1) / x^2 → q = 0, r = x + 1
        let dividend = poly(&[1, 1]);
        let divisor = poly(&[0, 0, 1]);
        let (q, r) = dividend.div_rem(&divisor);
        assert_eq!(q.coefficients, vec![fe(0)]);
        assert_eq!(r.coefficients, vec![fe(1), fe(1)]);
    }

    #[test]
    #[should_panic(expected = "0多項式")]
    fn div_rem_by_zero_polynomial_panics() {
        let dividend = poly(&[1, 1]);
        let divisor = poly(&[0]);
        let _ = dividend.div_rem(&divisor);
    }

    #[test]
    fn lagrange_interpolation_recovers_known_points() {
        // y_i = (i + 1)^2 mod 7 → [1, 4, 2]
        let y = vec![fe(1), fe(4), fe(2)];
        let p = Polynomial::lagrange_interpolation(&y);
        assert_eq!(p.evaluate(&fe(0)), fe(1));
        assert_eq!(p.evaluate(&fe(1)), fe(4));
        assert_eq!(p.evaluate(&fe(2)), fe(2));
    }

    #[test]
    fn lagrange_interpolation_single_point_is_constant() {
        let y = vec![fe(5)];
        let p = Polynomial::lagrange_interpolation(&y);
        assert_eq!(p.evaluate(&fe(0)), fe(5));
        assert_eq!(p.evaluate(&fe(99)), fe(5));
    }

    #[test]
    fn scale_multiplies_each_coefficient() {
        // (1 + 2x).scale(3) = 3 + 6x
        let p = poly(&[1, 2]);
        let scaled = p.scale(&fe(3));
        assert_eq!(scaled.coefficients, vec![fe(3), fe(6)]);
    }

    #[test]
    fn display_formats_polynomial() {
        // 1 + 0x + 2x^2 → "2x^2 + 1"
        let p = poly(&[1, 0, 2]);
        assert_eq!(format!("{}", p), "2x^2 + 1");
    }

    #[test]
    fn display_zero_polynomial() {
        assert_eq!(format!("{}", poly(&[0])), "0");
    }

    #[test]
    fn iz_zero_returns_true_for_zero_polynomial() {
        assert!(poly(&[0]).is_zero());
        assert!(poly(&[0, 0, 0]).is_zero()); // new() にて [0] に正規化される
    }
    
    #[test]
    fn is_zero_returns_false_for_nonzero_polynomial() {
        assert!(!poly(&[1]).is_zero());     // 定数 1
        assert!(!poly(&[0, 1]).is_zero());  // x
        assert!(!poly(&[1, 2]).is_zero());  // 1 + 2x
    }
}
