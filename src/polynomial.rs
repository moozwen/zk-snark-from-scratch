// 使い方
// =====
// let p = BigInt::from(17);

// // A(x) = 5x + 1
// let poly_a = Polynomial::new(vec![
//     FieldElement::new(BigInt::from(1), p.clone()),
//     FieldElement::new(BigInt::from(5), p.clone())
// ]);

// // B(x) x - 1; 有限体17 では x + 16
// let poly_b = Polynomial::new(vec![
//     FieldElement::new(BigInt::from(16), p.clone()),
//     FieldElement::new(BigInt::from(1), p.clone())
// ]);

// println!("A(x) = {}", poly_a);
// println!("B(x) = {}", poly_b);

// let poly_mul = &poly_a * &poly_b;
// println!("A(x) * B(x) = {}", poly_mul);
// =====

use crate::field::FieldElement;
use num_bigint::BigInt;
use std::ops::{Add, Div, Mul, RemAssign, Sub};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Polynomial {
    // coefficients[i] が x^i の係数
    // Dense（密）表現 を採用
    pub coefficients: Vec<FieldElement>,
}

impl Polynomial {
    pub fn new(mut coefficients: Vec<FieldElement>) -> Self {
        // 高次の係数が 0 の場合、その項を取り除く（例: 1 + 2x + 0x^2 -> 1 + 2x
        // ただし すべての係数が 0 の場合は [0] を返す
        while coefficients.len() > 1 && coefficients.last().unwrap().value == BigInt::from(0) {
            coefficients.pop();
        }
        Polynomial { coefficients }
    }

    pub fn degree(&self) -> usize {
        if self.coefficients.is_empty() {
            return 0;
        }
        self.coefficients.len() - 1
    }

    // P(x) を計算する
    pub fn evaluate(&self, x: &FieldElement) -> FieldElement {
        // ホーナー法: a_n*x^n + ... + a_0 = (...((a_n*x + a_{n-1}*x + a_{n-2})...))
        let mut result = FieldElement::new(BigInt::from(0), x.p.clone());
        for coeff in self.coefficients.iter().rev() {
            result = &(&result * x) + coeff;
        }
        result
    }

    // 多項式の割り算（self / divisor）
    // 戻り値: (商, 余り)
    pub fn div(&self, divisor: &Polynomial) -> (Polynomial, Polynomial) {
        let dividend = self.trim();
        let divisor = divisor.trim();

        if divisor.coefficients.is_empty() {
            panic!("0除算はできません");
        }

        if dividend.coefficients.is_empty() {
            return (Polynomial::new(vec![]), Polynomial::new(vec![]));
        }

        let p = divisor.coefficients[0].p.clone();
        let zero_fe = FieldElement::new(BigInt::from(0), p.clone());

        // 商 (quotient) と 余り (remainder)
        let mut quotient = Polynomial::new(vec![zero_fe.clone(); dividend.coefficients.len()]); // 十分なサイズで初期化
        let mut remainder = dividend.clone();

        // 筆算のループ： 余りの次数が割る数の次数より大きい間続ける
        while remainder.coefficients.len() >= divisor.coefficients.len() {
            if remainder.coefficients.is_empty() {
                break;
            }

            // 1. 最高次数の項同士を割って、係数を決める
            // （例： 3x^3 / x^1 = 3x^2）
            let rem_degree = remainder.degree();
            let div_degree = divisor.degree();
            let diff_degree = rem_degree - div_degree;

            let lead_rem = remainder.coefficients.last().unwrap();
            let lead_div = divisor.coefficients.last().unwrap();

            // 係数 = rem の頭 / div の頭 = rem の頭 * (div の頭の逆数)
            let factor = lead_rem * &lead_div.inverse();

            // 2. 引くための多項式を作る（factor * x^diff_degree）
            // 例： [0, 0, factor] みたいな多項式を作る
            let mut term_coeffs = vec![zero_fe.clone(); diff_degree];
            term_coeffs.push(factor.clone());
            let term_poly = Polynomial::new(term_coeffs);

            // 商に足す
            quotient = &quotient + &term_poly.clone();

            // 3. 余りから引く： remainder -= term * divisor
            let sub_poly = &term_poly * &divisor.clone();

            // 引き算
            // 簡易的に -1倍 して足す
            let minus_one = &FieldElement::new(BigInt::from(0), p.clone())
                - &FieldElement::new(BigInt::from(1), p.clone());
            let sub_poly_neg = sub_poly.scale(minus_one);
            remainder = &remainder + &sub_poly_neg;

            remainder = remainder.trim(); // 0になった最高次数の項を消す
        }

        (quotient.trim(), remainder.trim())
    }

    // 商と余りを返す（Quotient, Remainder）
    pub fn div_rem(&self, divisor: &Polynomial) -> (Polynomial, Polynomial) {
        let p = self.coefficients[0].p.clone();

        // 0 で割ろうとした場合はパニック
        if divisor.coefficients.len() == 1 && divisor.coefficients[0].value == BigInt::from(0) {
            panic!("0多項式で割ることはできません");
        }

        // 被除数の次数が除数より引く場合、商は 0、余りは被除数自身
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
        while remainder.degree() >= divisor.degree()
            && !(remainder.coefficients.len() == 1
                && remainder.coefficients[0].value == BigInt::from(0))
        {
            let deg_r = remainder.degree();
            let deg_d = divisor.degree();

            // a. 最高次の項同士の割り算（有限体なので逆元をかける）
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

    // ラグランジュ補間
    // y_values: x=0, x=1, x=2, ... に対応する y座標のリスト
    pub fn lagrange_interpolation(y_values: &Vec<FieldElement>) -> Polynomial {
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
                if i == j {
                    continue;
                } // 自分自身はスキップ

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
            let denom_inv = denominator.inverse();
            let basis_poly = numerator.scale(denom_inv);

            // 高さをあわせて合計に足す： total += y_i * basis_poly
            let weighted_poly = basis_poly.scale(y_i.clone());
            total_poly = &total_poly + &weighted_poly;
        }

        total_poly
    }

    // スカラー倍（係数を全部 k 倍する）
    pub fn scale(&self, factor: FieldElement) -> Polynomial {
        let new_coeffs = self.coefficients.iter().map(|c| c * &factor).collect();
        Polynomial::new(new_coeffs)
    }

    // 係数がゼロの項を末尾から削除してきれいにする（正規化）
    pub fn trim(&self) -> Polynomial {
        let mut coeffs = self.coefficients.clone();
        if coeffs.is_empty() {
            return self.clone();
        }

        let zero = BigInt::from(0);

        // 末尾から0を探して消す
        while coeffs.len() > 1 {
            if let Some(last) = coeffs.last() {
                if last.value == zero {
                    coeffs.pop();
                } else {
                    break;
                }
            } else {
                break;
            }
        }

        Polynomial::new(coeffs)
    }
}

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

impl<'a, 'b> Sub<&'b Polynomial> for &'a Polynomial {
    type Output = Polynomial;

    fn sub(self, other: &'b Polynomial) -> Polynomial {
        let max_len = std::cmp::max(self.coefficients.len(), other.coefficients.len());
        let mut res_coeffs = Vec::with_capacity(max_len);
        let p = self.coefficients[0].p.clone();

        for i in 0..max_len {
            let zero = FieldElement::new(BigInt::from(0), p.clone());
            let a = self.coefficients.get(i).unwrap_or(&zero);
            let b = self.coefficients.get(i).unwrap_or(&zero);
            res_coeffs.push(a - b);
        }

        Polynomial::new(res_coeffs)
    }
}

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
