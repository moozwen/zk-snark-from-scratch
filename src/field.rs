use num_bigint::BigInt;
use std::fmt;
use std::ops::{Add, Div, Mul, Sub};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FieldElement {
    pub value: BigInt, // 値
    pub p: BigInt,     // 法となる素数
}

impl FieldElement {
    pub fn new(value: BigInt, p: BigInt) -> Self {
        // 値が 0 <= value < p の範囲に収まるように正規化
        // Rust の % は余りを求める演算子であるため、負数を割ると結果がマイナスになる。
        // そこで、「負の数」を「正の正解」に無理やり引き戻す。
        let normalized_value = ((value % &p) + &p) % &p;
        FieldElement {
            value: normalized_value,
            p,
        }
    }

    // 逆元 a^-1 mod p を求める
    pub fn inverse(&self) -> Self {
        // 値が 0 の場合は逆元が存在しないのでパニック
        let inv_value = self.value.modinv(&self.p).expect("0の逆元は存在しません");
        FieldElement::new(inv_value, self.p.clone())
    }

    // 割り算 a / b は a * (b^-1) と同じ
    pub fn div(&self, other: &Self) -> Self {
        self * &other.inverse()
    }

    // べき乗（繰り返し二乗法 Square and Multiply）
    pub fn pow(&self, exponent: BigInt) -> Self {
        let mut res = FieldElement::new(BigInt::from(1), self.p.clone());
        let mut base = self.clone();
        let mut exp = exponent;

        let zero = BigInt::from(0);
        let two = BigInt::from(2);

        while exp > zero {
            // 指数の最下位ビットが1（奇数）なら、現在の base を結果にかける
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

    // モジュラ平方根を計算する関数
    // p % 4 == 3 の場合のみ対応 (Tonelli-Shanks法は未実装)
    pub fn sqrt(&self) -> Option<Self> {
        // 1. 定数の準備
        let three = BigInt::from(3);
        let four = BigInt::from(4);
        let one = BigInt::from(1);

        // 2. 素数の型チェック（p % 4 == 3 か？）
        if &self.p % &four != three {
            panic!("この素数は p % 4 == 3 の形式ではありません。Tonelli-Shanks法が必要です。");
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

impl<'a, 'b> Add<&'b FieldElement> for &'a FieldElement {
    type Output = FieldElement;

    fn add(self, other: &'b FieldElement) -> FieldElement {
        assert_eq!(self.p, other.p, "異なる標数の体では計算できません");
        FieldElement::new(&self.value + &other.value, self.p.clone())
    }
}

impl<'a, 'b> Sub<&'b FieldElement> for &'a FieldElement {
    type Output = FieldElement;

    fn sub(self, other: &'b FieldElement) -> FieldElement {
        assert_eq!(self.p, other.p, "異なる標数の体では計算できません");
        FieldElement::new(&self.value - &other.value, self.p.clone())
    }
}

impl<'a, 'b> Mul<&'b FieldElement> for &'a FieldElement {
    type Output = FieldElement;

    fn mul(self, other: &'b FieldElement) -> FieldElement {
        assert_eq!(self.p, other.p, "異なる標数の体では計算できません");
        FieldElement::new(&self.value * &other.value, self.p.clone())
    }
}

impl<'a, 'b> Div<&'b FieldElement> for &'a FieldElement {
    type Output = FieldElement;

    fn div(self, other: &'b FieldElement) -> FieldElement {
        assert_eq!(self.p, other.p, "異なる標数の体では計算できません");
        // 有限体の割り算は a * (bの逆元)
        self * &other.inverse()
    }
}

impl fmt::Display for FieldElement {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // "value mod p" という形式で表示するルールを定義
        write!(f, "{} mod {}", self.value, self.p)
    }
}
