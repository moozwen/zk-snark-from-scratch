//! 自作 [`FieldElement`](crate::field::FieldElement) と arkworks の `Fr` の相互変換を提供する。
//!
//! Groth16 実装の Layer 2 と Layer 3 の境界。Layer 1〜2 は教育目的で自作した
//! 有限体・多項式ライブラリを使い、Layer 3 のペアリング演算は arkworks の
//! BN254 実装に依存する。本モジュールはその橋渡しを担う。
//!
//! ## 主要関数
//! - [`field_element_to_fr`]: 自作 → arkworks
//! - [`fr_to_field_element`]: arkworks → 自作
//! - [`polynomial_to_fr_vec`] / [`polys_to_fr_vecs`]: 多項式の係数ベクトルをまとめて変換

use ark_bn254::Fr;
use ark_ff::{BigInteger256, PrimeField};
use num_bigint::BigInt;
use num_bigint::Sign;

use crate::field::FieldElement;
use crate::polynomial::Polynomial;

/// FieldElement -> ark_bn254::Fr に変換する
pub fn field_element_to_fr(fe: &FieldElement) -> Fr {
    // 1. BigInt から リトルエンディアン のバイト列を取得
    // to_bytes_le() は (Sing, Vec<u8>) を返す
    // つまり bytes は [u8; 32] のような 8ビット（1バイト）ごとの配列
    let (_sign, bytes) = fe.value.to_bytes_le();

    // 2. バイト列を [u64; 4] に詰め替える
    // arkworks の BigInteger256 は 256ビット = 64ビット（8バイト） × 4 の配列
    // 0u64: 64ビットの符号なし整数型。すべて 0 で初期化された、u64 が4つ並んだ配列を作る
    let mut limbs = [0u64; 4]; // Limb: 慣習的に u64 の塊を指す
                               // BN254 という曲線のスカラー体 Fr の最大値（素数 r）は254ビットの長さをもつ
                               // 254ビットのデータを格納するには 256ビット（32バイト＝64ビット × 4）のメモリ容量があれば事足りる
                               // 33バイト以上のデータを Fr に押し込もうとすると、有限体の最大値を突き抜けて、数学的なオーバーフローを引き起こしてしまう
    for (i, chunk) in bytes.chunks(8).enumerate() {
        if i >= 4 {
            // 32バイト（64ビット）を超える数字は扱えない
            // 8バイトの塊が5周目＝33バイト目以降に入ったら、それ以降のデータは切り捨ててループを抜ける
            break;
        }
        // データが8バイトに満たない場合（端数が出た場合）、そのままでは u64 に変換できない
        // そこで、一時的な8バイトの空箱 buf を作り、手持ちのデータをコピーする
        // 足りない部分は ゼロ埋め される
        let mut buf = [0u8; 8];
        buf[..chunk.len()].copy_from_slice(chunk);
        limbs[i] = u64::from_le_bytes(buf);
    }

    // 3. BigInteger256 を作って Fr に変換する
    let big_int = BigInteger256::new(limbs);
    Fr::from_bigint(big_int).expect("Fr の範囲外の値です")
}

/// ark_bn254::Fr -> 自作 FieldElement に変換する
pub fn fr_to_field_element(fr: &Fr, p: &BigInt) -> FieldElement {
    // 1. Fr から BigInteger256 を取り出す
    let big_int: BigInteger256 = fr.into_bigint();

    // 2. [u64; 4] からバイト列に変換する
    let limbs = big_int.0;
    let mut bytes = Vec::with_capacity(32);
    for limb in &limbs {
        bytes.extend_from_slice(&limb.to_le_bytes());
    }

    // 3. バイト列から num_bigint::BigInt を作る
    let value = BigInt::from_bytes_le(Sign::Plus, &bytes);

    FieldElement::new(value, p.clone())
}

/// 自作 Polynomial の係数を Vec<Fr> に変換する
pub fn polynomial_to_fr_vec(poly: &Polynomial) -> Vec<Fr> {
    poly.coefficients
        .iter()
        .map(|coeff| field_element_to_fr(coeff))
        .collect()
}

/// QAP の多項式群（Vec<Polynomial>）をまとめて変換する
pub fn polys_to_fr_vecs(polys: &Vec<Polynomial>) -> Vec<Vec<Fr>> {
    polys.iter().map(|p| polynomial_to_fr_vec(p)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_small_value_roundtrip() {
        // BN254 の曲線位数（Fr の法）
        let p = BigInt::parse_bytes(
            b"21888242871839275222246405745257275088548364400416034343698204186575808495617",
            10,
        )
        .unwrap();

        // 小さい値でテスト
        let fe = FieldElement::new(BigInt::from(42), p.clone());
        let fr = field_element_to_fr(&fe);

        // arkworks 側で同じ値を使って比較
        let fr_expected = Fr::from(42u64);
        assert_eq!(fr, fr_expected);

        // 逆変換して元に戻ることを確認
        let fe_back = fr_to_field_element(&fr, &p);
        assert_eq!(fe_back.value, BigInt::from(42));
    }

    #[test]
    fn test_zero() {
        let p = BigInt::parse_bytes(
            b"21888242871839275222246405745257275088548364400416034343698204186575808495617",
            10,
        )
        .unwrap();

        let fe = FieldElement::new(BigInt::from(0), p.clone());
        let fr = field_element_to_fr(&fe);
        assert_eq!(fr, Fr::from(0u64));
    }

    #[test]
    fn test_large_value() {
        let p = BigInt::parse_bytes(
            b"21888242871839275222246405745257275088548364400416034343698204186575808495617",
            10,
        )
        .unwrap();

        // p-1 は Fr で表現可能な最大値
        let val = &p - &BigInt::from(1);
        let fe = FieldElement::new(val, p.clone());
        let fr = field_element_to_fr(&fe);

        // 逆変換して一致するか
        let fe_back = fr_to_field_element(&fr, &p);
        assert_eq!(fe.value, fe_back.value);
    }
}
