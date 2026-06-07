//! 自作 [`FieldElement`] と arkworks の `Fr` の相互変換を提供する。
//!
//! Groth16 実装の Layer 2 と Layer 3 の境界。Layer 1〜2 は教育目的で自作した
//! 有限体・多項式ライブラリを使い、Layer 3 のペアリング演算は arkworks の
//! BN254 実装に依存する。本モジュールはその橋渡しを担う。
//!
//! ## 主要関数
//! - [`field_element_to_fr`]: 自作 → arkworks
//! - [`polynomial_to_fr_vec`] / [`polys_to_fr_vecs`]: 多項式の係数ベクトルをまとめて変換

use ark_bn254::Fr;
use ark_ff::{BigInteger256, PrimeField};

use crate::field::FieldElement;
use crate::polynomial::Polynomial;

/// 自作 [`FieldElement`] を arkworks の `Fr` に変換する。
///
/// 変換は次の段を踏む:
/// `BigInt` → リトルエンディアンのバイト列 → `[u64; 4]`（limb）→
/// `BigInteger256` → `Fr`。`Fr` の法（BN254 のスカラー体位数）を超える値は
/// `from_bigint` が `None` を返すため panic する。境界は片方向で、proof は
/// G1/G2 点・ペアリング結果は `GT` なので `Fr` への逆変換は不要。
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
    Fr::from_bigint(big_int).expect("FieldElement value exceeds Fr modulus")
}

/// 自作 [`Polynomial`] の係数を `Vec<Fr>` に変換する。
///
/// QAP の多項式を SRS 上で評価（prover の内積計算）するための前処理。
/// 各係数を [`field_element_to_fr`] で変換し、`coefficients[i]` の並び
/// （`x^i` の係数）をそのまま保つ。
pub fn polynomial_to_fr_vec(poly: &Polynomial) -> Vec<Fr> {
    poly.coefficients.iter().map(field_element_to_fr).collect()
}

/// QAP の多項式群（`a_polys` / `b_polys` / `c_polys` など）をまとめて変換する。
///
/// 各 [`Polynomial`] に [`polynomial_to_fr_vec`] を適用するだけ。prover が
/// witness で重み付けする前の、変数ごとの係数ベクトル列を得るのに使う。
pub fn polys_to_fr_vecs(polys: &[Polynomial]) -> Vec<Vec<Fr>> {
    polys.iter().map(polynomial_to_fr_vec).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigInt;

    /// BN254 の曲線位数（Fr の法）を返すテスト用ヘルパ。
    fn bn254_modulus() -> BigInt {
        BigInt::parse_bytes(
            b"21888242871839275222246405745257275088548364400416034343698204186575808495617",
            10,
        )
        .unwrap()
    }

    #[test]
    fn test_small_value() {
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

        // p-1 は Fr 上の ‐1 に対応する（変換が panic せず正しいことを確認）
        assert_eq!(fr, -Fr::from(1u64));
    }

    #[test]
    fn test_polynomial_to_fr_vec() {
        let p = bn254_modulus();

        // 1  2x  3x^2 の係数ベクトル
        let poly = Polynomial::new(vec![
            FieldElement::new(1, p.clone()),
            FieldElement::new(2, p.clone()),
            FieldElement::new(3, p.clone()),
        ]);

        let frs = polynomial_to_fr_vec(&poly);

        assert_eq!(frs, vec![Fr::from(1u64), Fr::from(2u64), Fr::from(3u64)]);
    }

    #[test]
    fn test_polynomial_to_fr_vec_zero() {
        let p = bn254_modulus();

        // ゼロ多項式（Polynomial::new は単一ゼロ係数を保持する）
        let poly = Polynomial::new(vec![FieldElement::new(0, p.clone())]);

        let frs = polynomial_to_fr_vec(&poly);

        assert_eq!(frs, vec![Fr::from(0u64)]);
    }

    #[test]
    fn test_polys_to_fr_vecs_matches_individual() {
        let p = bn254_modulus();

        // 2 本の多項式: 1  2x と 5  7x  11x^2
        let poly0 = Polynomial::new(vec![
            FieldElement::new(1, p.clone()),
            FieldElement::new(2, p.clone()),
        ]);
        let poly1 = Polynomial::new(vec![
            FieldElement::new(5, p.clone()),
            FieldElement::new(7, p.clone()),
            FieldElement::new(11, p.clone()),
        ]);

        let polys = vec![poly0.clone(), poly1.clone()];
        let batch = polys_to_fr_vecs(&polys);

        // 一括変換が個別変換の結果と一致すること
        assert_eq!(batch.len(), 2);
        assert_eq!(batch[0], polynomial_to_fr_vec(&poly0));
        assert_eq!(batch[1], polynomial_to_fr_vec(&poly1));
    }
}
