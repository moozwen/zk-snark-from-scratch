//! 証明をペアリング等式で検証する Verifier。
//!
//! Groth16 実装の Layer 3（プロトコル）。BN254 のペアリング `e: G1 × G2 → GT` を
//! 用いて、証明の正しさを定数時間で検証する。
//!
//! ## 主要関数
//! - [`verify_simple`]: simple QAP 版の検証 `e(A, B) = e(C, G2)`
//!
//! ## 注意
//! 現状は **simple 版** の検証等式。v0.6 で Groth16 本式
//! `e(A, B) = e(α, β) · e(Σ aᵢ, γ) · e(C, δ)` に置き換え予定。

use ark_bn254::{Bn254, G2Projective};
use ark_ec::{pairing::Pairing, CurveGroup, PrimeGroup};

use crate::prover::Proof;

/// シンプル版の検証
/// `[A]_1 * [B]_2 == [C]_1 * G_2` をペアリングで確認する
pub fn verify_simple(proof: &Proof) -> bool {
    let g2 = G2Projective::generator();

    // ペアリングには Affine 座標が必要
    let a_affine = proof.a_g1.into_affine();
    let b_affine = proof.b_g2.into_affine();
    let c_affine = proof.c_g1.into_affine();
    let g2_affine = g2.into_affine();

    // 左辺: e([A]_1, [B]_2)
    let lhs = Bn254::pairing(a_affine, b_affine);

    // 右辺: e([C]_1, G_2)
    let rhs = Bn254::pairing(c_affine, g2_affine);

    lhs == rhs
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bn254::{Fr, G1Projective};

    /// スカラー a, b, c から proof (a·G1, b·G2, c·G1) を組むヘルパ。
    ///
    /// 検証等式 `e(A, B) = e(C, G2)` は双線形性より `e(G1, G2)^{a·b} = e(G1, G2)^c`、
    /// つまり `a·b == c (mod r)` のとき成立する。QAP 全パイプラインを通さずに
    /// 検証ロジック単体を確認できる。
    fn make_proof(a: u64, b: u64, c: u64) -> Proof {
        let g1 = G1Projective::generator();
        let g2 = G2Projective::generator();
        Proof {
            a_g1: g1 * Fr::from(a),
            b_g2: g2 * Fr::from(b),
            c_g1: g1 * Fr::from(c),
        }
    }

    #[test]
    fn test_verify_accepts_valid_proof() {
        // a·b = 3·5 = 15 = c → accept
        let proof = make_proof(3, 5, 15);
        assert!(verify_simple(&proof));
    }

    #[test]
    fn test_verify_rejects_shifted_a() {
        // A を G1 generator 分ずらす（実質 a: 3 → 4）。4·5 = 20 ≠ 15 → reject
        let mut proof = make_proof(3, 5, 15);
        proof.a_g1 += G1Projective::generator();
        assert!(!verify_simple(&proof));
    }

    #[test]
    fn test_verify_rejects_doubled_c() {
        // C を 2 倍にする（実質 c: 15 → 30）。3·5 = 15 ≠ 30 → reject
        let mut proof = make_proof(3, 5, 15);
        proof.c_g1 = proof.c_g1 + proof.c_g1;
        assert!(!verify_simple(&proof));
    }
}
