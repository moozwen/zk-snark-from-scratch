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
/// [A]_1 * [B]_2 == [C]_1 * G_2 をペアリングで確認する
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
