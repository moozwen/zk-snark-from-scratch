//! 証明をペアリング等式で検証する Verifier。
//!
//! Groth16 実装の Layer 3（プロトコル）。BN254 のペアリング `e: G1 × G2 → GT` を
//! 用いて、証明の正しさを定数時間で検証する。
//!
//! ## 主要関数
//! - [`verify`]: Groth16 の検証 `e(A,B) = e(α,β)·e(vk_x,γ)·e(C,δ)`
//!
//! 双線形性 `e(aP, bQ) = e(P, Q)^{ab}` により、4 つのペアリングの等式で
//! QAP の充足を τ を知らずに確認する。

use ark_bn254::{Bn254, Fr};
use ark_ec::{pairing::Pairing, CurveGroup};

use crate::prover::Groth16Proof;
use crate::setup::VerifyingKey;

/// 本式 Groth16 の検証。`e(A,B) == e(α,β)·e(vk_x,γ)·e(C,δ)` をペアリングで確認する。
///
/// `vk`: trusted setup で生成した [`VerifyingKey`]
/// `public_inputs`: 公開入力 `a_1..a_ℓ`（CS_ONE の `a_0=1` は含めない。長さ ℓ）
/// `proof`: prover が生成した [`Groth16Proof`]
///
/// まず公開入力を vk の `IC` で線形結合して `vk_x = IC_0 + Σ a_i·IC_i` を作る
/// （`IC_0` は `a_0=1` の項）。これは検証者だけが持つ公開項を集約した G1 点。
/// 次に4つのペアリングで等式を確認する: 左辺 `e(A,B)` 1 本、右辺は
/// `e(α,β)`（型の整合）・`e(vk_x,γ)`（公開項）・`e(C,δ)`（秘密項 + h·t）の積。
/// γ/δ で割って焼き込んだ public/private 項が、ここで γ/δ とのペアリングにより
/// 「元の値」に戻り、A·B との一致が QAP の充足を意味する。
///
/// arkworks のペアリング API は Affine 座標を要求するため射影座標を変換する。
/// また `PairingOutput` は加法群表現なので、右辺の積は `+` で合成する。
///
/// # Panics
/// `public_inputs.len() != vk.ic.len() - 1`（公開入力数と vk の不整合）のとき panic。
pub fn verify(vk: &VerifyingKey, public_inputs: &[Fr], proof: &Groth16Proof) -> bool {
    assert_eq!(
        public_inputs.len(),
        vk.ic.len() - 1,
        "public_inputs length must equal ℓ (vk.ic.len() - 1)"
    );

    // vk_x = IC_0 + Σ_{i=1..ℓ} a_i·IC_i
    let mut vk_x = vk.ic[0];
    for (i, input) in public_inputs.iter().enumerate() {
        vk_x += vk.ic[i + 1] * input;
    }

    // ペアリングには Affine 座標が必要
    let a_aff = proof.a.into_affine();
    let b_aff = proof.b.into_affine();
    let c_aff = proof.c.into_affine();
    let alpha_aff = vk.alpha_g1.into_affine();
    let beta_aff = vk.beta_g2.into_affine();
    let gamma_aff = vk.gamma_g2.into_affine();
    let delta_aff = vk.delta_g2.into_affine();
    let vk_x_aff = vk_x.into_affine();

    // e(A,B) == e(α,β) · e(vk_x,γ) · e(C,δ)（PairingOutput は加法群なので積は +）
    let lhs = Bn254::pairing(a_aff, b_aff);
    let rhs = Bn254::pairing(alpha_aff, beta_aff)
        + Bn254::pairing(vk_x_aff, gamma_aff)
        + Bn254::pairing(c_aff, delta_aff);

    lhs == rhs
}
