//! 証明をペアリング等式で検証する Verifier。
//!
//! Groth16 実装の Layer 3（プロトコル）。BN254 のペアリング `e: G1 × G2 → GT` を
//! 用いて、証明の正しさを定数時間で検証する。
//!
//! ## 主要関数
//! - [`verify_simple`]: simple QAP 版の検証 `e(A, B) = e(C, G2)`
//! - [`verify`]: 本式 Groth16 の検証 `e(A,B) = e(α,β)·e(vk_x,γ)·e(C,δ)`
//!
//! ペアリング `e: G1 × G2 → GT` の双線形性 `e(aP, bQ) = e(P, Q)^{ab}` を使う。
//! 指数の上で `A(τ)·B(τ)` と `C(τ)·1` を比較することで、QAP の関係
//! `A(τ)·B(τ) = W(τ) + h(τ)·t(τ)` を点のまま（τ を知らずに）確認できる。
//!
//! ## 注意
//! simple 版（[`verify_simple`]）と本式（[`verify`]）が併存している。
//! 本式が E2E で通った後、simple 版は Phase 6b で削除予定。

use ark_bn254::{Bn254, Fr, G2Projective};
use ark_ec::{pairing::Pairing, CurveGroup, PrimeGroup};

use crate::prover::{Groth16Proof, Proof};
use crate::setup::VerifyingKey;

/// simple 版の検証。`e([A]_1, [B]_2) == e([C]_1, G2)` をペアリングで確認する。
///
/// 双線形性より左辺は `e(G1, G2)^{A(τ)·B(τ)}`、右辺は `e(G1, G2)^{C(τ)}`。
/// 両者が一致するのは指数 `A(τ)·B(τ) = C(τ)`、すなわち QAP が満たされるとき。
/// `C(τ) = W(τ) + h(τ)·t(τ)` なので、これは `A·B - W = h·t` の点上での確認に当たる。
///
/// arkworks のペアリング API が `G1Affine` / `G2Affine` を要求するため、
/// proof の射影座標を [`CurveGroup::into_affine`] で変換してから渡す。
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
        "public_inputs length must equal l (vk.ic.len() - 1"
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
