//! Groth16 の proving key / verifying key を生成する trusted setup。
//!
//! Groth16 実装の Layer 3（プロトコル）。秘密値 α, β, γ, δ, τ（toxic waste）から
//! proving key / verifying key の点群を計算し、生成後は toxic waste を破棄する。
//!
//! ## 主要型
//! - [`ToxicWaste`]: Groth16 の秘密値 α, β, γ, δ, τ
//! - [`QapFr`]: Fr 変換済みの QAP 多項式（setup/prover が τ 評価に使う）
//! - [`ProvingKey`] / [`VerifyingKey`]: proving key / verifying key の 2 本立て
//!
//! ## 主要関数
//! - [`generate_groth16_keys`]: QAP と toxic waste から pk/vk を生成

use ark_bn254::{Fr, G1Projective, G2Projective};
use ark_ec::PrimeGroup;
use ark_ff::Field;

/// 本式 Groth16 の toxic waste (trusted setupt の秘密値)。
///
/// setupt でこれらから pk/vk を構成し、生成後は破棄する前提。
/// いずれか 1 つでも漏れると、偽の証明を作れてしまう。
/// `gamma` / `delta` は逆元を取るため 0 であってはならない。
/// （本番では MPC ceremony で生成するが、本実装では単一生成 & 破棄前提とする。
pub struct ToxicWaste {
    pub alpha: Fr,
    pub beta: Fr,
    pub gamma: Fr,
    pub delta: Fr,
    pub tau: Fr,
}

/// QAP 多項式を `Fr` 係数に変換した薄いラッパ。
///
/// `a_polys[i]` / `b_polys[i]` / `c_polys[i]` が変数 `i` の `u_i / v_i / w_i` の
/// 昇順系数列（`x^j` の係数が index `j`）。
/// [`crate::adapter::polys_to_fr_vecs`] の出力を 3 本まとめて持つ。
/// setup と 本式 prover の双方が τ の評価入力に使う。
pub struct QapFr {
    pub a_polys: Vec<Vec<Fr>>,
    pub b_polys: Vec<Vec<Fr>>,
    pub c_polys: Vec<Vec<Fr>>,
}

/// 本式 Groth16 の proving key。prover が証明 `(A, B, C)` を作るのに必要な点群。
pub struct ProvingKey {
    /// `[α]_1`
    pub alpha_g1: G1Projective,
    /// `[β]_1`（C 計算で `r·B_1` の B_1 に β を混ぜる用）
    pub beta_g1: G1Projective,
    /// `[δ]_1`（C の `−r·s·δ` 項用）
    pub delta_g1: G1Projective,
    /// `[β]_2`
    pub beta_g2: G2Projective,
    /// `[δ]_2`
    pub delta_g2: G2Projective,
    /// `{[τ^i]_1}_{i=0..n-1}`。`Σ a_i·u_i(τ)` を G1 上で評価する用。長さ `n`。
    pub tau_g1: Vec<G1Projective>,
    /// `{[τ^i]_2}_{i=0..n-1}`。`Σ a_i·v_i(τ)` を G2 上で評価する用。長さ `n`。
    pub tau_g2: Vec<G2Projective>,
    /// private wire 項 `{ [ (β·u_i(τ) + α·v_i(τ) + w_i(τ)) / δ ]_1}_{i=ℓ+1..m-1}`。
    /// 長さ `m − num_public`。
    pub private_query: Vec<G1Projective>,
    /// h 項 `{ [ τ^i·t(τ) / δ ]_1 }_{i=0..n-2}`。長さ `n−1`。
    pub h_query: Vec<G1Projective>,
}

/// 本式 Groth16 の verifying key。verifier がペアリング等式を確認するのに必要な点群。
pub struct VerifyingKey {
    /// `[α]_1`
    pub alpha_g1: G1Projective,
    /// `[β]_2`
    pub beta_g2: G2Projective,
    /// `[γ]_2`
    pub gamma_g2: G2Projective,
    /// `[δ]_2`
    pub delta_g2: G2Projective,
    /// public wire 項 `IC = { [ (β·u_i(τ) + α·v_i(τ) + w_i(τ)) / γ ]_1}_{i=0..ℓ}`。
    /// 長さ `num_public`（= ℓ+1）。`IC[0]` は CS_ONE（`a_0 = 1`）の項。
    pub ic: Vec<G1Projective>,
}

/// 昇順係数の多項式を τ で評価する（Horner 法）。
///
/// `coeffs[j]` が `x^j` の係数。空ベクトルは 0 を返す。
fn eval_poly(coeffs: &[Fr], tau: Fr) -> Fr {
    let mut acc = Fr::from(0u64);
    for c in coeffs.iter().rev() {
        acc = acc * tau + c;
    }
    acc
}

/// 本式 Groth16 の proving key / verifying key を生成する trusted setup。
///
/// `qap_fr`: Fr 変換済みの QAP（変数ごとの `u_i / v_i / w_i` 係数）
/// `num_constraints`: 制約数 `n`（= QAP の補間点数）
/// `num_public`: public 変数の数 ℓ+1（CS_ONE を含む。R1CS の `num_public_variables`）
/// `toxic`: α, β, γ, δ, τ
///
/// 各 wire の `(β·u_i + α·v_i + w_i)(τ)` を計算し、public は `/γ`（→ vk の `IC`）、
/// private は `/δ`（→ pk の `private_query`）でスケールして G1 点に焼き込む。
/// この γ/δ 分離が public 入力項と h·t 項の混ぜ替えによる偽造を防ぐ（Groth16 の肝）。
///
/// # Panics
/// `num_constraints == 0`、`γ == 0`、`δ == 0`（いずれも逆元 / QAP が成立しない）
/// のとき panic する。
pub fn generate_groth16_keys(
    qap_fr: &QapFr,
    num_constraints: usize,
    num_public: usize,
    toxic: &ToxicWaste,
) -> (ProvingKey, VerifyingKey) {
    assert!(
        num_constraints >= 1,
        "Groth16 setup requires at least one constraint"
    );
    assert!(toxic.gamma != Fr::from(0u64), "gamma must be nonzero");
    assert!(toxic.delta != Fr::from(0u64), "delta must be nonzero");
    let gamma_inv = toxic.gamma.inverse().unwrap();
    let delta_inv = toxic.delta.inverse().unwrap();

    let g1 = G1Projective::generator();
    let g2 = G2Projective::generator();
    let n = num_constraints;
    let m = qap_fr.a_polys.len(); // 変数の数
    let tau = toxic.tau;

    // 1. [tau^i] の点列（i = 0..n-1）。
    // 冪のスカラーは使い回す
    let mut tau_powers = Vec::with_capacity(n);
    let mut current = Fr::from(1u64);
    for _ in 0..n {
        tau_powers.push(current);
        current *= tau;
    }
    let tau_g1: Vec<G1Projective> = tau_powers.iter().map(|tp| g1 * tp).collect();
    let tau_g2: Vec<G2Projective> = tau_powers.iter().map(|tp| g2 * tp).collect();

    // 2. t(tau) = Π_{j=0}^{n-1} (τ - j)。補間点 0..n-1 に対応するターゲット多項式。
    let mut t_tau = Fr::from(1u64);
    for j in 0..n {
        t_tau *= tau - Fr::from(j as u64);
    }

    // 3. 各 wire の (β·u_i + α·v_i + w_i)(τ) を public/private に振り分け。
    let mut ic = Vec::with_capacity(num_public);
    let mut private_query = Vec::with_capacity(m.saturating_sub(num_public));
    for i in 0..m {
        let u = eval_poly(&qap_fr.a_polys[i], tau);
        let v = eval_poly(&qap_fr.b_polys[i], tau);
        let w = eval_poly(&qap_fr.c_polys[i], tau);
        let combo = toxic.beta * u + toxic.alpha * v + w;
        if i < num_public {
            ic.push(g1 * (combo * gamma_inv)); // public: /gamma
        } else {
            private_query.push(g1 * (combo * delta_inv)); // private: /delta
        }
    }

    // 4. 項 [ τ^i·t(τ)/δ ]_1（i = 0..n-2）。h(x) の次数は高々 n-2。
    let h_query: Vec<G1Projective> = (0..n.saturating_sub(1))
        .map(|i| g1 * (tau_powers[i] * t_tau * delta_inv))
        .collect();

    let pk = ProvingKey {
        alpha_g1: g1 * toxic.alpha,
        beta_g1: g1 * toxic.beta,
        delta_g1: g1 * toxic.delta,
        beta_g2: g2 * toxic.beta,
        delta_g2: g2 * toxic.delta,
        tau_g1,
        tau_g2,
        private_query,
        h_query,
    };
    let vk = VerifyingKey {
        alpha_g1: g1 * toxic.alpha,
        beta_g2: g2 * toxic.beta,
        gamma_g2: g2 * toxic.gamma,
        delta_g2: g2 * toxic.delta,
        ic,
    };
    (pk, vk)
}

#[cfg(test)]
mod tests {
    use std::vec;

    use super::*;
    use ark_ff::Field; // tau.pow() のため

    fn sample_toxic() -> ToxicWaste {
        ToxicWaste {
            alpha: Fr::from(2u64),
            beta: Fr::from(3u64),
            gamma: Fr::from(5u64),
            delta: Fr::from(7u64),
            tau: Fr::from(11u64),
        }
    }

    /// m=3 変数・n=2 制約を想定した手組みQAP（構造確認用）
    /// wire0 = CS_ONE, wire1 = public 入力, wire2 = private
    fn sample_qap_fr() -> QapFr {
        QapFr {
            a_polys: vec![
                vec![Fr::from(1u64)],                 // wire0: u_0 = 1 -> u_0(tau) = 1
                vec![Fr::from(0u64), Fr::from(1u64)], // wire1: u_1 = 0 + x -> u_1(tau) = tau
                vec![Fr::from(2u64)],                 // wire2: u_2 = 2 -> u_2(tau) = 2
            ],
            b_polys: vec![
                vec![Fr::from(0u64)], // wire0: v_0 = 0
                vec![Fr::from(1u64)], // wire1: v_1 = 1
                vec![Fr::from(3u64)], // wire2: v_2 = 3
            ],
            c_polys: vec![
                vec![Fr::from(0u64)],                 // wire0: w_0 = 0
                vec![Fr::from(4u64)],                 // wire1: w_1 = 4
                vec![Fr::from(1u64), Fr::from(1u64)], // wire2: w_2 = 1 + x -> w_2(tau) = 1 + tau
            ],
        }
    }

    #[test]
    fn groth16_keys_have_expected_lengths() {
        let toxic = sample_toxic();
        let qap = sample_qap_fr();
        let n = 2;
        let num_public = 2; // CS_ONE + public 入力 1 つ（l = 1）
        let (pk, vk) = generate_groth16_keys(&qap, n, num_public, &toxic);

        let m = qap.a_polys.len(); // 3
        assert_eq!(vk.ic.len(), num_public); // l+1 = 2
        assert_eq!(pk.private_query.len(), m - num_public); // m-l-1 = 1

        assert_eq!(pk.h_query.len(), n - 1); // 1
        assert_eq!(pk.tau_g1.len(), n); // 2
        assert_eq!(pk.tau_g2.len(), n); // 2
    }

    #[test]
    fn groth16_keys_match_toxic_scalars() {
        let toxic = sample_toxic();
        let qap = sample_qap_fr();
        let (pk, vk) = generate_groth16_keys(&qap, 2, 2, &toxic);

        let g1 = G1Projective::generator();
        let g2 = G2Projective::generator();

        // toxic スカラーと点の対応を抜き取り確認
        assert_eq!(pk.alpha_g1, g1 * toxic.alpha);
        assert_eq!(vk.alpha_g1, g1 * toxic.alpha);
        assert_eq!(pk.beta_g2, g2 * toxic.beta);
        assert_eq!(vk.beta_g2, g2 * toxic.beta);
        assert_eq!(vk.gamma_g2, g2 * toxic.gamma);
        assert_eq!(vk.delta_g2, g2 * toxic.delta);
        // τ 冪列: index 1 は τ·G
        assert_eq!(pk.tau_g1[1], g1 * toxic.tau);
        assert_eq!(pk.tau_g2[1], g2 * toxic.tau);
    }

    #[test]
    fn ic_zero_matches_manual_combo() {
        let toxic = sample_toxic();
        let qap = sample_qap_fr();
        let (_pk, vk) = generate_groth16_keys(&qap, 2, 2, &toxic);

        // wire0: u_0(τ)=1, v_0(τ)=0, w_0(τ)=0
        // combo_0 = β·1 + α·0 + 0 = β、IC_0 = (combo_0/γ)·G1
        let combo0 = toxic.beta * Fr::from(1u64) + toxic.alpha * Fr::from(0u64) + Fr::from(0u64);
        let gamma_inv = toxic.gamma.inverse().unwrap();
        let expected = G1Projective::generator() * (combo0 * gamma_inv);
        assert_eq!(vk.ic[0], expected);
    }

    #[test]
    fn private_query_matches_manual_combo() {
        let toxic = sample_toxic();
        let qap = sample_qap_fr();
        let (pk, _vk) = generate_groth16_keys(&qap, 2, 2, &toxic);

        // wire2 (private): u=2, v=3, w=1+τ=12 → combo = β·2 + α·3 + 12 = 24
        // private_query[0] = (combo/δ)·G1
        let combo2 = toxic.beta * Fr::from(2u64)
            + toxic.alpha * Fr::from(3u64)
            + (Fr::from(1u64) + toxic.tau);
        let delta_inv = toxic.delta.inverse().unwrap();
        let expected = G1Projective::generator() * (combo2 * delta_inv);
        assert_eq!(pk.private_query[0], expected);
    }

    #[test]
    #[should_panic(expected = "delta must be nonzero")]
    fn groth16_keys_delta_zero_panics() {
        let mut toxic = sample_toxic();
        toxic.delta = Fr::from(0u64);
        let qap = sample_qap_fr();
        let _ = generate_groth16_keys(&qap, 2, 2, &toxic);
    }

    #[test]
    #[should_panic(expected = "gamma must be nonzero")]
    fn groth16_keys_gamma_zero_panics() {
        let mut toxic = sample_toxic();
        toxic.gamma = Fr::from(0u64);
        let qap = sample_qap_fr();
        let _ = generate_groth16_keys(&qap, 2, 2, &toxic);
    }
}
