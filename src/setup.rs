//! Groth16 の Structured Reference String (SRS) を生成する trusted setup。
//!
//! Groth16 実装の Layer 3（プロトコル）。秘密のスカラー τ から
//! `[1, τ, τ², ...]·G1`、同 G2、`τⁱ·t(τ)·G1` の点列を計算する。
//!
//! ## 主要型
//! - [`Srs`]: simple 版の G1 / G2 / ht 点ベクトル
//! - [`ToxicWate`]: 本式 Groth16 の秘密値 α, β, γ, δ, τ
//! - [`QapFr`]: Fr 変換済みの QAP 多項式（setup/prover が τ 評価に使う）
//! - [`ProvingKey`] / [`VerifyingKey`]: 本式 Groth16 の 2 本立て鍵
//!
//! ## 主要関数
//! - [`generate_srs`]: τ と制約数 n から simple 版 SRS を生成
//! - [`generate_groth16_keys`]: QAP と toxic waste から pk/vk を生成
//!
//! ## 注意
//! simple 版（[`generate_srs`]）と本式（[`generate_groth16_keys`]）が併存している。
//! 本式 prover/verifier が E2E を通過した後、simple 版は削除予定。

use ark_bn254::{Fr, G1Projective, G2Projective};
use ark_ec::PrimeGroup;
use ark_ff::Field;

/// 構造化文字列（Structured Reference String）。
///
/// simple QAP 版の SRS。秘密のスカラー τ の冪を楕円曲線点に焼き込んだ
/// 「多項式を τ で評価するためのデータベース」。prover は多項式の係数と
/// この点列の内積で `f(τ)·G` を τ を知らずに計算できる。
pub struct Srs {
    /// `[G1, τ·G1, τ²·G1, ..., τ^(2n-2)·G1]`。長さ `2n-1`。
    /// `A(x)·B(x)` が次数 `2n-2` まで届くため `2n-1` 個の冪が要る。
    pub g1_points: Vec<G1Projective>,
    /// `[G2, τ·G2, ..., τ^(2n-2)·G2]`。長さ `2n-1`。B(x) を G2 上で評価する用。
    pub g2_points: Vec<G2Projective>,
    /// `[t(τ)·G1, τ·t(τ)·G1, ..., τ^(n-2)·t(τ)·G1]`。長さ `n-1`。
    /// h(x)（次数高々 `n-2`）に `t(τ)` を掛けた項を評価する用。
    pub ht_points: Vec<G1Projective>,
}

/// SRS を生成する trusted setup。
///
/// `tau`: 秘密のスカラー。本来は setup ceremony で生成し、SRS を作ったら破棄する
/// （現状はデモのため [`main`](crate) で固定値を渡している）。`tau` が漏れると
/// 偽の証明を作れてしまう。
///
/// `num_constraints`: R1CS の制約数 `n` ＝ QAP の補間点の数。SRS 長さは
/// `g1_points`/`g2_points` が `2n-1`、`ht_points` が `n-1` になる。
///
/// 計算量は O(n) のスカラー倍が支配的。
///
/// # Panics
///
/// `num_constraints == 0`（precondition `>= 1` 違反）のとき panic する。
/// R1CS から QAP を作る上で制約 0 件は意味を持たず、`2n-1` が usize の
/// underflow を起こすため、明示的に弾く。
pub fn generate_srs(tau: Fr, num_constraints: usize) -> Srs {
    assert!(num_constraints >= 1, "SRS requires at least one constraint");
    let g1 = G1Projective::generator(); // G1 を取得する
    let g2 = G2Projective::generator(); // G2 を取得する
    let n = num_constraints;

    // 1. tau のべき乗を事前に計算する: [1, tau, tau^2, ..., tau^(2n-2)]
    // R1CS制約数が n のとき、ラグランジュ補間で得られる各多項式 u_i(x), v_i(x), w_i(x) の次数は最大 n-1（n個の点を補間するから）。
    // A(x)=Sigma a_i u_i(x) と B(x)=Sigma a_i v_i(x) はどちらも次数が n-1 なので、A(x)*B(x) の次数は最大で (n-1) + (n-1) = 2(n-1) = 2n-2 。
    // 次数が 2n-2 の多項式は x^0 から x^(2n-2) までの 2n-2 + 1 = 2n-1 個の項を持つ。
    // SRSには各項に対する楕円曲線点が必要なので、2n-1 個のべき乗が必要。
    let num_powers = 2 * n - 1; // tau^0 から tau^(2n-2) まで、計 2n-1 個
    let mut tau_powers = Vec::with_capacity(num_powers);
    let mut current = Fr::from(1u64);
    for _ in 0..num_powers {
        tau_powers.push(current);
        current *= tau;
    }

    // 2. G1 用 SRS: [G1, tau G1, tau^2 G1, ...]
    let g1_points: Vec<G1Projective> = tau_powers.iter().map(|tp| g1 * tp).collect();

    // 3. G2 用 SRS: [G2, tau G2, tau^2 G2, ...]
    let g2_points: Vec<G2Projective> = tau_powers.iter().map(|tp| g2 * tp).collect();

    // 4. ターゲット多項式 t(tau) = (tau-1)(tau-2)...(tau-(n-1)) を計算
    // QAP の補間点 0, 1, ..., n-1 に対応する
    // セットアップは tau を知っているのでスカラーとして直接計算できる
    let mut t_tau = Fr::from(1u64);
    for i in 0..n {
        t_tau *= tau - Fr::from(i as u64);
    }

    // 5. h(tau)t(tau) 用 SRS: [t(tau)G1, tau t(tau)G1, ..., tau^(n-2)t(tau)G1]
    // h(x) の次数は最大 n-2 なので、n-1 個の点が必要
    let ht_points: Vec<G1Projective> = (0..n.saturating_sub(1))
        .map(|i| g1 * (tau_powers[i] * t_tau))
        .collect();

    Srs {
        g1_points,
        g2_points,
        ht_points,
    }
}

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
/// [`adapter::polys_to_fr_vecs`] の出力を 3 本まとめて持つ。
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
    use super::*;
    use ark_ff::Field; // tau.pow() のため

    #[test]
    fn test_srs_lengths() {
        let tau = Fr::from(123u64);
        let n = 3; // 制約数

        let srs = generate_srs(tau, n);

        // G1, G2: 2n-1 = 5個
        assert_eq!(srs.g1_points.len(), 5);
        assert_eq!(srs.g2_points.len(), 5);

        // ht: n-1 = 2個
        assert_eq!(srs.ht_points.len(), 2);
    }

    #[test]
    fn test_srs_g1_first_is_generator() {
        let tau = Fr::from(42u64);
        let srs = generate_srs(tau, 3);

        // 最初の点は tau^0 * G1 = G1 そのもの
        assert_eq!(srs.g1_points[0], G1Projective::generator());
    }

    #[test]
    fn test_srs_g1_second_is_tau_times_g() {
        let tau = Fr::from(42u64);
        let srs = generate_srs(tau, 3);

        // 2番目の点は tau^1 * G1
        let expected = G1Projective::generator() * tau;
        assert_eq!(srs.g1_points[1], expected);
    }

    #[test]
    fn test_srs_g1_all_powers() {
        // 全 i について g1_points[i] == tau^i * G1 を確認する
        let tau = Fr::from(7u64);
        let n = 3;
        let srs = generate_srs(tau, n);

        // 2n-1 = 5 点すべて
        assert_eq!(srs.g1_points.len(), 5);
        let g1 = G1Projective::generator();
        for (i, point) in srs.g1_points.iter().enumerate() {
            let expected = g1 * tau.pow([i as u64]);
            assert_eq!(*point, expected, "g1_points[{i}] が tau^{i}*G1 と不一致");
        }
    }

    #[test]
    fn test_srs_g2_power() {
        // G2 側も同じ規則 g2_points[i] == tau^i * G2 に従う（i=2 で代表確認）
        let tau = Fr::from(7u64);
        let srs = generate_srs(tau, 3);

        let expected = G2Projective::generator() * tau.pow([2u64]);
        assert_eq!(srs.g2_points[2], expected);
    }

    #[test]
    fn test_srs_ht_first_is_t_tau() {
        // n=3 のとき t(tau) = (tau-0)(tau-1)(tau-2)
        // ht_points[0] = tau^0 * t(tau) * G1 = t(tau) * G1
        let tau = Fr::from(7u64);
        let srs = generate_srs(tau, 3);

        let t_tau = (tau - Fr::from(0u64)) * (tau - Fr::from(1u64)) * (tau - Fr::from(2u64));
        let expected = G1Projective::generator() * t_tau;
        assert_eq!(srs.ht_points[0], expected);
    }

    #[test]
    fn test_srs_n_equals_one() {
        // n=1 の境界: 2n-1=1 点、ht は n-1=0 点で空
        let tau = Fr::from(7u64);
        let srs = generate_srs(tau, 1);

        assert_eq!(srs.g1_points.len(), 1);
        assert_eq!(srs.g2_points.len(), 1);
        assert!(srs.ht_points.is_empty());
    }

    #[test]
    #[should_panic(expected = "SRS requires at least one constraint")]
    fn test_srs_n_equals_zero_panics() {
        // n=0 は precondition 違反（コミット4 で追加した assert）
        let tau = Fr::from(7u64);
        let _ = generate_srs(tau, 0);
    }
}
