//! Groth16 の Structured Reference String (SRS) を生成する trusted setup。
//!
//! Groth16 実装の Layer 3（プロトコル）。秘密のスカラー τ から
//! `[1, τ, τ², ...]·G1`、同 G2、`τⁱ·t(τ)·G1` の点列を計算する。
//!
//! ## 主要型
//! - [`Srs`]: G1 / G2 / ht の楕円曲線点ベクトル
//!
//! ## 主要関数
//! - [`generate_srs`]: τ と制約数 n から SRS を生成
//!
//! ## 注意
//! 現状は **simple QAP 版** の SRS（α, β, γ, δ なし）。
//! v0.5 で Groth16 本式の SRS に置き換え予定。

use ark_bn254::{Fr, G1Projective, G2Projective};
use ark_ec::{CurveGroup, PrimeGroup};
use ark_ff::Field;

/// 構造化文字列（Structured Reference String）
pub struct Srs {
    pub g1_points: Vec<G1Projective>, // [G1, tau G1, tau^2 G1, ...]
    pub g2_points: Vec<G2Projective>, // [G2, tau G2, tau^2 G2, ...]
    pub ht_points: Vec<G1Projective>, // [t(tau)G1, tau t(tau)G1, ...]
}

/// SRS を生成する
///
/// tau: 秘密のスカラー。セレモニー後に破棄する
/// num_constraints: R1CS の制約数＝QAPの補間点の数
pub fn generate_srs(tau: Fr, num_constraints: usize) -> Srs {
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

    // 4. ターゲット多項式 t(tau) = (tau-1)(tau-2)...(tau-n) を計算
    // セットアップは tau を知っているのでスカラーとして直接計算できる
    let mut t_tau = Fr::from(1u64);
    // 修正前: for i in 1..=n
    // 修正後: 自作 QAP の補間点 0, 1, 2, ..., n-1 に合わせる
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

#[cfg(test)]
mod tests {
    use super::*;

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
}
