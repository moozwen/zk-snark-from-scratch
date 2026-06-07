//! QAP + SRS から証明 `(A, B, C)` を生成する Prover。
//!
//! Groth16 実装の Layer 3（プロトコル）。Witness と QAP 多項式の係数を
//! SRS 上で評価し、楕円曲線点として証明を構築する。
//!
//! ## 主要型
//! - [`Proof`]: `(A ∈ G1, B ∈ G2, C ∈ G1)` の三点組
//!
//! ## 主要関数
//! - [`prove_simple`]: simple QAP 版の証明生成
//!
//! ## 注意
//! 現状は **simple 版** （α, β, γ, δ なし、zero-knowledge 化なし）。
//! v0.5 で Groth16 本式の Prover（ランダム化 r, s 含む）に置き換え予定。

use ark_bn254::{Fr, G1Projective, G2Projective};

use crate::setup::Srs;

/// simple 版の証明（楕円曲線上の3点）。
///
/// 各点は QAP 多項式を秘密の τ で評価した値を曲線上に焼き込んだもの:
/// `a_g1 = A(τ)·G1`、`b_g2 = B(τ)·G2`、`c_g1 = (W(τ) + h(τ)·t(τ))·G1`。
/// verifier はこの3点だけでペアリング等式を確認する。
#[derive(Debug)]
pub struct Proof {
    /// `[A]_1 = A(τ)·G1`。witness で重み付けした `u_i` の合成多項式の評価。
    pub a_g1: G1Projective,
    /// `[B]_2 = B(τ)·G2`。同じく `v_i` の合成多項式の評価（G2 上）。
    pub b_g2: G2Projective,
    /// `[C]_1 = (W(τ) + h(τ)·t(τ))·G1`。`w_i` の合成に h·t の項を足したもの。
    pub c_g1: G1Projective,
}

/// 多項式の係数ベクトルを SRS 上で評価し、`f(τ)·G` を点として得る。
///
/// `f(τ)·G = coeffs[0]·G + coeffs[1]·(τ·G) + coeffs[2]·(τ²·G) + ...`
/// すなわち係数ベクトルと SRS 点列の内積。SRS には τ の冪が点として
/// 入っているので、τ 自体を知らなくても `f(τ)·G` を計算できる。
///
/// `evaluate_on_g2` は G2 版で中身は同形。教育コードとして G1/G2 を
/// 並べて読めるよう、あえてジェネリック化せずコピーで残している。
fn evaluate_on_g1(coeffs: &[Fr], srs_g1: &[G1Projective]) -> G1Projective {
    let mut result = G1Projective::default(); // 無限遠点（単位元）
    for (i, coeff) in coeffs.iter().enumerate() {
        result += srs_g1[i] * coeff;
    }
    result
}

/// [`evaluate_on_g1`] の G2 版。係数ベクトルと SRS（G2 点列）の内積で
/// `f(τ)·G2` を計算する。
fn evaluate_on_g2(coeffs: &[Fr], srs_g2: &[G2Projective]) -> G2Projective {
    let mut result = G2Projective::default();
    for (i, coeff) in coeffs.iter().enumerate() {
        result += srs_g2[i] * coeff;
    }
    result
}

/// simple 版の証明を生成する（α, β なし、zero-knowledge 化なし）。
///
/// `u_polys`: 各変数の `u_i(x)` の係数（A 行列由来）
/// `v_polys`: 各変数の `v_i(x)` の係数（B 行列由来）
/// `w_polys`: 各変数の `w_i(x)` の係数（C 行列由来）
/// `witness`: ウィットネス `[1, x, v1, v2, y, ...]` を `Fr` に変換したもの
/// `h_coeffs`: `h(x)` の係数を `Fr` に変換したもの
/// `srs`: trusted setup で生成した [`Srs`]
///
/// witness で各変数の多項式を重み付け合成して `A(x)/B(x)/W(x)` を作り、
/// [`evaluate_on_g1`] / [`evaluate_on_g2`] で SRS 上の点に評価する。
/// `[C]_1` には h·t の項を足す。計算量は O(num_vars · num_constraints)。
///
/// 各 `*_polys[i]` は QAP の Lagrange 補間結果で、末尾のゼロが除去されるため
/// 長さは `0..=num_constraints` の範囲でばらつく（ある行列に登場しない変数は
/// ゼロ多項式 = 長さ 1）。合成係数 `a/b/c_coeffs` は長さ `num_constraints` で確保し
/// 各多項式は存在する係数だけ加算する（補間の次数は高々 `num_constraints - 1`）。
///
/// 補間点数 = 多項式の係数の数 = R1CS 制約数 の3つが一致する設計で、
/// `num_constraints` は `u_polys[0]`（CS_ONE 列、全制約に現れるため最長）の
/// 長さから取っている。
pub fn prove_simple(
    u_polys: &[Vec<Fr>],
    v_polys: &[Vec<Fr>],
    w_polys: &[Vec<Fr>],
    witness: &[Fr],
    h_coeffs: &[Fr],
    srs: &Srs,
) -> Proof {
    // 1. 合成多項式 A(x) = Sigma a_i * u_i(x) の係数を計算する
    // 各 u_i(x) の係数に witness[i] を掛けて足し合わせる
    let num_constraints = u_polys[0].len(); // 補間点数 = R1CS 制約数（u_polys[0] = CS_ONE 列は全制約に出るので最長）
    let mut a_coeffs = vec![Fr::from(0u64); num_constraints];
    let mut b_coeffs = vec![Fr::from(0u64); num_constraints];
    let mut c_coeffs = vec![Fr::from(0u64); num_constraints];

    for (i, w_val) in witness.iter().enumerate() {
        // 各多項式は存在する係数だけ加算する（短い列も index out of bounds にならない）
        for (j, coeff) in u_polys[i].iter().enumerate() {
            a_coeffs[j] += *coeff * w_val;
        }
        for (j, coeff) in v_polys[i].iter().enumerate() {
            b_coeffs[j] += *coeff * w_val;
        }
        for (j, coeff) in w_polys[i].iter().enumerate() {
            c_coeffs[j] += *coeff * w_val;
        }
    }

    // 2. SRS 上で評価して楕円曲線点にする
    let a_g1 = evaluate_on_g1(&a_coeffs, &srs.g1_points);
    let b_g2 = evaluate_on_g2(&b_coeffs, &srs.g2_points);
    let c_g1_w = evaluate_on_g1(&c_coeffs, &srs.g1_points);

    // 3. h(tau)t(tau) を ht_points (SRS) との内積で計算する
    let ht_g1 = evaluate_on_g1(h_coeffs, &srs.ht_points);

    // 4. [C]_1 = W(tau) + h(tau) t(tau)
    let c_g1 = c_g1_w + ht_g1;

    Proof { a_g1, b_g2, c_g1 }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ec::PrimeGroup; // generator() のため
    use crate::adapter::{field_element_to_fr, polynomial_to_fr_vec, polys_to_fr_vecs};
    use crate::field::FieldElement;
    use crate::polynomial::Polynomial;
    use crate::qap::Qap;
    use crate::r1cs::ConstraintSystem;
    use crate::setup::generate_srs;
    use crate::verifier::verify_simple;
    use num_bigint::BigInt;

    /// x^3 + 5 の回路を作って、QAP -> SRS -> Prove -> Verify を通す
    #[test]
    fn test_simple_qap_end_to_end() {
        // === BN254 の曲線位数を p として使う ===
        let p = BigInt::parse_bytes(
            b"21888242871839275222246405745257275088548364400416034343698204186575808495617",
            10,
        )
        .unwrap();

        // === 1. R1CS を構築する（x^3 + 5 の回路）===
        let mut cs = ConstraintSystem::new();
        let one = FieldElement::new(BigInt::from(1), p.clone());
        cs.init_one(one);

        // 入力 x = 3
        let x = cs.alloc_variable();
        cs.assign(x, FieldElement::new(BigInt::from(3), p.clone()));

        // v1 = x * x
        let v1 = cs.mul(x, x);

        // v2 = v1 * x
        let v2 = cs.mul(v1, x);

        // y = v2 + 5
        let five = FieldElement::new(BigInt::from(5), p.clone());
        let _y = cs.add_const(v2, five);

        // === 2. QAP に変換 ===
        let qap = Qap::from_r1cs(&cs);

        // === 3. 多項式を Fr に変換 ===
        let u_polys = polys_to_fr_vecs(&qap.a_polys);
        let v_polys = polys_to_fr_vecs(&qap.b_polys);
        let w_polys = polys_to_fr_vecs(&qap.c_polys);

        // === 4. ウィットネスを Fr に変換 ===
        let witness_fe = cs.generate_witness();
        let witness: Vec<Fr> = witness_fe.iter().map(field_element_to_fr).collect();

        // === 5. h(x) を計算する ===
        // A(x) * B(x) - C(x) を Z(x) で割った商が h(x)
        let num_constraints = cs.constraints.len();

        // 合成多項式を自作の Polynomial として計算
        let zero_fe = FieldElement::new(BigInt::from(0), p.clone());
        let one_fe = FieldElement::new(BigInt::from(1), p.clone());

        let mut a_poly = Polynomial::new(vec![zero_fe.clone()]);
        let mut b_poly = Polynomial::new(vec![zero_fe.clone()]);
        let mut c_poly = Polynomial::new(vec![zero_fe.clone()]);

        for (i, w_val) in witness_fe.iter().enumerate() {
            let scaled_a = qap.a_polys[i].scale(w_val);
            a_poly = &a_poly + &scaled_a;

            let scaled_b = qap.b_polys[i].scale(w_val);
            b_poly = &b_poly + &scaled_b;

            let scaled_c = qap.c_polys[i].scale(w_val);
            c_poly = &c_poly + &scaled_c;
        }

        // P(x) = A(x) * B(x) - C(x)
        let ab = &a_poly * &b_poly;
        let minus_one = &zero_fe - &one_fe;
        let neg_c = c_poly.scale(&minus_one);
        let p_poly = &ab + &neg_c;

        // Z(x) = (x-0)(x-1)(x-2)...
        // 注意: 自作 QAP の補間点は 0, 1, 2, ... 始まり
        let mut z_poly = Polynomial::new(vec![FieldElement::new(BigInt::from(1), p.clone())]);
        for i in 0..num_constraints {
            let i_fe = FieldElement::new(BigInt::from(i), p.clone());
            let neg_i = &zero_fe - &i_fe;
            let term = Polynomial::new(vec![neg_i, one_fe.clone()]);
            z_poly = &z_poly * &term;
        }

        // h(x) = P(x) / Z(x)
        let (h_poly, remainder) = p_poly.div_rem(&z_poly);

        // 余りがゼロであることを確認
        assert!(remainder.is_zero(), "P(x) が Z(x) で割り切れません");

        // h(x) の係数を Fr に変換
        let h_coeffs = polynomial_to_fr_vec(&h_poly);

        // === 6. SRS を生成 ===
        let tau = Fr::from(777u64); // テスト用の固定値
        let srs = generate_srs(tau, num_constraints);

        // === 7. 証明を生成 ===
        let proof = prove_simple(&u_polys, &v_polys, &w_polys, &witness, &h_coeffs, &srs);

        // === 8. 検証 ===
        assert!(verify_simple(&proof), "検証に失敗しました");
    }

    #[test]
    fn test_evaluate_on_g1_linear() {
        // 係数 [1, 2] を SRS [G, τG] で評価 → (1 + 2τ)·G = G + 2·(τG)
        let g1 = G1Projective::generator();
        let tau = Fr::from(5u64);
        let srs_g1 = vec![g1, g1 * tau];
        let coeffs = vec![Fr::from(1u64), Fr::from(2u64)];

        let result = evaluate_on_g1(&coeffs, &srs_g1);
        let expected = g1 + (g1 * tau) * Fr::from(2u64);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_evaluate_on_g1_all_zero() {
        // 係数がすべて 0 → 無限遠点（単位元）
        let g1 = G1Projective::generator();
        let tau = Fr::from(5u64);
        let srs_g1 = vec![g1, g1 * tau, g1 * tau * tau];
        let coeffs = vec![Fr::from(0u64); 3];

        let result = evaluate_on_g1(&coeffs, &srs_g1);
        assert_eq!(result, G1Projective::default());
    }

    #[test]
    fn test_evaluate_on_g2_linear() {
        // G2 側も同じ内積規則: 係数 [1, 2] を SRS [G2, τG2] で評価 → (1 + 2τ)·G2
        let g2 = G2Projective::generator();
        let tau = Fr::from(5u64);
        let srs_g2 = vec![g2, g2 * tau];
        let coeffs = vec![Fr::from(1u64), Fr::from(2u64)];

        let result = evaluate_on_g2(&coeffs, &srs_g2);
        let expected = g2 + (g2 * tau) * Fr::from(2u64);
        assert_eq!(result, expected);
    }
}
