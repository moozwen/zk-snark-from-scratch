//! QAP + proving key から証明 `(A, B, C)` を生成する Prover。
//!
//! Groth16 実装の Layer 3（プロトコル）。Witness と QAP 多項式の係数を
//! proving key の点群上で評価し、楕円曲線点として証明を構築する。
//!
//! ## 主要型
//! - [`Groth16Proof`]: Groth16 の証明 `(A, B, C)`（ランダム化込み）
//!
//! ## 主要関数
//! - [`prove`]: Groth16 の証明生成（ランダム r, s 込み）

use ark_bn254::{Fr, G1Projective, G2Projective};

use crate::setup::{ProvingKey, QapFr};

/// 本式 Groth16 の証明（楕円曲線上の 3 点）。
///
/// 各点 `(A, B, C)` に α/β/γ/δ とランダム値 r, s が織り込まれており、
/// 検証は 4 ペアリングの等式で行う。
/// r, s により同じ witness でも proof が毎回変わる（zero-knowledge）。
#[derive(Debug)]
pub struct Groth16Proof {
    /// `[A]_1 = [ α + Σ_i a_i·u_i(τ) + r·δ ]_1`
    pub a: G1Projective,
    /// `[B]_2 = [ β + Σ_i a_i·v_i(τ) + s·δ ]_2`
    pub b: G2Projective,
    /// `[C]_1`。private wire 合成 + h·t/δ に `s·A + r·B_1 − r·s·δ` を足したもの。
    pub c: G1Projective,
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

/// 本式 Groth16 の証明を生成する（ランダム化 r, s 込み）。
///
/// `pk`: trusted setup で生成した [`ProvingKey`]
/// `qap_fr`: Fr 変換済みの QAP（変数ごとの `u_i / v_i / w_i` 係数）
/// `witness`: `[1, 公開入力..., 秘密/中間...]` を `Fr` 変換したもの
/// `h_coeffs`: `h(x) = (A·B − C)/Z` の係数を `Fr` 変換したもの
/// `r`, `s`: zero-knowledge 用のランダム値（毎回ランダムに引くことで同じ witness でも
///   proof が変わる。テスト再現性のため引数で受ける）
///
/// public/private の境界 ℓ+1 は `witness.len() − pk.private_query.len()` から導く。
/// `A = α + Σ a_i·u_i(τ) + r·δ`、`B = β + Σ a_i·v_i(τ) + s·δ`（C 計算用に G1 版 B_1 も作る）。
/// `C` は private wire 合成 `Σ_{i>ℓ} a_i·(β u_i+α v_i+w_i)(τ)/δ` に `h(τ)t(τ)/δ` を足し、
/// さらに `s·A + r·B_1 − r·s·δ` を加える（最後の `−r·s·δ` 項が r,s の二重計上を相殺）。
/// 計算量は O(num_vars · num_constraints)。
pub fn prove(
    pk: &ProvingKey,
    qap_fr: &QapFr,
    witness: &[Fr],
    h_coeffs: &[Fr],
    r: Fr,
    s: Fr,
) -> Groth16Proof {
    // 1. 全 wire の合成多項式 A(x)=Σ a_i u_i(x), B(x)=Σ a_i v_i(x) の係数を計算する。
    //    補間次数は高々 n-1 なので係数列の長さは n（= pk.tau_g1.len()）で確保する。
    let n = pk.tau_g1.len();
    let mut a_coeffs = vec![Fr::from(0u64); n];
    let mut b_coeffs = vec![Fr::from(0u64); n];
    for (i, w_val) in witness.iter().enumerate() {
        for (j, coeff) in qap_fr.a_polys[i].iter().enumerate() {
            a_coeffs[j] += *coeff * w_val;
        }
        for (j, coeff) in qap_fr.b_polys[i].iter().enumerate() {
            b_coeffs[j] += *coeff * w_val;
        }
    }

    // 2. τ 冪 SRS 上で評価。B は C 計算用に G1 でも作る（B_1）。
    let sum_a_g1 = evaluate_on_g1(&a_coeffs, &pk.tau_g1); // Σ a_i u_i(τ)·G1
    let sum_b_g2 = evaluate_on_g2(&b_coeffs, &pk.tau_g2); // Σ a_i v_i(τ)·G2
    let sum_b_g1 = evaluate_on_g1(&b_coeffs, &pk.tau_g1); // Σ a_i v_i(τ)·G1

    // 3. A = α + Σ a_i u_i(τ) + r·δ、B = β + Σ a_i v_i(τ) + s·δ、B_1 は B の G1 版。
    let a = pk.alpha_g1 + sum_a_g1 + pk.delta_g1 * r;
    let b = pk.beta_g2 + sum_b_g2 + pk.delta_g2 * s;
    let b1 = pk.beta_g1 + sum_b_g1 + pk.delta_g1 * s;

    // 4. C の /δ 項: private wire 合成 Σ_{i>ℓ} a_i·[(βu+αv+w)/δ]_1 + h(τ)t(τ)/δ。
    //    private_query[j] は wire i = num_public + j に対応する。
    let num_public = witness.len() - pk.private_query.len();
    let mut c_div_delta = G1Projective::default();
    for (j, point) in pk.private_query.iter().enumerate() {
        c_div_delta += *point * witness[num_public + j];
    }
    c_div_delta += evaluate_on_g1(h_coeffs, &pk.h_query); // h(τ)t(τ)/δ

    // 5. C = (上記)/δ + s·A + r·B_1 − r·s·δ。
    let c = c_div_delta + a * s + b1 * r - pk.delta_g1 * (r * s);

    Groth16Proof { a, b, c }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::adapter::{field_element_to_fr, polynomial_to_fr_vec, polys_to_fr_vecs};
    use crate::field::FieldElement;
    use crate::polynomial::Polynomial;
    use crate::qap::Qap;
    use crate::r1cs::{ConstraintSystem, LinearCombination, CS_ONE};
    use crate::setup::{generate_groth16_keys, ToxicWaste, VerifyingKey};
    use crate::verifier::verify;
    use ark_ec::PrimeGroup; // generator() のため
    use num_bigint::BigInt;

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

    // === 本式 Groth16 の E2E / soundness / ZK テスト ===

    /// x³ + 5 = y（y を public 出力）の本式 Groth16 一式。x = 3, y = 32。
    struct Groth16Fixture {
        pk: ProvingKey,
        vk: VerifyingKey,
        qap_fr: QapFr,
        witness: Vec<Fr>,       // [CS_ONE, y, x, v1, v2]
        h_coeffs: Vec<Fr>,
        public_inputs: Vec<Fr>, // [y]
    }

    fn build_x3_plus5_fixture() -> Groth16Fixture {
        let p = BigInt::parse_bytes(
            b"21888242871839275222246405745257275088548364400416034343698204186575808495617",
            10,
        )
        .unwrap();
        let fe = |v: u64| FieldElement::new(BigInt::from(v), p.clone());
        let one = fe(1);
        let zero = fe(0);

        // === R1CS（public: CS_ONE, y / private: x, v1, v2）===
        let mut cs = ConstraintSystem::new();
        cs.init_one(one.clone());
        let y = cs.alloc_public_input(); // public 出力を前方固め
        cs.assign(y, fe(32)); // 3³ + 5 = 32
        let x = cs.alloc_variable(); // private 入力
        cs.assign(x, fe(3));
        let v1 = cs.mul(x, x); // 9
        let v2 = cs.mul(v1, x); // 27
        // 制約: (v2 + 5)·1 = y
        let mut lc_a = LinearCombination::new();
        lc_a.add_term(v2, one.clone());
        lc_a.add_term(CS_ONE, fe(5));
        let mut lc_b = LinearCombination::new();
        lc_b.add_term(CS_ONE, one.clone());
        let mut lc_c = LinearCombination::new();
        lc_c.add_term(y, one.clone());
        cs.enforce(lc_a, lc_b, lc_c);

        let num_constraints = cs.constraints.len(); // 3
        let num_public = cs.num_public_variables; // 2 (CS_ONE + y)

        // === QAP → Fr ===
        let qap = Qap::from_r1cs(&cs);
        let qap_fr = QapFr {
            a_polys: polys_to_fr_vecs(&qap.a_polys),
            b_polys: polys_to_fr_vecs(&qap.b_polys),
            c_polys: polys_to_fr_vecs(&qap.c_polys),
        };

        // === witness ===
        let witness_fe = cs.generate_witness();
        let witness: Vec<Fr> = witness_fe.iter().map(field_element_to_fr).collect();

        // === h(x) = (A·B − C) / Z ===
        let mut a_poly = Polynomial::new(vec![zero.clone()]);
        let mut b_poly = Polynomial::new(vec![zero.clone()]);
        let mut c_poly = Polynomial::new(vec![zero.clone()]);
        for (i, w_val) in witness_fe.iter().enumerate() {
            a_poly = &a_poly + &qap.a_polys[i].scale(w_val);
            b_poly = &b_poly + &qap.b_polys[i].scale(w_val);
            c_poly = &c_poly + &qap.c_polys[i].scale(w_val);
        }
        let ab = &a_poly * &b_poly;
        let minus_one = &zero - &one;
        let neg_c = c_poly.scale(&minus_one);
        let p_poly = &ab + &neg_c;
        let mut z_poly = Polynomial::new(vec![one.clone()]);
        for i in 0..num_constraints {
            let neg_i = &zero - &fe(i as u64);
            z_poly = &z_poly * &Polynomial::new(vec![neg_i, one.clone()]);
        }
        let (h_poly, remainder) = p_poly.div_rem(&z_poly);
        assert!(remainder.is_zero(), "P(x) が Z(x) で割り切れません");
        let h_coeffs = polynomial_to_fr_vec(&h_poly);

        // === 鍵生成 ===
        let toxic = ToxicWaste {
            alpha: Fr::from(11u64),
            beta: Fr::from(13u64),
            gamma: Fr::from(17u64),
            delta: Fr::from(19u64),
            tau: Fr::from(23u64),
        };
        let (pk, vk) = generate_groth16_keys(&qap_fr, num_constraints, num_public, &toxic);

        // public_inputs = a_1..a_ℓ = [y]
        let public_inputs = vec![witness[y.0]];

        Groth16Fixture {
            pk,
            vk,
            qap_fr,
            witness,
            h_coeffs,
            public_inputs,
        }
    }

    #[test]
    fn test_groth16_end_to_end_accepts() {
        let f = build_x3_plus5_fixture();
        let proof = prove(
            &f.pk,
            &f.qap_fr,
            &f.witness,
            &f.h_coeffs,
            Fr::from(5u64),
            Fr::from(7u64),
        );
        assert!(
            verify(&f.vk, &f.public_inputs, &proof),
            "正当な証明が reject された"
        );
    }

    #[test]
    fn test_groth16_accepts_for_different_randomness() {
        // r, s を変えても同じ witness の proof は accept される（ZK ランダム化の健全性）
        let f = build_x3_plus5_fixture();
        let proof1 = prove(
            &f.pk,
            &f.qap_fr,
            &f.witness,
            &f.h_coeffs,
            Fr::from(1u64),
            Fr::from(2u64),
        );
        let proof2 = prove(
            &f.pk,
            &f.qap_fr,
            &f.witness,
            &f.h_coeffs,
            Fr::from(99u64),
            Fr::from(100u64),
        );
        assert!(verify(&f.vk, &f.public_inputs, &proof1));
        assert!(verify(&f.vk, &f.public_inputs, &proof2));
        // proof が毎回変化していること（zero-knowledge）
        assert_ne!(proof1.a, proof2.a);
    }

    #[test]
    fn test_groth16_rejects_wrong_public_input() {
        let f = build_x3_plus5_fixture();
        let proof = prove(
            &f.pk,
            &f.qap_fr,
            &f.witness,
            &f.h_coeffs,
            Fr::from(5u64),
            Fr::from(7u64),
        );
        // y = 32 のはずを 33 と主張 → reject
        assert!(!verify(&f.vk, &[Fr::from(33u64)], &proof));
    }

    #[test]
    fn test_groth16_rejects_tampered_witness() {
        // 正しい h と不整合な witness（x を改ざん）で証明 → reject
        let f = build_x3_plus5_fixture();
        let mut bad = f.witness.clone();
        bad[2] += Fr::from(1u64); // index 2 = x
        let proof = prove(
            &f.pk,
            &f.qap_fr,
            &bad,
            &f.h_coeffs,
            Fr::from(5u64),
            Fr::from(7u64),
        );
        assert!(!verify(&f.vk, &f.public_inputs, &proof));
    }

    #[test]
    fn test_groth16_rejects_tampered_proof() {
        let f = build_x3_plus5_fixture();
        // A を細工 → reject
        let mut proof = prove(
            &f.pk,
            &f.qap_fr,
            &f.witness,
            &f.h_coeffs,
            Fr::from(5u64),
            Fr::from(7u64),
        );
        proof.a += G1Projective::generator();
        assert!(!verify(&f.vk, &f.public_inputs, &proof));
        // C を細工 → reject
        let mut proof2 = prove(
            &f.pk,
            &f.qap_fr,
            &f.witness,
            &f.h_coeffs,
            Fr::from(5u64),
            Fr::from(7u64),
        );
        proof2.c += G1Projective::generator();
        assert!(!verify(&f.vk, &f.public_inputs, &proof2));
    }
}
