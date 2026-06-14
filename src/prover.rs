//! QAP + SRS から証明 `(A, B, C)` を生成する Prover。
//!
//! Groth16 実装の Layer 3（プロトコル）。Witness と QAP 多項式の係数を
//! SRS 上で評価し、楕円曲線点として証明を構築する。
//!
//! ## 主要型
//! - [`Proof`]: simple 版の `(A ∈ G1, B ∈ G2, C ∈ G1)` 三点組
//! - [`Groth16Proof`]: 本式の `(A, B, C)` （ランダム化込み）
//!
//! ## 主要関数
//! - [`prove_simple`]: simple QAP 版の証明生成
//! - [`prove`]: 本式 Groth16 の証明生成（ランダム r, s 込み）
//!
//! ## 注意
//! simple 版（[`prove_simple`]）と本式（[`prove`]）が併存している。
//! 本式が E2E で通った後、simple 版は Phase 6b で削除予定。

use ark_bn254::{Fr, G1Projective, G2Projective};

use crate::setup::{ProvingKey, QapFr, Srs};

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

/// 本式 Groth16 の証明（楕円曲線上の 3 点）。
///
/// simple 版 [`Proof`] と点の構成は同じ `(A, B, C)` だが、各点に α/β/γ/δ と
/// ランダム値 r, s が織り込まれており、検証は 4 ペアリングの等式で行う。
/// r, s により同じ witness でも proof が毎回変わる（zero-knowledge）。
#[derive(Debug)]
#[allow(dead_code)] // Phase 6b の verifier / main 配線（B4）で使う。それまで dead_code を許可
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

/// 本式 Groth16 の証明を生成する（ランダム化 r, s 込み）。
///
/// `pk`: trusted setup で生成した [`ProvingKey`]
/// `qap_fr`: Fr 変換済みの QAP（変数ごとの `u_i / v_i / w_i` 係数）
/// `witness`: `[1, 公開入力..., 秘密/中間...]` を `Fr` 変換したもの
/// `h_coeffs`: `h(x) = (A·B − C)/Z` の係数を `Fr` 変換したもの
/// `r`, `s`: zero-knowledge 用のランダム値（テスト再現性のため引数で受ける）
///
///
/// public/private の境界 ℓ+1 は `witness.len() − pk.private_query.len()` から導く。
/// A, B(G2), B_1(G1) を組み、C は「private wire 合成 + h·t/δ」に
/// `s·A + r·B_1 − r·s·δ` を足して構成する。計算量は O(num_vars · num_constraints)。
#[allow(dead_code)] // B4 で main に配線したら外す
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
    use crate::setup::{generate_groth16_keys, generate_srs, ToxicWaste, VerifyingKey};
    use crate::verifier::{verify, verify_simple};
    use ark_ec::PrimeGroup; // generator() のため
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
