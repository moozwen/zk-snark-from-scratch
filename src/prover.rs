use ark_bn254::{Fr, G1Projective, G2Projective};
use ark_ec::CurveGroup;
use ark_ff::Field;

use crate::setup::Srs;

/// シンプル版の証明（楕円曲線上の3点）
pub struct Proof {
    pub a_g1: G1Projective, // [A]_1
    pub b_g2: G2Projective, // [B]_2
    pub c_g1: G1Projective, // [C]_1
}

/// 多項式の係数ベクトルを SRS 上で評価する
/// f(tau) = coeffs[0]*G + coeffs[1]* tau G + coeffs[2] * tau^2 G + ...
/// つまり SRS との内積を計算する
fn evaluate_on_g1(coeffs: &[Fr], srs_g1: &[G1Projective]) -> G1Projective {
    let mut result = G1Projective::default(); // 無限遠点（単位元）
    for (i, coeff) in coeffs.iter().enumerate() {
        result += srs_g1[i] * coeff;
    }
    result
}

fn evaluate_on_g2(coeffs: &[Fr], srs_g2: &[G2Projective]) -> G2Projective {
    let mut result = G2Projective::default();
    for (i, coeff) in coeffs.iter().enumerate() {
        result += srs_g2[i] * coeff;
    }
    result
}

/// シンプル版の証明を生成する（alpha, beta なし）
///
/// u_polys: 各変数の u_i(x) の係数（L行列由来）
/// v_polys: 各変数の v_i(x) の係数（R行列由来）
/// w_polys: 各変数の w_i(x) の係数（O行列由来）
/// witness: ウィットネス [1, x, v1, v2, y, ...] を Fr に変換したもの
/// h_coeffs: h(x) の係数を Fr に変換したもの
/// srs: Trusted Setup で生成した SRS
pub fn prove_simple(
    u_polys: &Vec<Vec<Fr>>,
    v_polys: &Vec<Vec<Fr>>,
    w_polys: &Vec<Vec<Fr>>,
    witness: &Vec<Fr>,
    h_coeffs: &Vec<Fr>,
    srs: &Srs,
) -> Proof {
    // 1. 合成多項式 A(x) = Sigma a_i * u_i(x) の係数を計算する
    // 各 u_i(x) の係数に witness[i] を掛けて足し合わせる
    let num_constraints = u_polys[0].len(); // 多項式の係数の数
    let mut a_coeffs = vec![Fr::from(0u64); num_constraints];
    let mut b_coeffs = vec![Fr::from(0u64); num_constraints];
    let mut c_coeffs = vec![Fr::from(0u64); num_constraints];

    for (i, w_val) in witness.iter().enumerate() {
        for j in 0..num_constraints {
            // u_polys[i] が i番目の変数の多項式の係数
            // 係数が足りない場合は 0 として扱う
            if j < u_polys[i].len() {
                a_coeffs[j] += u_polys[i][j] * w_val;
            }
            if j < v_polys[i].len() {
                b_coeffs[j] += v_polys[i][j] * w_val;
            }
            if j < w_polys[i].len() {
                c_coeffs[j] += w_polys[i][j] * w_val;
            }
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
    use crate::adapter::{field_element_to_fr, polynomial_to_fr_vec, polys_to_fr_vecs};
    use crate::field::FieldElement;
    use crate::polynomial::Polynomial;
    use crate::qap::Qap;
    use crate::r1cs::{ConstraintSystem, CS_ONE};
    use crate::setup::generate_srs;
    use crate::verifier::verify_simple;
    use num_bigint::{BigInt, BigUint};

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
        let y = cs.add_const(v2, five);

        // === 2. QAP に変換 ===
        let qap = Qap::from_r1cs(&cs);

        // === 3. 多項式を Fr に変換 ===
        let u_polys = polys_to_fr_vecs(&qap.a_polys);
        let v_polys = polys_to_fr_vecs(&qap.b_polys);
        let w_polys = polys_to_fr_vecs(&qap.c_polys);

        // === 4. ウィットネスを Fr に変換 ===
        let witness_fe = cs.generate_witness();
        let witness: Vec<Fr> = witness_fe.iter().map(|w| field_element_to_fr(w)).collect();

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
            let scaled_a = qap.a_polys[i].scale(w_val.clone());
            a_poly = &a_poly + &scaled_a;

            let scaled_b = qap.b_polys[i].scale(w_val.clone());
            b_poly = &b_poly + &scaled_b;

            let scaled_c = qap.c_polys[i].scale(w_val.clone());
            c_poly = &c_poly + &scaled_c;
        }

        // P(x) = A(x) * B(x) - C(x)
        let ab = &a_poly * &b_poly;
        let minus_one = &zero_fe - &one_fe;
        let neg_c = c_poly.scale(minus_one);
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
        let is_divisible = remainder
            .coefficients
            .iter()
            .all(|c| c.value == BigInt::from(0));
        assert!(is_divisible, "P(x) が Z(x) で割り切れません");

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
}
