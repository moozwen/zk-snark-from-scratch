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
