mod field;
mod polynomial;
mod qap;
mod r1cs;
mod adapter;
mod setup;

use field::FieldElement;
use num_bigint::BigInt;

use crate::{
    polynomial::Polynomial,
    qap::Qap,
    r1cs::{ConstraintSystem, LinearCombination},
};

fn main() {
    let p = BigInt::from(17);
    let mut cs = ConstraintSystem::new();

    // 0. 定数 CS_ONE の初期化
    // これを忘れると Index 0 が None になり panic する
    let one = FieldElement::new(BigInt::from(1), p.clone());
    cs.init_one(one);

    // 入力 x = 3
    let x = cs.alloc_variable();
    cs.assign(x, FieldElement::new(BigInt::from(3), p.clone()));

    // v1 = x * x
    let v1 = cs.mul(x, x);

    // v2 = v1 * x
    let v2 = cs.mul(v1, x);

    // 値計算： y = v2 + 5
    let five = FieldElement::new(BigInt::from(5), p.clone());
    let y = cs.add_const(v2, five);

    // 検証
    println!("制約数: {}", cs.constraints.len()); // mul 2回 + add 1回 = 3つになるはず
    let mut witness = cs.generate_witness();
    println!("計算結果 y = {}", witness[y.0]);

    if is_satisfied(&cs, &witness) {
        println!("x^3 + 5 = y (抽象化版) 成功！");
    }

    // === 不正な Witness に変更 ===
    // witness[4] = FieldElement::new(BigInt::from(999), p.clone());
    // println!("改ざん後の y: {:?}", witness[4].value);
    // === ここまで ===

    // QAP 変換
    let qap = Qap::from_r1cs(&cs);

    println!("変数の数：{}", qap.a_polys.len());
    println!("QAP A多項式の数: {}", qap.a_polys.len());
    println!("A_poly[0] (ONE): {:?}", qap.a_polys[0]);
    println!("A_poly[1] (x): {:?}", qap.a_polys[1]);
    println!("A_poly[2] (v1): {:?}", qap.a_polys[2]);
    println!("A_poly[3] (v2): {:?}", qap.a_polys[3]);
    println!("A_poly[4] (y): {:?}", qap.a_polys[4]);

    // 1. 合成多項式 A(x)、B(x)、C(x) を作る
    // 公式: A(x) = sum( witness[i] * A_poly_i(x) )
    // つまり、Witness の値で重み付けして足し合わせる
    let mut a_x = Polynomial::new(vec![]);
    let mut b_x = Polynomial::new(vec![]);
    let mut c_x = Polynomial::new(vec![]);

    // Witness は [ONE, x, v1, v2, y] の順に並んでいる
    for (i, w_val) in witness.iter().enumerate() {
        // Aの合成
        let poly_a = &qap.a_polys[i];
        let scaled_a = poly_a.scale(w_val.clone());
        a_x = &a_x + &scaled_a;

        // Bの合成
        let poly_b = &qap.b_polys[i];
        let scaled_b = poly_b.scale(w_val.clone());
        b_x = &b_x + &scaled_b;

        // Cの合成
        let poly_c = &qap.c_polys[i];
        let scaled_c = poly_c.scale(w_val.clone());
        c_x = &c_x + &scaled_c;
    }

    println!("A(x) 次数: {}", a_x.degree());

    // 2. P(x) = A(x) * B(x) - C(x)
    let prod_ab = &a_x * &b_x;

    // 引き算: P = AB + (-1 * C)
    let minus_one = &FieldElement::new(BigInt::from(0), p.clone())
        - &FieldElement::new(BigInt::from(1), p.clone());
    let neg_c = c_x.scale(minus_one);

    let p_x = &prod_ab + &neg_c;
    println!("P(x) 計算完了. 次数: {}", p_x.degree());

    // 3. ターゲット多項式 Z(x) を作る
    // Z(x) = (x - 0)(x - 1)...(x - (制約数 - 1))
    // x=0, 1, 2 で必ず0になる多項式
    let num_constraints = cs.constraints.len();
    let mut z_x = Polynomial::new(vec![FieldElement::new(BigInt::from(1), p.clone())]); // 初期値1

    let one_fe = FieldElement::new(BigInt::from(1), p.clone());
    let zero_fe = FieldElement::new(BigInt::from(0), p.clone());

    for i in 0..num_constraints {
        // (x - i) を作る -> [-i, 1]
        let i_fe = FieldElement::new(BigInt::from(i), p.clone());
        let neg_i = &zero_fe - &i_fe;

        let term = Polynomial::new(vec![neg_i, one_fe.clone()]);
        z_x = &z_x * &term;
    }

    println!("Z(x) 計算完了. 次数: {}", z_x.degree());

    // 4. 割り算: H(x) = P(x) / Z(x)
    let (h_x, remainder) = p_x.div(&z_x);

    println!("H(x) 次数: {}", h_x.degree());
    println!("割り算の余り (次数): {}", remainder.degree());

    // 余りがゼロ（係数が空 OR すべて0）なら証明成功
    let is_valid_proof = remainder
        .coefficients
        .iter()
        .all(|c| c.value == BigInt::from(0));

    if is_valid_proof {
        println!("🎉 大勝利！ H(x) が割り切れました。");
        println!("これにて『計算が正しいこと』の数学的証明が完成です。");
    } else {
        println!("💀 失敗... 余りが出てしまいました。Witnessか回路が間違っています。");
        println!("余り: {:?}", remainder);
    }
}

// 指定した Witness が、 ConstraintSystem のすべての制約を満たしているかチェックする
fn is_satisfied(cs: &ConstraintSystem, witness: &Vec<FieldElement>) -> bool {
    for constraint in &cs.constraints {
        let a_val = evaluate_lc(&constraint.a, witness);
        let b_val = evaluate_lc(&constraint.b, witness);
        let c_val = evaluate_lc(&constraint.c, witness);

        // A * B == C かどうかを判定
        if &(&a_val * &b_val) != &c_val {
            return false;
        }
    }
    true
}

// LinearCombination（線形結合）に Witness を代入して値を計算する
fn evaluate_lc(lc: &LinearCombination, witness: &Vec<FieldElement>) -> FieldElement {
    let p = witness[0].p.clone();

    // 1. 合計値を 0 で初期化
    let mut total = FieldElement::new(BigInt::from(0), p.clone());

    // 2. LC に含まれる「項（term）」を一つずつ取り出す
    for (var, coeff) in &lc.terms {
        // 3. var.0 (インデックス) を使って、witness ベクトルから実際の値を取り出す
        let val = &witness[var.0];
        // 4. (係数 × 実際の値) を計算する
        let product = coeff * val;
        // 5. これを合計に足していく
        total = &total + &product;
    }
    total
}
