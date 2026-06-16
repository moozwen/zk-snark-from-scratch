mod adapter;
mod field;
mod polynomial;
mod prover;
mod qap;
mod r1cs;
mod setup;
mod verifier;

use field::FieldElement;
use num_bigint::BigInt;

use ark_bn254::Fr;

use crate::{
    adapter::{field_element_to_fr, polynomial_to_fr_vec, polys_to_fr_vecs},
    polynomial::Polynomial,
    prover::prove,
    qap::Qap,
    r1cs::{ConstraintSystem, LinearCombination, CS_ONE},
    setup::{generate_groth16_keys, QapFr, ToxicWaste},
    verifier::verify,
};

fn main() {
    println!("=== zk-snark-from-scratch: x^3 + 5 Groth16 proof demo ===\n");

    // BN254 のスカラー体 位数
    let p = BigInt::parse_bytes(
        b"21888242871839275222246405745257275088548364400416034343698204186575808495617",
        10,
    )
    .unwrap();
    let fe = |v: u64| FieldElement::new(BigInt::from(v), p.clone());

    // Step 1: R1CS (y = x^3 + 5 with x = 3 は秘密入力 / y = 32 公開出力)
    println!("Step 1: Building R1CS for y = x^3 + 5 (x = 3 private, y = 32 public)...");
    let mut cs = ConstraintSystem::new();
    cs.init_one(fe(1));
    let y = cs.alloc_public_input(); // 公開出力 y を前方に固める
    cs.assign(y, fe(32));
    let x = cs.alloc_variable(); // 秘密入力 x
    cs.assign(x, fe(3));
    let v1 = cs.mul(x, x);
    let v2 = cs.mul(v1, x);
    // 制約: (v2 + 5) · 1 = y
    let mut lc_a = LinearCombination::new();
    lc_a.add_term(v2, fe(1));
    lc_a.add_term(CS_ONE, fe(5));
    let mut lc_b = LinearCombination::new();
    lc_b.add_term(CS_ONE, fe(1));
    let mut lc_c = LinearCombination::new();
    lc_c.add_term(y, fe(1));
    cs.enforce(lc_a, lc_b, lc_c);

    let num_constraints = cs.constraints.len();
    let num_public = cs.num_public_variables;
    println!(
        "  {} constraints, {} variables ({} public incl. CS_ONE)",
        num_constraints, cs.next_var_index, num_public
    );

    // Step 2: R1CS -> QAP -> Fr
    println!("\nStep 2: Converting R1CS -> QAP...");
    let qap = Qap::from_r1cs(&cs);
    let qap_fr = QapFr {
        a_polys: polys_to_fr_vecs(&qap.a_polys),
        b_polys: polys_to_fr_vecs(&qap.b_polys),
        c_polys: polys_to_fr_vecs(&qap.c_polys),
    };

    // Step 3: h(x) = (A(x)*B(x) - C(x)) / Z(x)
    println!("\nStep 3: Computing h(x)...");
    let witness_fe = cs.generate_witness();
    let h_poly = compute_h_poly(&qap, &witness_fe, num_constraints, &p);
    println!("  h(x) degree: {}", h_poly.degree());

    // Step 4: Trusted setup（本式 pk/vk。デモ用に toxic waste は固定値、本番は破棄）
    println!("\nStep 4: Generating proving/verifying keys (demo toxic waste)...");
    let toxic = ToxicWaste {
        alpha: Fr::from(11u64),
        beta: Fr::from(13u64),
        gamma: Fr::from(17u64),
        delta: Fr::from(19u64),
        tau: Fr::from(23u64),
    };
    let (pk, vk) = generate_groth16_keys(&qap_fr, num_constraints, num_public, &toxic);

    // Step 5: Prove（r, s はデモ用固定。本番では毎回ランダムに引く = zero-knowledge）
    println!("\nStep 5: Generating proof (r, s fixed for demo)...");
    let witness: Vec<Fr> = witness_fe.iter().map(field_element_to_fr).collect();
    let h_coeffs = polynomial_to_fr_vec(&h_poly);
    let proof = prove(
        &pk,
        &qap_fr,
        &witness,
        &h_coeffs,
        Fr::from(5u64),
        Fr::from(7u64),
    );

    // Step 6: Verify
    println!("\nStep 6: Verifying proof...");
    let public_inputs = vec![witness[y.0]];
    if verify(&vk, &public_inputs, &proof) {
        println!("  OK! Proof verified");
    } else {
        println!("  NG..Proof rejected");
    }
}

/// h(x) = (A(x)*B(x) - C(x)) / Z(x) を計算する
fn compute_h_poly(
    qap: &Qap,
    witness: &[FieldElement],
    num_constraints: usize,
    p: &BigInt,
) -> Polynomial {
    let zero = FieldElement::new(0, p.clone());
    let one = FieldElement::new(1, p.clone());

    // A(x), B(x), C(x) = sum_i witness[i] * poly_i(x)
    let mut a = Polynomial::new(vec![zero.clone()]);
    let mut b = Polynomial::new(vec![zero.clone()]);
    let mut c = Polynomial::new(vec![zero.clone()]);
    for (i, w) in witness.iter().enumerate() {
        a = &a + &qap.a_polys[i].scale(w);
        b = &b + &qap.b_polys[i].scale(w);
        c = &c + &qap.c_polys[i].scale(w);
    }

    // P(x) = A(x)*B(x) - C(x)
    let minus_one = &zero - &one;
    let p_poly = &(&a * &b) + &c.scale(&minus_one);

    // Z(x) = (x - 0)(x - 1)...(x - (n - 1))
    let mut z_poly = Polynomial::new(vec![one.clone()]);
    for i in 0..num_constraints {
        let neg_i = &zero - &FieldElement::new(i, p.clone());
        z_poly = &z_poly * &Polynomial::new(vec![neg_i, one.clone()]);
    }

    // h(x) = P(x) / Z(x)
    let (h, remainder) = p_poly.div_rem(&z_poly);
    assert!(remainder.is_zero(), "P(x) is not divisible by Z(x)");
    h
}
