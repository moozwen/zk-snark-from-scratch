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
    prover::prove_simple,
    qap::Qap,
    r1cs::ConstraintSystem,
    setup::generate_srs,
    verifier::verify_simple,
};

fn main() {
    println!("=== zk-snark-from-scratch: x^3 + 5 proof demo ===\n");

    // BN254 のスカラー体 位数
    let p = BigInt::parse_bytes(
        b"21888242871839275222246405745257275088548364400416034343698204186575808495617",
        10,
    )
    .unwrap();

    // Step 1: R1CS (y = x^3 + 5 with x = 3, y = 32)
    println!("Step 1: Building R1CS for y = x^3 + 5 (x = 3)...");
    let mut cs = ConstraintSystem::new();
    cs.init_one(FieldElement::new(BigInt::from(1), p.clone()));

    let x = cs.alloc_variable();
    cs.assign(x, FieldElement::new(BigInt::from(3), p.clone()));
    let v1 = cs.mul(x, x);
    let v2 = cs.mul(v1, x);
    let _y = cs.add_const(v2, FieldElement::new(BigInt::from(5), p.clone()));

    println!(
        "  {} constraints, {} variables",
        cs.constraints.len(),
        cs.next_var_index
    );

    // Step 2: R1CS -> QAP
    println!("\nStep 2: Converting R1CS -> QAP...");
    let qap = Qap::from_r1cs(&cs);

    // Step 3: h(x) = (A(x)*B(x) - C(x)) / Z(x)
    println!("\nStep 3: Computing h(x)...");
    let witness_fe = cs.generate_witness();
    let h_poly = compute_h_poly(&qap, &witness_fe, cs.constraints.len(), &p);
    println!("  h(x) degree: {}", h_poly.degree());

    // Step 4: SRS (デモ用 tau = 777)
    println!("\nStep 4: Generating SRS (tau = 777, demo only)...");
    let tau = Fr::from(777u64);
    let srs = generate_srs(tau, cs.constraints.len());

    // Step 5: Prove
    println!("\nStep 5: Generating proof...");
    let u_polys = polys_to_fr_vecs(&qap.a_polys);
    let v_polys = polys_to_fr_vecs(&qap.b_polys);
    let w_polys = polys_to_fr_vecs(&qap.c_polys);
    let witness: Vec<Fr> = witness_fe.iter().map(|w| field_element_to_fr(w)).collect();
    let h_coeffs = polynomial_to_fr_vec(&h_poly);
    let proof = prove_simple(&u_polys, &v_polys, &w_polys, &witness, &h_coeffs, &srs);

    // Step 6: Verify
    println!("\nStep 6: Verifying proof...");
    if verify_simple(&proof) {
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
    let zero = FieldElement::new(BigInt::from(0), p.clone());
    let one = FieldElement::new(BigInt::from(1), p.clone());

    // A(x), B(x), C(x) = sum_i witness[i] * poly_i(x)
    let mut a = Polynomial::new(vec![zero.clone()]);
    let mut b = Polynomial::new(vec![zero.clone()]);
    let mut c = Polynomial::new(vec![zero.clone()]);
    for (i, w) in witness.iter().enumerate() {
        a = &a + &qap.a_polys[i].scale(w.clone());
        b = &b + &qap.b_polys[i].scale(w.clone());
        c = &c + &qap.c_polys[i].scale(w.clone());
    }

    // P(x) = A(x)*B(x) - C(x)
    let p_poly = &(&a * &b) + &c.scale(&zero - &one);

    // Z(x) = (x - 0)(x - 1)...(x - (n - 1))
    let mut z_poly = Polynomial::new(vec![one.clone()]);
    for i in 0..num_constraints {
        let neg_i = &zero - &FieldElement::new(BigInt::from(i), p.clone());
        z_poly = &z_poly * &Polynomial::new(vec![neg_i, one.clone()]);
    }

    // h(x) = P(x) / Z(x)
    let (h, remainder) = p_poly.div_rem(&z_poly);
    assert!(
        remainder
            .coefficients
            .iter()
            .all(|c| c.value == BigInt::from(0)),
        "P(x) is not divisible by Z(x)"
    );
    h
}
