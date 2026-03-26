use ark_bn254::{Bn254, G1Projective, G2Projective};
use ark_ec::{pairing::Pairing, CurveGroup, PrimeGroup};

use crate::prover::Proof;

/// シンプル版の検証
/// [A]_1 * [B]_2 == [C]_1 * G_2 をペアリングで確認する
pub fn verify_simple(proof: &Proof) -> bool {
    let g2 = G2Projective::generator();

    // ペアリングには Affine 座標が必要
    let a_affine = proof.a_g1.into_affine();
    let b_affine = proof.b_g2.into_affine();
    let c_affine = proof.c_g1.into_affine();
    let g2_affine = g2.into_affine();

    // 左辺: e([A]_1, [B]_2)
    let lhs = Bn254::pairing(a_affine, b_affine);

    // 右辺: e([C]_1, G_2)
    let rhs = Bn254::pairing(c_affine, g2_affine);

    lhs == rhs
}
