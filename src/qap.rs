//! R1CS から Quadratic Arithmetic Program (QAP) への変換を提供する。
//!
//! Groth16 実装の Layer 2（回路の表現）。R1CS の制約系をラグランジュ補間で
//! 多項式表現に変換することで、ペアリングベースの効率的な検証を可能にする。
//!
//! ## 主要型
//! - [`Qap`]: 各変数 i に対する `u_i(x), v_i(x), w_i(x)` の組
//!   （`a_polys`, `b_polys`, `c_polys`）
//!
//! ## 補間点
//! 制約 i 番目を `x = i` の点に対応させる（補間点列は 0, 1, ..., n-1）。

use num_bigint::BigInt;

use crate::field::FieldElement;
use crate::polynomial::Polynomial;
use crate::r1cs::ConstraintSystem;

/// R1CS から変換した Quadratic Arithmetic Program (QAP)。
///
/// 各変数 `i` について 3 本の多項式 `(a_i(x), b_i(x), c_i(x))` を保持する。
/// 制約 `j` 番目を `x = j` に対応させ、ラグランジュ補間で構成する。
///
/// 「証明者が知っている Witness `s = (s_0, s_1, ...)` で全制約が満たされる」
/// ことは、`(Σ s_i a_i(x)) · (Σ s_i b_i(x)) − (Σ s_i c_i(x))` が
/// 各補間点 `x = 0, 1, ..., n-1` で 0 になる、と言い換えられる。
/// これは `Z(x) = Π (x − j)` で割り切れることと同値で、Groth16 はこの性質を
/// ペアリングで検証する。
#[derive(Debug, Clone)]
pub struct Qap {
    /// A 行列由来の多項式列。`a_polys[i]` は変数 `i` の A 列を補間したもの。
    /// `i = 0` は定数 1（[`CS_ONE`](crate::r1cs::CS_ONE)）に対応する。
    pub a_polys: Vec<Polynomial>,
    /// B 行列由来の多項式列。インデックス規約は `a_polys` と同じ。
    pub b_polys: Vec<Polynomial>,
    /// C 行列由来の多項式列。インデックス規約は `a_polys` と同じ。
    pub c_polys: Vec<Polynomial>,
}

impl Qap {
    /// 制約系から QAP を構築する。
    ///
    /// 各 (行列, 変数) ペアの列を「制約 index → 係数」の点列とみなし、
    /// 補間点列 `0, 1, ..., num_constraints − 1` でラグランジュ補間する。
    ///
    /// 計算量は `O(num_vars · num_constraints^2)`
    /// （変数ごとに `O(num_constraints^2)` の補間を 3 行列分）。
    /// 制約系は `init_one` 済みであることが前提（法 `p` を取り出すため
    /// `assignments[0]` を参照する）。
    pub fn from_r1cs(cs: &ConstraintSystem) -> Self {
        let num_vars = cs.next_var_index;
        let num_constraints = cs.constraints.len();
        let p = cs
            .assignments
            .first()
            .expect("CS未初期化")
            .as_ref()
            .unwrap()
            .p
            .clone();

        // 指定行列の各変数列を Lagrange 補間で多項式化する
        let interpolate_column = |matrix: Matrix| -> Vec<Polynomial> {
            (0..num_vars)
                .map(|i| {
                    let points = extract_column(cs, i, matrix);
                    let dense = to_dense_vector(points, num_constraints, &p);
                    Polynomial::lagrange_interpolation(&dense)
                })
                .collect()
        };

        Qap {
            a_polys: interpolate_column(Matrix::A),
            b_polys: interpolate_column(Matrix::B),
            c_polys: interpolate_column(Matrix::C),
        }
    }
}

/// スパースな点列 `[(row, value), ...]` を、長さ `num_constraints` の
/// 密ベクトルに展開する（欠けた行は 0 で埋める）。
///
/// `Polynomial::lagrange_interpolation` が補間点 `x = 0, 1, ..., n-1` の
/// `y` 値列を要求するための前処理。
fn to_dense_vector(
    sparse_points: Vec<(usize, FieldElement)>,
    num_constraints: usize,
    p: &BigInt,
) -> Vec<FieldElement> {
    let zero = FieldElement::new(0, p.clone());
    let mut dense = vec![zero; num_constraints];

    for (row_idx, val) in sparse_points {
        if row_idx < num_constraints {
            dense[row_idx] = val;
        }
    }
    dense
}

/// `extract_column` が見る行列を指定するセレクタ。
#[derive(Clone, Copy)]
enum Matrix {
    A,
    B,
    C,
}

/// 指定行列の指定変数列に出てくる係数を、`(制約 index, 係数)` のスパース列で返す。
///
/// 同一制約内に同じ変数が複数項として登録されている場合、それぞれ別エントリで返す
/// （`LinearCombination::add_term` がマージしない仕様に対応）。
fn extract_column(
    cs: &ConstraintSystem,
    var_idx: usize,
    matrix: Matrix,
) -> Vec<(usize, FieldElement)> {
    let mut points = Vec::new();

    for (i, constraint) in cs.constraints.iter().enumerate() {
        let lc = match matrix {
            Matrix::A => &constraint.a,
            Matrix::B => &constraint.b,
            Matrix::C => &constraint.c,
        };

        for (var, coeff) in &lc.terms {
            if var.0 == var_idx {
                // (x座標: 制約式, y座標: 係数)
                points.push((i, coeff.clone()));
            }
        }
        // エントリなければ 0 だが、スパース表現として詰めない（to_dense_vector で 0 埋め）
    }
    points
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::r1cs::{ConstraintSystem, CS_ONE};

    const P: i64 = 7;

    fn fe(v: i64) -> FieldElement {
        FieldElement::new(v, P)
    }

    // 1 制約の最小回路: x * x = y, x = 3
    fn build_x_squared_cs() -> ConstraintSystem {
        let mut cs = ConstraintSystem::new();
        cs.init_one(fe(1));
        let x = cs.alloc_variable();
        cs.assign(x, fe(3));
        let _y = cs.mul(x, x);
        cs
    }

    #[test]
    fn from_r1cs_single_constraint_recovers_coefficients() {
        // 制約 0: (x) * (x) = (y)
        // -> A 列: CS_ONE=0, x=1, y=0
        // -> B 列: CS_ONE=0, x=1, y=0
        // -> C 列: CS_ONE=0, x=0, y=1
        let cs = build_x_squared_cs();
        let qap = Qap::from_r1cs(&cs);
        let pt = fe(0);

        assert_eq!(qap.a_polys[0].evaluate(&pt), fe(0));
        assert_eq!(qap.a_polys[1].evaluate(&pt), fe(1));
        assert_eq!(qap.a_polys[2].evaluate(&pt), fe(0));

        assert_eq!(qap.b_polys[0].evaluate(&pt), fe(0));
        assert_eq!(qap.b_polys[1].evaluate(&pt), fe(1));
        assert_eq!(qap.b_polys[2].evaluate(&pt), fe(0));

        assert_eq!(qap.c_polys[0].evaluate(&pt), fe(0));
        assert_eq!(qap.c_polys[1].evaluate(&pt), fe(0));
        assert_eq!(qap.c_polys[2].evaluate(&pt), fe(1));
    }

    #[test]
    fn from_r1cs_two_constraints_recovers_coefficients_at_each_point() {
        // 制約 0: (x) * (x) = (v1)   → A[x]=1, B[x]=1, C[v1]=1
        // 制約 1: (v1) * (x) = (v2)  → A[v1]=1, B[x]=1, C[v2]=1
        let mut cs = ConstraintSystem::new();
        cs.init_one(fe(1));
        let x = cs.alloc_variable();
        cs.assign(x, fe(2));
        let v1 = cs.mul(x, x);
        let _v2 = cs.mul(v1, x);

        let qap = Qap::from_r1cs(&cs);
        // num_vars = CS_ONE + x + v1 + v2 = 4
        assert_eq!(qap.a_polys.len(), 4);

        let p0 = fe(0);
        let p1 = fe(1);

        // 制約 0 側
        assert_eq!(qap.a_polys[1].evaluate(&p0), fe(1)); // x
        assert_eq!(qap.b_polys[1].evaluate(&p0), fe(1)); // x
        assert_eq!(qap.c_polys[2].evaluate(&p0), fe(1)); // v1

        // 制約 1 側
        assert_eq!(qap.a_polys[2].evaluate(&p1), fe(1)); // v1
        assert_eq!(qap.b_polys[1].evaluate(&p1), fe(1)); // x
        assert_eq!(qap.c_polys[3].evaluate(&p1), fe(1)); // v2

        // 反対側の点では 0 になっていること（スパース列の埋め）
        assert_eq!(qap.a_polys[1].evaluate(&p1), fe(0)); // x は制約 1 の A に出ない
        assert_eq!(qap.c_polys[2].evaluate(&p1), fe(0)); // v1 は制約 1 の C に出ない
    }

    #[test]
    fn from_r1cs_sparse_columns_pad_with_zero() {
        // 制約 0: (x + 2·1) * 1 = z      ← add_const → A 側に CS_ONE が出る
        // 制約 1: (z) * (x) = w          ← mul       ← A 側に CS_ONE は出ない
        // CS_ONE 列 A は [2, 0] を補間するはず
        let mut cs = ConstraintSystem::new();
        cs.init_one(fe(1));
        let x = cs.alloc_variable();
        cs.assign(x, fe(3));
        let z = cs.add_const(x, fe(2));
        let _w = cs.mul(z, x);

        let qap = Qap::from_r1cs(&cs);
        assert_eq!(qap.a_polys[CS_ONE.0].evaluate(&fe(0)), fe(2));
        assert_eq!(qap.a_polys[CS_ONE.0].evaluate(&fe(1)), fe(0));
    }
}