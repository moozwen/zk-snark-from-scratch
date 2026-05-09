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

// R1CS を QAP に変換するための構造体
#[derive(Debug, Clone)]
pub struct Qap {
    // 各変数ごとに 補間された多項式を持つ
    // index 0: 定数1 の多項式
    // index 1: 定数x の多項式...
    pub a_polys: Vec<Polynomial>, // A行列由来のリスト（index 0 は定数1用、index 1 は変数x用...）
    pub b_polys: Vec<Polynomial>, // B行列由来
    pub c_polys: Vec<Polynomial>, // C行列由来
}

impl Qap {
    // R1CS から QAP を生成するメイン関数
    pub fn from_r1cs(cs: &ConstraintSystem) -> Self {
        let num_vars = cs.next_var_index; // 変数の総数（列の数）
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

// ヘルパー関数： extract_column で取得したスパースな点データを、
// ラグランジュ補間に渡せるように「0埋めされた密なベクトル」に変換する
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

// どの行列（A / B / C）の列を抜き出すかを指定するセレクタ
#[derive(Clone, Copy)]
enum Matrix {
    A,
    B,
    C,
}

// 行列の「ある列（変数 index）」の係数をすべて抜き出すヘルパー関数
// 戻り値： [(制約番号, 係数), (制約番号, 係数), ...]
fn extract_column(
    cs: &ConstraintSystem,
    var_idx: usize,
    matrix: Matrix,
) -> Vec<(usize, FieldElement)> {
    let mut points = Vec::new();

    for (i, constraint) in cs.constraints.iter().enumerate() {
        // A, B, C どの行列（LinearCombination）を見るか？
        let lc = match matrix {
            Matrix::A => &constraint.a,
            Matrix::B => &constraint.b,
            Matrix::C => &constraint.c,
        };

        // その制約式の中に ターゲット変数の係数はあるか？
        for (var, coeff) in &lc.terms {
            if var.0 == var_idx {
                // (x座標: 制約式, y座標: 係数)
                points.push((i, coeff.clone()));
            }
        }
        // なければ 0 だが、スパース表現なのでリストに入れない（補間時に処理する）
    }
    points
}
