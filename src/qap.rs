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

        let mut a_polys = Vec::new();
        let mut b_polys = Vec::new();
        let mut c_polys = Vec::new();

        // すべての変数（列）についてループ
        for i in 0..num_vars {
            // 1. Matrix A の i番目の列を抜き出して多項式化
            let points_a = extract_column(cs, i, 'A');

            // y座標だけのリストにする（x座標は 0,1,2... と決まっているため、interpolation側で処理される想定）
            // ※ lagrange_interpolation の実装に合わせて、(x,y) を渡すか y だけ渡すか確認してください。
            //   前回の実装では `y_values: &Vec<FieldElement>` (yだけ) でしたね。
            //   ただし、extract_column はスパース（0を飛ばす）なデータを返すので、
            //   ここで「密なベクトル（0埋め）」に変換する必要があります。
            let dense_points_a = to_dense_vector(points_a, cs.constraints.len(), cs);
            a_polys.push(Polynomial::lagrange_interpolation(&dense_points_a));

            // 2. Matrix B
            let points_b = extract_column(cs, i, 'B');
            let dense_points_b = to_dense_vector(points_b, cs.constraints.len(), cs);
            b_polys.push(Polynomial::lagrange_interpolation(&dense_points_b));

            // 3. Matrix C
            let points_c = extract_column(cs, i, 'C');
            let dense_points_c = to_dense_vector(points_c, cs.constraints.len(), cs);
            c_polys.push(Polynomial::lagrange_interpolation(&dense_points_c));
        }

        Qap {
            a_polys,
            b_polys,
            c_polys,
        }
    }
}

// ヘルパー関数： extract_column で取得したスパースな点データを、
// ラグランジュ補間に渡せるように「0埋めされた密なベクトル」に変換する
fn to_dense_vector(
    sparse_points: Vec<(usize, FieldElement)>,
    num_constraints: usize,
    cs: &ConstraintSystem, // ゼロ生成用に p を取得するために必要
) -> Vec<FieldElement> {
    let p = if !sparse_points.is_empty() {
        sparse_points[0].1.p.clone()
    } else {
        cs.assignments
            .get(0)
            .expect("CS未初期化")
            .as_ref()
            .unwrap()
            .p
            .clone()
    };

    let zero = FieldElement::new(BigInt::from(0), p.clone());
    let mut dense = vec![zero; num_constraints];

    for (row_idx, val) in sparse_points {
        if row_idx < num_constraints {
            dense[row_idx] = val;
        }
    }
    dense
}

// 行列の「ある列（変数 index）」の係数をすべて抜き出すヘルパー関数
// 戻り値： [(制約番号, 係数), (制約番号, 係数), ...]
fn extract_column(
    cs: &ConstraintSystem,
    var_idx: usize,
    matrix_selector: char, // 'A', 'B', or 'C'
) -> Vec<(usize, FieldElement)> {
    let mut points = Vec::new();

    for (i, constraint) in cs.constraints.iter().enumerate() {
        // A, B, C どの行列（LinearCombination）を見るか？
        let lc = match matrix_selector {
            'A' => &constraint.a,
            'B' => &constraint.b,
            'C' => &constraint.c,
            _ => panic!("Invalid matrix selector"),
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
