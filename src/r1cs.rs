//! Rank-1 Constraint System (R1CS) の表現と回路構築 API を提供する。
//!
//! Groth16 実装の Layer 2（回路の表現）。算術回路を制約系
//! `A·B = C` の形式で表現し、Witness を生成する。
//!
//! ## 主要型
//! - [`ConstraintSystem`][]: 制約と変数代入を保持する回路全体
//! - [`Variable`][]: 変数（インデックス）。[`CS_ONE`] は定数 1 を表す予約変数
//! - [`LinearCombination`][]: 変数の線形結合
//! - [`Constraint`][]: 単一の `A·B = C` 制約
//!
//! ## 回路構築 API
//! - [`ConstraintSystem::mul`][]: 掛け算ゲート
//! - [`ConstraintSystem::add`][]: 足し算ゲート
//! - [`ConstraintSystem::add_const`][]: 定数加算ゲート

use crate::field::FieldElement;

// R1CS において変数は「インデックス」
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct Variable(pub usize);

// インデックス 0 は常に 1 を表す
pub const CS_ONE: Variable = Variable(0);

#[derive(Clone, Debug)]
pub struct LinearCombination {
    // (変数のインデックス, その係数) のリスト
    // 3x + 2y + 5 は [(Variable(1), 3), (Variable(2), 2), (Variable(0), 5)] となる
    pub terms: Vec<(Variable, FieldElement)>,
}

impl LinearCombination {
    pub fn new() -> Self {
        Self { terms: Vec::new() }
    }

    pub fn add_term(&mut self, var: Variable, coeff: FieldElement) {
        self.terms.push((var, coeff));
    }
}

impl Default for LinearCombination {
    fn default() -> Self {
        Self::new()
    }
}

#[derive(Clone, Debug)]
pub struct Constraint {
    pub a: LinearCombination,
    pub b: LinearCombination,
    pub c: LinearCombination,
}

pub struct ConstraintSystem {
    pub next_var_index: usize,
    pub constraints: Vec<Constraint>,
    // 各変数の値を保持するリスト
    pub assignments: Vec<Option<FieldElement>>,
}

impl ConstraintSystem {
    pub fn new() -> Self {
        // インデックス 0 は定数 1 のために予約済みなので 1 から開始
        Self {
            next_var_index: 0,
            constraints: Vec::new(),
            // alloc 時に None を埋める方針
            assignments: Vec::new(),
        }
    }

    pub fn assign(&mut self, var: Variable, value: FieldElement) {
        if var.0 < self.assignments.len() {
            self.assignments[var.0] = Some(value);
        } else {
            panic!("variable {} is out of bounds; alloc it first", var.0);
        }
    }

    // 定数1（Index 0）を初期化するための専用メソッド
    // ※ p（素数）が必要なので、外部から呼んでもらう
    pub fn init_one(&mut self, one: FieldElement) {
        // Index 0 がまだなければ作る
        if self.assignments.is_empty() {
            self.alloc_variable(); // Index 0 を確保
        }
        self.assign(CS_ONE, one);
    }

    // 記録された値から Witness ベクトルを生成する
    pub fn generate_witness(&self) -> Vec<FieldElement> {
        self.assignments
            .iter()
            .map(|val| {
                val.as_ref()
                    .expect("witness contains an unassigned variable")
                    .clone()
            })
            .collect()
    }

    // 新しい変数を「発行」する
    pub fn alloc_variable(&mut self) -> Variable {
        let var = Variable(self.next_var_index);
        self.next_var_index += 1;
        self.assignments.push(None);
        var
    }

    // 回路に新しい制約（A * B = C）を追加する
    pub fn enforce(&mut self, a: LinearCombination, b: LinearCombination, c: LinearCombination) {
        self.constraints.push(Constraint { a, b, c });
    }

    // 掛け算ゲート（a * b = c）を作成し、結果の変数 c を返す
    pub fn mul(&mut self, a: Variable, b: Variable) -> Variable {
        // 1. 結果用の変数 c を確保
        let c = self.alloc_variable();

        // 2. 値の計算（Witness 生成）
        // a と b の値を読み出して掛け算し、c に代入する
        let val_a = self.assignments[a.0]
            .as_ref()
            .expect("variable a is unassigned");
        let val_b = self.assignments[b.0]
            .as_ref()
            .expect("variable b is unassigned");
        let val_c = val_a * val_b;
        self.assign(c, val_c);

        // 3. 制約の追加（a * b = c）
        let mut lc_a = LinearCombination::new();
        lc_a.add_term(a, self.one());
        let mut lc_b = LinearCombination::new();
        lc_b.add_term(b, self.one());
        let mut lc_c = LinearCombination::new();
        lc_c.add_term(c, self.one());

        self.enforce(lc_a, lc_b, lc_c);

        c
    }

    // ヘルパー関数： 係数 1 のFieldElement を返す
    fn one(&self) -> FieldElement {
        // assignments[0] (CS_ONE) から p を取得して 1 を作る
        let p = self
            .assignments
            .first()
            .expect("constraint system not initialized; call init_one() first")
            .as_ref()
            .expect("CS_ONE is unassigned")
            .p
            .clone();
        FieldElement::new(1, p)
    }

    // 足し算ゲート： (a + b) * 1 = c
    pub fn add(&mut self, a: Variable, b: Variable) -> Variable {
        // 1. 結果用変数 c を確保
        let c = self.alloc_variable();

        // 2. 値の計算（Witness 生成）
        let val_a = self.assignments[a.0]
            .as_ref()
            .expect("variable a is unassigned");
        let val_b = self.assignments[b.0]
            .as_ref()
            .expect("variable b is unassigned");
        self.assign(c, val_a + val_b);

        // 3. 制約： (a + b) * 1 = c
        // A: a + b
        let mut lc_a = LinearCombination::new();
        lc_a.add_term(a, self.one());
        lc_a.add_term(b, self.one());

        // B: 1
        let mut lc_b = LinearCombination::new();
        lc_b.add_term(CS_ONE, self.one());

        // C: c
        let mut lc_c = LinearCombination::new();
        lc_c.add_term(c, self.one());

        self.enforce(lc_a, lc_b, lc_c);

        c
    }

    // 定数の足し算： (a + const) * 1 = c
    pub fn add_const(&mut self, a: Variable, constant: FieldElement) -> Variable {
        let c = self.alloc_variable();

        // 2. 値の計算
        let val_a = self.assignments[a.0]
            .as_ref()
            .expect("variable a is unassigned");
        self.assign(c, val_a + &constant); // 定数を足す

        // 3. 制約： (a + (1 * const)) * 1 = c
        // A: a * 1 + 1 * const
        let mut lc_a = LinearCombination::new();
        lc_a.add_term(a, self.one());
        lc_a.add_term(CS_ONE, constant);

        // B: 1
        let mut lc_b = LinearCombination::new();
        lc_b.add_term(CS_ONE, self.one());

        // C: c
        let mut lc_c = LinearCombination::new();
        lc_c.add_term(c, self.one());

        self.enforce(lc_a, lc_b, lc_c);

        c
    }
}

impl Default for ConstraintSystem {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const P: i64 = 7;

    fn fe(v: i64) -> FieldElement {
        FieldElement::new(v, P)
    }

    #[test]
    fn alloc_variable_assigns_sequential_indices() {
        let mut cs = ConstraintSystem::new();
        let v0 = cs.alloc_variable();
        let v1 = cs.alloc_variable();
        let v2 = cs.alloc_variable();
        assert_eq!(v0, Variable(0));
        assert_eq!(v1, Variable(1));
        assert_eq!(v2, Variable(2));
        assert_eq!(cs.next_var_index, 3);
        assert_eq!(cs.assignments.len(), 3);
        assert!(cs.assignments.iter().all(|a| a.is_none()));
    }

    #[test]
    fn init_one_sets_cs_one_to_one() {
        let mut cs = ConstraintSystem::new();
        cs.init_one(fe(1));
        assert_eq!(cs.assignments[CS_ONE.0], Some(fe(1)));
        assert_eq!(cs.next_var_index, 1);
    }

    #[test]
    #[should_panic(expected = "out of bounds")]
    fn assign_out_of_bounds_variable_panics() {
        let mut cs = ConstraintSystem::new();
        cs.assign(Variable(5), fe(3));
    }

    #[test]
    fn generate_witness_returns_assigned_values() {
        let mut cs = ConstraintSystem::new();
        cs.init_one(fe(1));
        let a = cs.alloc_variable();
        let b = cs.alloc_variable();
        cs.assign(a, fe(2));
        cs.assign(b, fe(3));
        assert_eq!(cs.generate_witness(), vec![fe(1), fe(2), fe(3)]);
    }

    #[test]
    #[should_panic(expected = "unassigned")]
    fn generate_witness_panics_on_unassigned() {
        let mut cs = ConstraintSystem::new();
        cs.init_one(fe(1));
        let _ = cs.alloc_variable(); // 未 assign のまま
        cs.generate_witness();
    }

    #[test]
    fn mul_computes_value_and_adds_constraint() {
        let mut cs = ConstraintSystem::new();
        cs.init_one(fe(1));
        let a = cs.alloc_variable();
        let b = cs.alloc_variable();
        cs.assign(a, fe(2));
        cs.assign(b, fe(3));

        let c = cs.mul(a, b);

        // 2 * 3 ≡ 6 (mod 7)
        assert_eq!(cs.assignments[c.0], Some(fe(6)));
        assert_eq!(cs.constraints.len(), 1);

        // 制約形: (a) * (b) = (c)
        let con = &cs.constraints[0];
        assert_eq!(con.a.terms, vec![(a, fe(1))]);
        assert_eq!(con.b.terms, vec![(b, fe(1))]);
        assert_eq!(con.c.terms, vec![(c, fe(1))]);
    }

    #[test]
    fn add_computes_value_and_adds_constraint() {
        let mut cs = ConstraintSystem::new();
        cs.init_one(fe(1));
        let a = cs.alloc_variable();
        let b = cs.alloc_variable();
        cs.assign(a, fe(5));
        cs.assign(b, fe(4));

        let c = cs.add(a, b);

        // 5 + 4 = 9 ≡ 2 (mod 7)
        assert_eq!(cs.assignments[c.0], Some(fe(2)));
        assert_eq!(cs.constraints.len(), 1);

        // 制約形: (a + b) * 1 = c
        let con = &cs.constraints[0];
        assert_eq!(con.a.terms, vec![(a, fe(1)), (b, fe(1))]);
        assert_eq!(con.b.terms, vec![(CS_ONE, fe(1))]);
        assert_eq!(con.c.terms, vec![(c, fe(1))]);
    }

    #[test]
    fn add_const_computes_value_and_adds_constraint() {
        let mut cs = ConstraintSystem::new();
        cs.init_one(fe(1));
        let a = cs.alloc_variable();
        cs.assign(a, fe(3));

        let c = cs.add_const(a, fe(5));

        // 3 + 5 = 8 ≡ 1 (mod 7)
        assert_eq!(cs.assignments[c.0], Some(fe(1)));
        assert_eq!(cs.constraints.len(), 1);

        // 制約形: (a + 1·k) * 1 = c
        let con = &cs.constraints[0];
        assert_eq!(con.a.terms, vec![(a, fe(1)), (CS_ONE, fe(5))]);
        assert_eq!(con.b.terms, vec![(CS_ONE, fe(1))]);
        assert_eq!(con.c.terms, vec![(c, fe(1))]);
    }

    #[test]
    fn linear_combination_add_term_allows_duplicates() {
        let mut lc = LinearCombination::new();
        lc.add_term(Variable(1), fe(2));
        lc.add_term(Variable(1), fe(3));
        // マージせずに 2 件のまま保持される
        assert_eq!(lc.terms.len(), 2);
        assert_eq!(lc.terms[0], (Variable(1), fe(2)));
        assert_eq!(lc.terms[1], (Variable(1), fe(3)));
    }
}