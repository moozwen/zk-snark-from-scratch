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

/// 制約系内の変数を識別するインデックス。
///
/// `Variable(0)` は定数 1 に予約済み（[`CS_ONE`]）。通常の変数は
/// [`ConstraintSystem::alloc_variable`] で発行される。
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct Variable(pub usize);

/// 定数 1 を表す予約変数。`assignments[0]` に値 1 が入っていることが前提。
///
/// 制約系を作った直後に [`ConstraintSystem::init_one`] を呼んで初期化する。
pub const CS_ONE: Variable = Variable(0);

/// 変数の線形結合 `Σ c_i · x_i` を表す。
///
/// 例： `3x + 2y + 5` は
/// `[(Variable(1), 3), (Variable(2), 2), (CS_ONE, 5)]` という並び。
///
/// 同一変数の重複は許され、自動マージは行わない。
#[derive(Clone, Debug)]
pub struct LinearCombination {
    // (変数のインデックス, その係数) のリスト
    // 3x + 2y + 5 は [(Variable(1), 3), (Variable(2), 2), (Variable(0), 5)] となる
    pub terms: Vec<(Variable, FieldElement)>,
}

impl LinearCombination {
    /// 空の線形結合を生成する。
    pub fn new() -> Self {
        Self { terms: Vec::new() }
    }

    /// 項 `coeff · var` を末尾に追加する。
    ///
    /// 既存の同変数項とはマージせず、別エントリとして保持する。
    pub fn add_term(&mut self, var: Variable, coeff: FieldElement) {
        self.terms.push((var, coeff));
    }
}

impl Default for LinearCombination {
    fn default() -> Self {
        Self::new()
    }
}

/// 単一の R1CS 制約 `A · B = C` を表す。
///
/// `A`, `B`, `C` はそれぞれ Witness ベクトルとの内積によりスカラー値となり、
/// 「左辺どうしの積が右辺と一致するか」が証明対象になる。
#[derive(Clone, Debug)]
pub struct Constraint {
    pub a: LinearCombination,
    pub b: LinearCombination,
    pub c: LinearCombination,
}

/// 算術回路全体を保持する制約系。
///
/// 制約のリストと、各変数の現在値（Witness 候補）を持つ。
/// `mul` / `add` / `add_const` を使う前に [`init_one`](Self::init_one) を呼んで
/// [`CS_ONE`] を初期化する必要がある（内部で係数 1 を作るときに参照するため）。
///
/// ## 変数レイアウト
///
/// 変数は常に `[CS_ONE, 公開入力, ..., 中間/秘密変数, ...]` の順に並ぶ。
/// 先頭 `num_public_variables` 個（[`CS_ONE`] を含む）が public で、
/// この境界 l+1 が Groth16 の検証等式で公開入力 `a_0, ..., a_l` を切り出すのに使われる。
/// 公開入力は [`alloc_public_input`](Self::alloc_public_input) で、秘密/中間変数は
/// [`alloc_variable`](Self::alloc_variable) で発行する。
/// ｐublic は private より前に確保しなければならない（前方に固める Groth16 の慣習に従う）。
///
/// # 例
///
/// ```text
/// let mut cs = ConstraintSystem::new();
/// cs.init_one(FieldElement::new(1, 7));
/// let x = cs.alloc_variable();
/// cs.assign(x, FieldElement::new(3, 7));
/// let y = cs.mul(x, x); // y = x^2
/// ```
pub struct ConstraintSystem {
    pub next_var_index: usize,
    pub constraints: Vec<Constraint>,
    // 各変数の値を保持するリスト
    pub assignments: Vec<Option<FieldElement>>,
    // 先頭から数えた public 変数の数（CS_ONE 含む = l+1）。
    // 不変条件: public 変数は常にインデックス 0..num_public_variables
    pub num_public_variables: usize,
}

impl ConstraintSystem {
    /// 空の制約系を生成する。
    ///
    /// 変数も制約もまだ存在しない状態。最初に [`init_one`](Self::init_one) を
    /// 呼んで [`CS_ONE`] を確保してから使う。
    pub fn new() -> Self {
        Self {
            next_var_index: 0,
            constraints: Vec::new(),
            // alloc 時に None を埋める方針
            assignments: Vec::new(),
            // init_one で CS_ONE を public として 1 に設定する
            num_public_variables: 0,
        }
    }

    /// 変数 `var` に値 `value` を代入する。
    ///
    /// `var` が [`alloc_variable`](Self::alloc_variable) 未発行のときは panic する。
    pub fn assign(&mut self, var: Variable, value: FieldElement) {
        if var.0 < self.assignments.len() {
            self.assignments[var.0] = Some(value);
        } else {
            panic!("variable {} is out of bounds; alloc it first", var.0);
        }
    }

    /// 定数 1 を保持する [`CS_ONE`] を初期化する。
    ///
    /// 内部で `Variable(0)` を確保して `one` を代入する。制約系を作った直後、
    /// 他の `alloc_variable` を呼ぶ前に一度だけ呼び出すこと。
    /// `FieldElement` から法 `p` を取得するため、外部から渡してもらう設計。
    pub fn init_one(&mut self, one: FieldElement) {
        // Index 0 がまだなければ作る
        if self.assignments.is_empty() {
            self.alloc_variable(); // Index 0 を確保
        }
        self.assign(CS_ONE, one);
        // CS_ONE (a_0 = 1) は常に public。これが public 領域の起点になる。
        self.num_public_variables = 1;
    }

    /// 全変数の現在値を Witness ベクトルとして取り出す。
    ///
    /// 未代入の変数（`None`）が残っていれば panic する。
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

    /// 新しい変数を発行し、その [`Variable`] ハンドルを返す。
    ///
    /// 値は未代入（`None`）状態で確保される。`assign` で値を入れる必要がある。
    pub fn alloc_variable(&mut self) -> Variable {
        let var = Variable(self.next_var_index);
        self.next_var_index += 1;
        self.assignments.push(None);
        var
    }

    /// 公開入力変数を 1 つ発行し、その [`Variable`] ハンドルを返す。
    ///
    /// 公開変数はインデックス前方（[`CS_ONE`] の直後）に固める必要があるため、
    /// [`alloc_variable`](Self::alloc_variable) で秘密/中間変数を 1 つでも確保した後に
    /// 呼ぶと panic する。
    /// [`init_one`](Self::init_one) 済みであることも前提とする。
    ///
    /// 値は未代入（`None`）状態で確保される。`assign` で値を入れる必要がある。
    pub fn alloc_public_input(&mut self) -> Variable {
        assert!(
            self.num_public_variables >= 1,
            "call init_one() before alloc_public_input()"
        );
        // public 領域は 0..num_public_variables。
        // private を alloc 済みだと next_var_index がこれを追い越すので、前方固めが崩れる。
        assert_eq!(
            self.next_var_index, self.num_public_variables,
            "alloc_public_input() must be called before any private alloc_variable()"
        );

        let var = self.alloc_variable();
        self.num_public_variables += 1;
        var
    }

    /// 制約 `A · B = C` を制約系に直接追加する。
    ///
    /// 通常は `mul` / `add` / `add_const` 経由で間接的に呼ばれる。
    pub fn enforce(&mut self, a: LinearCombination, b: LinearCombination, c: LinearCombination) {
        self.constraints.push(Constraint { a, b, c });
    }

    /// 掛け算ゲートを追加する。
    ///
    /// 新変数 `c` を確保して `c = a * b` を計算し、制約 `(a) · (b) = (c)` を追加する。
    /// 戻り値は `c`。
    pub fn mul(&mut self, a: Variable, b: Variable) -> Variable {
        let c = self.alloc_variable();

        // 値の計算（Witness 生成）
        let val_a = self.assignments[a.0]
            .as_ref()
            .expect("variable a is unassigned");
        let val_b = self.assignments[b.0]
            .as_ref()
            .expect("variable b is unassigned");
        let val_c = val_a * val_b;
        self.assign(c, val_c);

        // 制約: (a) * (b) = (c)
        let mut lc_a = LinearCombination::new();
        lc_a.add_term(a, self.one());
        let mut lc_b = LinearCombination::new();
        lc_b.add_term(b, self.one());
        let mut lc_c = LinearCombination::new();
        lc_c.add_term(c, self.one());

        self.enforce(lc_a, lc_b, lc_c);

        c
    }

    /// 法 `p` のもとでの `FieldElement` 1 を返す。
    ///
    /// `assignments[0]` ([`CS_ONE`]) から法を取り出すため、`init_one` 済み前提。
    fn one(&self) -> FieldElement {
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

    /// 足し算ゲートを追加する。
    ///
    /// 新変数 `c` を確保して `c = a + b` を計算し、
    /// 制約 `(a + b) · 1 = (c)` を追加する。戻り値は `c`。
    ///
    /// 現在は unit test からのみ呼ばれる。
    /// Phase 5 以降の回路で使われ始めたら attribute を外す。
    #[allow(dead_code)]
    pub fn add(&mut self, a: Variable, b: Variable) -> Variable {
        let c = self.alloc_variable();

        // 値の計算
        let val_a = self.assignments[a.0]
            .as_ref()
            .expect("variable a is unassigned");
        let val_b = self.assignments[b.0]
            .as_ref()
            .expect("variable b is unassigned");
        self.assign(c, val_a + val_b);

        // 制約： (a + b) * 1 = c
        let mut lc_a = LinearCombination::new();
        lc_a.add_term(a, self.one());
        lc_a.add_term(b, self.one());

        let mut lc_b = LinearCombination::new();
        lc_b.add_term(CS_ONE, self.one());

        let mut lc_c = LinearCombination::new();
        lc_c.add_term(c, self.one());

        self.enforce(lc_a, lc_b, lc_c);

        c
    }

    /// 定数加算ゲートを追加する。
    ///
    /// 新変数 `c` を確保して `c = a + k` を計算し、
    /// 制約 `(a + k · 1) · 1 = (c)` を追加する。戻り値は `c`。
    pub fn add_const(&mut self, a: Variable, constant: FieldElement) -> Variable {
        let c = self.alloc_variable();

        // 値の計算
        let val_a = self.assignments[a.0]
            .as_ref()
            .expect("variable a is unassigned");
        self.assign(c, val_a + &constant);

        // 制約： (a + 1 * k) * 1 = c
        let mut lc_a = LinearCombination::new();
        lc_a.add_term(a, self.one());
        lc_a.add_term(CS_ONE, constant);

        let mut lc_b = LinearCombination::new();
        lc_b.add_term(CS_ONE, self.one());

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

    #[test]
    fn alloc_public_input_increments_public_count() {
        let mut cs = ConstraintSystem::new();
        cs.init_one(fe(1));

        // init_one 時点では CS_ONE のみが public
        assert_eq!(cs.num_public_variables, 1);

        let p0 = cs.alloc_public_input();
        let p1 = cs.alloc_public_input();

        // CS_ONE の直後（index 1, 2）に公開入力が固めて配置されること
        assert_eq!(p0, Variable(1));
        assert_eq!(p1, Variable(2));
        assert_eq!(cs.num_public_variables, 3); // CS_ONE + 2

        // 以降の秘密/中間変数は public 領域を侵さない
        let _priv = cs.alloc_variable();
        assert_eq!(cs.num_public_variables, 3);
        assert_eq!(cs.next_var_index, 4);
    }

    #[test]
    fn circuit_without_public_inputs_has_one_public_var() {
        // public を一切使わない回路では CS_ONE だけが public (l = 0)
        let mut cs = ConstraintSystem::new();
        cs.init_one(fe(1));
        let a = cs.alloc_variable();
        let b = cs.alloc_variable();
        cs.assign(a, fe(2));
        cs.assign(b, fe(3));
        let _c = cs.mul(a, b);
        assert_eq!(cs.num_public_variables, 1);
    }

    #[test]
    #[should_panic(expected = "before any private")]
    fn alloc_public_iput_after_private_panics() {
        let mut cs = ConstraintSystem::new();
        cs.init_one(fe(1));
        let _priv = cs.alloc_variable();
        cs.alloc_public_input();
    }

    #[test]
    #[should_panic(expected = "init_one")]
    fn alloc_public_input_before_init_one_panics() {
        let mut cs = ConstraintSystem::new();
        cs.alloc_public_input();
    }
}
