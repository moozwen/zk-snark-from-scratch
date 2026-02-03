use num_bigint::BigInt;

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
            panic!("存在しない変数に代入しようとしました");
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
                val.clone()
                    .expect("未定義の変数があります！値をassignしてください")
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
            .clone()
            .expect("変数aの値が未設定です");
        let val_b = self.assignments[b.0]
            .clone()
            .expect("変数bの値が未設定です");
        let val_c = &val_a * &val_b;
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
        let p = self.assignments[0].as_ref().unwrap().p.clone();
        FieldElement::new(BigInt::from(1), p)
    }

    // 足し算ゲート： (a + b) * 1 = c
    pub fn add(&mut self, a: Variable, b: Variable) -> Variable {
        // 1. 結果用変数 c を確保
        let c = self.alloc_variable();

        // 2. 値の計算（Witness 生成）
        let val_a = self.assignments[a.0].as_ref().expect("a is missing");
        let val_b = self.assignments[b.0].as_ref().expect("b is missing");
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
        let val_a = self.assignments[a.0].as_ref().expect("a is missing");
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
