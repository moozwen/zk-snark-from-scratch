# zk-snark-from-scratch

> Rust で書かれた Groth16 の教育用スクラッチ実装。理論と実装のギャップを埋めるためのリファレンス。

![status](https://img.shields.io/badge/status-WIP-orange) ![version](https://img.shields.io/badge/version-v0.6-blue) ![license](https://img.shields.io/badge/license-MIT-green)

## これは何か / What is this

**ZKP（特に Groth16）を "コードを読みながら" 理解したいソフトウェアエンジニア向けの、Rust 製教材＋リファレンス実装です。**

### 想定読者 / Target reader

- 有限体・楕円曲線の基礎は知っている、またはこれから学びたい
- 既存 ZK ライブラリ（arkworks, snarkjs 等）を触ったが中身がブラックボックスに感じる
- 日本語での解説を求めている

### 提供するもの / What this provides

- Groth16 Prover/Verifier の Rust 実装
- R1CS → QAP 変換、SRS 生成、証明生成・検証までの一貫した流れ
- 数理背景の解説ドキュメント（`docs/`）
- 基本的なサーキット例とベンチマーク

### 提供しないもの / Non-goals

- 本番用の最適化実装 — **arkworks / bellman を使ってください**
- 複雑なサーキット DSL（Circom のような）
- Groth16 以外の証明系（Plonk, Halo2, STARKs 等）
- 新しい暗号学的貢献

---

## ステータスとロードマップ / Status & Roadmap

**v0.6 (Groth16 Verifier)** まで完了。本式 Groth16 の prover / verifier が E2E で動作する（`cargo run` で `x^3 + 5` のデモを実行可能）。次は **v0.7 (サンプル回路)**。

| Version | 内容 | 状態 |
|---|---|---|
| v0.3 | 有限体、楕円曲線、ペアリング基盤 | ✅ |
| v0.4 | R1CS → QAP 変換 | ✅ |
| v0.5 | Groth16 プルーバー（pk + ランダム化 r, s による zero-knowledge） | ✅ |
| v0.6 | Groth16 検証者（4 ペアリングの検証等式） | ✅ |
| v0.7 | サンプルサーキット 3 種（sudoku, range proof, hash preimage） | ⏳ |
| v0.8 | ベンチマーク（arkworks / snarkjs との比較） | ⏳ |
| v0.9 | ドキュメント整備 | ⏳ |
| v1.0 | 公開リリース・Zenn 書籍完結 | ⏳ |

進捗は [GitHub Issues](https://github.com/moozwen/zk-snark-from-scratch/issues) と [Zenn スクラップ](https://zenn.dev/nacekimura/scraps/6b876e07ee3f18)で追えます。

---

## クイックスタート / Quick Start

```bash
git clone https://github.com/<user>/zk-snark-from-scratch
cd zk-snark-from-scratch
cargo run
cargo test
```

### 最小サンプル（E2E: R1CS → Prove → Verify）

```rust
// 例: y = x^3 + 5 を満たす x を知っていることを、x を明かさず証明する
//   1. R1CS 構築（y を公開入力、x を秘密入力）
//   2. R1CS -> QAP -> h(x)
//   3. generate_groth16_keys で pk/vk を生成
//   4. prove(pk, ...) -> Groth16Proof
//   5. verify(vk, public_inputs, proof) -> true
```

動作する E2E デモは [`src/main.rs`](src/main.rs)（`cargo run`）を、各段の単体・結合テストは
各モジュール内のインライン `#[cfg(test)]`（`cargo test`）を参照してください。

---

## アーキテクチャ / Architecture

Groth16 を 5 つのレイヤーに分けて実装しています。

| Layer | 担当モジュール | 内容 |
|---|---|---|
| 1. 数学的基盤 | `field.rs`, `polynomial.rs` | 有限体、多項式演算 |
| 2. 回路の表現 | `r1cs.rs`, `qap.rs`, `adapter.rs` | R1CS, QAP 変換, witness 生成 |
| 3. プロトコル | `setup.rs`, `prover.rs`, `verifier.rs` | 鍵生成 (pk/vk), 証明, 検証 |
| 4. 実装的側面 | (WIP) | シリアライゼーション, 定時間演算 |
| 5. 応用 | (WIP) `examples/` | サンプル回路, Web デモ |

---

## 依存関係の方針 / A note on dependencies

教育目的のため、できる限り自作していますが、**楕円曲線とペアリングの実装は現状 [arkworks](https://github.com/arkworks-rs) (`ark-bn254`, `ark-ec`, `ark-ff`) に依存しています**。

- **理由**: BN254 上のペアリングをゼロから正しく実装するには数千行と数ヶ月かかるため、プロジェクト主題（Groth16 プロトコル本体）に集中するためのトレードオフ
- **将来**: v1.x 以降、ペアリング層も自作に置き換える可能性あり

---

## ドキュメント / Documentation

3 層構造を予定しています（執筆中）。

- **Layer 1** — この README（5 分で全体像）
- **Layer 2** — [`docs/getting_started.md`](docs/getting_started.md) (30 分で ZKP と本実装の概要)
- **Layer 3** — [`docs/theory.md`](docs/theory.md), [`docs/design.md`](docs/design.md) (数時間かけて数理と設計判断)

---

## 参考文献 / References

- Groth, J. (2016). *On the Size of Pairing-based Non-interactive Arguments*.
- [Moonmath Manual](https://leastauthority.com/community-matters/moonmath-manual/)
- [The RareSkills Book of Zero Knowledge](https://rareskills.io/zk-book)
- [arkworks](https://github.com/arkworks-rs) — reference implementation
- [Zenn 連載記事](#) (公開予定)

---

## License

MIT. See [LICENSE](LICENSE).
