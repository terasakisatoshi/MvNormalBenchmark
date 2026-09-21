# Rust 実装の乱数生成仕様

## 一様乱数生成器

`StandardRng` は依存クレートを使わない `xorshift64` 型の疑似乱数生成器です。

状態 `s` は次の順で更新します。

```text
s = s xor (s << 13)
s = s xor (s >> 7)
s = s xor (s << 17)
```

更新後の値から上位53ビットを取り出し、半単位を加えて `(0, 1)` の倍精度一様乱数に変換します。

```text
u = (f64(next_u64() >> 11) + 0.5) / 2^53
```

シード 0 は xorshift の固定点になるため、非ゼロの初期状態へ置き換えます。
この RNG はベンチマーク用であり、暗号用途には適しません。
共通ベンチマークのサンプル生成シードは `0x5EED2021` です。

## 標準正規乱数

`NormalAlgorithm::MarsagliaPolar` は Marsaglia の polar 法で2個ずつ生成します。

1. `u` と `v` を区間 `(-1, 1)` から生成する。
2. `s = u^2 + v^2` を計算する。
3. `s <= 0` または `s >= 1` の場合は棄却して再試行する。
4. `a = sqrt(-2 log(s) / s)` を計算する。
5. `u*a` と `v*a` を標準正規乱数として返す。

三角関数を呼ばないため、元の Box--Muller 直接法よりスカラー生成の負荷を抑えます。

`fill_standard_normals` は受理したペアを直接スライスへ書き込みます。
奇数長のスライスを `standard_normal` で埋める場合だけ、余剰の1値を `Option<f64>` に保持します。

`NormalAlgorithm::Ziggurat` は256層のZiggurat法を使います。
層境界から `k`、`w`、`f` のテーブルを一度だけ生成し、矩形内部の候補を整数比較だけで高速に受理します。
矩形外の候補は密度比較で再判定し、層0では指数分布を使って正規分布の裾を生成します。

```rust
let mut rng = StandardRng::with_algorithm(42, NormalAlgorithm::Ziggurat);
```

## MvNormal の変換

Cholesky 因子は Rust の連続配列に行優先で格納します。

```text
cholesky[row * dimension + column]
```

`sample_inplace` は標準正規ベクトルを `output` に生成した後、行を大きい番号から処理します。
行 `row` の計算前には `output[..=row]` がまだ必要な標準正規値なので、追加の作業ベクトルなしで
`output[row] = mean[row] + Σ cholesky[row,column] * output[column]` を計算します。
行優先レイアウトに合わせて、各行の係数を連続領域から走査します。

`sample_into` は `scratch` と `output` を分ける互換APIとして残しています。

## rand_distr 経路

`--normal rand-distr` は `rand_distr` クレートの `StandardNormal` を使います。
`StandardNormal` は ZIGNOR 版の Ziggurat 法で実装されています。

```rust
use rand::SeedableRng;
use rand::rngs::StdRng;
use rand_distr::{Distribution, StandardNormal};

let mut rng = StdRng::seed_from_u64(0x5EED2021);
let normal = StandardNormal;
distribution.sample_inplace_with(|| normal.sample(&mut rng), &mut output)?;
```

`MvNormal::sample_inplace_with` は外部の標準正規生成器をクロージャで受け取り、
`mean + L z` の変換だけを共通化します。`StandardNormal` は `StdRng` 独自の乱数
列を使うため、共通 xorshift64 を使う `polar`/`ziggurat` とは乱数列が異なります。

## statrs 経路

`--normal statrs` は `statrs` クレートの `MultivariateNormal` を使います。
`statrs::distribution::MultivariateNormal::new(mean, cov)` は共分散の対称性と
正定値性を検証して Cholesky 因子を保持し、`Distribution` 実装の `sample` は
`L z + μ` を計算します。`Distributions.jl` の `MvNormal` に最も近い crate
実装です。

```rust
use statrs::distribution::MultivariateNormal;

let distribution = MultivariateNormal::new(mean, covariance_flat)?;
let sample = distribution.sample(&mut rng);
```

平均は `Vec<f64>`、共分散は行優先に平坦化した `Vec<f64>` を渡します。共分散は
対称なので列優先の nalgebra 内部表現でも同じ並びになります。`z` は `StdRng` を
乱源とする `Normal(0, 1)` から生成されます。

## 再現性

同じシードと同じビルド条件では、Rust 実装内で同じ乱数列を再生成できます。
`polar`/`ziggurat` は共通 xorshift64 を使うため、乱数消費順序を揃えれば他言語と
checksum が一致します。`rand-distr` は独自の乱数列を使うため一致しません。
他言語とは RNG の状態更新、棄却回数、浮動小数点ライブラリが異なるため、
チェックサム一致は要求しません。
