# Rust MvNormal benchmark

Rust 2024、標準ライブラリのみで実装した多変量正規分布のサンプラーです。
共分散行列 `Σ` を構築時に Cholesky 分解して `Σ = L L'` とし、各サンプルを

```text
x = μ + L z
```

で生成します。`z` は `--normal` で選択した標準正規乱数アルゴリズムで生成します。

## ビルド

```sh
cargo build --release
```

## ベンチマーク実行

```sh
cargo run --release -- --dim 8 --samples 100 --repeats 2
cargo run --release -- --dim 8 --samples 100 --repeats 2 --normal ziggurat
```

成功時は次のCSV形式の1行だけを標準出力へ出します。

```text
rust-polar,<dim>,<samples>,<repeats>,<setup_sec>,<avg_sample_sec>,<min_sample_sec>,<checksum>
```

`--normal ziggurat` を指定した場合は、先頭列が `rust-ziggurat` になります。

`avg_sample_sec` と `min_sample_sec` は、それぞれ1回の `samples` バッチにかかった時間の平均・最小です。CLIは出力バッファだけを再利用する `sample_inplace` を使います。
`--normal polar`（既定値）と `--normal ziggurat` で標準正規乱数アルゴリズムを選択できます。

## ライブラリAPI

```rust
use mvnormal_benchmark::{MvNormal, StandardRng};

let distribution = MvNormal::new(
    vec![0.0, 1.0],
    vec![vec![1.0, 0.2], vec![0.2, 2.0]],
)?;
let mut rng = StandardRng::new(42);

let sample = distribution.sample(&mut rng);

let mut output = vec![0.0; distribution.dimension()];
distribution.sample_inplace(&mut rng, &mut output)?;
```

`MvNormal::new` は、平均との次元整合性、共分散行列の正方形・有限値・対称性・正定値性を検証します。`StandardRng` は再現可能なベンチマーク用であり、暗号用途には適しません。
