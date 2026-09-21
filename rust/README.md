# Rust MvNormal benchmark

Rust 2024、標準ライブラリのみで実装した多変量正規分布のサンプラーです。
共分散行列 `Σ` を構築時に Cholesky 分解して `Σ = L L'` とし、各サンプルを

```text
x = μ + L z
```

で生成します。`z` は Marsaglia の polar 法による標準正規乱数です。

## ビルド

```sh
cargo build --release
```

## ベンチマーク実行

```sh
cargo run --release -- --dim 8 --samples 100 --repeats 2
```

成功時は次のCSV形式の1行だけを標準出力へ出します。

```text
rust,<dim>,<samples>,<repeats>,<setup_sec>,<avg_sample_sec>,<min_sample_sec>,<checksum>
```

`avg_sample_sec` と `min_sample_sec` は、それぞれ1回の `samples` バッチにかかった時間の平均・最小です。CLIは出力バッファだけを再利用する `sample_inplace` を使います。

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
