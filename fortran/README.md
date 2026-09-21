# Fortran MvNormal

`mvnormal.f90` は、平均ベクトル `μ` と共分散行列 `Σ` から多変量正規分布を生成する再利用可能な Fortran モジュールです。
生成時に `Σ` の正方形・対称性・正定値、および `μ` との次元整合性を検証します。
サンプリングは次の手順で行います。

1. `Σ = L L'` を Cholesky 分解する。
2. Marsaglia の polar 法で標準正規ベクトル `z` を生成する。
3. `x = μ + L z` を計算する。

## ビルド

gfortran の Fortran 2023 モードでビルドします。

```sh
gfortran -std=f2023 -O3 -o benchmark mvnormal.f90 benchmark.f90
```

## 実行

```sh
./benchmark --dim 32 --samples 10000 --repeats 3
```

`--dim=N`、`--samples=N`、`--repeats=N` の形式も使えます。出力は CSV 形式の1行です。

```text
fortran,<dim>,<samples>,<repeats>,<setup_sec>,<avg_sample_sec>,<min_sample_sec>,<checksum>
```

`avg_sample_sec` は各 repeat で `samples` 個を生成した全体時間の平均、`min_sample_sec` は各 repeat の全サンプル時間の最小値です。
ベンチマークは `sample_inplace` を使い、標準正規乱数と結果を `output` へ直接書き込みます。
`sample_into` は `scratch` と `output` を分ける互換APIとして残しています。
`normal_rng_t` はベンチマーク用の再現可能な乱数生成器で、サンプリング手続きに明示的に渡します。
