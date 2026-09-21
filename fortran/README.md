# Fortran MvNormal

`mvnormal.f90` は、平均ベクトル `μ` と共分散行列 `Σ` から多変量正規分布を生成する再利用可能な Fortran モジュールです。
生成時に `Σ` の正方形・対称性・正定値、および `μ` との次元整合性を検証します。
サンプリングは次の手順で行います。

1. `Σ = L L'` を Cholesky 分解する。
2. `--normal` で選択した Marsaglia の polar 法または Ziggurat 法で標準正規ベクトル `z` を生成する。
3. `x = μ + L z` を計算する。

## ビルド

gfortran の Fortran 2023 モードでビルドします。

```sh
gfortran -std=f2023 -O3 -o benchmark mvnormal.f90 benchmark.f90
```

## 実行

```sh
./benchmark --dim 32 --samples 10000 --repeats 3
./benchmark --dim 32 --samples 10000 --repeats 3 --normal ziggurat
```

`--dim=N`、`--samples=N`、`--repeats=N` の形式も使えます。出力は CSV 形式の1行です。
`--normal polar`（既定値）と `--normal ziggurat` で標準正規乱数アルゴリズムを選択できます。

```text
fortran-polar,<dim>,<samples>,<repeats>,<setup_sec>,<avg_sample_sec>,<min_sample_sec>,<checksum>
```

`--normal ziggurat` を指定した場合は、先頭列が `fortran-ziggurat` になります。

`avg_sample_sec` は各 repeat で `samples` 個を生成した全体時間の平均、`min_sample_sec` は各 repeat の全サンプル時間の最小値です。
ベンチマークは `sample_inplace` を使い、標準正規乱数と結果を `output` へ直接書き込みます。
`sample_into` は `scratch` と `output` を分ける互換APIとして残しています。
`normal_rng_t` はベンチマーク用の再現可能な乱数生成器で、サンプリング手続きに明示的に渡します。

`sample_inplace` は列優先配置に合わせて1列ずつ処理します。
Ziggurat の通常の受理経路は、インライン展開しやすい小さな関数に分離しています。
通常の `do` だけを使い、複数スレッドでの並列化は行いません。
詳細は [spec.md](spec.md)、全言語の測定結果は [比較レポート](../benchmark/REPORT.md) を参照してください。

## テスト

`fortran/` から次を実行します。

```sh
gfortran -std=f2023 -O0 -g -fcheck=all -o /tmp/test_mvnormal_debug \
  mvnormal.f90 test_mvnormal.f90
/tmp/test_mvnormal_debug
gfortran -std=f2023 -O3 -o /tmp/test_mvnormal_release \
  mvnormal.f90 test_mvnormal.f90
/tmp/test_mvnormal_release
```

既知の下三角行列から作った共分散行列で `μ + Lz` と照合し、奇数と偶数の次元、非連続配列、両方の乱数方式を検証します。
乱数列については、配列の一括生成、奇数個を含む分割生成、ペア生成が一致することを確認します。
