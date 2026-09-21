# MvNormal implementations and benchmarks

Julia の `Distributions.jl` が提供する `MvNormal` のサンプリング方法を、Julia（自前実装と公式実装）、C++23、Fortran 2023、Rust 2024、Python（NumPy/Numba）で実装したリポジトリです。

## TL;DR

実行結果と各言語の速度比較は [benchmark/REPORT.md](./benchmark/REPORT.md) を参照してください。

レポートを再生成するには、次を実行します。

```sh
./benchmark/run_and_report.sh
```

## アルゴリズム

平均ベクトル `μ` と正定値共分散行列 `Σ` に対して、構築時に

```text
Σ = L Lᵀ
```

を Cholesky 分解で求めます。サンプリング時は独立な標準正規ベクトル `z` を生成し、次を計算します。

```text
x = μ + L z
```

これは `Distributions.jl` の `MvNormal` が利用する基本的な whitening / unwhitening の考え方と同じです。

共通ベンチマークでは、入力も統一しています。0始まりの添字に対して
`μ[i] = 0.01 * i`、`Σ[i,j] = 0.25^abs(i-j)` を使います。

## 初回セットアップ

Julia の公式実装を使うための依存環境を作成します。

```sh
julia --project=./julia -e 'using Pkg; Pkg.instantiate()'
```

必要なツール:

- Julia（`julia/Project.toml` の環境）
- C++23 対応コンパイラ（例: `g++`）
- Fortran 2023 対応コンパイラ（例: GNU Fortran 16 以降）
- Rust 2024 対応の Cargo / rustc
- Python 3.12 以降と `uv`（NumPy/Numba を `python/pyproject.toml` で管理）

Python の依存環境は初回だけ次で作成します。

```sh
uv sync --project python
```

## 各言語の実行

すべてのベンチマークプログラムは、次のオプションを受け取ります。

```text
--dim N --samples N --repeats N
```

### Julia 自前実装

```sh
julia --project=./julia julia/benchmark.jl \
  --dim 32 --samples 10000 --repeats 3

# 標準正規乱数を Marsaglia polar 法または Ziggurat 法に固定する場合
julia --project=./julia julia/benchmark.jl \
  --dim 32 --samples 10000 --repeats 3 --normal polar
julia --project=./julia julia/benchmark.jl \
  --dim 32 --samples 10000 --repeats 3 --normal ziggurat
```

`--normal julia`（既定値）は Julia/Random の標準実装を使います。`polar` と
`ziggurat` は他言語版と比較するための自前 xorshift64 実装です。

### Julia / Distributions.jl 公式実装

```sh
julia --project=./julia julia/benchmark_distributions.jl \
  --dim 32 --samples 10000 --repeats 3
```

### C++23

```sh
g++ -std=c++23 -O3 -DNDEBUG \
  -I cxx cxx/benchmark.cpp -o /tmp/mvnormal_cxx
/tmp/mvnormal_cxx --dim 32 --samples 10000 --repeats 3
```

### Fortran 2023

```sh
mkdir -p /tmp/mvnormal-fortran-mod
gfortran -std=f2023 -O3 \
  -J /tmp/mvnormal-fortran-mod \
  fortran/mvnormal.f90 fortran/benchmark.f90 \
  -o /tmp/mvnormal_fortran
/tmp/mvnormal_fortran --dim 32 --samples 10000 --repeats 3
```

### Rust 2024

```sh
cargo test --manifest-path rust/Cargo.toml
cargo run --release --manifest-path rust/Cargo.toml -- \
  --dim 32 --samples 10000 --repeats 3
cargo run --release --manifest-path rust/Cargo.toml -- \
  --dim 32 --samples 10000 --repeats 3 --normal rand-distr
cargo run --release --manifest-path rust/Cargo.toml -- \
  --dim 32 --samples 10000 --repeats 3 --normal statrs
```

`--normal rand-distr` は `rand_distr` クレートの ZIGNOR 版 Ziggurat 法を使う
追加経路です。`--normal statrs` は `statrs::distribution::MultivariateNormal`
（nalgebra の Cholesky と `Distributions.jl` の `MvNormal` に相当する `L z + μ`）
を使います。どちらも `StdRng` の独自乱数列を使うため、先頭列はそれぞれ
`rust-rand-distr`、`rust-statrs` になります。

### Python（NumPy/Numba）

```sh
uv run --project python python/test_mvnormal.py
uv run --project python python/benchmark.py \
  --dim 32 --samples 10000 --repeats 3 --normal polar
uv run --project python python/benchmark.py \
  --dim 32 --samples 10000 --repeats 3 --normal ziggurat
```

## 全言語のベンチマーク

`benchmark/run.sh` が各言語をビルドして同じ条件で実行し、CSVに集約します。

```sh
./benchmark/run.sh benchmark/results.csv
```

条件は環境変数で変更できます。

```sh
DIM=128 SAMPLES=100000 REPEATS=10 ./benchmark/run.sh
```

ベンチマーク結果から Markdown の比較レポートを生成します。

```sh
python3 benchmark/generate_report.py \
  benchmark/results.csv benchmark/REPORT.md
```

実行とレポート生成を一度に行う場合:

```sh
./benchmark/run_and_report.sh
```

生成されるレポートには、セットアップ時間、サンプリング時間、最速実装に対する相対速度、checksum が含まれます。`julia-polar`、`cxx-polar`、`python-ziggurat` などは乱数アルゴリズムを表す別行です。入力と自前 RNG のシードは統一していますが、Cholesky 分解・浮動小数点演算・加算順序が異なるため、checksum は言語間で一致する必要はありません。Rust の `rust-rand-distr`／`rust-statrs` と Julia の `julia`／`julia-distributions` は各言語の標準 RNG を使う別系統です。

共通 runner が計測する主な行は次のとおりです。

| 行 | 内容 |
|---|---|
| `julia` / `julia-polar` / `julia-ziggurat` | Julia 標準 RNG / 自前 xorshift64 polar / Ziggurat |
| `julia-distributions` | 公式 `Distributions.jl` の `MvNormal` |
| `cxx-polar` / `cxx-ziggurat` | C++23 自前 xorshift64 |
| `fortran-polar` / `fortran-ziggurat` | Fortran 2023 自前 xorshift64 |
| `rust-polar` / `rust-ziggurat` | Rust 2024 自前 xorshift64 |
| `rust-rand-distr` | `rand_distr` の Zignor Ziggurat + `StdRng` |
| `rust-statrs` | `statrs::MultivariateNormal` + `StdRng` |
| `python-polar` / `python-ziggurat` | NumPy/Numba 自前 xorshift64（Numba JIT） |

Julia で複数サンプルを allocation-free に生成する場合は、出力行列と作業行列を事前確保して `sample!(rng, d, out, scratch)` を使います。

## ベンチマークの解釈

- `setup_sec`: `μ` と `Σ` から分布を構築し、Cholesky 因子を作る時間
- `avg_sample_sec`: 1 repeat 分の全サンプル生成時間の平均
- `min_sample_sec`: 1 repeat 分の全サンプル生成時間の最小値
- `checksum`: 最適化でサンプリング処理が除去されないための検査値

サンプリング時間は分布の構築後、JIT を使う処理では warmup 後に測定します。コンパイル時間は測定対象に含めません。
