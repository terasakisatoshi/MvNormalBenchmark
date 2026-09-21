# AGENTS.md

このリポジトリは、Julia の `Distributions.jl` にある `MvNormal` のサンプリングアルゴリズムを複数言語で実装し、同じ条件で速度比較するためのものです。

## 実装方針

各実装は次のアルゴリズムを共有します。

1. 正定値共分散行列 `Σ` を Cholesky 分解して `Σ = L Lᵀ` を得る。
2. 独立な標準正規ベクトル `z` を生成する。
3. `x = μ + L z` を計算する。

言語ごとの実装は次のディレクトリに置きます。

- `julia/`: 自前実装と公式 `Distributions.jl` のベンチマーク
- `cxx/`: C++23 実装
- `fortran/`: Fortran 2023 実装
- `rust/`: Rust 2024 実装
- `benchmark/`: 共通runner、CSV、Markdownレポート生成

共分散行列は正方形・対称・有限・正定値でなければなりません。次元不一致や不正な共分散行列は、各言語で明示的にエラーにしてください。

## CodeGraph

リポジトリルートに `.codegraph/` が存在する場合は、コードの場所や呼び出し関係を調べる前に CodeGraph を使います。

<!-- CODEGRAPH_START -->
## CodeGraph

In repositories indexed by CodeGraph (a `.codegraph/` directory exists at the repo root), reach for it BEFORE grep/find or reading files when you need to understand or locate code:

- **MCP tool** (when available): `codegraph_explore` answers most code questions in one call — the relevant symbols' verbatim source plus the call paths between them, including dynamic-dispatch hops grep can't follow. Name a file or symbol in the query to read its current line-numbered source. If it's listed but deferred, load it by name via tool search.
- **Shell** (always works): `codegraph explore "<symbol names or question>"` prints the same output.
<!-- CODEGRAPH_END -->

## 標準ワークフロー

### 1. Julia 環境を準備する

初回のみ、リポジトリルートで実行します。

```sh
julia --project=./julia -e 'using Pkg; Pkg.instantiate()'
```

`julia/Project.toml` と `julia/Manifest.toml` は公式 `Distributions.jl` の依存環境です。公式実装を使う場合は、必ずこの project 環境を指定します。

### 2. 各言語をコンパイル・テストする

個別に確認する場合は次を使います。

```sh
# Julia
julia --project=./julia julia/benchmark.jl \
  --dim 8 --samples 100 --repeats 2
julia --project=./julia julia/benchmark_distributions.jl \
  --dim 8 --samples 100 --repeats 2

# C++23
g++ -std=c++23 -O3 -DNDEBUG -I cxx \
  cxx/benchmark.cpp -o /tmp/mvnormal_cxx
/tmp/mvnormal_cxx --dim 8 --samples 100 --repeats 2

# Fortran 2023
mkdir -p /tmp/mvnormal-fortran-mod
gfortran -std=f2023 -O3 -J /tmp/mvnormal-fortran-mod \
  fortran/mvnormal.f90 fortran/benchmark.f90 \
  -o /tmp/mvnormal_fortran
/tmp/mvnormal_fortran --dim 8 --samples 100 --repeats 2

# Rust 2024
cargo test --manifest-path rust/Cargo.toml
cargo run --release --manifest-path rust/Cargo.toml -- \
  --dim 8 --samples 100 --repeats 2
```

各ベンチマークは次のCSV 1行を標準出力に出します。

```text
language,dim,samples,repeats,setup_sec,avg_sample_sec,min_sample_sec,checksum
```

### 3. 共通条件でベンチマークする

全言語をビルドして実行し、`benchmark/results.csv` に集約します。

```sh
./benchmark/run.sh benchmark/results.csv
```

条件を変更する場合:

```sh
DIM=128 SAMPLES=100000 REPEATS=10 \
  ./benchmark/run.sh benchmark/results.csv
```

`benchmark/run.sh` は Julia project を自動指定します。`Distributions.jl` が利用できれば、公式 `MvNormal` の行も `julia-distributions` として追加します。

### 4. Markdown レポートを生成する

CSVを手編集せず、次のスクリプトで比較表を生成します。

```sh
python3 benchmark/generate_report.py \
  benchmark/results.csv benchmark/REPORT.md
```

ベンチマークとレポート生成をまとめて行う場合:

```sh
./benchmark/run_and_report.sh \
  benchmark/results.csv benchmark/REPORT.md
```

レポートは `avg_sample_sec` の昇順で並び、セットアップ時間、平均・最小サンプリング時間、最速実装に対する相対速度を表示します。

## ベンチマーク比較時の注意

- サンプリング時間は、分布構築とJIT warmupの後に測定します。
- `avg_sample_sec` と `min_sample_sec` は、1 repeat 内の全サンプル生成時間です。
- 各言語で乱数生成器が異なるため、`checksum` は実行確認用であり、言語間で一致する必要はありません。
- コンパイル時間は比較対象に含めません。
- `benchmark/.build/`、`rust/target/`、Fortranの `.mod` と実行ファイルは生成物です。ベンチマーク結果は再生成できるため、手作業で修正しません。

## 変更後の確認

アルゴリズムやベンチマークを変更したら、少なくとも次を実行します。

```sh
cargo test --manifest-path rust/Cargo.toml
DIM=8 SAMPLES=200 REPEATS=2 ./benchmark/run.sh /tmp/mvnormal-results.csv
python3 benchmark/generate_report.py \
  /tmp/mvnormal-results.csv /tmp/mvnormal-report.md
```

出力CSVの全行に同じ `dim`、`samples`、`repeats` が入り、Markdown表に全実装が現れることを確認してください。
