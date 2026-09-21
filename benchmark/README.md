# MvNormal benchmark

`run.sh` は、4言語の実装を同じ条件でビルド・実行し、結果をCSVに集約します。Julia で `Distributions.jl` が利用可能な場合は、公式 `MvNormal` の結果も `julia-distributions` として追加します。

Julia の依存環境は `julia/Project.toml` で管理しています。初回だけ次を実行してください。

```sh
julia --project=./julia -e 'using Pkg; Pkg.instantiate()'
```

```sh
./benchmark/run.sh
```

出力先を指定する場合:

```sh
./benchmark/run.sh /tmp/mvnormal-results.csv
```

条件は環境変数で変更できます。

```sh
DIM=128 SAMPLES=100000 REPEATS=10 ./benchmark/run.sh
```

CSVからMarkdownの速度比較表を生成する場合:

```sh
python3 benchmark/generate_report.py \
  benchmark/results.csv benchmark/REPORT.md
```

実行とレポート生成をまとめる場合:

```sh
./benchmark/run_and_report.sh
```

CSVの列は次の通りです。

- `setup_sec`: 平均・共分散からCholesky因子を作る時間
- `avg_sample_sec`: サンプル生成時間の反復平均
- `min_sample_sec`: サンプル生成時間の最小値
- `checksum`: 最適化でサンプリング処理が除去されないための検査値

サンプリング時間は、分布の構築とJIT warmupの後に測定します。言語ごとに乱数生成器は異なるため、サンプル列そのものではなく、アルゴリズムの実行時間を比較します。C++・Fortran・Rustのビルド成果物は `benchmark/.build/` に置かれます。
