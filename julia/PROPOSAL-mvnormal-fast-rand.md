# 提案書: `Distributions.jl` の `MvNormal` サンプリング高速化

対象: `julia/extern/Distributions.jl`（v0.25.131, commit `aa53c319`）
関連: `julia/mvnormal.jl`（自前実装）, `benchmark/results.csv`

## 1. 目的

`Distributions.jl` の `MvNormal` の per-sample サンプリングを、本リポジトリの
自前実装（`julia/mvnormal.jl`）と同等の速度まで引き上げる。数理アルゴリズム・
乱数列は変えない。

## 2. 現状

`benchmark/results.csv`（dim=64, samples=20000, repeats=5）:

| language | avg_sample_sec |
|---|---|
| `julia`（自前 per-sample） | 0.00888 |
| `julia-distributions` | 0.01774 |
| `julia-batch`（自前 batch） | 0.00631 |
| `julia-distributions-batch` | 0.00575 |

per-sample では自前実装が約2倍速く、batch では公式実装のほうが速い。

## 3. 原因分析

対象コード: `src/multivariate/mvnormal.jl:257-263`

```julia
function _rand!(rng::AbstractRNG, d::MvNormal, x::AbstractVecOrMat)
    unwhiten!(d.Σ, randn!(rng, x))
    x .+= d.μ
    return x
end
```

`Profile` で調べると、`randn!`（自前実装と共通）に次いで `unwhiten!` が重い。
理由は次の2点。

1. **Cholesky 因子のアクセスが非連続**
   `PDMats.unwhiten!(a::PDMat, x)` は `lmul!(chol_lower(chol), x)` を呼ぶ。
   LAPACK は上三角で保持するため、PDMat の `chol.L` は
   `LowerTriangular{Float64, Adjoint{Float64, Matrix{Float64}}}`
   になる。`Adjoint`/`LowerTriangular` 経由の要素アクセスは列優先の素の
   `Matrix` より遅い。自前実装は `L = tril(Matrix(factor.L))` を素の
   `Matrix` として保持しているため速い。

2. **`x .+= d.μ` が分離パス**
   三角形積のあと全要素をもう一度走査する。自前実装は逆順列累積で
   `μ + Lz` を1パスで計算する。

計測（dim=64, 20000 samples, `Xoshiro`, `BLAS.set_num_threads(1)`）:

| 経路 | avg_sample_sec |
|---|---|
| 現状 `_rand!`（`unwhiten!` + 加算） | 0.0163 |
| 自前 `sample!` | 0.0093 |

## 4. 提案

`_rand!` の直後（`src/multivariate/mvnormal.jl:263` の後）に、`FullNormal` かつ
strided vector 用の特化メソッドを追加する。`parent(cholesky(d.Σ).U)` で
**素の列優先 `Matrix`** を取り出し、`μ + Lz` を逆順列累積＋`muladd` で1パス
計算する。

```julia
function _rand!(rng::AbstractRNG,
                d::FullNormal{<:Real,<:PDMat{<:Real,<:StridedMatrix}},
                x::StridedVector{Float64})
    randn!(rng, x)
    U = parent(cholesky(d.Σ).U)   # 素の Matrix（列優先）: chol.L == U'
    m = d.μ
    n = length(x)
    @inbounds for r in n:-1:1     # 逆順なので x[c] はまだ z[c]
        acc = m[r]
        for c in 1:r              # U[c,r] は列 r の連続アクセス
            acc = muladd(U[c, r], x[c], acc)
        end
        x[r] = acc
    end
    return x
end
```

### 数理的根拠

`Σ = U'U` かつ `chol.L == chol.U'` なので、`L = U'`。
`x = μ + Lz` の第 `r` 成分は

```math
x_r = \mu_r + \sum_{c \le r} U_{c,r}\, z_c
```

`r` を降順に処理すると、内側の `c` ループ時点で `x[c]`（`c < r`）はまだ
`z[c]` のままなので、追加バッファなしで書ける。`U[c,r]` は列 `r` の連続
アクセスになる。

### 実測（提案コードを `_rand!` に差し替えた場合）

| 経路 | avg_sample_sec |
|---|---|
| 現状 `_rand!` | 0.0163 |
| 提案 `_rand!` | **0.0094** |
| 自前 `sample!` | 0.0093 |

提案コードと自前実装の出力差は `max abs diff = 8.9e-16`（同一乱数列・同一数理、
浮動小数点演算順序の違いのみ）。

## 5. 実装上の注意

- **`parent` が必須**。`d.Σ.chol.U`（`UpperTriangular`）や `chol.L`
  （`Adjoint`）のまま `@inbounds` で回すと三角形の分岐が残り、逆に遅くなる。
  実測では 0.0248s で、現状より悪化した。
- **型制約でフォールバック**。`PDMat{<:Real,<:StridedMatrix}` に限定し、
  `PDiagMat`/`ScalMat`、`Float32`、非 strided 配列などは既存の汎用
  `_rand!` に委譲する。
- **batch（行列）経路は変更しない**。既に `randn!(rng, 行列)` + BLAS `trmm` で
  動いており、`julia-distributions-batch` は自前 `julia-batch` より速い。
  特化メソッドは `StridedVector{Float64}` のみを対象にする。
- `ZeroMeanFullNormal`（`Zeros` の平均）は今回は対象外。`m[r]` は `Zeros`
  でも動くが、必要になったら別メソッドで対応する。

## 6. 検証手順

```sh
# クローンを差し替え（Manifest を汚したくなければ専用環境推奨）
julia --project=./julia -e 'using Pkg; Pkg.develop(path="./julia/extern/Distributions.jl")'

julia --project=./julia julia/benchmark_distributions.jl \
  --dim 64 --samples 20000 --repeats 5
```

確認事項:

- `avg_sample_sec` が約4割改善すること。
- 乱数列（`randn!` の呼び出し）は変えていないため、同一シードでの
  サンプル列が提案前後で浮動小数点誤差の範囲に収まること。
- `cargo test` 等の他言語テストに影響がないこと（Julia のみの変更）。

## 7. 代替案

1. **PDMats 側の `unwhiten!` を改善**
   `chol_lower` が返す `LowerTriangular`/`Adjoint` のインデックスを
   素の `Matrix` 経由に変える。ただし本リポジトリには PDMats がなく、
   上流2パッケージにまたがる変更になるため見送り。
2. **`MvNormal` 構築時に密な下三角 `L` をキャッシュ**
   最も速いが、`PDMat` の構造変更が必要で上流互換の観点から侵襲的。
   今回の提案は `parent` だけで同等性能が出る。
3. **`unwhiten!` + `x .+= μ` をブロードキャスト融合のみ**
   `@.` 等で1パス化。実測では Adjoint の非連続アクセスが残るため効果は限定的。

## 8. 未解決点 / 判断してほしい点

- `parent(cholesky(d.Σ).U)` が常に素の `StridedMatrix` を返すという前提を、
  型制約 `PDMat{<:Real,<:StridedMatrix}` だけで十分に保証できるか。
- `Float32` や一般 `Real` にも対応するか（`muladd` の型、`randn!` の経路）。
- 上流への PR を想定するか、リポジトリ内のベンチマーク専用パッチに留めるか。
- バッチ経路にも同様の特化を入れて `julia-batch` を超えられるか
  （現状は不要と判断）。
