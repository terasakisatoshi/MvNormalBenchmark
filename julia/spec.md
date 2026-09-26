# Julia 実装の乱数生成仕様

## 対象

このディレクトリには、次の2種類のサンプリング経路があります。

- `mvnormal.jl`: `MvNormal` の自前実装
- `benchmark_distributions.jl`: `Distributions.jl` が提供する `MvNormal`

どちらも標準正規ベクトル `z` を生成し、Cholesky 因子 `L` を使って
`x = μ + Lz` を計算します。

## 一様乱数生成器

ベンチマークでは `Random.seed!` で種を固定したうえで、`Random.default_rng()` を使います。

```julia
Random.seed!(0xc0ffee)
rng = Random.default_rng()
```

`default_rng()` はタスクローカルな乱数生成器を返します。
その実体は Xoshiro256++ 系の `Xoshiro` です。
暗号用途には適しません。

## 標準正規乱数

自前実装は `Random.randn!` を使って、出力バッファ全体を標準正規乱数で埋めます。

```julia
randn!(rng, out)
```

`Float64` の標準正規乱数は **Ziggurat 法**で生成します。
Julia の `Random` は Marsaglia と Tsang が提案した Ziggurat 法を実装しています。

### Ziggurat 法

標準正規密度 `f(x) = exp(-x²/2)` の下側の面積を、等しい面積を持つ 256 枚の水平な層に分割します。
各層は内側の矩形と、その外側のくさび形の領域からなります。
層の境界は事前計算済みの定数テーブルに格納されています。

- `wi`: 各層の幅
- `ki`: 各層の内側矩形を判定する整数しきい値
- `fi`: 各層境界での密度値

1 個の乱数は次の手順で生成します。

1. 一様乱数から 52 bit の整数 `r` を取り出します。
2. 最下位 bit を符号に使い、残りから層の番号 `idx = rabs & 0xFF` を選びます。
3. `x = rabs * wi[idx+1]` として層内の候補値を作ります。
4. `rabs < ki[idx+1]` なら、候補値は内側矩形の中にあるため、そのまま採用します。
5. 内側矩形の外側なら `randn_unlikely` に進み、くさび形領域の採否判定か末尾層の処理を行います。

手順 4 で採用される割合は約 99.3% です。
多くの乱数は整数乱数 1 個とテーブル参照だけで決まり、これが Ziggurat 法の高速性の理由です。

内側矩形の外側に落ちた場合は、`idx` の値で分岐します。

- `idx == 0` の層は正規分布の裾 (`|x| >= 3.654...`) に対応し、指数分布に従う候補を生成して採否を判定します。
- それ以外の層では、くさび形領域の下側の三角形から候補を生成し、密度 `exp(-x²/2)` と比較して採否を決めます。
- 棄却された場合は `randn` を再帰呼び出しして最初からやり直します。

裾のしきい値は定数 `ziggurat_nor_r = 3.6541528853610088` です。
この値を超える確率は小さく、末尾層の処理に入る頻度は低いです。

### 型による違い

`randn(rng, Float64)` が上記の Ziggurat 法を使います。
`Float16` と `Float32` は `Float64` の Ziggurat 値を受け取ってから型変換します。
そのため `Float16`、`Float32`、`Float64` はすべて Ziggurat 法を経由します。

`BigFloat` など `BitFloatType` に含まれない浮動小数点型は、`Random` のフォールバック実装である Box--Muller 変換 (Marsaglia の極座標法) を使います。

ベンチマークが使う `default_rng()` (`TaskLocalRNG`) と `Vector{Float64}` の組み合わせでは、`randn!` は長さ 7 以上でバルク経路の `_randn` を呼び出します。
短いベクトルではスカラーの `randn` を 1 要素ずつ呼びますが、どちらの経路も同じ Ziggurat 法です。

内部方式は Julia `Random` の実装に依存するため、Julia のバージョンによって変わる可能性があります。

`Distributions.jl` のベンチマークも `rand!(rng, d, out)` を呼び出します。
その経路で使われる標準正規乱数生成は `Distributions.jl` と Julia `Random` の実装に委譲されます。

## 自前RNGによる比較

`mvnormal.jl` には、他言語版と同じ `xorshift64` の状態更新を使う `NormalRNG` もあります。
`MarsagliaPolarRNG(seed)` は Marsaglia の polar 法を使い、`ZigguratRNG(seed)` は256層のZiggurat法を使います。
共通ベンチマークでは、全言語で標準正規乱数のシードを `0x5EED2021`
（10進数 `1592598561`）に統一しています。

```julia
polar_rng = MarsagliaPolarRNG(0x5EED2021)
ziggurat_rng = ZigguratRNG(0x5EED2021)
sample!(polar_rng, d, out)
sample!(ziggurat_rng, d, out)
```

これらはJulia標準の `Random.default_rng()` とは別の比較用RNGです。
同じシードと同じアルゴリズムを各言語で使っても、浮動小数点演算順序や数学関数の実装が異なる場合、ビット単位の一致は保証されません。

複数サンプルを繰り返し生成する場合は、出力行列を再利用できます。

```julia
out = Matrix{Float64}(undef, dimension(d), nsamples)
sample!(ziggurat_rng, d, out)
```

この経路は標準正規乱数を `out` に直接生成し、`lmul!` で `L * out` を
in-place に書き込み、反復ごとの行列確保を避けます。
`L` は `LowerTriangular` として渡すため、BLAS の三角行列積 (`trmm`) が使われます。

## MvNormal の変換

自前実装では、生成した `z` を出力バッファに置き、列 `j` を大きい番号から順に処理します。
各列の処理前には `out[j]` がまだ元の `z[j]` であるため、それをローカル変数に退避してから `out[j]` を `μ[j]` に置き換えます。

Julia の行列は列優先なので、内側の `i` を `j:n` として `L[i, j]` と `out[i]` を連続アクセスします。
この逆向き列累積により、追加の正規乱数用ベクトルを確保せずに `out = μ + Lz` を計算します。
内側の `i` ループは 2 要素ずつ展開し、`muladd` で積和演算として実行します。
これにより依存関係のない更新を連続して発行できます。

## ベンチマーク上の注意

`julia/benchmark.jl` は自前実装を測定し、`julia/benchmark_distributions.jl` は `Distributions.jl` の公式実装を測定します。

両方とも JIT コンパイル後のサンプリング経路を測定するため、測定前に同じサンプル数のウォームアップを実行します。

`--batch` を付けると、行列版の `sample!(rng, d, out)` を使って1 repeat 分のサンプルをまとめて生成します。
標準正規乱数を `out` に直接生成し、`lmul!` で `L * out` を in-place に計算するため、乱数列は per-sample 経路と同じです。
out-of-place の `mul!` と別の scratch 行列を使う場合に比べ、in-place の三角行列積のほうがわずかに速いです。
`lmul!` は BLAS を使うため、比較を単一スレッドに揃える目的で `BLAS.set_num_threads(1)` を呼びます。
バッチ経路の行名は `julia-batch`、`julia-polar-batch`、`julia-ziggurat-batch` です。
