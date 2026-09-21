# C++ 実装の乱数生成仕様

## 一様乱数生成器

ベンチマークでは標準ライブラリの `std::mt19937_64` を使います。

```cpp
std::mt19937_64 rng(0x4d764e6e6f726dULL);
```

これは 64 ビット出力の Mersenne Twister 系疑似乱数生成器です。
暗号用途には適しません。

## 標準正規乱数

標準正規乱数は `std::normal_distribution<double>(0.0, 1.0)` から取得します。

```cpp
value = standard_normal_(rng);
```

C++ 標準は `std::normal_distribution` の内部変換アルゴリズムを規定していません。
したがって、Box--Muller 法や Ziggurat 法などの特定方式をこの仕様では仮定しません。
使用する標準ライブラリの実装やバージョンによって、内部方式と速度が変わる可能性があります。

`standard_normal_` は `MvNormal` オブジェクトに保持して、サンプルごとに再構築しません。
標準ライブラリ実装が余剰の正規乱数を内部に保持する場合、その状態もサンプル間で再利用されます。

## MvNormal の変換

Cholesky 分解で得た `L` は、C++ の連続配列に行優先で格納します。

```text
cholesky_[row * dimension + column]
```

サンプリングでは `z` を生成した後、行 `row` ごとに
`x[row] = μ[row] + Σ L[row,column] z[column]` を計算します。

`sample_inplace` は `z` を出力バッファへ直接生成し、行を大きい番号から処理します。
行 `row` の計算前には `output[column]`（`column <= row`）がまだ `z[column]` なので、追加の作業ベクトルなしで変換できます。
行優先配列に合わせて、各行の係数を連続領域から読み出します。

## 再現性

乱数シードと標準ライブラリが同じでも、`std::normal_distribution` の仕様が内部方式を固定していないため、異なる標準ライブラリ間で標本列の一致は保証されません。
