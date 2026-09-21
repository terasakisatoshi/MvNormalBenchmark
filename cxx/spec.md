# C++ 実装の乱数生成仕様

## ベンチマーク用の乱数生成器

ベンチマークでは、言語間でアルゴリズムを比較できるように、`NormalRng` が
64ビット xorshift64 と `NormalAlgorithm` の組み合わせを使います。
シード `0` は内部で固定の非ゼロ状態に置き換えます。暗号用途には適しません。

`NormalRng` は、`NormalAlgorithm::MarsagliaPolar` と
`NormalAlgorithm::Ziggurat` を選択できます。
両方式とも `xorshift64` から一様乱数を生成します。
Marsaglia polar 法は2個の正規乱数をペアで生成し、Ziggurat法は256層のテーブルを使います。
一様乱数は xorshift64 出力の上位53ビットを使って
`(double(next_u64() >> 11) + 0.5) / 2^53` として求めます。
共通ベンチマークのサンプル生成シードは `0x5EED2021` です。

```cpp
mvnormal::NormalRng rng(42, mvnormal::NormalAlgorithm::Ziggurat);
distribution.sample_inplace(rng, output);
```

`MvNormal` の汎用インターフェースには `std::normal_distribution<double>` を渡すこともできます。
C++ 標準はその内部変換アルゴリズムを規定していないため、標準ライブラリ経路の方式や速度は
実装・バージョンに依存します。`NormalRng` を渡す経路が、他言語版と比較するベンチマーク経路です。

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

`NormalRng` についても、全言語で同じ標本列を得るには、xorshift64 のビット演算、シード、
浮動小数点演算、Ziggurat テーブル、棄却判定まで一致させる必要があります。
標準ライブラリの `std::normal_distribution` 経路では、標本列の一致は保証されません。
