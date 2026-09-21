# Python 実装の乱数生成仕様

## 一様乱数生成器

`NormalRng` は依存クレートを使わない `xorshift64` 型の疑似乱数生成器を
NumPy の `uint64` スカラー配列 `state` に保持します。

状態 `s` は次の順で更新します。

```text
s = s xor (s << 13)
s = s xor (s >> 7)
s = s xor (s << 17)
```

更新後の値から上位53ビットを取り出し、半単位を加えて `(0, 1)` の倍精度一様乱数に変換します。

```text
u = (f64(next_u64() >> 11) + 0.5) / 2^53
```

シード 0 は xorshift の固定点になるため `0x9E3779B97F4A7C15` へ置き換えます。
共通ベンチマークのサンプル生成シードは `0x5EED2021` です。

状態は長さ1の NumPy 配列として持ち、Numba でコンパイルした関数から破壊的に
更新できるようにしています。

## 標準正規乱数

`NormalAlgorithm.MARSAGLIA_POLAR` は Marsaglia の polar 法で2個ずつ生成します。

1. `u` と `v` を区間 `(-1, 1)` から生成する。
2. `s = u^2 + v^2` を計算する。
3. `s <= 0` または `s >= 1` の場合は棄却して再試行する。
4. `a = sqrt(-2 log(s) / s)` を計算する。
5. `u*a` と `v*a` を標準正規乱数として返す。

`NormalAlgorithm.ZIGGURAT` は256層の Ziggurat 法を使います。テーブル
`k`、`w`、`f` は初回に一度だけ生成し、矩形内部の候補を整数比較で受理します。
矩形外の候補は密度比較で再判定し、層0では指数分布を使って裾を生成します。

polar 法で生じた余剰の1値は `spare`/`has_spare` 配列に保持し、次のサンプルへ
持ち越します。C++/Fortran/Rust と同じ順序で乱数を消費するため、`dim=64`,
`--normal polar` の checksum は他言語実装と一致します。

## Numba によるコンパイル

`_next_u64`、`_polar_pair`、`_ziggurat`、`_fill_standard_normals`、
`_sample_matrix`、`_sample_checksum` を `numba.njit` でコンパイルします。ヘルパー
関数は `inline="always"` で呼び出し側へ展開します。Cholesky 分解は
`numpy.linalg.cholesky` を使います。

## 変換

Cholesky 因子は NumPy の行優先2次元配列に格納します。

```text
lower[row, column]
```

`_sample_matrix` / `_sample_checksum` は標準正規ベクトルを `z` に生成した後、
行を大きい番号から処理します。行 `row` の計算前には `z[..=row]` がまだ必要な
標準正規値なので、追加の作業ベクトルなしで

```text
z[row] = mean[row] + Σ lower[row, column] * z[column]
```

を計算できます (cxx/rust と同じ順序)。

## 再現性

同じシードとビルド条件では、Python 実装内で同じ乱数列を再生成できます。
他言語とはコンパイラや数学関数の実装が異なるため、checksum 一致は要求しません。
ただし共通の xorshift64 を使う polar/ziggurat 経路は、他実装と同じ乱数消費順序に
なるよう設計しています。
