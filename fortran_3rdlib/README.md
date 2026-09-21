# fortran_3rdlib

サードパーティの Fortran ライブラリで多変量正規分布 `MvNormal` を扱えるかを、
fpm (Fortran Package Manager) で試すための実験ディレクトリです。

## 分かったこと

- fpm registry に完成品の `MvNormal` パッケージは無い。
- `stdlib` は一変量 `rvs_normal` と Cholesky (`chol`) を持つので、
  `x = μ + L z` を自分で組めば MvNormal を実装できる（本ディレクトリの `stdlib_mvn`）。
- `forlab` には `randn` と `chol` はあるが MvNormal は無い。

つまり「完成したサードパーティ MvNormal」は存在せず、
stdlib の部品を土台に組む形になります。

## 構成

- `src/mvn_common.f90`: 共通入力（`μ[i]=0.01*i`, `Σ[i,j]=0.25^|i-j|`、0始まり相当）と checksum
- `src/stdlib_mvn.f90`: `stdlib` の `rvs_normal` + `chol` による実装（再利用できる `mvn_plan_t`）
- `app/main.f90`: 1サンプル生成して表示
- `app/benchmark.f90`: 共通 runner と同じCSV形式で `fortran-stdlib` 行を出力
- `test/check.f90`: 10万サンプルの標本平均・標本共分散が `μ`, `Σ` に一致するかを確認

## ビルドと実行

```sh
fpm run
fpm test
```

`stdlib` は fpm registry から取得します（初回ビルド時のみネット接続が必要）。
`benchmark/run.sh` は fpm がある場合に `fortran-stdlib` 行も自動で追加します。

CSV 行だけを直接確認する場合:

```sh
fpm build --profile release
./build/*/app/fortran_3rdlib_benchmark --dim 64 --samples 20000 --repeats 5
```

## 注意

- `stdlib` の分布乱数は `stdlib_random` の RNG を使うため、他の言語実装の
  xorshift64 とビット単位では一致しません。checksum は実行確認用です。
