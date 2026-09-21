# C++23 MvNormal

`mvnormal.hpp` は、平均ベクトル `μ` と共分散行列 `Σ` を受け取り、

```text
Σ = L * L'
x = μ + L * z
```

で多変量正規分布をサンプリングする、標準ライブラリだけの C++23 実装です。
コンストラクタで共分散行列の正方形・対称性・正定値性と、平均との次元整合性を検証します。

## ビルド

```sh
g++ -std=c++23 -O3 -Wall -Wextra -pedantic benchmark.cpp -o benchmark
```

## 実行

```sh
./benchmark --dim 32 --samples 10000 --repeats 5
./benchmark --dim 32 --samples 10000 --repeats 5 --normal ziggurat
```

ベンチマークは次の1行だけを標準出力へ出します。

```text
cxx-polar,<dim>,<samples>,<repeats>,<setup_sec>,<avg_sample_sec>,<min_sample_sec>,<checksum>
```

`--normal ziggurat` を指定した場合は、先頭列が `cxx-ziggurat` になります。

`sample_inplace` は、標準正規乱数を `output` に直接生成し、行優先レイアウトに合わせて後ろの行から変換します。
そのため、サンプルごとの一時ベクトルを使いません。
`sample_into` は、`scratch` と `output` を分けて使う互換APIとして残しています。
`--normal polar`（既定値）と `--normal ziggurat` でベンチマークの乱数生成法を選択できます。

```cpp
#include "mvnormal.hpp"

#include <random>
#include <vector>

std::vector<double> mean{1.0, -1.0};
std::vector<double> covariance{
    2.0, 0.5,
    0.5, 1.0,
};
mvnormal::MvNormal distribution(mean, covariance, 2);

std::mt19937_64 rng(42);
std::vector<double> output(2);
distribution.sample_inplace(rng, output);
```
