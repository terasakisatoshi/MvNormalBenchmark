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
```

ベンチマークは次の1行だけを標準出力へ出します。

```text
cxx,<dim>,<samples>,<repeats>,<setup_sec>,<avg_sample_sec>,<min_sample_sec>,<checksum>
```

`sample_into` は、標準正規乱数を入れる `scratch` と結果を入れる `output` を呼び出し側で事前確保し、サンプルごとの一時ベクトル確保を避けます。

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
std::vector<double> scratch(2);
std::vector<double> output(2);
distribution.sample_into(rng, scratch, output);
```
