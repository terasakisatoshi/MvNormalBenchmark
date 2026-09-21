# Python MvNormal benchmark

NumPy と Numba を使った多変量正規分布のサンプラーです。共分散行列 `Σ` を
構築時に Cholesky 分解して `Σ = L L'` とし、各サンプルを

```text
x = μ + L z
```

で生成します。`z` は他言語と同じ xorshift64 を共通シード `0x5EED2021` で
初期化し、Marsaglia の polar 法または256層の Ziggurat 法で生成します。

サンプル生成のホットループは Numba でコンパイルするため、C++/Rust 実装と
同じ for ループ構造のままネイティブ速度で実行されます。

## 依存と環境

`python/pyproject.toml` で NumPy と Numba を管理します。初回だけ次を実行します。

```sh
uv sync --project python
```

Numba の JIT コンパイル時間はベンチマークの測定対象外です。ウォームアップで
コンパイルしてから計測します。

## ベンチマーク実行

```sh
uv run --project python python/benchmark.py \
  --dim 8 --samples 100 --repeats 2 --normal polar
uv run --project python python/benchmark.py \
  --dim 8 --samples 100 --repeats 2 --normal ziggurat
```

成功時は次のCSV形式の1行だけを標準出力へ出します。

```text
python-polar,<dim>,<samples>,<repeats>,<setup_sec>,<avg_sample_sec>,<min_sample_sec>,<checksum>
```

`avg_sample_sec` と `min_sample_sec` は、それぞれ1回の `samples` バッチに
かかった時間の平均・最小です。

## テスト

```sh
uv run --project python python/test_mvnormal.py
```

## ライブラリAPI

```python
import numpy as np
from mvnormal import MvNormal, NormalRng, NormalAlgorithm

distribution = MvNormal(
    np.array([0.0, 1.0]),
    np.array([[1.0, 0.2], [0.2, 2.0]]),
)
rng = NormalRng(42, NormalAlgorithm.MARSAGLIA_POLAR)

sample = distribution.sample(rng)
samples = distribution.sample_matrix(rng, 1000)  # shape (1000, dimension)
```

`MvNormal` は、平均との次元整合性、共分散行列の正方形・有限値・対称性・
正定値性を検証します。`NormalRng` は再現可能なベンチマーク用であり、暗号用途
には適しません。
