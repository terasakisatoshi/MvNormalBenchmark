"""Multivariate normal sampling with a shared xorshift64 random stream.

The sampling algorithm matches the other language implementations in this
repository:

1. Factorize the positive definite covariance matrix as ``Sigma = L * L'``.
2. Draw a vector ``z`` of independent standard normal values.
3. Return ``mean + L * z``.

Two standard normal algorithms are provided: Marsaglia's polar method and a
256-layer Ziggurat method.  Both use the same xorshift64 generator so the
results are reproducible for a fixed seed.

The hot loops are compiled with Numba.  Without Numba the per-sample Python
loop would dominate the runtime; Numba keeps the code close to the C++/Rust
implementations while remaining a single Python source file.
"""

from __future__ import annotations

import math
from enum import IntEnum
from functools import lru_cache

import numpy as np
from numba import njit

ZIGGURAT_LAYERS = 256
ZIGGURAT_R = 3.6541528853610088
ZIGGURAT_V = 0.004928673233974658
TWO_POW_63 = 9.223372036854775808e18
INV_TWO_POW_53 = 1.0 / 9007199254740992.0
ZERO_SEED_REPLACEMENT = 0x9E3779B97F4A7C15
UINT64_MASK = 0xFFFFFFFFFFFFFFFF

NORMAL_POLAR = 1
NORMAL_ZIGGURAT = 2


class NormalAlgorithm(IntEnum):
    """Selects the standard normal algorithm for :class:`NormalRng`."""

    MARSAGLIA_POLAR = NORMAL_POLAR
    ZIGGURAT = NORMAL_ZIGGURAT

    def label(self) -> str:
        return "polar" if self is NormalAlgorithm.MARSAGLIA_POLAR else "ziggurat"


class MvNormalError(ValueError):
    """Raised when a covariance matrix or its dimensions are invalid."""


@lru_cache(maxsize=1)
def _ziggurat_tables() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build the ``(w, f, k)`` tables of the 256-layer Ziggurat method."""

    x = np.zeros(ZIGGURAT_LAYERS)
    k = np.zeros(ZIGGURAT_LAYERS, dtype=np.uint64)
    w = np.zeros(ZIGGURAT_LAYERS)
    f = np.zeros(ZIGGURAT_LAYERS)

    tail_density = math.exp(-0.5 * ZIGGURAT_R * ZIGGURAT_R)
    q = ZIGGURAT_V / tail_density

    x[ZIGGURAT_LAYERS - 1] = ZIGGURAT_R
    f[0] = 1.0
    f[ZIGGURAT_LAYERS - 1] = tail_density
    w[0] = q / TWO_POW_63
    w[ZIGGURAT_LAYERS - 1] = ZIGGURAT_R / TWO_POW_63
    k[0] = np.uint64(int(ZIGGURAT_R / q * TWO_POW_63))
    k[1] = np.uint64(0)

    for index in range(ZIGGURAT_LAYERS - 2, 0, -1):
        x[index] = math.sqrt(
            -2.0 * math.log(ZIGGURAT_V / x[index + 1] + f[index + 1])
        )
        f[index] = math.exp(-0.5 * x[index] * x[index])
        w[index] = x[index] / TWO_POW_63
        k[index + 1] = np.uint64(int(x[index] / x[index + 1] * TWO_POW_63))

    return w, f, k


@njit(cache=True, inline="always")
def _next_u64(state: np.ndarray) -> np.uint64:
    value = state[0]
    value = value ^ (value << np.uint64(13))
    value = value ^ (value >> np.uint64(7))
    value = value ^ (value << np.uint64(17))
    state[0] = value
    return value


@njit(cache=True, inline="always")
def _uniform_open01(state: np.ndarray) -> np.float64:
    raw = _next_u64(state) >> np.uint64(11)
    return (np.float64(raw) + 0.5) * INV_TWO_POW_53


@njit(cache=True, inline="always")
def _polar_pair(state: np.ndarray) -> tuple[np.float64, np.float64]:
    while True:
        u = 2.0 * _uniform_open01(state) - 1.0
        v = 2.0 * _uniform_open01(state) - 1.0
        radius_squared = u * u + v * v
        if radius_squared > 0.0 and radius_squared < 1.0:
            scale = math.sqrt(-2.0 * math.log(radius_squared) / radius_squared)
            return u * scale, v * scale


@njit(cache=True, inline="always")
def _ziggurat(
    state: np.ndarray, w: np.ndarray, f: np.ndarray, k: np.ndarray
) -> np.float64:
    while True:
        bits = _next_u64(state)
        index = np.int64(bits & np.uint64(0xFF))
        if (bits >> np.uint64(63)) == np.uint64(0):
            sign = 1.0
        else:
            sign = -1.0
        magnitude = bits & np.uint64(0x7FFFFFFFFFFFFFFF)
        x = np.float64(magnitude) * w[index]
        if magnitude < k[index]:
            return sign * x

        if index == 0:
            while True:
                tail_x = -math.log(_uniform_open01(state)) / ZIGGURAT_R
                tail_y = -math.log(_uniform_open01(state))
                if 2.0 * tail_y >= tail_x * tail_x:
                    return sign * (ZIGGURAT_R + tail_x)

        y = f[index] + _uniform_open01(state) * (f[index - 1] - f[index])
        if y < math.exp(-0.5 * x * x):
            return sign * x


@njit(cache=True, inline="always")
def _fill_standard_normals(
    state: np.ndarray,
    output: np.ndarray,
    algorithm: int,
    spare: np.ndarray,
    has_spare: np.ndarray,
    w: np.ndarray,
    f: np.ndarray,
    k: np.ndarray,
) -> None:
    length = output.size
    if length == 0:
        return

    if algorithm == NORMAL_ZIGGURAT:
        for index in range(length):
            output[index] = _ziggurat(state, w, f, k)
        return

    index = 0
    if has_spare[0] != 0:
        output[0] = spare[0]
        has_spare[0] = 0
        index = 1

    while index + 1 < length:
        first, second = _polar_pair(state)
        output[index] = first
        output[index + 1] = second
        index += 2

    if index < length:
        first, second = _polar_pair(state)
        output[index] = first
        spare[0] = second
        has_spare[0] = 1


@njit(cache=True)
def _sample_matrix(
    state: np.ndarray,
    mean: np.ndarray,
    lower: np.ndarray,
    output: np.ndarray,
    algorithm: int,
    spare: np.ndarray,
    has_spare: np.ndarray,
    w: np.ndarray,
    f: np.ndarray,
    k: np.ndarray,
) -> float:
    """Fill ``output`` with ``nsamples`` draws and return their total sum."""

    nsamples = output.shape[0]
    dimension = mean.size
    total = 0.0
    z = np.empty(dimension)
    for sample in range(nsamples):
        _fill_standard_normals(state, z, algorithm, spare, has_spare, w, f, k)

        # Generate z in place, then transform rows from bottom to top so the
        # still-needed z values remain intact (same order as cxx and rust).
        for row in range(dimension - 1, -1, -1):
            value = mean[row]
            for column in range(row + 1):
                value += lower[row, column] * z[column]
            z[row] = value

        for index in range(dimension):
            output[sample, index] = z[index]
            total += z[index]
    return total


@njit(cache=True)
def _sample_checksum(
    state: np.ndarray,
    mean: np.ndarray,
    lower: np.ndarray,
    nsamples: int,
    algorithm: int,
    spare: np.ndarray,
    has_spare: np.ndarray,
    w: np.ndarray,
    f: np.ndarray,
    k: np.ndarray,
) -> float:
    """Same as :func:`_sample_matrix` but without materializing the samples."""

    dimension = mean.size
    total = 0.0
    z = np.empty(dimension)
    for _ in range(nsamples):
        _fill_standard_normals(state, z, algorithm, spare, has_spare, w, f, k)
        for row in range(dimension - 1, -1, -1):
            value = mean[row]
            for column in range(row + 1):
                value += lower[row, column] * z[column]
            z[row] = value
        for index in range(dimension):
            total += z[index]
    return total


class NormalRng:
    """Deterministic xorshift64 generator with polar/Ziggurat normals."""

    def __init__(
        self, seed: int, algorithm: NormalAlgorithm = NormalAlgorithm.MARSAGLIA_POLAR
    ) -> None:
        algorithm = NormalAlgorithm(algorithm)
        value = int(seed) & UINT64_MASK
        if value == 0:
            value = ZERO_SEED_REPLACEMENT

        self.algorithm = algorithm
        self._state = np.array([value], dtype=np.uint64)
        self._spare = np.zeros(1, dtype=np.float64)
        self._has_spare = np.zeros(1, dtype=np.int64)
        if algorithm is NormalAlgorithm.ZIGGURAT:
            self._w, self._f, self._k = _ziggurat_tables()
        else:
            self._w = np.zeros(0, dtype=np.float64)
            self._f = np.zeros(0, dtype=np.float64)
            self._k = np.zeros(0, dtype=np.uint64)

    def fill_standard_normals(self, output: np.ndarray) -> None:
        output = np.ascontiguousarray(output, dtype=np.float64)
        _fill_standard_normals(
            self._state,
            output.reshape(-1),
            int(self.algorithm),
            self._spare,
            self._has_spare,
            self._w,
            self._f,
            self._k,
        )


class MvNormal:
    """Multivariate normal distribution with a cached lower Cholesky factor."""

    def __init__(self, mean, covariance) -> None:
        mean = np.asarray(mean, dtype=np.float64)
        covariance = np.asarray(covariance, dtype=np.float64)

        if mean.ndim != 1:
            raise MvNormalError("mean must be one-dimensional")
        if covariance.ndim != 2:
            raise MvNormalError("covariance matrix must be two-dimensional")

        dimension = mean.size
        if dimension == 0:
            raise MvNormalError("distribution dimension must be positive")
        if covariance.shape[0] != covariance.shape[1]:
            raise MvNormalError("covariance matrix must be square")
        if covariance.shape[0] != dimension:
            raise MvNormalError(
                "mean and covariance dimensions must agree"
            )
        if not np.all(np.isfinite(mean)):
            raise MvNormalError("mean must contain only finite values")
        if not np.all(np.isfinite(covariance)):
            raise MvNormalError("covariance matrix must contain only finite values")

        scale = max(1.0, float(np.max(np.abs(covariance))))
        if not np.allclose(covariance, covariance.T, rtol=0.0, atol=1e-12 * scale):
            raise MvNormalError("covariance matrix must be symmetric")

        try:
            lower = np.linalg.cholesky(covariance)
        except np.linalg.LinAlgError as error:
            raise MvNormalError(
                "covariance matrix must be positive definite"
            ) from error

        self._mean = np.ascontiguousarray(mean)
        self._lower = np.ascontiguousarray(lower)

    @property
    def dimension(self) -> int:
        return self._mean.size

    def mean(self) -> np.ndarray:
        return self._mean.copy()

    def covariance(self) -> np.ndarray:
        return self._lower @ self._lower.T

    def cholesky_factor(self) -> np.ndarray:
        return self._lower.copy()

    def sample_matrix(self, rng: NormalRng, nsamples: int) -> np.ndarray:
        if nsamples < 0:
            raise ValueError("number of samples must be non-negative")
        output = np.empty((nsamples, self.dimension), dtype=np.float64)
        _sample_matrix(
            rng._state,
            self._mean,
            self._lower,
            output,
            int(rng.algorithm),
            rng._spare,
            rng._has_spare,
            rng._w,
            rng._f,
            rng._k,
        )
        return output

    def sample_checksum(self, rng: NormalRng, nsamples: int) -> float:
        if nsamples < 0:
            raise ValueError("number of samples must be non-negative")
        return _sample_checksum(
            rng._state,
            self._mean,
            self._lower,
            int(nsamples),
            int(rng.algorithm),
            rng._spare,
            rng._has_spare,
            rng._w,
            rng._f,
            rng._k,
        )

    def sample(self, rng: NormalRng) -> np.ndarray:
        return self.sample_matrix(rng, 1)[0]
