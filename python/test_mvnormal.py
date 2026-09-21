"""Tests for the Python MvNormal implementation.

Run with the project environment:

    uv run --project python python/test_mvnormal.py
"""

from __future__ import annotations

import unittest

import numpy as np

from mvnormal import MvNormal, MvNormalError, NormalAlgorithm, NormalRng

ALGORITHMS = (NormalAlgorithm.MARSAGLIA_POLAR, NormalAlgorithm.ZIGGURAT)


def identity_problem(dimension: int) -> MvNormal:
    return MvNormal(np.zeros(dimension), np.eye(dimension))


class MvNormalConstructionTest(unittest.TestCase):
    def test_reconstructs_covariance(self) -> None:
        covariance = np.array([[2.0, 0.5], [0.5, 1.0]])
        distribution = MvNormal(np.array([1.0, -1.0]), covariance)
        self.assertEqual(distribution.dimension, 2)
        np.testing.assert_allclose(distribution.covariance(), covariance, atol=1e-12)
        factor = distribution.cholesky_factor()
        np.testing.assert_allclose(factor, np.tril(factor), atol=0.0)

    def test_rejects_invalid_covariance(self) -> None:
        with self.assertRaises(MvNormalError):
            MvNormal(np.zeros(2), np.zeros((2, 3)))
        with self.assertRaises(MvNormalError):
            MvNormal(np.zeros(2), np.ones((3, 3)))
        with self.assertRaises(MvNormalError):
            MvNormal(np.zeros(2), np.array([[1.0, 0.2], [0.3, 1.0]]))
        with self.assertRaises(MvNormalError):
            MvNormal(np.zeros(2), np.array([[1.0, 2.0], [2.0, 1.0]]))
        with self.assertRaises(MvNormalError):
            MvNormal(np.array([0.0, np.nan]), np.eye(2))
        with self.assertRaises(MvNormalError):
            MvNormal(np.zeros(0), np.zeros((0, 0)))


class SamplingTest(unittest.TestCase):
    def test_reproducible(self) -> None:
        for algorithm in ALGORITHMS:
            distribution = identity_problem(4)
            first = distribution.sample_matrix(NormalRng(42, algorithm), 1000)
            second = distribution.sample_matrix(NormalRng(42, algorithm), 1000)
            np.testing.assert_array_equal(first, second)

    def test_sample_shape_and_finiteness(self) -> None:
        for algorithm in ALGORITHMS:
            distribution = identity_problem(3)
            samples = distribution.sample_matrix(NormalRng(7, algorithm), 129)
            self.assertEqual(samples.shape, (129, 3))
            self.assertTrue(np.all(np.isfinite(samples)))

    def test_polar_odd_dimension_reuses_spare(self) -> None:
        distribution = identity_problem(5)
        samples = distribution.sample_matrix(
            NormalRng(11, NormalAlgorithm.MARSAGLIA_POLAR), 257
        )
        self.assertTrue(np.all(np.isfinite(samples)))
        again = distribution.sample_matrix(
            NormalRng(11, NormalAlgorithm.MARSAGLIA_POLAR), 257
        )
        np.testing.assert_array_equal(samples, again)

    def test_standard_normal_moments(self) -> None:
        for algorithm in ALGORITHMS:
            distribution = identity_problem(2)
            samples = distribution.sample_matrix(NormalRng(2021, algorithm), 200_000)
            self.assertLess(abs(float(samples.mean())), 0.02)
            self.assertLess(abs(float(samples.var()) - 1.0), 0.03)

    def test_sample_checksum_matches_matrix_total(self) -> None:
        for algorithm in ALGORITHMS:
            distribution = identity_problem(4)
            matrix_total = float(
                distribution.sample_matrix(NormalRng(99, algorithm), 500).sum()
            )
            checksum_total = distribution.sample_checksum(
                NormalRng(99, algorithm), 500
            )
            self.assertAlmostEqual(matrix_total, checksum_total, places=9)


if __name__ == "__main__":
    unittest.main()
