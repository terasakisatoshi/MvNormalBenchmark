//! A small, dependency-free multivariate normal sampler.
//!
//! Sampling follows the usual Cholesky construction used by
//! `Distributions.jl`'s `MvNormal`: for `Sigma = L * L'` and a vector `z`
//! of independent standard normal values, return `mu + L * z`.

use std::error::Error;
use std::fmt;

/// Errors returned while constructing a multivariate normal distribution.
#[derive(Debug, Clone, PartialEq)]
pub enum MvNormalError {
    /// The mean vector and covariance matrix do not have matching dimensions.
    DimensionMismatch { mean: usize, covariance_rows: usize },
    /// A mean entry is NaN or infinite.
    NonFiniteMean { index: usize, value: f64 },
    /// A covariance row has a different length from the covariance dimension.
    NonSquare {
        row: usize,
        expected: usize,
        actual: usize,
    },
    /// The covariance matrix is not symmetric within the numerical tolerance.
    NotSymmetric {
        row: usize,
        column: usize,
        left: f64,
        right: f64,
    },
    /// A covariance entry is NaN or infinite.
    NonFinite {
        row: usize,
        column: usize,
        value: f64,
    },
    /// The covariance matrix is not positive definite.
    NotPositiveDefinite { index: usize, value: f64 },
}

impl fmt::Display for MvNormalError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DimensionMismatch {
                mean,
                covariance_rows,
            } => write!(
                f,
                "mean has length {mean}, but covariance has {covariance_rows} rows"
            ),
            Self::NonFiniteMean { index, value } => {
                write!(f, "mean entry {index} is not finite: {value}")
            }
            Self::NonSquare {
                row,
                expected,
                actual,
            } => write!(
                f,
                "covariance row {row} has length {actual}, expected {expected}"
            ),
            Self::NotSymmetric {
                row,
                column,
                left,
                right,
            } => write!(
                f,
                "covariance is not symmetric at ({row}, {column}): {left} != {right}"
            ),
            Self::NonFinite { row, column, value } => write!(
                f,
                "covariance entry ({row}, {column}) is not finite: {value}"
            ),
            Self::NotPositiveDefinite { index, value } => write!(
                f,
                "covariance is not positive definite at diagonal {index}: {value}"
            ),
        }
    }
}

impl Error for MvNormalError {}

/// Errors returned when a sampling buffer has the wrong size.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SampleError {
    /// `scratch` or `output` does not have one entry per dimension.
    BufferLength {
        name: &'static str,
        expected: usize,
        actual: usize,
    },
}

impl fmt::Display for SampleError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::BufferLength {
                name,
                expected,
                actual,
            } => {
                write!(f, "{name} has length {actual}, expected {expected}")
            }
        }
    }
}

impl Error for SampleError {}

/// A multivariate normal distribution with a cached lower Cholesky factor.
///
/// `covariance` is a row-major square matrix. Its entries are validated for
/// finiteness, symmetry, and positive definiteness during construction.
#[derive(Debug, Clone)]
pub struct MvNormal {
    mean: Vec<f64>,
    // Lower triangular L, stored row-major. Entries above the diagonal are 0.
    cholesky: Vec<f64>,
    dimension: usize,
}

impl MvNormal {
    /// Construct `N(mean, covariance)` and cache `L` such that `Sigma = L*L'`.
    pub fn new(mean: Vec<f64>, covariance: Vec<Vec<f64>>) -> Result<Self, MvNormalError> {
        let dimension = covariance.len();
        if mean.len() != dimension {
            return Err(MvNormalError::DimensionMismatch {
                mean: mean.len(),
                covariance_rows: dimension,
            });
        }
        if let Some((index, &value)) = mean
            .iter()
            .enumerate()
            .find(|(_, value)| !value.is_finite())
        {
            return Err(MvNormalError::NonFiniteMean { index, value });
        }

        for (row, values) in covariance.iter().enumerate() {
            if values.len() != dimension {
                return Err(MvNormalError::NonSquare {
                    row,
                    expected: dimension,
                    actual: values.len(),
                });
            }
            for (column, &value) in values.iter().enumerate() {
                if !value.is_finite() {
                    return Err(MvNormalError::NonFinite { row, column, value });
                }
            }
        }

        // A relative tolerance avoids rejecting harmless floating-point noise
        // while still detecting a genuinely asymmetric covariance matrix.
        const SYMMETRY_TOLERANCE: f64 = 1.0e-12;
        for row in 0..dimension {
            for column in (row + 1)..dimension {
                let left = covariance[row][column];
                let right = covariance[column][row];
                let scale = 1.0_f64.max(left.abs()).max(right.abs());
                if (left - right).abs() > SYMMETRY_TOLERANCE * scale {
                    return Err(MvNormalError::NotSymmetric {
                        row,
                        column,
                        left,
                        right,
                    });
                }
            }
        }

        let mut cholesky = vec![0.0; dimension * dimension];
        for row in 0..dimension {
            for column in 0..=row {
                let mut value = covariance[row][column];
                for index in 0..column {
                    value -=
                        cholesky[row * dimension + index] * cholesky[column * dimension + index];
                }
                if row == column {
                    if !value.is_finite() || value <= 0.0 {
                        return Err(MvNormalError::NotPositiveDefinite { index: row, value });
                    }
                    cholesky[row * dimension + column] = value.sqrt();
                } else {
                    cholesky[row * dimension + column] =
                        value / cholesky[column * dimension + column];
                }
            }
        }

        Ok(Self {
            mean,
            cholesky,
            dimension,
        })
    }

    /// Number of dimensions.
    pub fn dimension(&self) -> usize {
        self.dimension
    }

    /// The mean vector.
    pub fn mean(&self) -> &[f64] {
        &self.mean
    }

    /// Draw one sample, allocating a result vector and a temporary normal vector.
    pub fn sample(&self, rng: &mut StandardRng) -> Vec<f64> {
        let mut scratch = vec![0.0; self.dimension];
        let mut output = vec![0.0; self.dimension];
        self.sample_into(rng, &mut scratch, &mut output)
            .expect("internal sample buffers have the correct length");
        output
    }

    /// Draw one sample while reusing caller-provided scratch and output buffers.
    pub fn sample_into(
        &self,
        rng: &mut StandardRng,
        scratch: &mut [f64],
        output: &mut [f64],
    ) -> Result<(), SampleError> {
        if scratch.len() != self.dimension {
            return Err(SampleError::BufferLength {
                name: "scratch",
                expected: self.dimension,
                actual: scratch.len(),
            });
        }
        if output.len() != self.dimension {
            return Err(SampleError::BufferLength {
                name: "output",
                expected: self.dimension,
                actual: output.len(),
            });
        }

        for value in scratch.iter_mut() {
            *value = rng.standard_normal();
        }
        for row in 0..self.dimension {
            let mut value = self.mean[row];
            for column in 0..=row {
                value += self.cholesky[row * self.dimension + column] * scratch[column];
            }
            output[row] = value;
        }
        Ok(())
    }
}

/// A deterministic, dependency-free pseudo-random generator.
///
/// It is intended to make examples and benchmarks reproducible, not for
/// cryptographic use. `standard_normal` uses the Box--Muller transform.
#[derive(Debug, Clone)]
pub struct StandardRng {
    state: u64,
    spare_normal: Option<f64>,
}

impl StandardRng {
    /// Create a generator from a non-zero or zero seed.
    pub fn new(seed: u64) -> Self {
        Self {
            state: seed,
            spare_normal: None,
        }
    }

    fn next_u64(&mut self) -> u64 {
        // SplitMix64: compact and adequate for this numerical benchmark.
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut value = self.state;
        value = (value ^ (value >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        value = (value ^ (value >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        value ^ (value >> 31)
    }

    fn uniform_open01(&mut self) -> f64 {
        // The added half-unit keeps the Box--Muller logarithm away from zero.
        ((self.next_u64() >> 11) as f64 + 0.5) * (1.0 / 9_007_199_254_740_992.0)
    }

    /// Generate one N(0, 1) value using Box--Muller and one-value caching.
    pub fn standard_normal(&mut self) -> f64 {
        if let Some(value) = self.spare_normal.take() {
            return value;
        }
        let radius = (-2.0 * self.uniform_open01().ln()).sqrt();
        let angle = 2.0 * std::f64::consts::PI * self.uniform_open01();
        let (sin, cos) = angle.sin_cos();
        self.spare_normal = Some(radius * sin);
        radius * cos
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cholesky_sampling_is_reproducible() {
        let distribution =
            MvNormal::new(vec![1.0, -1.0], vec![vec![2.0, 0.5], vec![0.5, 1.0]]).unwrap();
        let mut first_rng = StandardRng::new(42);
        let mut second_rng = StandardRng::new(42);
        assert_eq!(
            distribution.sample(&mut first_rng),
            distribution.sample(&mut second_rng)
        );
    }

    #[test]
    fn rejects_invalid_covariance() {
        assert!(MvNormal::new(vec![0.0, 0.0], vec![vec![1.0]]).is_err());
        assert!(MvNormal::new(vec![0.0, 0.0], vec![vec![1.0, 0.2], vec![0.3, 1.0],]).is_err());
        assert!(MvNormal::new(vec![0.0, 0.0], vec![vec![1.0, 2.0], vec![2.0, 1.0],]).is_err());
    }

    #[test]
    fn sample_into_reuses_buffers() {
        let distribution =
            MvNormal::new(vec![0.0, 0.0], vec![vec![1.0, 0.0], vec![0.0, 1.0]]).unwrap();
        let mut rng = StandardRng::new(7);
        let mut scratch = [0.0; 2];
        let mut output = [0.0; 2];
        distribution
            .sample_into(&mut rng, &mut scratch, &mut output)
            .unwrap();
        assert!(output.iter().all(|value| value.is_finite()));
    }
}
