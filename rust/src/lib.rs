//! A small, dependency-free multivariate normal sampler.
//!
//! Sampling follows the usual Cholesky construction used by
//! `Distributions.jl`'s `MvNormal`: for `Sigma = L * L'` and a vector `z`
//! of independent standard normal values, return `mu + L * z`.

use std::error::Error;
use std::fmt;
use std::sync::OnceLock;

const ZIGGURAT_LAYERS: usize = 256;
const ZIGGURAT_R: f64 = 3.654_152_885_361_008_8;
const ZIGGURAT_V: f64 = 0.004_928_673_233_974_658;
const TWO_POW_63: f64 = 9.223_372_036_854_775_808e18;
const SIGNED_MAGNITUDE_MASK: u64 = 0x7fff_ffff_ffff_ffff;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NormalAlgorithm {
    MarsagliaPolar,
    Ziggurat,
}

impl NormalAlgorithm {
    pub fn label(self) -> &'static str {
        match self {
            Self::MarsagliaPolar => "polar",
            Self::Ziggurat => "ziggurat",
        }
    }
}

struct ZigguratTables {
    k: [u64; ZIGGURAT_LAYERS],
    w: [f64; ZIGGURAT_LAYERS],
    f: [f64; ZIGGURAT_LAYERS],
}

impl ZigguratTables {
    fn new() -> Self {
        let mut x = [0.0; ZIGGURAT_LAYERS];
        let mut k = [0; ZIGGURAT_LAYERS];
        let mut w = [0.0; ZIGGURAT_LAYERS];
        let mut f = [0.0; ZIGGURAT_LAYERS];
        let tail_density = (-0.5 * ZIGGURAT_R * ZIGGURAT_R).exp();
        let q = ZIGGURAT_V / tail_density;

        x[ZIGGURAT_LAYERS - 1] = ZIGGURAT_R;
        f[0] = 1.0;
        f[ZIGGURAT_LAYERS - 1] = tail_density;
        w[0] = q / TWO_POW_63;
        w[ZIGGURAT_LAYERS - 1] = ZIGGURAT_R / TWO_POW_63;
        k[0] = (ZIGGURAT_R / q * TWO_POW_63) as u64;
        k[1] = 0;

        for index in (1..ZIGGURAT_LAYERS - 1).rev() {
            x[index] = (-2.0 * (ZIGGURAT_V / x[index + 1] + f[index + 1]).ln()).sqrt();
            f[index] = (-0.5 * x[index] * x[index]).exp();
            w[index] = x[index] / TWO_POW_63;
            k[index + 1] = (x[index] / x[index + 1] * TWO_POW_63) as u64;
        }

        Self { k, w, f }
    }
}

fn ziggurat_tables() -> &'static ZigguratTables {
    static TABLES: OnceLock<ZigguratTables> = OnceLock::new();
    TABLES.get_or_init(ZigguratTables::new)
}

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

    /// Draw one sample, allocating a result vector.
    pub fn sample(&self, rng: &mut StandardRng) -> Vec<f64> {
        let mut output = vec![0.0; self.dimension];
        self.sample_inplace(rng, &mut output)
            .expect("internal sample buffers have the correct length");
        output
    }

    /// Draw one sample directly into an output buffer.
    ///
    /// Standard normals are generated in `output`, then rows are processed
    /// from bottom to top so the still-needed values remain untouched.
    pub fn sample_inplace(
        &self,
        rng: &mut StandardRng,
        output: &mut [f64],
    ) -> Result<(), SampleError> {
        if output.len() != self.dimension {
            return Err(SampleError::BufferLength {
                name: "output",
                expected: self.dimension,
                actual: output.len(),
            });
        }

        rng.fill_standard_normals(output);
        for row in (0..self.dimension).rev() {
            let row_start = row * self.dimension;
            let row_end = row_start + row + 1;
            let row_factor = &self.cholesky[row_start..row_end];
            let mut value = self.mean[row];
            for (factor, normal) in row_factor.iter().zip(&output[..=row]) {
                value += factor * normal;
            }
            output[row] = value;
        }
        Ok(())
    }

    /// Draw one sample using an external standard-normal generator.
    ///
    /// This lets the benchmark plug in `rand_distr`'s Ziggurat-based standard
    /// normal while keeping the same `mean + L * z` transformation and the
    /// bottom-up in-place update used by the other samplers.
    pub fn sample_inplace_with<G>(
        &self,
        mut generate_normal: G,
        output: &mut [f64],
    ) -> Result<(), SampleError>
    where
        G: FnMut() -> f64,
    {
        if output.len() != self.dimension {
            return Err(SampleError::BufferLength {
                name: "output",
                expected: self.dimension,
                actual: output.len(),
            });
        }

        for value in output.iter_mut() {
            *value = generate_normal();
        }

        for row in (0..self.dimension).rev() {
            let row_start = row * self.dimension;
            let mut value = self.mean[row];
            for column in 0..=row {
                value += self.cholesky[row_start + column] * output[column];
            }
            output[row] = value;
        }
        Ok(())
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

        rng.fill_standard_normals(scratch);
        for (row, (mean, output_value)) in self.mean.iter().zip(output.iter_mut()).enumerate() {
            let row_start = row * self.dimension;
            let row_end = row_start + row + 1;
            let row_factor = &self.cholesky[row_start..row_end];
            let mut value = *mean;
            for (factor, normal) in row_factor.iter().zip(&scratch[..=row]) {
                value += factor * normal;
            }
            *output_value = value;
        }
        Ok(())
    }
}

/// A deterministic, dependency-free pseudo-random generator.
///
/// It is intended to make examples and benchmarks reproducible, not for
/// cryptographic use. Standard normal values use the Marsaglia polar method.
#[derive(Debug, Clone)]
pub struct StandardRng {
    state: u64,
    spare_normal: Option<f64>,
    algorithm: NormalAlgorithm,
}

impl StandardRng {
    /// Create a generator from a non-zero or zero seed.
    pub fn new(seed: u64) -> Self {
        Self::with_algorithm(seed, NormalAlgorithm::MarsagliaPolar)
    }

    pub fn with_algorithm(seed: u64, algorithm: NormalAlgorithm) -> Self {
        if algorithm == NormalAlgorithm::Ziggurat {
            let _ = ziggurat_tables();
        }
        Self {
            state: if seed == 0 {
                0x9E37_79B9_7F4A_7C15
            } else {
                seed
            },
            spare_normal: None,
            algorithm,
        }
    }

    pub fn algorithm(&self) -> NormalAlgorithm {
        self.algorithm
    }

    #[inline]
    fn next_u64(&mut self) -> u64 {
        // Xorshift64: a compact, fast generator for this numerical benchmark.
        let mut value = self.state;
        value ^= value << 13;
        value ^= value >> 7;
        value ^= value << 17;
        self.state = value;
        value
    }

    #[inline]
    fn uniform_open01(&mut self) -> f64 {
        // The added half-unit keeps the polar transform away from exact bounds.
        ((self.next_u64() >> 11) as f64 + 0.5) * (1.0 / 9_007_199_254_740_992.0)
    }

    #[inline]
    fn standard_normal_pair(&mut self) -> (f64, f64) {
        loop {
            let u = 2.0 * self.uniform_open01() - 1.0;
            let v = 2.0 * self.uniform_open01() - 1.0;
            let radius_squared = u * u + v * v;
            if radius_squared > 0.0 && radius_squared < 1.0 {
                let scale = (-2.0 * radius_squared.ln() / radius_squared).sqrt();
                return (u * scale, v * scale);
            }
        }
    }

    #[inline]
    fn standard_normal_ziggurat(&mut self) -> f64 {
        let tables = ziggurat_tables();
        loop {
            let bits = self.next_u64();
            let index = (bits & (ZIGGURAT_LAYERS as u64 - 1)) as usize;
            let sign = if bits & (1_u64 << 63) == 0 { 1.0 } else { -1.0 };
            let magnitude = bits & SIGNED_MAGNITUDE_MASK;
            let x = magnitude as f64 * tables.w[index];
            if magnitude < tables.k[index] {
                return sign * x;
            }

            if index == 0 {
                loop {
                    let tail_x = -self.uniform_open01().ln() / ZIGGURAT_R;
                    let tail_y = -self.uniform_open01().ln();
                    if 2.0 * tail_y >= tail_x * tail_x {
                        return sign * (ZIGGURAT_R + tail_x);
                    }
                }
            }

            let y =
                tables.f[index] + self.uniform_open01() * (tables.f[index - 1] - tables.f[index]);
            if y < (-0.5 * x * x).exp() {
                return sign * x;
            }
        }
    }

    /// Fill a slice with independent N(0, 1) values.
    #[inline]
    pub fn fill_standard_normals(&mut self, output: &mut [f64]) {
        if self.algorithm == NormalAlgorithm::Ziggurat {
            for value in output.iter_mut() {
                *value = self.standard_normal_ziggurat();
            }
            return;
        }

        let mut index = 0;
        if let Some(value) = self.spare_normal.take() {
            if let Some(first) = output.first_mut() {
                *first = value;
                index = 1;
            }
        }

        while index + 1 < output.len() {
            let (first, second) = self.standard_normal_pair();
            output[index] = first;
            output[index + 1] = second;
            index += 2;
        }
        if index < output.len() {
            output[index] = self.standard_normal();
        }
    }

    /// Generate one N(0, 1) value using the Marsaglia polar method.
    #[inline]
    pub fn standard_normal(&mut self) -> f64 {
        if self.algorithm == NormalAlgorithm::Ziggurat {
            return self.standard_normal_ziggurat();
        }
        if let Some(value) = self.spare_normal.take() {
            return value;
        }
        let (first, second) = self.standard_normal_pair();
        self.spare_normal = Some(second);
        first
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

    #[test]
    fn sample_inplace_matches_buffered_sampling() {
        let distribution = MvNormal::new(
            vec![0.5, -1.0, 2.0],
            vec![
                vec![2.0, 0.1, 0.2],
                vec![0.1, 1.5, 0.0],
                vec![0.2, 0.0, 1.2],
            ],
        )
        .unwrap();
        let mut buffered_rng = StandardRng::new(123);
        let mut inplace_rng = StandardRng::new(123);
        let mut scratch = vec![0.0; 3];
        let mut buffered_output = vec![0.0; 3];
        let mut inplace_output = vec![0.0; 3];

        distribution
            .sample_into(&mut buffered_rng, &mut scratch, &mut buffered_output)
            .unwrap();
        distribution
            .sample_inplace(&mut inplace_rng, &mut inplace_output)
            .unwrap();

        assert_eq!(buffered_output, inplace_output);
    }

    #[test]
    fn sample_inplace_with_external_generator() {
        let distribution = MvNormal::new(
            vec![0.5, -1.0],
            vec![vec![1.0, 0.0], vec![0.0, 1.0]],
        )
        .unwrap();
        let mut output = [0.0; 2];
        distribution
            .sample_inplace_with(|| 0.0, &mut output)
            .unwrap();
        assert_eq!(output, [0.5, -1.0]);
    }

    #[test]
    fn ziggurat_sampling_is_finite() {
        let mut rng = StandardRng::with_algorithm(123, NormalAlgorithm::Ziggurat);
        let mut values = [0.0; 257];
        rng.fill_standard_normals(&mut values);
        assert!(values.iter().all(|value| value.is_finite()));
    }
}
