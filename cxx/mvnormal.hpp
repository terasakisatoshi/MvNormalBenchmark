#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mvnormal {

// Multivariate normal distribution parameterized by a mean vector and a
// covariance matrix.  The covariance matrix is factorized once at
// construction time, then samples are generated as mu + L * z.
class MvNormal {
public:
    using value_type = double;

    MvNormal(std::vector<double> mean,
             std::vector<std::vector<double>> covariance,
             double symmetry_tolerance = 1e-12)
        : mean_(std::move(mean)) {
        if (!std::isfinite(symmetry_tolerance) || symmetry_tolerance < 0.0) {
            throw std::invalid_argument(
                "symmetry_tolerance must be finite and non-negative");
        }
        if (covariance.size() != mean_.size()) {
            throw std::invalid_argument(
                "covariance dimension must match mean dimension");
        }
        if (mean_.empty()) {
            throw std::invalid_argument("distribution dimension must be positive");
        }

        const std::size_t dimension = mean_.size();
        for (const double value : mean_) {
            if (!std::isfinite(value)) {
                throw std::invalid_argument("mean must contain finite values");
            }
        }

        std::vector<double> covariance_flat;
        covariance_flat.reserve(checked_square_size(dimension));
        for (std::size_t row = 0; row < dimension; ++row) {
            if (covariance[row].size() != dimension) {
                throw std::invalid_argument(
                    "covariance matrix must be square");
            }
            covariance_flat.insert(covariance_flat.end(), covariance[row].begin(),
                                   covariance[row].end());
        }
        initialize_from_flat(covariance_flat, dimension, symmetry_tolerance);
    }

    // The flat covariance is row-major: covariance[row * dimension + column].
    MvNormal(std::vector<double> mean,
             std::vector<double> covariance,
             std::size_t dimension,
             double symmetry_tolerance = 1e-12)
        : mean_(std::move(mean)) {
        if (dimension != mean_.size()) {
            throw std::invalid_argument(
                "covariance dimension must match mean dimension");
        }
        if (dimension == 0) {
            throw std::invalid_argument("distribution dimension must be positive");
        }
        initialize_from_flat(covariance, dimension, symmetry_tolerance);
    }

    [[nodiscard]] std::size_t dimension() const noexcept { return mean_.size(); }
    [[nodiscard]] const std::vector<double>& mean() const noexcept { return mean_; }
    [[nodiscard]] const std::vector<double>& cholesky_factor() const noexcept {
        return cholesky_factor_;
    }

    template <class RandomNumberGenerator>
    [[nodiscard]] std::vector<double> sample(RandomNumberGenerator& rng) const {
        std::vector<double> scratch(dimension());
        std::vector<double> output(dimension());
        sample_into(rng, scratch, output);
        return output;
    }

    // Fills scratch with standard normal variates and output with the sample.
    // Both vectors must be pre-sized to dimension(); their storage is reused.
    template <class RandomNumberGenerator>
    void sample_into(RandomNumberGenerator& rng,
                     std::vector<double>& scratch,
                     std::vector<double>& output) const {
        const std::size_t dimension_value = dimension();
        if (scratch.size() != dimension_value || output.size() != dimension_value) {
            throw std::invalid_argument(
                "scratch and output must have size equal to distribution dimension");
        }
        if (&scratch == &output) {
            throw std::invalid_argument("scratch and output must be distinct vectors");
        }

        std::normal_distribution<double> standard_normal(0.0, 1.0);
        for (double& value : scratch) {
            value = standard_normal(rng);
        }

        for (std::size_t row = 0; row < dimension_value; ++row) {
            double value = mean_[row];
            for (std::size_t column = 0; column <= row; ++column) {
                value += cholesky_factor_[row * dimension_value + column] *
                         scratch[column];
            }
            output[row] = value;
        }
    }

private:
    [[nodiscard]] static std::size_t checked_square_size(std::size_t dimension) {
        if (dimension != 0 &&
            dimension > std::numeric_limits<std::size_t>::max() / dimension) {
            throw std::invalid_argument("distribution dimension is too large");
        }
        return dimension * dimension;
    }

    void initialize_from_flat(const std::vector<double>& covariance,
                              std::size_t dimension,
                              double symmetry_tolerance) {
        if (!std::isfinite(symmetry_tolerance) || symmetry_tolerance < 0.0) {
            throw std::invalid_argument(
                "symmetry_tolerance must be finite and non-negative");
        }
        const std::size_t matrix_size = checked_square_size(dimension);
        if (covariance.size() != matrix_size) {
            throw std::invalid_argument(
                "covariance matrix must be square and match its dimension");
        }
        for (const double value : mean_) {
            if (!std::isfinite(value)) {
                throw std::invalid_argument("mean must contain finite values");
            }
        }
        for (const double value : covariance) {
            if (!std::isfinite(value)) {
                throw std::invalid_argument(
                    "covariance must contain finite values");
            }
        }

        for (std::size_t row = 0; row < dimension; ++row) {
            for (std::size_t column = row + 1; column < dimension; ++column) {
                const double a = covariance[row * dimension + column];
                const double b = covariance[column * dimension + row];
                const double scale = std::max({1.0, std::abs(a), std::abs(b)});
                if (std::abs(a - b) > symmetry_tolerance * scale) {
                    throw std::invalid_argument(
                        "covariance matrix must be symmetric");
                }
            }
        }

        cholesky_factor_.assign(matrix_size, 0.0);
        for (std::size_t row = 0; row < dimension; ++row) {
            for (std::size_t column = 0; column <= row; ++column) {
                double value = covariance[row * dimension + column];
                for (std::size_t k = 0; k < column; ++k) {
                    value -= cholesky_factor_[row * dimension + k] *
                              cholesky_factor_[column * dimension + k];
                }

                if (row == column) {
                    if (!(value > 0.0) || !std::isfinite(value)) {
                        throw std::invalid_argument(
                            "covariance matrix must be positive definite");
                    }
                    cholesky_factor_[row * dimension + column] = std::sqrt(value);
                } else {
                    const double pivot =
                        cholesky_factor_[column * dimension + column];
                    if (!(pivot > 0.0) || !std::isfinite(pivot)) {
                        throw std::invalid_argument(
                            "covariance matrix must be positive definite");
                    }
                    cholesky_factor_[row * dimension + column] = value / pivot;
                }
            }
        }
    }

    std::vector<double> mean_;
    std::vector<double> cholesky_factor_;
};

}  // namespace mvnormal
