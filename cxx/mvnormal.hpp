#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mvnormal {

enum class NormalAlgorithm {
    MarsagliaPolar,
    Ziggurat,
};

class NormalRng {
public:
    explicit NormalRng(std::uint64_t seed,
                       NormalAlgorithm algorithm = NormalAlgorithm::MarsagliaPolar)
        : state_(seed == 0 ? 0x9E37'79B9'7F4A'7C15ULL : seed),
          algorithm_(algorithm) {
        if (algorithm_ == NormalAlgorithm::Ziggurat) {
            (void)tables();
        }
    }

    [[nodiscard]] NormalAlgorithm algorithm() const noexcept { return algorithm_; }

    double standard_normal() {
        if (algorithm_ == NormalAlgorithm::Ziggurat) {
            return standard_normal_ziggurat();
        }
        if (spare_normal_.has_value()) {
            const double value = *spare_normal_;
            spare_normal_.reset();
            return value;
        }
        const auto [first, second] = polar_pair();
        spare_normal_ = second;
        return first;
    }

private:
    static constexpr std::size_t layers = 256;
    static constexpr double ziggurat_r = 3.6541528853610088;
    static constexpr double ziggurat_v = 0.004928673233974658;
    static constexpr double two_pow_63 = 9.223372036854775808e18;
    static constexpr std::uint64_t magnitude_mask = 0x7FFF'FFFF'FFFF'FFFFULL;

    struct Tables {
        std::array<std::uint64_t, layers> k{};
        std::array<double, layers> w{};
        std::array<double, layers> f{};

        Tables() {
            std::array<double, layers> x{};
            const double tail_density = std::exp(-0.5 * ziggurat_r * ziggurat_r);
            const double q = ziggurat_v / tail_density;
            x[layers - 1] = ziggurat_r;
            f[0] = 1.0;
            f[layers - 1] = tail_density;
            w[0] = q / two_pow_63;
            w[layers - 1] = ziggurat_r / two_pow_63;
            k[0] = static_cast<std::uint64_t>(ziggurat_r / q * two_pow_63);
            k[1] = 0;
            for (std::size_t index = layers - 2; index > 0; --index) {
                x[index] = std::sqrt(
                    -2.0 * std::log(ziggurat_v / x[index + 1] + f[index + 1]));
                f[index] = std::exp(-0.5 * x[index] * x[index]);
                w[index] = x[index] / two_pow_63;
                k[index + 1] = static_cast<std::uint64_t>(
                    x[index] / x[index + 1] * two_pow_63);
            }
        }
    };

    static const Tables& tables() {
        static const Tables value;
        return value;
    }

    std::uint64_t next_u64() {
        std::uint64_t value = state_;
        value ^= value << 13;
        value ^= value >> 7;
        value ^= value << 17;
        state_ = value;
        return value;
    }

    double uniform_open01() {
        return (static_cast<double>(next_u64() >> 11) + 0.5) /
               9.007199254740992e15;
    }

    std::pair<double, double> polar_pair() {
        for (;;) {
            const double u = 2.0 * uniform_open01() - 1.0;
            const double v = 2.0 * uniform_open01() - 1.0;
            const double radius_squared = u * u + v * v;
            if (radius_squared > 0.0 && radius_squared < 1.0) {
                const double scale =
                    std::sqrt(-2.0 * std::log(radius_squared) / radius_squared);
                return {u * scale, v * scale};
            }
        }
    }

    double standard_normal_ziggurat() {
        const Tables& table = tables();
        for (;;) {
            const std::uint64_t bits = next_u64();
            const std::size_t index = bits & (layers - 1);
            const double sign = (bits & (std::uint64_t{1} << 63)) == 0 ? 1.0 : -1.0;
            const std::uint64_t magnitude = bits & magnitude_mask;
            const double x = static_cast<double>(magnitude) * table.w[index];
            if (magnitude < table.k[index]) {
                return sign * x;
            }
            if (index == 0) {
                for (;;) {
                    const double tail_x = -std::log(uniform_open01()) / ziggurat_r;
                    const double tail_y = -std::log(uniform_open01());
                    if (2.0 * tail_y >= tail_x * tail_x) {
                        return sign * (ziggurat_r + tail_x);
                    }
                }
            }
            const double y = table.f[index] +
                             uniform_open01() * (table.f[index - 1] - table.f[index]);
            if (y < std::exp(-0.5 * x * x)) {
                return sign * x;
            }
        }
    }

    std::uint64_t state_;
    NormalAlgorithm algorithm_;
    std::optional<double> spare_normal_;
};

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
        std::vector<double> output(dimension());
        sample_inplace(rng, output);
        return output;
    }

    [[nodiscard]] std::vector<double> sample(NormalRng& rng) const {
        std::vector<double> output(dimension());
        sample_inplace(rng, output);
        return output;
    }

    // Generate z directly in output and transform it in place.  Rows are
    // processed from bottom to top so the still-needed z values remain intact.
    // The row-major layout makes each inner coefficient access contiguous.
    template <class RandomNumberGenerator>
    void sample_inplace(RandomNumberGenerator& rng,
                        std::vector<double>& output) const {
        const std::size_t dimension_value = dimension();
        if (output.size() != dimension_value) {
            throw std::invalid_argument(
                "output must have size equal to distribution dimension");
        }

        sample_inplace_generated(
            [&]() { return standard_normal_(rng); }, output);
    }

    void sample_inplace(NormalRng& rng, std::vector<double>& output) const {
        if (output.size() != dimension()) {
            throw std::invalid_argument(
                "output must have size equal to distribution dimension");
        }
        sample_inplace_generated([&]() { return rng.standard_normal(); }, output);
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

        for (double& value : scratch) {
            value = standard_normal_(rng);
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
    template <class GenerateNormal>
    void sample_inplace_generated(GenerateNormal&& generate_normal,
                                  std::vector<double>& output) const {
        const std::size_t dimension_value = dimension();
        for (double& value : output) {
            value = generate_normal();
        }

        for (std::size_t row = dimension_value; row > 0; --row) {
            const std::size_t row_index = row - 1;
            double value = mean_[row_index];
            const std::size_t row_start = row_index * dimension_value;
            for (std::size_t column = 0; column <= row_index; ++column) {
                value += cholesky_factor_[row_start + column] * output[column];
            }
            output[row_index] = value;
        }
    }

    // Reuse the distribution object, including any cached normal value.
    mutable std::normal_distribution<double> standard_normal_{0.0, 1.0};

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
