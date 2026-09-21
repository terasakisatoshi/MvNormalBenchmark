#include "mvnormal.hpp"

#include <chrono>
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <vector>

namespace {

constexpr std::uint64_t comparison_normal_seed = 0x5EED2021ULL;

struct Arguments {
    std::size_t dimension = 0;
    std::size_t samples = 0;
    std::size_t repeats = 0;
    mvnormal::NormalAlgorithm normal_algorithm =
        mvnormal::NormalAlgorithm::MarsagliaPolar;
};

[[nodiscard]] const char* normal_algorithm_label(mvnormal::NormalAlgorithm algorithm) {
    return algorithm == mvnormal::NormalAlgorithm::MarsagliaPolar ? "polar" : "ziggurat";
}

[[nodiscard]] std::size_t parse_positive_size(const std::string& text,
                                              const char* option) {
    if (text.empty() || text.front() == '-') {
        throw std::invalid_argument(std::string(option) + " must be positive");
    }
    std::size_t consumed = 0;
    const unsigned long long parsed = std::stoull(text, &consumed);
    if (consumed != text.size() || parsed == 0 ||
        parsed > std::numeric_limits<std::size_t>::max()) {
        throw std::invalid_argument(std::string(option) + " must be positive");
    }
    return static_cast<std::size_t>(parsed);
}

[[nodiscard]] Arguments parse_arguments(int argc, char** argv) {
    Arguments arguments;
    for (int index = 1; index < argc; index += 2) {
        if (index + 1 >= argc) {
            throw std::invalid_argument("each option requires a value");
        }
        const std::string option(argv[index]);
        const std::string value(argv[index + 1]);
        if (option == "--dim") {
            arguments.dimension = parse_positive_size(value, "--dim");
        } else if (option == "--samples") {
            arguments.samples = parse_positive_size(value, "--samples");
        } else if (option == "--repeats") {
            arguments.repeats = parse_positive_size(value, "--repeats");
        } else if (option == "--normal") {
            if (value == "polar") {
                arguments.normal_algorithm = mvnormal::NormalAlgorithm::MarsagliaPolar;
            } else if (value == "ziggurat") {
                arguments.normal_algorithm = mvnormal::NormalAlgorithm::Ziggurat;
            } else {
                throw std::invalid_argument(
                    "--normal must be polar or ziggurat");
            }
        } else {
            throw std::invalid_argument("unknown option: " + option);
        }
    }
    if (arguments.dimension == 0 || arguments.samples == 0 ||
        arguments.repeats == 0) {
        throw std::invalid_argument(
            "required options: --dim, --samples, and --repeats");
    }
    return arguments;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Arguments arguments = parse_arguments(argc, argv);

        const auto setup_start = std::chrono::steady_clock::now();
        std::vector<double> mean(arguments.dimension);
        std::vector<double> covariance(arguments.dimension * arguments.dimension,
                                       0.0);
        // Shared benchmark case: mu[i] = 0.01*i and
        // Sigma[i,j] = 0.25^abs(i-j). This Toeplitz matrix is SPD.
        for (std::size_t row = 0; row < arguments.dimension; ++row) {
            mean[row] = 0.01 * static_cast<double>(row);
            for (std::size_t column = 0; column < arguments.dimension; ++column) {
                const std::size_t distance = row > column ? row - column : column - row;
                const double value =
                    std::ldexp(1.0, -2 * static_cast<int>(distance));
                covariance[row * arguments.dimension + column] = value;
            }
        }
        const mvnormal::MvNormal distribution(
            mean, covariance, arguments.dimension);
        const auto setup_end = std::chrono::steady_clock::now();

        mvnormal::NormalRng rng(comparison_normal_seed, arguments.normal_algorithm);
        std::vector<double> output(arguments.dimension);
        std::vector<double> sample_times;
        sample_times.reserve(arguments.repeats);
        double checksum = 0.0;

        for (std::size_t repeat = 0; repeat < arguments.repeats; ++repeat) {
            const auto sample_start = std::chrono::steady_clock::now();
            for (std::size_t sample = 0; sample < arguments.samples; ++sample) {
                distribution.sample_inplace(rng, output);
                for (const double value : output) {
                    checksum += value;
                }
            }
            const auto sample_end = std::chrono::steady_clock::now();
            sample_times.push_back(
                std::chrono::duration<double>(sample_end - sample_start).count());
        }

        double total_sample_seconds = 0.0;
        double minimum_sample_seconds = std::numeric_limits<double>::infinity();
        for (const double seconds : sample_times) {
            total_sample_seconds += seconds;
            minimum_sample_seconds = std::min(minimum_sample_seconds, seconds);
        }
        const double average_sample_seconds =
            total_sample_seconds / static_cast<double>(sample_times.size());
        const double setup_seconds =
            std::chrono::duration<double>(setup_end - setup_start).count();

        std::cout << std::setprecision(17) << "cxx-"
                  << normal_algorithm_label(arguments.normal_algorithm) << ','
                  << arguments.dimension << ','
                  << arguments.samples << ',' << arguments.repeats << ','
                  << setup_seconds << ',' << average_sample_seconds << ','
                  << minimum_sample_seconds << ',' << checksum << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}
