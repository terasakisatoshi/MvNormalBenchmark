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

struct Arguments {
    std::size_t dimension = 0;
    std::size_t samples = 0;
    std::size_t repeats = 0;
};

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
        std::vector<double> mean(arguments.dimension, 0.0);
        std::vector<double> covariance(arguments.dimension * arguments.dimension,
                                       0.0);
        // A deterministic, well-conditioned SPD covariance for comparable runs.
        for (std::size_t row = 0; row < arguments.dimension; ++row) {
            for (std::size_t column = 0; column <= row; ++column) {
                const double value =
                    (row == column) ? 1.0 + 0.01 * static_cast<double>(row + 1)
                                    : 0.01;
                covariance[row * arguments.dimension + column] = value;
                covariance[column * arguments.dimension + row] = value;
            }
        }
        const mvnormal::MvNormal distribution(
            mean, covariance, arguments.dimension);
        const auto setup_end = std::chrono::steady_clock::now();

        std::mt19937_64 rng(0x4d764e6e6f726dULL);
        std::vector<double> scratch(arguments.dimension);
        std::vector<double> output(arguments.dimension);
        std::vector<double> sample_times;
        sample_times.reserve(arguments.repeats);
        double checksum = 0.0;

        for (std::size_t repeat = 0; repeat < arguments.repeats; ++repeat) {
            const auto sample_start = std::chrono::steady_clock::now();
            for (std::size_t sample = 0; sample < arguments.samples; ++sample) {
                distribution.sample_into(rng, scratch, output);
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

        std::cout << std::setprecision(17) << "cxx," << arguments.dimension << ','
                  << arguments.samples << ',' << arguments.repeats << ','
                  << setup_seconds << ',' << average_sample_seconds << ','
                  << minimum_sample_seconds << ',' << checksum << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}
