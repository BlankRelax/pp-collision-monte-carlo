#pragma once
// Thin wrapper around a 64-bit Mersenne Twister so every component draws from
// one reproducible stream (same seed -> identical output).

#include <cstdint>
#include <random>

namespace ppmc {

class Random {
public:
    explicit Random(std::uint64_t seed) : engine_(seed) {}

    // Uniform in [0, 1).
    double uniform() { return uni_(engine_); }
    // Uniform in [a, b).
    double uniform(double a, double b) { return a + (b - a) * uniform(); }
    // Gaussian with the given mean and standard deviation.
    double gauss(double mean, double sigma) { return std::normal_distribution<double>(mean, sigma)(engine_); }

private:
    std::mt19937_64 engine_;
    std::uniform_real_distribution<double> uni_{0.0, 1.0};
};

} // namespace ppmc
