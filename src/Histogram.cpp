#include "ppmc/Histogram.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <utility>

namespace ppmc {

Histogram::Histogram(std::string name, std::string xLabel, int nBins, double lo, double hi)
    : name_(std::move(name)), xLabel_(std::move(xLabel)), n_(nBins), lo_(lo), hi_(hi),
      width_((hi - lo) / nBins), sumW_(nBins, 0.0), sumW2_(nBins, 0.0) {}

void Histogram::fill(double x, double w) {
    if (x < lo_) { under_ += w; return; }
    if (x >= hi_) { over_ += w; return; }
    const int k = std::min(n_ - 1, static_cast<int>((x - lo_) / width_));
    sumW_[k] += w;
    sumW2_[k] += w * w;
}

void Histogram::scale(double f) {
    for (int k = 0; k < n_; ++k) {
        sumW_[k] *= f;
        sumW2_[k] *= f * f;
    }
    under_ *= f;
    over_ *= f;
}

double Histogram::integral() const {
    double total = 0;
    for (double w : sumW_) total += w;
    return total;
}

void Histogram::writeCSV(const std::string& path, const std::string& yLabel) const {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot open " + path);
    out << "# " << name_ << "\n# x: " << xLabel_ << "\n# y: " << yLabel << "\n";
    out << "x_low,x_high,value,error\n";
    out << std::setprecision(8);
    for (int k = 0; k < n_; ++k)
        out << lo_ + k * width_ << "," << lo_ + (k + 1) * width_ << "," << sumW_[k] / width_ << ","
            << std::sqrt(sumW2_[k]) / width_ << "\n";
}

void Histogram::printASCII(std::ostream& os, const std::string& yLabel, int width) const {
    double maxVal = 0;
    for (double w : sumW_) maxVal = std::max(maxVal, w / width_);
    os << "\n  " << name_ << "   [y = " << yLabel << "]\n";
    os << "  " << std::string(width + 30, '-') << "\n";
    for (int k = 0; k < n_; ++k) {
        const double val = sumW_[k] / width_;
        const int len = maxVal > 0 ? static_cast<int>(std::lround(width * val / maxVal)) : 0;
        os << "  " << std::setw(9) << std::fixed << std::setprecision(2) << lo_ + k * width_ << " | "
           << std::left << std::setw(width) << std::string(len, '#') << std::right << " "
           << std::scientific << std::setprecision(3) << val << "\n";
    }
    os << std::defaultfloat;
    os << "  " << std::string(width + 30, '-') << "\n";
    os << "  x = " << xLabel_ << "\n";
}

} // namespace ppmc
