#pragma once
// A weighted 1D histogram that can print itself as text and save itself as CSV.

#include <ostream>
#include <string>
#include <vector>

namespace ppmc {

class Histogram {
public:
    Histogram(std::string name, std::string xLabel, int nBins, double lo, double hi);

    // Adds weight w to the bin containing x (values outside [lo, hi) are counted as
    // underflow/overflow and do not appear in the plot).
    void fill(double x, double w);
    // Multiplies every bin by f (used to divide by the number of generated events).
    void scale(double f);

    // Sum of all in-range bin contents.
    double integral() const;
    const std::string& name() const { return name_; }

    // Writes columns x_low, x_high, value, error where value is the content divided by
    // the bin width, i.e. a differential cross section d(sigma)/dx.
    void writeCSV(const std::string& path, const std::string& yLabel) const;
    // Prints a horizontal bar chart of the differential distribution.
    void printASCII(std::ostream& os, const std::string& yLabel, int width = 50) const;

private:
    std::string name_, xLabel_;
    int n_;
    double lo_, hi_, width_;
    std::vector<double> sumW_, sumW2_;
    double under_ = 0, over_ = 0;
};

} // namespace ppmc
