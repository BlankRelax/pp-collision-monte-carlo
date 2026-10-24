#pragma once
// Toy parton distribution functions (PDFs) for the proton.
//
// A PDF f_i(x) is the number density of parton species i carrying a fraction x of
// the proton's momentum. Real PDFs (CT18, NNPDF, ...) are fitted to data and evolve
// with the probing scale Q. This toy uses simple fixed-scale shapes
//
//     x f(x) = N x^a (1 - x)^b
//
// whose normalisations N are fixed by the proton's quantum numbers:
//   * valence sum rules:  integral (u - ubar) dx = 2,  integral (d - dbar) dx = 1
//   * momentum sum rule:  sum_i integral x f_i(x) dx = 1

#include <ostream>

namespace ppmc {

class ToyPDF {
public:
    ToyPDF();

    // Number density f_i(x) for a PDG parton id (+-1 d, +-2 u, +-3 s, 21 g).
    // Returns 0 for anything else (no charm or bottom in this toy proton).
    double f(int id, double x) const;
    // Momentum density x f_i(x).
    double xf(int id, double x) const;

    // Numerically integrates the PDFs and prints the sum rules as a self-check.
    void printSumRules(std::ostream& os) const;
    // Numerical integral of f_i (number) or x f_i (momentum) over 0 < x < 1.
    double integrate(int id, bool momentum) const;

private:
    double uValence(double x) const; // x u_v(x)
    double dValence(double x) const; // x d_v(x)
    double sea(double x) const;      // x S(x), the same for ubar, dbar, s, sbar
    double gluon(double x) const;    // x g(x)

    double nUv_, nDv_, nSea_, nG_;
};

} // namespace ppmc
