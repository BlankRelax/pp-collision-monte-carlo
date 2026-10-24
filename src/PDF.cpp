#include "ppmc/PDF.h"

#include "ppmc/Constants.h"

#include <cmath>
#include <cstdlib>
#include <iomanip>

namespace ppmc {

namespace {

// Shape exponents for x f(x) = N x^a (1-x)^b, chosen to roughly resemble real
// fitted PDFs at a scale Q ~ 100 GeV (where W, Z and jet production happen).
constexpr double UV_A = 0.5, UV_B = 3.0;
constexpr double DV_A = 0.5, DV_B = 4.0;
constexpr double SEA_A = -0.3, SEA_B = 8.0;
constexpr double G_A = -0.3, G_B = 5.0;

// Momentum carried by EACH of the six sea distributions (u, d, s quarks and antiquarks).
constexpr double SEA_MOMENTUM_EACH = 0.035;

// Euler beta function B(a, b) = integral_0^1 x^(a-1) (1-x)^(b-1) dx.
double betaFn(double a, double b) { return std::exp(std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b)); }

} // namespace

ToyPDF::ToyPDF() {
    // integral u_v dx = N B(a, b+1) where u_v = N x^(a-1) (1-x)^b.
    nUv_ = 2.0 / betaFn(UV_A, UV_B + 1);
    nDv_ = 1.0 / betaFn(DV_A, DV_B + 1);
    // integral x S dx = N B(a+1, b+1).
    nSea_ = SEA_MOMENTUM_EACH / betaFn(SEA_A + 1, SEA_B + 1);
    // The gluon takes whatever momentum the quarks leave over.
    const double valenceMomentum = nUv_ * betaFn(UV_A + 1, UV_B + 1) + nDv_ * betaFn(DV_A + 1, DV_B + 1);
    const double gluonMomentum = 1.0 - valenceMomentum - 6 * SEA_MOMENTUM_EACH;
    nG_ = gluonMomentum / betaFn(G_A + 1, G_B + 1);
}

double ToyPDF::uValence(double x) const { return nUv_ * std::pow(x, UV_A) * std::pow(1 - x, UV_B); }
double ToyPDF::dValence(double x) const { return nDv_ * std::pow(x, DV_A) * std::pow(1 - x, DV_B); }
double ToyPDF::sea(double x) const { return nSea_ * std::pow(x, SEA_A) * std::pow(1 - x, SEA_B); }
double ToyPDF::gluon(double x) const { return nG_ * std::pow(x, G_A) * std::pow(1 - x, G_B); }

double ToyPDF::xf(int id, double x) const {
    if (x <= 0 || x >= 1) return 0;
    switch (id) {
    case pdg::u: return uValence(x) + sea(x);
    case pdg::d: return dValence(x) + sea(x);
    case -pdg::u:
    case -pdg::d:
    case pdg::s:
    case -pdg::s: return sea(x);
    case pdg::gluon: return gluon(x);
    default: return 0;
    }
}

double ToyPDF::f(int id, double x) const { return x > 0 ? xf(id, x) / x : 0; }

double ToyPDF::integrate(int id, bool momentum) const {
    // Substitute x = e^t so the small-x region is sampled densely; then
    // dx = x dt and the integrand becomes x f(x) (number) or x^2 f(x) (momentum).
    const double tMin = std::log(1e-12);
    const int n = 200000; // even, for Simpson's rule
    const double h = -tMin / n;
    double sum = 0;
    for (int k = 0; k <= n; ++k) {
        const double x = std::exp(tMin + k * h);
        const double val = momentum ? x * xf(id, x) : xf(id, x);
        const double coeff = (k == 0 || k == n) ? 1 : (k % 2 ? 4 : 2);
        sum += coeff * val;
    }
    return sum * h / 3;
}

void ToyPDF::printSumRules(std::ostream& os) const {
    const auto num = [&](int id) { return integrate(id, false); };
    const auto mom = [&](int id) { return integrate(id, true); };
    const int ids[] = {pdg::u, pdg::d, -pdg::u, -pdg::d, pdg::s, -pdg::s, pdg::gluon};
    const char* names[] = {"u", "d", "ubar", "dbar", "s", "sbar", "g"};

    os << std::fixed << std::setprecision(4);
    os << "Toy proton PDF sum rules (numerical integration)\n";
    os << "  integral (u - ubar) dx = " << num(pdg::u) - num(-pdg::u) << "   (expected 2)\n";
    os << "  integral (d - dbar) dx = " << num(pdg::d) - num(-pdg::d) << "   (expected 1)\n";
    os << "  integral (s - sbar) dx = " << num(pdg::s) - num(-pdg::s) << "   (expected 0)\n";
    os << "\n  Momentum fraction carried by each parton species:\n";
    double total = 0;
    for (int i = 0; i < 7; ++i) {
        const double m = mom(ids[i]);
        total += m;
        os << "    " << std::setw(5) << names[i] << " : " << m << "\n";
    }
    os << "    total : " << total << "   (expected 1)\n";
}

} // namespace ppmc
