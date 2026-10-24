#include "ppmc/QCDDijet.h"

#include "ppmc/Constants.h"

#include <cmath>
#include <cstdlib>
#include <vector>

namespace ppmc {

using ME = QCDDijet::MatrixElements;

// Leading-order 2->2 QCD matrix elements, e.g. Ellis, Stirling & Webber,
// "QCD and Collider Physics", Table 7.1.
double ME::qqPrime(double s, double t, double u) { return 4.0 / 9.0 * (s * s + u * u) / (t * t); }
double ME::qqSame(double s, double t, double u) {
    return 4.0 / 9.0 * ((s * s + u * u) / (t * t) + (s * s + t * t) / (u * u)) - 8.0 / 27.0 * s * s / (u * t);
}
double ME::qqbarToQpQpbar(double s, double t, double u) { return 4.0 / 9.0 * (t * t + u * u) / (s * s); }
double ME::qqbarSame(double s, double t, double u) {
    return 4.0 / 9.0 * ((s * s + u * u) / (t * t) + (t * t + u * u) / (s * s)) - 8.0 / 27.0 * u * u / (s * t);
}
double ME::qqbarToGG(double s, double t, double u) {
    return 32.0 / 27.0 * (t * t + u * u) / (t * u) - 8.0 / 3.0 * (t * t + u * u) / (s * s);
}
double ME::ggToQQbar(double s, double t, double u) {
    return 1.0 / 6.0 * (t * t + u * u) / (t * u) - 3.0 / 8.0 * (t * t + u * u) / (s * s);
}
double ME::qgToQG(double s, double t, double u) {
    return -4.0 / 9.0 * (s * s + u * u) / (s * u) + (u * u + s * s) / (t * t);
}
double ME::ggToGG(double s, double t, double u) {
    return 9.0 / 2.0 * (3.0 - t * u / (s * s) - s * u / (t * t) - s * t / (u * u));
}

QCDDijet::QCDDijet(const ToyPDF& pdf, DijetConfig cfg) : pdf_(pdf), cfg_(cfg) {}

double QCDDijet::alphaS(double q) {
    const int nf = 5;
    const double b0 = (33.0 - 2.0 * nf) / (12.0 * PI);
    return ALPHA_S_MZ / (1.0 + ALPHA_S_MZ * b0 * std::log(q * q / (MZ * MZ)));
}

double QCDDijet::generate(Random& rng, Event& ev) {
    ev.clear();
    const double s = cfg_.sqrtS * cfg_.sqrtS;

    // 1) Transverse momentum. The cross section falls roughly like 1/pT^4, so sample
    //    1/pT^2 uniformly: this puts points where the cross section actually is.
    const double pTMax = cfg_.sqrtS / 2;
    const double a = 1.0 / (cfg_.pTMin * cfg_.pTMin), b = 1.0 / (pTMax * pTMax);
    const double invPT2 = a - rng.uniform() * (a - b);
    const double pT2 = 1.0 / invPT2;
    const double pT = std::sqrt(pT2);
    const double jacPT2 = pT2 * pT2 * (a - b); // d(pT^2) = pT^4 (a - b) dr

    // 2) Rapidities of the two outgoing partons, uniform in [-yMax, yMax].
    const double y3 = rng.uniform(-cfg_.yMax, cfg_.yMax);
    const double y4 = rng.uniform(-cfg_.yMax, cfg_.yMax);
    const double jacY = (2 * cfg_.yMax) * (2 * cfg_.yMax);

    // 3) Momentum conservation fixes the incoming momentum fractions.
    const double x1 = pT / cfg_.sqrtS * (std::exp(y3) + std::exp(y4));
    const double x2 = pT / cfg_.sqrtS * (std::exp(-y3) + std::exp(-y4));
    if (x1 >= 1 || x2 >= 1) return 0; // not enough energy in the protons

    const double sh = x1 * x2 * s;
    const double th = -pT2 * (1 + std::exp(y4 - y3));
    const double uh = -pT2 * (1 + std::exp(y3 - y4));

    // d(sigma)/(dy3 dy4 dpT^2) = sum_ij x1 f_i(x1) x2 f_j(x2) * dsigma_hat/dt,
    // with dsigma_hat/dt = pi alpha_s^2 / s^2 * (|M|^2 / g^4).
    const double as = alphaS(pT);
    const double norm = PI * as * as / (sh * sh);

    // 4) Enumerate every (initial pair -> final pair) channel with its weight.
    //    Parton "3" is the one at rapidity y3. For identical outgoing partons, summing
    //    over all (y3, y4) counts each configuration twice, hence the factors 1/2.
    struct Channel { int in1, in2, out3, out4; double w; const char* label; };
    std::vector<Channel> chans;
    chans.reserve(64);
    const int partons[] = {-3, -2, -1, 1, 2, 3, pdg::gluon};
    const int nf = cfg_.nfFinal;

    for (int i : partons) {
        const double fi = pdf_.xf(i, x1);
        if (fi <= 0) continue;
        for (int j : partons) {
            const double fj = pdf_.xf(j, x2);
            if (fj <= 0) continue;
            const double lum = fi * fj * norm;
            const bool gi = i == pdg::gluon, gj = j == pdg::gluon;

            if (gi && gj) {
                chans.push_back({i, j, 21, 21, lum * 0.5 * ME::ggToGG(sh, th, uh), "gg -> gg"});
                const double wq = lum * ME::ggToQQbar(sh, th, uh); // per produced flavour
                for (int q = 1; q <= nf; ++q) chans.push_back({i, j, q, -q, wq, "gg -> q qbar"});
            } else if (!gi && gj) {
                chans.push_back({i, j, i, 21, lum * ME::qgToQG(sh, th, uh), "qg -> qg"});
            } else if (gi && !gj) {
                // g(p1) q(p2) -> q(p3) g(p4): same as qg->qg with p1 <-> p2, i.e. t <-> u.
                chans.push_back({i, j, j, 21, lum * ME::qgToQG(sh, uh, th), "qg -> qg"});
            } else if (i == j) {
                chans.push_back({i, j, i, j, lum * 0.5 * ME::qqSame(sh, th, uh), "qq -> qq"});
            } else if (i == -j) {
                chans.push_back({i, j, i, j, lum * ME::qqbarSame(sh, th, uh), "q qbar -> q qbar"});
                const double wq = lum * ME::qqbarToQpQpbar(sh, th, uh);
                const int sign = i > 0 ? 1 : -1;
                for (int q = 1; q <= nf; ++q)
                    if (q != std::abs(i)) chans.push_back({i, j, sign * q, -sign * q, wq, "q qbar -> q' qbar'"});
                chans.push_back({i, j, 21, 21, lum * 0.5 * ME::qqbarToGG(sh, th, uh), "q qbar -> gg"});
            } else {
                chans.push_back({i, j, i, j, lum * ME::qqPrime(sh, th, uh), "qq' -> qq'"});
            }
        }
    }

    double total = 0;
    for (const auto& c : chans) total += c.w;
    const double weight = total * jacPT2 * jacY * GEV2_TO_PB;
    if (!(weight > 0)) return 0;

    // 5) Choose one channel for this event, proportional to its contribution.
    double r = rng.uniform() * total;
    const Channel* chosen = &chans.back();
    for (const auto& c : chans) {
        if (r < c.w) { chosen = &c; break; }
        r -= c.w;
    }

    // 6) Four-momenta. The outgoing partons are back-to-back in the transverse plane.
    const double phi = rng.uniform(0, 2 * PI);
    const double cx = pT * std::cos(phi), cy = pT * std::sin(phi);
    const double halfSqrtS = cfg_.sqrtS / 2;
    ev.particles.push_back({chosen->in1, -1, Vec4(x1 * halfSqrtS, 0, 0, x1 * halfSqrtS)});
    ev.particles.push_back({chosen->in2, -1, Vec4(x2 * halfSqrtS, 0, 0, -x2 * halfSqrtS)});
    ev.particles.push_back({chosen->out3, 1, Vec4(pT * std::cosh(y3), cx, cy, pT * std::sinh(y3))});
    ev.particles.push_back({chosen->out4, 1, Vec4(pT * std::cosh(y4), -cx, -cy, pT * std::sinh(y4))});

    ev.weight = weight;
    ev.channel = chosen->label;
    ev.x1 = x1;
    ev.x2 = x2;
    ev.scale = pT;
    return weight;
}

} // namespace ppmc
