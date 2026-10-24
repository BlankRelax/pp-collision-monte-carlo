#include "ppmc/DrellYan.h"

#include "ppmc/Constants.h"

#include <cmath>
#include <complex>

namespace ppmc {

namespace {

struct Charges {
    double q;  // electric charge in units of e
    double t3; // third component of weak isospin of the left-handed field
};

Charges chargesOf(int id) {
    switch (std::abs(id)) {
    case pdg::u:
    case pdg::c: return {2.0 / 3.0, 0.5};
    case pdg::d:
    case pdg::s:
    case pdg::b: return {-1.0 / 3.0, -0.5};
    default: return {-1.0, -0.5}; // charged leptons
    }
}

} // namespace

DrellYan::DrellYan(const ToyPDF& pdf, DrellYanConfig cfg) : pdf_(pdf), cfg_(cfg) {}

double DrellYan::dSigmaDCosTheta(int quarkId, double shat, double cosTheta) {
    const double s2 = SIN2_THETA_W, c2 = 1 - s2;
    const Charges q = chargesOf(quarkId), l = chargesOf(pdg::muon);

    // Chiral Z couplings g_L = T3 - Q sin^2, g_R = -Q sin^2 (in units of e / (sin cos)).
    const double gq[2] = {q.t3 - q.q * s2, -q.q * s2}; // [L, R]
    const double gl[2] = {l.t3 - l.q * s2, -l.q * s2};

    // Z propagator relative to the photon one, with a fixed-width Breit-Wigner.
    const std::complex<double> chi =
        shat / std::complex<double>(shat - MZ * MZ, MZ * GAMMA_Z) / (s2 * c2);

    // Helicity amplitudes A_ij: photon exchange + Z exchange, for quark helicity i
    // and lepton helicity j. Same-helicity pairs prefer (1+cos)^2, opposite (1-cos)^2.
    double same = 0, opposite = 0;
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j) {
            const double a2 = std::norm(q.q * l.q + gq[i] * gl[j] * chi);
            (i == j ? same : opposite) += a2;
        }

    const double c = cosTheta;
    return PI * ALPHA_EM * ALPHA_EM / (8 * shat * N_COLORS) *
           (same * (1 + c) * (1 + c) + opposite * (1 - c) * (1 - c));
}

double DrellYan::generate(Random& rng, Event& ev) {
    ev.clear();
    const double s = cfg_.sqrtS * cfg_.sqrtS;

    // 1) Dimuon mass^2. Sample through a Breit-Wigner mapping so most points land
    //    near the Z peak, where the cross section is largest:
    //       M^2 = MZ^2 + MZ*GZ*tan(rho),  rho uniform.
    const double mg = MZ * GAMMA_Z;
    const double rhoMin = std::atan((cfg_.massMin * cfg_.massMin - MZ * MZ) / mg);
    const double rhoMax = std::atan((cfg_.massMax * cfg_.massMax - MZ * MZ) / mg);
    const double rho = rng.uniform(rhoMin, rhoMax);
    const double m2 = MZ * MZ + mg * std::tan(rho);
    const double jacM2 = (rhoMax - rhoMin) * ((m2 - MZ * MZ) * (m2 - MZ * MZ) + mg * mg) / mg;

    // 2) Rapidity of the dimuon system, uniform over the kinematically allowed range.
    //    tau = x1 x2 = M^2/s and y = 1/2 ln(x1/x2), so dx1 dx2 = dtau dy.
    const double tau = m2 / s;
    const double yMax = -0.5 * std::log(tau);
    const double y = rng.uniform(-yMax, yMax);
    const double x1 = std::sqrt(tau) * std::exp(y);
    const double x2 = std::sqrt(tau) * std::exp(-y);

    // 3) Decay angle of the mu- relative to the +z axis, in the dimuon rest frame.
    const double cosTheta = rng.uniform(-1, 1);
    const double phi = rng.uniform(0, 2 * PI);

    // 4) Sum over the ways to get a quark-antiquark pair: the quark can come from
    //    either proton. If it comes from proton 2 (moving along -z), the angle between
    //    quark and mu- is pi - theta.
    struct Option { int id1, id2; double w; };
    Option opts[6];
    int nOpt = 0;
    double total = 0;
    for (int q = 1; q <= 3; ++q) {
        const double wA = pdf_.f(q, x1) * pdf_.f(-q, x2) * dSigmaDCosTheta(q, m2, cosTheta);
        const double wB = pdf_.f(-q, x1) * pdf_.f(q, x2) * dSigmaDCosTheta(q, m2, -cosTheta);
        opts[nOpt++] = {q, -q, wA};
        opts[nOpt++] = {-q, q, wB};
        total += wA + wB;
    }

    // Weight = (integrand) x (volume of each sampled variable), converted to pb.
    const double weight = total * (jacM2 / s) * (2 * yMax) * 2.0 * GEV2_TO_PB;
    if (!(weight > 0)) return 0;

    // Pick which flavour configuration this event "is", proportional to its share.
    double r = rng.uniform() * total;
    const Option* chosen = &opts[nOpt - 1];
    for (int k = 0; k < nOpt; ++k) {
        if (r < opts[k].w) { chosen = &opts[k]; break; }
        r -= opts[k].w;
    }

    // 5) Build four-momenta. Muons back-to-back in the dimuon rest frame...
    const double m = std::sqrt(m2);
    const double sinTheta = std::sqrt(1 - cosTheta * cosTheta);
    const Vec4 muMinusRest(m / 2, m / 2 * sinTheta * std::cos(phi), m / 2 * sinTheta * std::sin(phi), m / 2 * cosTheta);
    const Vec4 muPlusRest(m / 2, -muMinusRest.px, -muMinusRest.py, -muMinusRest.pz);

    // ...then boosted to the lab, where the dimuon system has rapidity y and a small
    // transverse momentum from the toy kT kick (a stand-in for initial-state radiation).
    const double kx = cfg_.kTSigma > 0 ? rng.gauss(0, cfg_.kTSigma) : 0;
    const double ky = cfg_.kTSigma > 0 ? rng.gauss(0, cfg_.kTSigma) : 0;
    const double mT = std::sqrt(m2 + kx * kx + ky * ky);
    const Vec4 dimuon(mT * std::cosh(y), kx, ky, mT * std::sinh(y));
    double bx, by, bz;
    dimuon.betaVector(bx, by, bz);

    const double halfSqrtS = cfg_.sqrtS / 2;
    ev.particles.push_back({chosen->id1, -1, Vec4(x1 * halfSqrtS, 0, 0, x1 * halfSqrtS)});
    ev.particles.push_back({chosen->id2, -1, Vec4(x2 * halfSqrtS, 0, 0, -x2 * halfSqrtS)});
    ev.particles.push_back({pdg::muon, 1, muMinusRest.boosted(bx, by, bz)});
    ev.particles.push_back({-pdg::muon, 1, muPlusRest.boosted(bx, by, bz)});

    ev.weight = weight;
    ev.channel = std::abs(chosen->id1) == pdg::u ? "u ubar -> gamma*/Z -> mu mu"
               : std::abs(chosen->id1) == pdg::d ? "d dbar -> gamma*/Z -> mu mu"
                                                 : "s sbar -> gamma*/Z -> mu mu";
    ev.x1 = x1;
    ev.x2 = x2;
    ev.scale = m;
    // Angle between incoming quark and mu-, the quantity the matrix element depends on.
    ev.cosThetaStar = chosen->id1 > 0 ? cosTheta : -cosTheta;
    return weight;
}

} // namespace ppmc
