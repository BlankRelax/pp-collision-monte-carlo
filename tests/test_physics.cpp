// Physics self-tests. Each check compares the code against something we know
// independently: conservation laws, sum rules, or a textbook formula.

#include "ppmc/Constants.h"
#include "ppmc/DrellYan.h"
#include "ppmc/PDF.h"
#include "ppmc/QCDDijet.h"
#include "ppmc/Random.h"
#include "ppmc/Vec4.h"

#include <cmath>
#include <iostream>
#include <string>

using namespace ppmc;

static int failures = 0;

static void check(bool ok, const std::string& what) {
    std::cout << (ok ? "  [pass] " : "  [FAIL] ") << what << "\n";
    if (!ok) ++failures;
}

static bool close(double a, double b, double relTol) { return std::abs(a - b) <= relTol * std::max(std::abs(a), std::abs(b)); }

static bool close(const Vec4& a, const Vec4& b, double absTol) {
    return std::abs(a.e - b.e) < absTol && std::abs(a.px - b.px) < absTol && std::abs(a.py - b.py) < absTol &&
           std::abs(a.pz - b.pz) < absTol;
}

static void testLorentz() {
    std::cout << "Lorentz kinematics\n";
    const Vec4 p(50, 10, -20, 30);
    const Vec4 q = p.boosted(0.3, -0.4, 0.5);
    check(close(p.m2(), q.m2(), 1e-12), "a boost leaves the invariant mass unchanged");
    check(close(q.boosted(-0.3, 0.4, -0.5), p, 1e-9), "boosting by -beta undoes a boost by beta");
    // Pure z boosts compose by adding rapidities.
    const Vec4 r = p.boosted(0, 0, std::tanh(0.7));
    check(close(r.rapidity(), p.rapidity() + 0.7, 1e-12), "a z boost shifts rapidity by a constant");
    const Vec4 back = r.boosted(0, 0, -std::tanh(0.7));
    check(close(back, p, 1e-9), "boosting forward then back returns the original vector");
}

static void testPDF() {
    std::cout << "Parton distribution sum rules\n";
    ToyPDF pdf;
    const double uVal = pdf.integrate(pdg::u, false) - pdf.integrate(-pdg::u, false);
    const double dVal = pdf.integrate(pdg::d, false) - pdf.integrate(-pdg::d, false);
    double momentum = 0;
    for (int id : {1, 2, 3, -1, -2, -3, 21}) momentum += pdf.integrate(id, true);
    check(close(uVal, 2.0, 1e-3), "proton contains 2 net up quarks");
    check(close(dVal, 1.0, 1e-3), "proton contains 1 net down quark");
    check(close(momentum, 1.0, 1e-3), "partons carry 100% of the proton momentum");
    check(pdf.f(pdg::c, 0.1) == 0 && pdf.f(pdg::gluon, 1.0) == 0, "no charm, and nothing at x = 1");
}

static void testDrellYanPeak() {
    std::cout << "Drell-Yan partonic cross section\n";
    // At sqrt(shat) = MZ the Z dominates and the cross section must equal the textbook
    // Breit-Wigner peak value  sigma = 12 pi Gamma(Z->qq, one colour) Gamma(Z->mumu) / (Nc MZ^2 GZ^2).
    const double s2 = SIN2_THETA_W, c2 = 1 - s2;
    auto width = [&](double q, double t3) {
        const double gl = t3 - q * s2, gr = -q * s2;
        return ALPHA_EM * MZ * (gl * gl + gr * gr) / (6 * s2 * c2);
    };
    const double gU = width(2.0 / 3.0, 0.5), gMu = width(-1, -0.5);
    const double expected = 12 * PI * gU * gMu / (N_COLORS * MZ * MZ * GAMMA_Z * GAMMA_Z);
    // Integrate d(sigma)/d(cos theta) over cos theta with Simpson's rule.
    const int n = 200;
    double sigma = 0;
    for (int k = 0; k <= n; ++k) {
        const double c = -1 + 2.0 * k / n;
        sigma += (k == 0 || k == n ? 1 : (k % 2 ? 4 : 2)) * DrellYan::dSigmaDCosTheta(pdg::u, MZ * MZ, c);
    }
    sigma *= (2.0 / n) / 3;
    check(close(sigma, expected, 0.01), "u ubar -> Z -> mu mu peak matches Breit-Wigner formula (1%)");
    // Far below the Z, the photon dominates: sigma -> 4 pi alpha^2 Q^2 / (3 Nc s).
    const double sh = 4.0; // 2 GeV
    double low = 0;
    for (int k = 0; k <= n; ++k) {
        const double c = -1 + 2.0 * k / n;
        low += (k == 0 || k == n ? 1 : (k % 2 ? 4 : 2)) * DrellYan::dSigmaDCosTheta(pdg::u, sh, c);
    }
    low *= (2.0 / n) / 3;
    const double qed = 4 * PI * ALPHA_EM * ALPHA_EM * (4.0 / 9.0) / (3 * N_COLORS * sh);
    check(close(low, qed, 0.01), "far below the Z it reduces to pure QED (1%)");
    // Above the Z, gamma/Z interference makes the mu- prefer the quark direction.
    const double sh2 = 120.0 * 120.0;
    check(DrellYan::dSigmaDCosTheta(pdg::u, sh2, 0.8) > DrellYan::dSigmaDCosTheta(pdg::u, sh2, -0.8),
          "positive forward-backward asymmetry above the Z peak");
}

static void testEvents() {
    std::cout << "Generated events\n";
    ToyPDF pdf;
    Random rng(7);

    DrellYan dy(pdf, DrellYanConfig{13000, 60, 120, 0.0});
    DijetConfig jc;
    QCDDijet jets(pdf, jc);

    bool dyOk = true, jetOk = true, massOk = true;
    double sumDY = 0;
    const int n = 20000;
    for (int i = 0; i < n; ++i) {
        Event ev;
        const double w = dy.generate(rng, ev);
        sumDY += w;
        if (w > 0) {
            dyOk &= close(ev.sumMomentum(-1), ev.sumMomentum(1), 1e-6 * 13000);
            const double m = (ev.particles[2].p + ev.particles[3].p).m();
            massOk &= m > 60 - 1e-6 && m < 120 + 1e-6;
        }
        Event ej;
        if (jets.generate(rng, ej) > 0) jetOk &= close(ej.sumMomentum(-1), ej.sumMomentum(1), 1e-6 * 13000);
    }
    check(dyOk, "Drell-Yan: four-momentum conserved in every event");
    check(massOk, "Drell-Yan: dimuon mass inside the requested window");
    check(jetOk, "Dijets: four-momentum conserved in every event");
    const double sigmaNb = sumDY / n / 1e3;
    std::cout << "    (toy sigma(pp -> Z -> mu mu, 60-120 GeV, 13 TeV) = " << sigmaNb << " nb; LHC measures ~1.9 nb)\n";
    check(sigmaNb > 0.5 && sigmaNb < 5, "Drell-Yan cross section has the right order of magnitude");
}

static void testQCD() {
    std::cout << "QCD ingredients\n";
    check(close(QCDDijet::alphaS(MZ), ALPHA_S_MZ, 1e-12), "alpha_s(MZ) = 0.118");
    check(QCDDijet::alphaS(1000) < QCDDijet::alphaS(100), "alpha_s decreases at high energy (asymptotic freedom)");
    // gg -> gg must be symmetric under swapping t and u (identical gluons).
    const double s = 1000, t = -300, u = -700;
    using ME = QCDDijet::MatrixElements;
    check(close(ME::ggToGG(s, t, u), ME::ggToGG(s, u, t), 1e-12), "gg->gg symmetric under t <-> u");
    check(close(ME::qqSame(s, t, u), ME::qqSame(s, u, t), 1e-12), "qq->qq symmetric under t <-> u");
}

int main() {
    testLorentz();
    testPDF();
    testDrellYanPeak();
    testEvents();
    testQCD();
    std::cout << (failures ? "\nSOME TESTS FAILED\n" : "\nAll tests passed.\n");
    return failures ? 1 : 0;
}
