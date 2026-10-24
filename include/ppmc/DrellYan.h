#pragma once
// Drell-Yan production:  p p -> (q qbar -> gamma*/Z -> mu- mu+) + X
//
// A quark from one proton annihilates with an antiquark from the other into a
// virtual photon or a Z boson, which decays to a muon pair. The full gamma*/Z
// interference is included, so the forward-backward asymmetry comes out right.

#include "ppmc/PDF.h"
#include "ppmc/Process.h"

namespace ppmc {

struct DrellYanConfig {
    double sqrtS = 13000;  // proton-proton centre-of-mass energy (GeV)
    double massMin = 60;   // generated dimuon mass window (GeV)
    double massMax = 120;
    double kTSigma = 2.0;  // width of the toy Gaussian transverse kick given to the muon pair (GeV)
};

class DrellYan : public Process {
public:
    DrellYan(const ToyPDF& pdf, DrellYanConfig cfg);

    std::string name() const override { return "Drell-Yan  pp -> gamma*/Z -> mu+ mu-"; }
    double generate(Random& rng, Event& event) override;

    // Partonic differential cross section d(sigma)/d(cos theta) for q qbar -> mu- mu+
    // in GeV^-2, colour-averaged. `quarkId` is the (positive) quark flavour, `shat`
    // the partonic energy squared, and theta the angle between the incoming quark and
    // the outgoing mu-.
    static double dSigmaDCosTheta(int quarkId, double shat, double cosTheta);

private:
    const ToyPDF& pdf_;
    DrellYanConfig cfg_;
};

} // namespace ppmc
