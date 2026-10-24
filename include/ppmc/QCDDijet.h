#pragma once
// QCD dijet production:  p p -> (parton parton -> parton parton) + X
//
// Every leading-order 2 -> 2 scattering of quarks and gluons is included
// (gg->gg, qg->qg, qq->qq, qqbar->gg, gg->qqbar, ...). At high energy these are by
// far the most common hard collisions at a hadron collider; each outgoing parton
// would become a jet in a real detector.

#include "ppmc/PDF.h"
#include "ppmc/Process.h"

namespace ppmc {

struct DijetConfig {
    double sqrtS = 13000; // proton-proton centre-of-mass energy (GeV)
    double pTMin = 50;    // minimum transverse momentum of the outgoing partons (GeV)
    double yMax = 4.5;    // maximum |rapidity| of each outgoing parton
    int nfFinal = 5;      // number of (massless) quark flavours that can be produced
};

class QCDDijet : public Process {
public:
    QCDDijet(const ToyPDF& pdf, DijetConfig cfg);

    std::string name() const override { return "QCD dijets  pp -> 2 partons"; }
    double generate(Random& rng, Event& event) override;

    // One-loop running strong coupling alpha_s(Q).
    static double alphaS(double q);

    // Spin- and colour-averaged |M|^2 / g_s^4 for each subprocess, as functions of the
    // Mandelstam variables s, t = (p1 - p3)^2, u = (p1 - p4)^2.
    struct MatrixElements {
        static double qqPrime(double s, double t, double u);        // q q'    -> q q'
        static double qqSame(double s, double t, double u);         // q q     -> q q
        static double qqbarToQpQpbar(double s, double t, double u); // q qbar  -> q' qbar'
        static double qqbarSame(double s, double t, double u);      // q qbar  -> q qbar
        static double qqbarToGG(double s, double t, double u);      // q qbar  -> g g
        static double ggToQQbar(double s, double t, double u);      // g g     -> q qbar
        static double qgToQG(double s, double t, double u);         // q g     -> q g
        static double ggToGG(double s, double t, double u);         // g g     -> g g
    };

private:
    const ToyPDF& pdf_;
    DijetConfig cfg_;
};

} // namespace ppmc
