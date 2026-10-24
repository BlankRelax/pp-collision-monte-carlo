#pragma once
// The event record: everything we know about one simulated collision.

#include "ppmc/Vec4.h"

#include <string>
#include <vector>

namespace ppmc {

struct Particle {
    int id = 0;     // PDG code (e.g. 2 = u quark, -13 = mu+, 21 = gluon)
    int status = 0; // -1 = incoming parton, +1 = final-state particle
    Vec4 p;         // four-momentum in the lab frame (GeV)
};

struct Event {
    double weight = 0;   // contribution to the cross section, in pb (see README)
    std::string channel; // which partonic subprocess produced it, e.g. "qqbar->gamma*/Z->mumu"
    double x1 = 0, x2 = 0; // momentum fractions of the two colliding partons
    double scale = 0;      // hard scale Q of the collision (GeV)
    double cosThetaStar = 0; // decay angle in the partonic rest frame (Drell-Yan only)
    std::vector<Particle> particles;

    void clear() { *this = Event{}; }

    Vec4 sumMomentum(int status) const {
        Vec4 total;
        for (const auto& prt : particles)
            if (prt.status == status) total += prt.p;
        return total;
    }
};

} // namespace ppmc
