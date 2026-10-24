#pragma once
// Interface shared by all hard-scattering processes.

#include "ppmc/Event.h"
#include "ppmc/Random.h"

#include <string>

namespace ppmc {

class Process {
public:
    virtual ~Process() = default;
    virtual std::string name() const = 0;

    // Draws one random phase-space point, fills `event`, and returns its weight in pb.
    //
    // The weights are importance-sampling weights: the average weight over many calls
    // converges to the total cross section, and the weights falling into a histogram
    // bin (divided by the number of calls) converge to the cross section in that bin.
    // A weight of 0 means the point fell outside the physical region.
    virtual double generate(Random& rng, Event& event) = 0;
};

} // namespace ppmc
