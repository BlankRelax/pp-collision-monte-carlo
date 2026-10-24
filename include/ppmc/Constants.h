#pragma once
// Physical constants and particle codes. Energies in GeV.

namespace ppmc {

constexpr double PI = 3.14159265358979323846;

// Conversion from natural units (GeV^-2) to picobarns: (hbar c)^2 = 0.3894 mb GeV^2.
constexpr double GEV2_TO_PB = 0.3893794e9;

// Electroweak parameters (tree-level, PDG-like values).
constexpr double MZ = 91.1876;          // Z boson mass
constexpr double GAMMA_Z = 2.4952;      // Z boson total width
constexpr double SIN2_THETA_W = 0.2312; // weak mixing angle, sin^2(theta_W)
constexpr double ALPHA_EM = 1.0 / 128.0; // QED coupling evaluated near the Z mass

// Strong coupling at the Z mass, run to other scales at one loop.
constexpr double ALPHA_S_MZ = 0.118;

constexpr int N_COLORS = 3;

// PDG Monte Carlo particle numbering scheme. Antiparticles carry a minus sign.
namespace pdg {
constexpr int d = 1, u = 2, s = 3, c = 4, b = 5;
constexpr int electron = 11, muon = 13;
constexpr int gluon = 21;
constexpr int proton = 2212;
} // namespace pdg

} // namespace ppmc
