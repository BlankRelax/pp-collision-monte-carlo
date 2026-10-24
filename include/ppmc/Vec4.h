#pragma once
// Minimal Lorentz four-vector in natural units (GeV), metric (+,-,-,-).

#include <cmath>
#include <ostream>

namespace ppmc {

struct Vec4 {
    double e = 0, px = 0, py = 0, pz = 0;

    Vec4() = default;
    Vec4(double e_, double px_, double py_, double pz_) : e(e_), px(px_), py(py_), pz(pz_) {}

    Vec4 operator+(const Vec4& o) const { return {e + o.e, px + o.px, py + o.py, pz + o.pz}; }
    Vec4 operator-(const Vec4& o) const { return {e - o.e, px - o.px, py - o.py, pz - o.pz}; }
    Vec4 operator*(double s) const { return {e * s, px * s, py * s, pz * s}; }
    Vec4& operator+=(const Vec4& o) { e += o.e; px += o.px; py += o.py; pz += o.pz; return *this; }

    // Invariant mass squared, E^2 - |p|^2.
    double m2() const { return e * e - px * px - py * py - pz * pz; }
    double m() const { const double v = m2(); return v > 0 ? std::sqrt(v) : 0.0; }
    double p() const { return std::sqrt(px * px + py * py + pz * pz); }
    // Momentum transverse to the beam (z) axis.
    double pT() const { return std::hypot(px, py); }
    double phi() const { return std::atan2(py, px); }
    // Rapidity y = 1/2 ln((E+pz)/(E-pz)); additive under boosts along z.
    double rapidity() const { return 0.5 * std::log((e + pz) / (e - pz)); }
    // Pseudorapidity eta = -ln tan(theta/2); equals y for massless particles.
    double eta() const { return std::asinh(pz / pT()); }

    // Velocity of a frame in which this four-vector is at rest.
    void betaVector(double& bx, double& by, double& bz) const { bx = px / e; by = py / e; bz = pz / e; }

    // Returns this vector as seen by an observer for whom the original frame moves
    // with velocity (bx, by, bz).
    Vec4 boosted(double bx, double by, double bz) const {
        const double b2 = bx * bx + by * by + bz * bz;
        if (b2 <= 0) return *this;
        const double gamma = 1.0 / std::sqrt(1.0 - b2);
        const double bp = bx * px + by * py + bz * pz;
        const double g2 = (gamma - 1.0) / b2;
        return {gamma * (e + bp),
                px + g2 * bp * bx + gamma * bx * e,
                py + g2 * bp * by + gamma * by * e,
                pz + g2 * bp * bz + gamma * bz * e};
    }
};

inline double dot(const Vec4& a, const Vec4& b) {
    return a.e * b.e - a.px * b.px - a.py * b.py - a.pz * b.pz;
}

inline std::ostream& operator<<(std::ostream& os, const Vec4& v) {
    return os << "(" << v.e << ", " << v.px << ", " << v.py << ", " << v.pz << ")";
}

} // namespace ppmc
