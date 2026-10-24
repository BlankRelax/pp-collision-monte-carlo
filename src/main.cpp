// ppmc: a toy Monte Carlo event generator for proton-proton collisions.
//
// Run `ppmc --help` for usage. See README.md for the physics.

#include "ppmc/Constants.h"
#include "ppmc/DrellYan.h"
#include "ppmc/Histogram.h"
#include "ppmc/PDF.h"
#include "ppmc/QCDDijet.h"
#include "ppmc/Random.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

using namespace ppmc;

namespace {

struct Options {
    std::string process = "dy";
    long events = 200000;
    unsigned long long seed = 12345;
    double sqrtS = 13000;
    double massMin = 60, massMax = 120, kT = 2.0;
    double pTMin = 50, yMax = 4.5;
    std::string outDir = "output";
    int show = 3;
    bool checkPdf = false;
    bool quiet = false;
};

void printUsage() {
    std::cout <<
R"(Usage: ppmc [options]

Simulates proton-proton collisions with the Monte Carlo method and prints the
cross section plus histograms of the produced particles.

General options
  --process NAME   dy     Drell-Yan: q qbar -> gamma*/Z -> mu+ mu-   (default)
                   dijet  QCD 2->2 scattering of quarks and gluons
  --events N       number of collisions to simulate            (default 200000)
  --seed S         random seed; same seed -> same results       (default 12345)
  --sqrts E        proton-proton collision energy in GeV        (default 13000)
  --out DIR        directory for the CSV histogram files        (default output)
  --show N         print the first N events in full             (default 3)
  --quiet          do not draw the text histograms
  --check-pdf      print the parton-distribution sum rules and exit
  --help           show this message

Drell-Yan options
  --mmin M         lowest dimuon mass generated, GeV             (default 60)
  --mmax M         highest dimuon mass generated, GeV            (default 120)
  --kt K           width of toy transverse kick to the Z, GeV    (default 2)

Dijet options
  --ptmin P        minimum parton transverse momentum, GeV       (default 50)
  --ymax Y         maximum |rapidity| of each parton             (default 4.5)

Examples
  ppmc                                  # 200k Drell-Yan events at 13 TeV
  ppmc --process dijet --ptmin 100      # jets above 100 GeV
  ppmc --sqrts 7000 --events 1000000    # LHC Run 1 energy, more statistics
)";
}

bool parseArgs(int argc, char** argv, Options& o) {
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto next = [&]() -> std::string {
            if (i + 1 >= argc) { std::cerr << "missing value for " << a << "\n"; std::exit(1); }
            return argv[++i];
        };
        if (a == "--help" || a == "-h") { printUsage(); std::exit(0); }
        else if (a == "--process") o.process = next();
        else if (a == "--events") o.events = std::stol(next());
        else if (a == "--seed") o.seed = std::stoull(next());
        else if (a == "--sqrts") o.sqrtS = std::stod(next());
        else if (a == "--mmin") o.massMin = std::stod(next());
        else if (a == "--mmax") o.massMax = std::stod(next());
        else if (a == "--kt") o.kT = std::stod(next());
        else if (a == "--ptmin") o.pTMin = std::stod(next());
        else if (a == "--ymax") o.yMax = std::stod(next());
        else if (a == "--out") o.outDir = next();
        else if (a == "--show") o.show = std::stoi(next());
        else if (a == "--quiet") o.quiet = true;
        else if (a == "--check-pdf") o.checkPdf = true;
        else { std::cerr << "unknown option " << a << " (try --help)\n"; return false; }
    }
    if (o.process != "dy" && o.process != "dijet") {
        std::cerr << "unknown process '" << o.process << "' (use dy or dijet)\n";
        return false;
    }
    if (o.events <= 0) { std::cerr << "--events must be positive\n"; return false; }
    if (o.massMin <= 0 || o.massMax <= o.massMin || o.massMax >= o.sqrtS) {
        std::cerr << "need 0 < --mmin < --mmax < --sqrts\n";
        return false;
    }
    if (o.pTMin <= 0 || 2 * o.pTMin >= o.sqrtS) { std::cerr << "need 0 < 2*--ptmin < --sqrts\n"; return false; }
    return true;
}

std::string particleName(int id) {
    static const std::map<int, std::string> names = {
        {1, "d"}, {2, "u"}, {3, "s"}, {4, "c"}, {5, "b"}, {-1, "dbar"}, {-2, "ubar"}, {-3, "sbar"},
        {-4, "cbar"}, {-5, "bbar"}, {21, "g"}, {13, "mu-"}, {-13, "mu+"}};
    const auto it = names.find(id);
    return it != names.end() ? it->second : std::to_string(id);
}

void printEvent(const Event& ev, long index) {
    std::cout << "\n  Event " << index << ":  " << ev.channel << "\n";
    std::cout << "    weight = " << std::scientific << std::setprecision(3) << ev.weight << " pb,  x1 = "
              << std::fixed << std::setprecision(5) << ev.x1 << ",  x2 = " << ev.x2
              << ",  Q = " << std::setprecision(2) << ev.scale << " GeV\n";
    std::cout << "    " << std::left << std::setw(6) << "name" << std::setw(9) << "status" << std::right
              << std::setw(11) << "E" << std::setw(10) << "px" << std::setw(10) << "py" << std::setw(11) << "pz"
              << std::setw(9) << "pT" << std::setw(8) << "y" << "\n";
    for (const auto& p : ev.particles) {
        std::cout << "    " << std::left << std::setw(6) << particleName(p.id) << std::setw(9)
                  << (p.status < 0 ? "in" : "out") << std::right << std::fixed << std::setprecision(2)
                  << std::setw(11) << p.p.e << std::setw(10) << p.p.px << std::setw(10) << p.p.py
                  << std::setw(11) << p.p.pz << std::setw(9) << p.p.pT();
        if (p.status > 0) std::cout << std::setw(8) << p.p.rapidity();
        std::cout << "\n";
    }
    const Vec4 diff = ev.sumMomentum(-1) - ev.sumMomentum(1);
    std::cout << "    momentum balance (in - out): " << std::setprecision(3) << diff << "\n";
    if (std::abs(diff.px) > 1e-6 || std::abs(diff.py) > 1e-6)
        std::cout << "    (non-zero because of the toy kT kick; the rest of the proton, not simulated, recoils)\n";
}

// Running mean and standard error of the event weights = the cross section estimate.
struct WeightStats {
    long n = 0, nonZero = 0;
    double sumW = 0, sumW2 = 0, maxW = 0;
    void add(double w) {
        ++n;
        if (w > 0) ++nonZero;
        sumW += w;
        sumW2 += w * w;
        if (w > maxW) maxW = w;
    }
    double mean() const { return sumW / n; }
    double error() const {
        const double m = mean();
        return std::sqrt(std::max(0.0, sumW2 / n - m * m) / n);
    }
};

// A named histogram plus the quantity it shows, so analyses stay declarative.
struct Plot {
    std::string file, yLabel;
    Histogram hist;
};

std::string formatXsec(double v, double e) {
    std::ostringstream ss;
    if (v > 1e6) ss << std::setprecision(4) << v / 1e6 << " +- " << std::setprecision(2) << e / 1e6 << " ub (microbarn)";
    else if (v > 1e3) ss << std::setprecision(4) << v / 1e3 << " +- " << std::setprecision(2) << e / 1e3 << " nb";
    else ss << std::setprecision(4) << v << " +- " << std::setprecision(2) << e << " pb";
    return ss.str();
}

} // namespace

int main(int argc, char** argv) {
    Options opt;
    if (!parseArgs(argc, argv, opt)) return 1;

    ToyPDF pdf;
    if (opt.checkPdf) {
        pdf.printSumRules(std::cout);
        return 0;
    }

    Random rng(opt.seed);
    std::unique_ptr<Process> process;
    std::vector<Plot> plots;
    std::function<void(const Event&)> analyse;

    // Extra per-process tallies, all in pb (weights summed, divided by N at the end).
    double fiducial = 0, forward = 0, backward = 0;
    std::map<std::string, double> channelXsec;

    if (opt.process == "dy") {
        process = std::make_unique<DrellYan>(pdf, DrellYanConfig{opt.sqrtS, opt.massMin, opt.massMax, opt.kT});
        plots = {
            {"dy_mass.csv", "dsigma/dM [pb/GeV]", Histogram("Dimuon invariant mass", "M(mu+ mu-) [GeV]", 30, opt.massMin, opt.massMax)},
            {"dy_rapidity.csv", "dsigma/dy [pb]", Histogram("Dimuon rapidity", "y(mu+ mu-)", 24, -6, 6)},
            {"dy_ptZ.csv", "dsigma/dpT [pb/GeV]", Histogram("Dimuon transverse momentum", "pT(mu+ mu-) [GeV]", 20, 0, 10)},
            {"dy_muon_pt.csv", "dsigma/dpT [pb/GeV]", Histogram("Muon transverse momentum (both muons)", "pT(mu) [GeV]", 30, 0, 75)},
            {"dy_muon_eta.csv", "dsigma/deta [pb]", Histogram("Muon pseudorapidity (both muons)", "eta(mu)", 32, -8, 8)},
            {"dy_costheta.csv", "dsigma/dcos [pb]", Histogram("Decay angle: quark vs mu- in dimuon rest frame", "cos(theta*)", 20, -1, 1)},
        };
        analyse = [&](const Event& ev) {
            const Vec4& muM = ev.particles[2].p;
            const Vec4& muP = ev.particles[3].p;
            const Vec4 z = muM + muP;
            const double w = ev.weight;
            plots[0].hist.fill(z.m(), w);
            plots[1].hist.fill(z.rapidity(), w);
            plots[2].hist.fill(z.pT(), w);
            plots[3].hist.fill(muM.pT(), w);
            plots[3].hist.fill(muP.pT(), w);
            plots[4].hist.fill(muM.eta(), w);
            plots[4].hist.fill(muP.eta(), w);
            plots[5].hist.fill(ev.cosThetaStar, w);
            (ev.cosThetaStar > 0 ? forward : backward) += w;
            // A typical LHC detector only sees muons that are central and not too soft.
            const auto seen = [](const Vec4& p) { return p.pT() > 25 && std::abs(p.eta()) < 2.5; };
            if (seen(muM) && seen(muP)) fiducial += w;
        };
    } else {
        process = std::make_unique<QCDDijet>(pdf, DijetConfig{opt.sqrtS, opt.pTMin, opt.yMax, 5});
        const double p = opt.pTMin;
        plots = {
            {"dijet_pt.csv", "dsigma/dpT [pb/GeV]", Histogram("Parton transverse momentum", "pT [GeV]", 35, p, 8 * p)},
            {"dijet_mass.csv", "dsigma/dM [pb/GeV]", Histogram("Dijet invariant mass", "M(jj) [GeV]", 40, 0, 20 * p)},
            {"dijet_rapidity.csv", "dsigma/dy [pb]", Histogram("Parton rapidity (both partons)", "y", 36, -opt.yMax, opt.yMax)},
            {"dijet_chi.csv", "dsigma/dchi [pb]", Histogram("Dijet angular variable chi = exp|y1 - y2|", "chi", 29, 1, 30)},
        };
        analyse = [&](const Event& ev) {
            const Vec4& j1 = ev.particles[2].p;
            const Vec4& j2 = ev.particles[3].p;
            const double w = ev.weight;
            plots[0].hist.fill(j1.pT(), w);
            plots[1].hist.fill((j1 + j2).m(), w);
            plots[2].hist.fill(j1.rapidity(), w);
            plots[2].hist.fill(j2.rapidity(), w);
            plots[3].hist.fill(std::exp(std::abs(j1.rapidity() - j2.rapidity())), w);
            channelXsec[ev.channel] += w;
        };
    }

    std::cout << "==============================================================\n";
    std::cout << " ppmc - toy Monte Carlo for proton-proton collisions\n";
    std::cout << "==============================================================\n";
    std::cout << " Process : " << process->name() << "\n";
    std::cout << " sqrt(s) : " << opt.sqrtS << " GeV\n";
    std::cout << " Events  : " << opt.events << "   (seed " << opt.seed << ")\n";
    if (opt.process == "dy")
        std::cout << " Cuts    : " << opt.massMin << " < M(mu mu) < " << opt.massMax << " GeV\n";
    else
        std::cout << " Cuts    : pT > " << opt.pTMin << " GeV, |y| < " << opt.yMax << "\n";

    // ----------------------------- the event loop -----------------------------
    WeightStats stats;
    Event ev;
    int shown = 0;
    const auto t0 = std::chrono::steady_clock::now();
    for (long i = 0; i < opt.events; ++i) {
        const double w = process->generate(rng, ev);
        stats.add(w);
        if (w <= 0) continue;
        analyse(ev);
        if (shown < opt.show) {
            if (shown == 0) std::cout << "\n--- First few events (lab frame, GeV) ---";
            printEvent(ev, i);
            ++shown;
        }
    }
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    // Weights summed per bin / N = cross section per bin.
    const double invN = 1.0 / static_cast<double>(stats.n);
    for (auto& pl : plots) pl.hist.scale(invN);

    // --------------------------------- results ---------------------------------
    const double sigma = stats.mean(), err = stats.error();
    std::cout << "\n--- Results ---\n";
    std::cout << std::defaultfloat << std::setprecision(4);
    std::cout << " Total cross section     : " << formatXsec(sigma, err) << "\n";
    if (opt.process == "dy") {
        const double fid = fiducial * invN;
        std::cout << " Seen by a detector      : " << formatXsec(fid, err * fid / sigma)
                  << "   (both muons pT > 25 GeV, |eta| < 2.5)\n";
        std::cout << " Detector acceptance     : " << std::setprecision(3) << 100 * fid / sigma << " %\n";
        std::cout << " Forward-backward asym.  : A_FB = " << std::setprecision(4)
                  << (forward - backward) / (forward + backward) << "\n";
        std::cout << " Expected events in 140 fb^-1 (LHC Run 2): " << std::setprecision(3)
                  << std::scientific << fid * 140e3 << std::defaultfloat << "\n";
    } else {
        std::cout << " Share of the cross section by subprocess:\n";
        std::vector<std::pair<double, std::string>> sorted;
        for (const auto& [name, w] : channelXsec) sorted.push_back({w * invN, name});
        std::sort(sorted.rbegin(), sorted.rend());
        for (const auto& [x, name] : sorted)
            std::cout << "     " << std::left << std::setw(20) << name << std::right << std::fixed
                      << std::setprecision(1) << std::setw(6) << 100 * x / sigma << " %\n";
        std::cout << std::defaultfloat;
    }
    std::cout << std::setprecision(4);
    std::cout << " Points inside phase space: " << 100.0 * stats.nonZero / stats.n << " %\n";
    std::cout << " Weight efficiency <w>/max: " << sigma / stats.maxW
              << "   (fraction kept if you unweighted the events)\n";
    std::cout << " Run time                 : " << std::setprecision(3) << seconds << " s  ("
              << std::setprecision(3) << stats.n / seconds / 1e6 << " M events/s)\n";

    if (!opt.quiet)
        for (const auto& pl : plots) pl.hist.printASCII(std::cout, pl.yLabel);

    std::filesystem::create_directories(opt.outDir);
    for (const auto& pl : plots) pl.hist.writeCSV((std::filesystem::path(opt.outDir) / pl.file).string(), pl.yLabel);
    std::cout << "\nHistograms written to " << opt.outDir << "/ as CSV (plot them with scripts/plot.py).\n";
    return 0;
}
