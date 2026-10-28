// Full-evolution checks of CoulombAfterburner::Apply on simple systems with known solutions.
// Build and run with checks/coulomb/run.sh. Spectators are switched off through the configuration
// (CoulombUseSpectators = 0), so that the test events contain only the particles defined here.
//
// Units as in THERMINATOR: inside the afterburner lengths and times are in GeV^-1, momenta in GeV.
// Here positions and times are given in fm (fm/c) and converted with kHbarC.
//
//   A. free motion                          F. decay of a static source during Rutherford scattering
//   B. relativistic Rutherford scattering   G. decay of a moving source passing a probe
//   C. two charges released from rest       H. Platonic solids of charges released from rest
//   D. the pair of C boosted to beta = 0.9  I. a hexagon of charges, at rest and boosted to beta = 0.9
//   E. acceleration (radiation) term        J. random cluster of 20 charges vs instantaneous Coulomb N-body
//   K. settings of the afterburner from the configuration (events.ini, events_template)
//   L. a particle with NaN momentum does not affect the others
//   M. close encounters (impact parameters down to 0.01 fm) with sub-steps, fixed and moving source
//
// Besides the pass/fail checks, trajectories of the code and of the exact solutions are written as TGraphs
// to a ROOT file (first argument, default coulomb_checks.root), one directory per check. The second argument is
// the top directory of the repository, with events.ini and events_template for K.
//   <label>_code, <label>_exact, <label>_diff (= code - exact) as functions of t (fm/c),
//   and a few other graphs described where they are made.

// The radiation check (E) evaluates the force for a prescribed source trajectory, which needs the internal state
// of CoulombAfterburner: everything its header includes is included first, then the header itself with private
// members made public.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <utility>
#include <vector>

#include "AbstractAfterburner.h"
#include "AbstractEventSaver.h"
#include "Configurator.h"
#include "CoulombUnits.h"
#include "Event.h"
#include "Particle.h"
#include "ParticleDB.h"
#include "ParticleType.h"
#include "Parser.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TGlauberMC.h"
#include "TGraph.h"
#include "THGlobal.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TProfile.h"
#include "TRandom3.h"
#include "TVector3.h"

#define private public
#define protected public
#include "CoulombAfterburner.h"
#undef private
#undef protected

TString sModelINI;   // globals expected by the library code
TString sTimeStamp;

static const double kAlphaHbarC = 1.439964e-3; // GeV fm, alpha * hbar * c (CODATA)

static int sFailures = 0;

// Configuration of the afterburner in all runs: no spectators; with sub-steps of close pairs in M
static Configurator sChecksConfig;
static Configurator sSubstepConfig;
static Configurator *sRunConfig = &sChecksConfig;

static void Check(const char *aName, double aGot, double aExpected, double aTolerance, bool aRelative = true)
{
    double tDiff = std::fabs(aGot - aExpected);
    if (aRelative) tDiff /= std::fabs(aExpected);
    bool tOk = tDiff <= aTolerance;
    printf("%s  %-62s got %.8g, expected %.8g, %s diff %.2e (tolerance %.0e)\n",
           tOk ? "PASS" : "FAIL", aName, aGot, aExpected, aRelative ? "rel." : "abs.", tDiff, aTolerance);
    if (!tOk) ++sFailures;
}

static void Info(const char *aText) { printf("INFO  %s\n", aText); }

// ---------------------------------------------------------------------------------------------
// Graphs, written at the end to the output file in one directory per check

static std::vector<std::pair<std::string, TGraph *> > sGraphs;

static void AddGraph(const char *aDir, const std::string &aName, const std::string &aTitle,
                     const std::vector<double> &aX, const std::vector<double> &aY)
{
    TGraph *tGraph = new TGraph(aX.size(), aX.data(), aY.data());
    tGraph->SetName(aName.c_str());
    tGraph->SetTitle(aTitle.c_str());
    sGraphs.push_back(std::make_pair(std::string(aDir), tGraph));
}

// Adds <label>_code, <label>_exact and <label>_diff (code - exact) as functions of t
static void AddCodeExact(const char *aDir, const std::string &aLabel, const std::string &aDescription, const std::string &aYTitle,
                         const std::vector<double> &aT, const std::vector<double> &aCode, const std::vector<double> &aExact)
{
    std::vector<double> tDiff(aCode.size());
    for (size_t i = 0; i < aCode.size(); ++i) tDiff[i] = aCode[i] - aExact[i];
    AddGraph(aDir, aLabel + "_code", aDescription + " (code);t (fm/c);" + aYTitle, aT, aCode);
    AddGraph(aDir, aLabel + "_exact", aDescription + " (exact);t (fm/c);" + aYTitle, aT, aExact);
    AddGraph(aDir, aLabel + "_diff", aDescription + " (code - exact);t (fm/c);#Delta" + aYTitle, aT, tDiff);
}

static void WriteGraphs(const char *aFileName)
{
    TFile tFile(aFileName, "RECREATE");
    for (auto &tEntry : sGraphs) {
        TDirectory *tDir = tFile.GetDirectory(tEntry.first.c_str());
        if (!tDir) tDir = tFile.mkdir(tEntry.first.c_str());
        tDir->cd();
        tEntry.second->Write();
    }
    tFile.Close();
    for (auto &tEntry : sGraphs) delete tEntry.second;
    sGraphs.clear();
}

static double MaxAbs(const std::vector<double> &aV)
{
    double tMax = 0;
    for (double x : aV) tMax = std::max(tMax, std::fabs(x));
    return tMax;
}

// ---------------------------------------------------------------------------------------------
// Running the afterburner

struct PartInit {
    ParticleType *type;
    TVector3 x;           // fm
    TVector3 p;           // GeV
    double t = 0;         // creation time, fm/c
    int father = -1;      // index of the father in the list, -1 for primary particles
    bool decayed = false; // the particle decays (its daughters follow in the list)
};

struct RunResult {
    std::vector<TVector3> finalMomenta;                    // GeV
    std::vector<std::vector<TLorentzVector> > trajectory;  // (x, y, z, t) in fm, fm/c, at the time steps
    std::vector<std::vector<TLorentzVector> > fourMomenta; // (px, py, pz, E) in GeV, at the time steps
    bool restored;                                         // positions and times restored after Apply()
};

// Runs CoulombAfterburner::Apply on one event made of the given particles. The time steps start at the earliest
// creation time; daughters are created like in ParticleDecayer, with SetParticlePX(..., father).
static RunResult RunCoulomb(const char *aDir, int aSteps, double aStepFm, const std::vector<PartInit> &aInit, ParticleDB *aDB)
{
    mkdir(aDir, 0755);
    if (chdir(aDir) != 0) { perror(aDir); exit(2); }
    RunResult tResult;
    {
        CoulombAfterburner tCoulomb(aSteps, aStepFm, 0., nullptr, aDB);
        tCoulomb.ReadParameters(sRunConfig);
        tCoulomb.SetKeepTrajectories(true);
        tCoulomb.Init();
        Event tEvent; // resets the eid counter
        std::list<Particle> *tList = tEvent.GetParticleList();
        for (const PartInit &tInit : aInit) {
            Particle tParticle(tInit.type);
            double tE = std::sqrt(tParticle.mass * tParticle.mass + tInit.p.Mag2());
            Particle *tFather = tInit.father >= 0 ? &*std::next(tList->begin(), tInit.father) : nullptr;
            tParticle.SetParticlePX(tE, tInit.p.X(), tInit.p.Y(), tInit.p.Z(),
                                    tInit.t / kHbarC, tInit.x.X() / kHbarC, tInit.x.Y() / kHbarC, tInit.x.Z() / kHbarC, tFather);
            if (tInit.decayed) tParticle.SetDecayed();
            tList->push_back(tParticle);
        }
        tCoulomb.Apply(&tEvent);
        tResult.restored = true;
        size_t i = 0;
        for (Particle &tParticle : *tList) {
            tResult.finalMomenta.push_back(tParticle.VecP());
            TVector3 tX0 = (1 / kHbarC) * aInit[i].x;
            if (std::fabs(tParticle.T() - aInit[i].t / kHbarC) > 1e-6 * (1 + aInit[i].t / kHbarC)
                || (tParticle.VecX() - tX0).Mag() > 1e-6 * (1 + tX0.Mag())) tResult.restored = false;
            std::vector<TLorentzVector> tTrajectory;
            for (const TLorentzVector &tX : tCoulomb.GetTrajectory(tParticle.eid)) tTrajectory.push_back(kHbarC * tX);
            tResult.trajectory.push_back(tTrajectory);
            tResult.fourMomenta.push_back(tCoulomb.GetFourMomenta(tParticle.eid));
            ++i;
        }
        tCoulomb.Finalize();
    }
    if (chdir("..") != 0) { perror(".."); exit(2); }
    return tResult;
}

static ParticleType MakeType(const char *aName, int aPDG, double aMass, int aCharge)
{
    ParticleType tType;
    tType.SetName(aName);
    tType.SetPDGCode(aPDG);
    tType.SetMass(aMass);
    tType.SetBarionN(0);
    tType.SetStrangeN(0);
    tType.SetI3(aCharge); // GetCharge() = I3 + (B + S) / 2
    return tType;
}

// ---------------------------------------------------------------------------------------------
// Reference solutions

// Relativistic motion of a charge in the field of a fixed point charge at the origin, by RK4.
// aKappa = k q Q in GeV fm (> 0 repulsive); x in fm, p in GeV, t in fm/c.
struct State { TVector3 x, p; };
static State Derivative(const State &aS, double aMass, double aKappa)
{
    double tE = std::sqrt(aMass * aMass + aS.p.Mag2());
    double tR = aS.x.Mag();
    return State{(1 / tE) * aS.p, (aKappa / (tR * tR * tR)) * aS.x};
}
static State RK4Step(const State &aS, double aH, double aMass, double aKappa)
{
    auto add = [](const State &a, const State &b, double f) { return State{a.x + f * b.x, a.p + f * b.p}; };
    State k1 = Derivative(aS, aMass, aKappa);
    State k2 = Derivative(add(aS, k1, aH / 2), aMass, aKappa);
    State k3 = Derivative(add(aS, k2, aH / 2), aMass, aKappa);
    State k4 = Derivative(add(aS, k3, aH), aMass, aKappa);
    return State{aS.x + (aH / 6) * (k1.x + 2 * k2.x + 2 * k3.x + k4.x),
                 aS.p + (aH / 6) * (k1.p + 2 * k2.p + 2 * k3.p + k4.p)};
}

// Exact relativistic Coulomb orbit (Darwin 1913): u = 1/r satisfies u'' + g^2 u = -kappa E / L^2, g^2 = 1 - kappa^2/L^2,
// so u = A cos(g phi) - C, with the periapsis at phi = 0. Constructed from a start at distance r0 with momentum p0
// and impact parameter b.
struct ExactOrbit {
    double mMass, mE, mL, mG, mC, mA;
    ExactOrbit(double aMass, double aKappa, double aR0, double aP0, double aB) : mMass(aMass)
    {
        mE = std::sqrt(aMass * aMass + aP0 * aP0) + aKappa / aR0; // total energy
        mL = aB * aP0;                                             // angular momentum
        double tG2 = 1 - aKappa * aKappa / (mL * mL);
        mG = std::sqrt(tG2);
        mC = mE * aKappa / (mL * mL * tG2);
        mA = std::sqrt(mC * mC + (mE * mE - aMass * aMass) / (mL * mL * tG2));
    }
    double R(double aPhi) const { return 1 / (mA * std::cos(mG * aPhi) - mC); }
    double Phi(double aR) const { return std::acos((1 / aR + mC) / mA) / mG; } // >= 0, outgoing branch
    // Angle between the momenta at the incoming and the outgoing point at distance r0 (mirror points of the orbit)
    double DeflectionAtDistance(double aR0) const
    {
        double tU0 = 1 / aR0;
        double tPhi0 = Phi(aR0);
        double tDuDphi = -mA * mG * std::sin(mG * tPhi0);
        // p_r = -L u', p_phi = L u
        TVector3 tOut = (-tDuDphi) * TVector3(std::cos(tPhi0), std::sin(tPhi0), 0)
                        + tU0 * TVector3(-std::sin(tPhi0), std::cos(tPhi0), 0);
        TVector3 tIn(-tOut.X(), tOut.Y(), 0); // mirror image, time reversed
        return tIn.Angle(tOut);
    }
};

// Radial motion r'' = K / r^2 from rest at r0 (two charges with K = kappa/mu, or a symmetric configuration):
// the time to reach r = x r0 is t(x) = sqrt(r0^3 / (2K)) [sqrt(x(x-1)) + ln(sqrt(x) + sqrt(x-1))].
// The functions take (mu, kappa) with K = kappa / mu.
static double NonRelRepulsionX(double aMu, double aKappa, double aR0, double aT)
{
    auto tOfX = [&](double x) {
        return std::sqrt(aMu * aR0 * aR0 * aR0 / (2 * aKappa)) * (std::sqrt(x * (x - 1)) + std::log(std::sqrt(x) + std::sqrt(x - 1)));
    };
    double lo = 1, hi = 2;
    while (tOfX(hi) < aT) hi *= 2;
    for (int i = 0; i < 200; ++i) { double m = 0.5 * (lo + hi); (tOfX(m) < aT ? lo : hi) = m; }
    return 0.5 * (lo + hi);
}
static double NonRelRepulsionDistance(double aMu, double aKappa, double aR0, double aT)
{
    return aR0 * NonRelRepulsionX(aMu, aKappa, aR0, aT);
}
static double NonRelRepulsionMomentum(double aMu, double aKappa, double aR0, double aT)
{
    double x = NonRelRepulsionX(aMu, aKappa, aR0, aT);
    return std::sqrt(2 * aMu * aKappa / aR0 * (1 - 1 / x)); // mu * dr/dt
}

// Lienard-Wiechert force on a unit test charge at rest, kappa = k q Q, from a source at retarded distance vector R (fm)
// with velocity v and acceleration a (fm^-1) at the retarded time: velocity and acceleration parts, GeV/fm.
static std::pair<TVector3, TVector3> ExactLienardWiechert(double aKappa, const TVector3 &aR, const TVector3 &aV, const TVector3 &aA)
{
    double tR = aR.Mag();
    TVector3 tU = (1 / tR) * aR - aV;
    double tRdotU = aR.Dot(tU);
    double tPre = aKappa * tR / (tRdotU * tRdotU * tRdotU);
    TVector3 tVel = (tPre * (1 - aV.Mag2())) * tU;
    TVector3 tAcc = tPre * (aR.Dot(aA) * tU - tRdotU * aA);
    return std::make_pair(tVel, tAcc);
}

// Instantaneous Coulomb forces between N charges with relativistic kinematics, by RK4 (x in fm, p in GeV, t in fm/c)
struct NBody {
    std::vector<double> mass, charge;
    std::vector<TVector3> x, p;
    void Derivative(const std::vector<TVector3> &aX, const std::vector<TVector3> &aP,
                    std::vector<TVector3> &aDX, std::vector<TVector3> &aDP) const
    {
        size_t n = aX.size();
        aDX.assign(n, TVector3());
        aDP.assign(n, TVector3());
        for (size_t i = 0; i < n; ++i) {
            aDX[i] = (1 / std::sqrt(mass[i] * mass[i] + aP[i].Mag2())) * aP[i];
            for (size_t j = 0; j < n; ++j) {
                if (i == j) continue;
                TVector3 tD = aX[i] - aX[j];
                double tR = tD.Mag();
                aDP[i] += (kAlphaHbarC * charge[i] * charge[j] / (tR * tR * tR)) * tD;
            }
        }
    }
    void Step(double aH)
    {
        size_t n = x.size();
        std::vector<TVector3> k1x, k1p, k2x, k2p, k3x, k3p, k4x, k4p, tX(n), tP(n);
        Derivative(x, p, k1x, k1p);
        for (size_t i = 0; i < n; ++i) { tX[i] = x[i] + (aH / 2) * k1x[i]; tP[i] = p[i] + (aH / 2) * k1p[i]; }
        Derivative(tX, tP, k2x, k2p);
        for (size_t i = 0; i < n; ++i) { tX[i] = x[i] + (aH / 2) * k2x[i]; tP[i] = p[i] + (aH / 2) * k2p[i]; }
        Derivative(tX, tP, k3x, k3p);
        for (size_t i = 0; i < n; ++i) { tX[i] = x[i] + aH * k3x[i]; tP[i] = p[i] + aH * k3p[i]; }
        Derivative(tX, tP, k4x, k4p);
        for (size_t i = 0; i < n; ++i) {
            x[i] += (aH / 6) * (k1x[i] + 2 * k2x[i] + 2 * k3x[i] + k4x[i]);
            p[i] += (aH / 6) * (k1p[i] + 2 * k2p[i] + 2 * k3p[i] + k4p[i]);
        }
    }
};

// Coulomb force constant S of a configuration of unit charges at the given positions (unit circumradius):
// the force on the charge at aVertices[0] is S kappa / R^2 directed radially. Returns S, and the sine of the
// angle between that force and the radial direction (0 for a symmetric configuration).
static double ConfigurationS(const std::vector<TVector3> &aVertices, double &aSinAngle)
{
    TVector3 tF;
    for (size_t j = 1; j < aVertices.size(); ++j) {
        TVector3 tD = aVertices[0] - aVertices[j];
        tF += (1 / std::pow(tD.Mag(), 3)) * tD;
    }
    aSinAngle = std::sin(tF.Angle(aVertices[0]));
    return tF.Mag();
}

// ---------------------------------------------------------------------------------------------

// Asymptotic deflection angle of the exact relativistic Coulomb orbit (see ExactOrbit) of a particle with energy aE
// (at infinite distance) and impact parameter aB: the orbit sweeps 2 phi_inf between the asymptotes, u = 0 at
// cos(g phi_inf) = C / A. Positive away from a repulsive source.
static double AsymptoticDeflection(double aMass, double aKappa, double aE, double aB)
{
    double tP = std::sqrt(aE * aE - aMass * aMass), tL = aB * tP;
    double tG2 = 1 - aKappa * aKappa / (tL * tL), tG = std::sqrt(tG2);
    double tC = aE * aKappa / (tL * tL * tG2);
    double tA = std::sqrt(tC * tC + tP * tP / (tL * tL * tG2));
    return TMath::Pi() - 2 * std::acos(tC / tA) / tG;
}

int main(int argc, char **argv)
{
    const char *tOutputFile = argc > 1 ? argv[1] : "coulomb_checks.root";
    const char *tTopDir = argc > 2 ? argv[2] : nullptr;
    Parameter tNoSpectators{"CoulombUseSpectators", "0"};
    sChecksConfig.AddParameter(&tNoSpectators);
    Parameter tSubsteps{"CoulombSubsteps", "1"};
    sSubstepConfig.AddParameter(&tNoSpectators);
    sSubstepConfig.AddParameter(&tSubsteps);

    ParticleDB tDB;
    ParticleType tTypes[] = {
        MakeType("proton", 2212, 0.938272, 1),        // also needed by the constructor
        MakeType("neutron", 2112, 0.939565, 0),       // also needed by the constructor; neutral decay product
        MakeType("heavyZ79", 9079, 1.0e5, 79),        // practically fixed source
        MakeType("heavyZ79daughter", 9080, 1.0e5, 79),
        MakeType("slowQ1", 9001, 10.0, 1),            // heavy unit charge, slow motion
        MakeType("probeQ1", 9002, 1.0e5, 1),          // practically fixed probe
        MakeType("platonicQ1", 9003, 100.0, 1),
        MakeType("clusterQ1", 9004, 20.0, 1),
        MakeType("movingZ10", 9010, 1.0e5, 10),
        MakeType("movingZ10daughter", 9011, 1.0e5, 10),
        MakeType("probeQm1", 9005, 1.0e5, -1),        // practically fixed negative probe
    };
    for (ParticleType &tType : tTypes) tDB.AddParticleType(&tType);
    // pointers into the database are valid only after all types are added
    ParticleType *tProton = tDB.FindByPID(2212), *tNeutron = tDB.FindByPID(2112);
    ParticleType *tHeavy = tDB.FindByPID(9079), *tHeavyDaughter = tDB.FindByPID(9080);
    ParticleType *tSlow = tDB.FindByPID(9001), *tProbe = tDB.FindByPID(9002);
    ParticleType *tPlatonic = tDB.FindByPID(9003), *tClusterType = tDB.FindByPID(9004);
    ParticleType *tMoving = tDB.FindByPID(9010), *tMovingDaughter = tDB.FindByPID(9011);
    ParticleType *tProbeNegative = tDB.FindByPID(9005);
    char tBuf[300];

    // ----- A. No interaction: one charged and one neutral particle
    Info("A. one charged + one neutral particle: no force");
    {
        const char *tDir = "A_free";
        std::vector<PartInit> tInit = {{tProton, TVector3(1, 2, 3), TVector3(0.1, 0.2, 0.3)},
                                       {tNeutron, TVector3(-1, 0, 0), TVector3(0, 0, 0.5)}};
        RunResult tRun = RunCoulomb(tDir, 50, 0.1, tInit, &tDB);
        for (int i = 0; i < 2; ++i) {
            snprintf(tBuf, sizeof(tBuf), "A: |p_final - p_initial| of particle %d [GeV]", i);
            Check(tBuf, (tRun.finalMomenta[i] - tInit[i].p).Mag(), 0, 1e-6, false);
            // Exact: straight line x0 + v t
            TVector3 tV = (1 / std::sqrt(tInit[i].type->GetMass() * tInit[i].type->GetMass() + tInit[i].p.Mag2())) * tInit[i].p;
            std::vector<double> tT, tCode[3], tExact[3];
            for (const TLorentzVector &tX : tRun.trajectory[i]) {
                tT.push_back(tX.T());
                TVector3 tXExact = tInit[i].x + tX.T() * tV;
                for (int c = 0; c < 3; ++c) { tCode[c].push_back(tX[c]); tExact[c].push_back(tXExact[c]); }
            }
            const char *tAxis[3] = {"x", "y", "z"};
            for (int c = 0; c < 3; ++c) {
                std::string tLabel = std::string("p") + std::to_string(i) + "_" + tAxis[c];
                AddCodeExact(tDir, tLabel, std::string(i == 0 ? "charged" : "neutral") + " particle, free motion",
                             std::string(tAxis[c]) + " (fm)", tT, tCode[c], tExact[c]);
            }
        }
        Check("A: positions and times restored (1 = yes)", tRun.restored, 1, 0, false);
    }

    // ----- Setup of the Rutherford scattering, used in B (no decay) and F (decay of the source)
    const double tRuthM = tProton->GetMass(); // as used by the code (float)
    const double tRuthKappa = 79 * kAlphaHbarC;
    const double tRuthB = 4.0, tRuthZ0 = 40.0;
    const double tRuthP0 = std::sqrt(std::pow(tRuthM + 0.200, 2) - tRuthM * tRuthM);
    const double tRuthR0 = std::hypot(tRuthB, tRuthZ0);
    const ExactOrbit tRuthOrbit(tRuthM, tRuthKappa, tRuthR0, tRuthP0, tRuthB);
    const double tRuthStep = 0.1;
    // Scattering plane: horizontal z, vertical x; polar angle psi = atan2(x, z). The orbit is traversed clockwise
    // (L_y < 0) and starts at phi = -phi0, so psi = psi_peri - phi.
    const double tRuthPsiPeri = std::atan2(tRuthB, -tRuthZ0) - tRuthOrbit.Phi(tRuthR0);
    int tRuthSteps;
    State tRuthRefFinal;
    std::vector<State> tRuthRefAtSteps;
    {
        // Time at which the proton is back at distance r0, from a fine RK4 integration
        State tS{TVector3(tRuthB, 0, -tRuthZ0), TVector3(0, 0, tRuthP0)};
        double tH = 1e-3, t = 0, tPrevR = tRuthR0;
        bool tPastPeriapsis = false;
        double tTReturn = 0;
        while (true) {
            State tNext = RK4Step(tS, tH, tRuthM, tRuthKappa);
            double tR = tNext.x.Mag();
            if (tR > tPrevR) tPastPeriapsis = true;
            if (tPastPeriapsis && tR >= tRuthR0) { tTReturn = t + tH * (tRuthR0 - tPrevR) / (tR - tPrevR); break; }
            tS = tNext; tPrevR = tR; t += tH;
        }
        tRuthSteps = (int)std::lround(tTReturn / tRuthStep);
        snprintf(tBuf, sizeof(tBuf), "Rutherford setup: back at r0 = %.2f fm after %.2f fm/c, run %d steps of %.1f fm/c",
                 tRuthR0, tTReturn, tRuthSteps, tRuthStep);
        Info(tBuf);
        // RK4 reference at the time steps of the code (100 substeps per step) and at the final time
        State tRef{TVector3(tRuthB, 0, -tRuthZ0), TVector3(0, 0, tRuthP0)};
        for (int n = 0; n < tRuthSteps; ++n) {
            tRuthRefAtSteps.push_back(tRef);
            for (int k = 0; k < 100; ++k) tRef = RK4Step(tRef, tRuthStep / 100, tRuthM, tRuthKappa);
        }
        tRuthRefFinal = tRef;
    }
    // Distance of the proton from the exact orbit at the same polar angle, along a trajectory
    auto tDistanceFromOrbit = [&](const std::vector<TLorentzVector> &aTraj, std::vector<double> &aT, std::vector<double> &aDr) {
        for (const TLorentzVector &tX : aTraj) {
            double tPsi = std::atan2(tX.X(), tX.Z());
            aT.push_back(tX.T());
            aDr.push_back(std::hypot(tX.X(), tX.Z()) - tRuthOrbit.R(tRuthPsiPeri - tPsi));
        }
    };
    const double tThetaExact = tRuthOrbit.DeflectionAtDistance(tRuthR0);
    RunResult tRunB;

    // ----- B. Relativistic Rutherford scattering on a fixed Z = 79 charge
    Info("B. proton, T = 200 MeV, b = 4 fm, on a fixed Z = 79 charge: exact relativistic orbit");
    {
        const char *tDir = "B_rutherford";
        std::vector<PartInit> tInit = {{tProton, TVector3(tRuthB, 0, -tRuthZ0), TVector3(0, 0, tRuthP0)},
                                       {tHeavy, TVector3(0, 0, 0), TVector3(0, 0, 0)}};
        tRunB = RunCoulomb(tDir, tRuthSteps, tRuthStep, tInit, &tDB);
        double tThetaRef = TVector3(0, 0, 1).Angle(tRuthRefFinal.p);
        double tThetaCode = TVector3(0, 0, 1).Angle(tRunB.finalMomenta[0]);
        Check("B: reference RK4 deflection vs exact orbit [rad]", tThetaRef, tThetaExact, 1e-4);
        Check("B: code deflection vs exact orbit [rad]", tThetaCode, tThetaExact, 3e-5);
        Check("B: code |p| vs initial |p| (energy conservation) [GeV]", tRunB.finalMomenta[0].Mag(), tRuthP0, 3e-5);
        Check("B: code |p - p_RK4| / |p|", (tRunB.finalMomenta[0] - tRuthRefFinal.p).Mag() / tRuthP0, 0, 1e-5, false);
        Check("B: positions and times restored (1 = yes)", tRunB.restored, 1, 0, false);

        // Graphs: paths in the scattering plane, distance from the exact orbit, and x, z, p_x, p_z vs RK4
        std::vector<double> tPathZ, tPathX;
        double tPhi0 = tRuthOrbit.Phi(tRuthR0);
        for (int k = 0; k <= 1000; ++k) {
            double tPhi = -tPhi0 + 2 * tPhi0 * k / 1000;
            double tR = tRuthOrbit.R(tPhi), tPsi = tRuthPsiPeri - tPhi;
            tPathZ.push_back(tR * std::cos(tPsi));
            tPathX.push_back(tR * std::sin(tPsi));
        }
        AddGraph(tDir, "path_exact", "Exact relativistic Coulomb orbit;z (fm);x (fm)", tPathZ, tPathX);
        std::vector<double> tTimes, tCodeZ, tCodeX, tRefZ, tRefX, tCodePx, tRefPx, tCodePz, tRefPz, tDrT, tDrOrbit;
        const std::vector<TLorentzVector> &tTraj = tRunB.trajectory[0];
        const std::vector<TLorentzVector> &tMom = tRunB.fourMomenta[0];
        for (size_t n = 0; n < tTraj.size(); ++n) {
            tTimes.push_back(tTraj[n].T());
            tCodeZ.push_back(tTraj[n].Z());
            tCodeX.push_back(tTraj[n].X());
            tRefZ.push_back(tRuthRefAtSteps[n].x.Z());
            tRefX.push_back(tRuthRefAtSteps[n].x.X());
            tCodePx.push_back(tMom[n].Px());
            tRefPx.push_back(tRuthRefAtSteps[n].p.X());
            tCodePz.push_back(tMom[n].Pz());
            tRefPz.push_back(tRuthRefAtSteps[n].p.Z());
        }
        tDistanceFromOrbit(tTraj, tDrT, tDrOrbit);
        AddGraph(tDir, "path_code", "Trajectory (code);z (fm);x (fm)", tCodeZ, tCodeX);
        AddGraph(tDir, "path_rk4", "Trajectory (RK4 of the exact equation of motion);z (fm);x (fm)", tRefZ, tRefX);
        AddGraph(tDir, "dr_orbit", "r(code) - r(exact orbit) at the same polar angle;t (fm/c);#Deltar (fm)", tDrT, tDrOrbit);
        AddCodeExact(tDir, "x", "proton x, exact = RK4", "x (fm)", tTimes, tCodeX, tRefX);
        AddCodeExact(tDir, "z", "proton z, exact = RK4", "z (fm)", tTimes, tCodeZ, tRefZ);
        AddCodeExact(tDir, "px", "proton p_{x}, exact = RK4", "p_{x} (GeV)", tTimes, tCodePx, tRefPx);
        AddCodeExact(tDir, "pz", "proton p_{z}, exact = RK4", "p_{z} (GeV)", tTimes, tCodePz, tRefPz);
        Check("B: max |r(code) - r(exact orbit)| along the path [fm]", MaxAbs(tDrOrbit), 0, 1e-3, false);
    }

    // ----- C. Two equal charges released from rest (slow, heavy)
    Info("C. two charges q = 1, m = 10 GeV, released from rest 1 fm apart, 60 fm/c: non-relativistic solution");
    RunResult tRunRest;
    double tMSlow = tSlow->GetMass();
    double tKappa1 = kAlphaHbarC;
    int tStepsC = 600;
    double tStepC = 0.1;
    {
        const char *tDir = "C_rest";
        std::vector<PartInit> tInit = {{tSlow, TVector3(-0.5, 0, 0), TVector3()},
                                       {tSlow, TVector3(0.5, 0, 0), TVector3()}};
        tRunRest = RunCoulomb(tDir, tStepsC, tStepC, tInit, &tDB);
        const std::vector<TVector3> &tP = tRunRest.finalMomenta;
        double tPExact = NonRelRepulsionMomentum(tMSlow / 2, tKappa1, 1.0, tStepsC * tStepC);
        snprintf(tBuf, sizeof(tBuf), "C: final speed of each particle v = %.4f (relativistic corrections ~ v^2 = %.1e)",
                 tPExact / tMSlow, std::pow(tPExact / tMSlow, 2));
        Info(tBuf);
        Check("C: p_x of particle 2 vs non-relativistic solution [GeV]", tP[1].X(), tPExact, 5e-4);
        Check("C: p_x of particle 1 = -p_x of particle 2 [GeV]", tP[0].X(), -tP[1].X(), 1e-5);
        Check("C: transverse momenta |p_y| + |p_z| [GeV]",
              std::fabs(tP[0].Y()) + std::fabs(tP[0].Z()) + std::fabs(tP[1].Y()) + std::fabs(tP[1].Z()), 0, 1e-9, false);

        std::vector<double> tT, tX2, tX2Exact, tPx2, tPx2Exact, tPxSum;
        for (size_t n = 0; n < tRunRest.trajectory[1].size(); ++n) {
            double t = tRunRest.trajectory[1][n].T();
            tT.push_back(t);
            tX2.push_back(tRunRest.trajectory[1][n].X());
            tX2Exact.push_back(0.5 * NonRelRepulsionDistance(tMSlow / 2, tKappa1, 1.0, t));
            tPx2.push_back(tRunRest.fourMomenta[1][n].Px());
            tPx2Exact.push_back(NonRelRepulsionMomentum(tMSlow / 2, tKappa1, 1.0, t));
            tPxSum.push_back(tRunRest.fourMomenta[0][n].Px() + tRunRest.fourMomenta[1][n].Px());
        }
        AddCodeExact(tDir, "p2_x", "x of particle 2, exact = non-relativistic", "x (fm)", tT, tX2, tX2Exact);
        AddCodeExact(tDir, "p2_px", "p_{x} of particle 2, exact = non-relativistic", "p_{x} (GeV)", tT, tPx2, tPx2Exact);
        AddGraph(tDir, "px_sum", "p_{x}(1) + p_{x}(2) (code), exact = 0;t (fm/c);p_{x} (GeV)", tT, tPxSum);
    }

    // ----- D. The same pair moving with beta = 0.9 along z (Lorentz covariance)
    Info("D. pair of C boosted to beta = 0.9 along z, lab step = gamma x 0.1 fm/c: p_T(lab, t) = p_T(rest, t/gamma)");
    {
        const char *tDir = "D_boosted";
        double tBeta = 0.9, tGamma = 1 / std::sqrt(1 - tBeta * tBeta);
        double tPz = tGamma * tMSlow * tBeta;
        std::vector<PartInit> tInit = {{tSlow, TVector3(-0.5, 0, 0), TVector3(0, 0, tPz)},
                                       {tSlow, TVector3(0.5, 0, 0), TVector3(0, 0, tPz)}};
        RunResult tRun = RunCoulomb(tDir, tStepsC, tGamma * tStepC, tInit, &tDB);
        const std::vector<TVector3> &tP = tRun.finalMomenta;
        const std::vector<TVector3> &tPRest = tRunRest.finalMomenta;
        Check("D: p_x(lab) vs p_x(rest frame run C) [GeV]", tP[1].X(), tPRest[1].X(), 1e-5);
        Check("D: p_x(lab) vs non-relativistic solution in rest frame [GeV]", tP[1].X(),
              NonRelRepulsionMomentum(tMSlow / 2, tKappa1, 1.0, tStepsC * tStepC), 5e-4);
        // p_z(lab) = gamma beta E'(rest): the kinetic energy gained in the rest frame shows up in p_z. This is sensitive
        // to the acceleration term. The change, 4.6e-4 GeV, is small compared to p_z = 20.6 GeV stored as Float_t.
        double tDPzExpected = tGamma * tBeta * (std::sqrt(tMSlow * tMSlow + tPRest[1].Mag2()) - tMSlow);
        Check("D: change of p_z(lab) vs gamma beta x kinetic energy(rest) [GeV]", tP[1].Z() - tPz, tDPzExpected, 3e-3);
        Check("D: p_x of particle 1 = -p_x of particle 2 [GeV]", tP[0].X(), -tP[1].X(), 1e-5);

        // Lorentz transformation of the rest frame motion (t', x', z' = 0): t = gamma t', x = x'(t/gamma), z = beta t.
        // Step n of the lab run corresponds to step n of the rest frame run C (same proper time).
        std::vector<double> tT, tX2, tX2Exact, tX2RunC, tZ2, tZ2Exact, tPx2, tPx2RunC;
        for (size_t n = 0; n < tRun.trajectory[1].size(); ++n) {
            double t = tRun.trajectory[1][n].T();
            tT.push_back(t);
            tX2.push_back(tRun.trajectory[1][n].X());
            tX2Exact.push_back(0.5 * NonRelRepulsionDistance(tMSlow / 2, tKappa1, 1.0, t / tGamma));
            tX2RunC.push_back(tRunRest.trajectory[1][n].X());
            tZ2.push_back(tRun.trajectory[1][n].Z());
            tZ2Exact.push_back(tBeta * t);
            tPx2.push_back(tRun.fourMomenta[1][n].Px());
            tPx2RunC.push_back(tRunRest.fourMomenta[1][n].Px());
        }
        AddCodeExact(tDir, "p2_x", "x of particle 2 (lab), exact = non-relativistic x'(t/#gamma)", "x (fm)", tT, tX2, tX2Exact);
        AddCodeExact(tDir, "p2_x_covariance", "x of particle 2 (lab), exact = run C at t/#gamma", "x (fm)", tT, tX2, tX2RunC);
        AddCodeExact(tDir, "p2_z", "z of particle 2 (lab), exact = #beta t", "z (fm)", tT, tZ2, tZ2Exact);
        AddCodeExact(tDir, "p2_px_covariance", "p_{x} of particle 2 (lab), exact = run C at t/#gamma", "p_{x} (GeV)", tT, tPx2, tPx2RunC);
    }

    // ----- E. Acceleration (radiation) term in isolation
    // A charge in hyperbolic motion z(t) = sqrt(alpha^2 + t^2) along z (proper acceleration 1/alpha), recorded on the
    // time grid like in Apply() (accelerations as backward differences of velocities), seen by a test charge at rest
    // at x = 30 fm. The acceleration part of the code's force is F(recorded accelerations) - F(accelerations = 0);
    // the reference is the exact Lienard-Wiechert field with the exact retarded time.
    Info("E. acceleration term: charge in hyperbolic motion (alpha = 2 fm) seen from 30 fm, exact Lienard-Wiechert field");
    {
        const char *tDir = "E_radiation";
        mkdir(tDir, 0755);
        if (chdir(tDir) != 0) { perror(tDir); exit(2); }
        {
            const double tAlpha = 2.0, tT0 = -40.0, tStep = 0.1;
            const unsigned int tNPoints = 801;
            auto tSourceZ = [&](double t) { return std::sqrt(tAlpha * tAlpha + t * t); };
            auto tSourceV = [&](double t) { return t / std::sqrt(tAlpha * tAlpha + t * t); };
            auto tSourceA = [&](double t) { return tAlpha * tAlpha / std::pow(tAlpha * tAlpha + t * t, 1.5); };

            CoulombAfterburner tCoulomb(1, tStep, 0., nullptr, &tDB);
            tCoulomb.Init();
            // The internal state of the afterburner has units (CoulombUnits.h)
            using namespace CoulombUnits;
            tCoulomb.m_coordinates = new std::vector<CoulombAfterburner::SpaceTimePoint>[2];
            tCoulomb.m_velocities = new std::vector<Vec3<Speed> >[2];
            tCoulomb.m_accelerations = new std::vector<Vec3<Acceleration> >[2];
            tCoulomb.m_FirstStepExists = new unsigned int[2]{0, 0};
            tCoulomb.m_LastStepExists = new unsigned int[2]{tNPoints - 1, tNPoints - 1};
            tCoulomb.m_CreationTimes = new Time[2]{tT0 * fm_c, tT0 * fm_c};
            tCoulomb.m_DecayTimes = new Time[2]{std::numeric_limits<double>::infinity() * fm_c, std::numeric_limits<double>::infinity() * fm_c};
            tCoulomb.m_termsCoul->emplace_back(); tCoulomb.m_termsVelo->emplace_back(); tCoulomb.m_termsAcce->emplace_back();
            tCoulomb.m_termsMagn->emplace_back(); tCoulomb.m_Rs->emplace_back(); tCoulomb.m_us->emplace_back(); tCoulomb.m_vs->emplace_back();
            for (unsigned int n = 0; n < tNPoints; ++n) {
                double t = tT0 + n * tStep;
                tCoulomb.m_coordinates[1].push_back(CoulombAfterburner::SpaceTimePoint{t * fm_c, Vec3<Length>(0 * fm, 0 * fm, tSourceZ(t) * fm)});
                tCoulomb.m_velocities[1].push_back(Vec3<Speed>(0 * c, 0 * c, tSourceV(t) * c));
                tCoulomb.m_accelerations[1].push_back(n == 0 ? Vec3<Acceleration>()
                    : Vec3<Acceleration>((tCoulomb.m_velocities[1][n] - tCoulomb.m_velocities[1][n - 1]) / (tStep * fm_c)));
            }
            Particle::ZeroEID();
            std::list<Particle> tList;
            for (int i = 0; i < 2; ++i) {
                Particle tParticle(tProton);
                tParticle.SetParticlePX(tProton->GetMass(), 0, 0, 0, tT0 / kHbarC, 0, 0, tSourceZ(tT0) / kHbarC);
                tList.push_back(tParticle);
            }
            auto tObserver = tList.begin();
            int **tBestIs = new int *[2];
            for (int i = 0; i < 2; ++i) tBestIs[i] = new int[2]{-1, -1};

            const TVector3 tX(30, 0, 0); // fm
            std::vector<double> tTObs, tAccX[2], tAccZ[2], tVelX[2], tVelZ[2], tTRet;
            for (double t = 15; t <= 45 + 1e-9; t += 0.25) {
                tObserver->SetTXYZ(t / kHbarC, tX.X() / kHbarC, tX.Y() / kHbarC, tX.Z() / kHbarC);
                auto tInGeVPerFm = [](const Vec3<Force> &aF) {
                    return TVector3(aF.x.numerical_value_in(GeV / fm), aF.y.numerical_value_in(GeV / fm), aF.z.numerical_value_in(GeV / fm));
                };
                TVector3 tFull = tInGeVPerFm(tCoulomb.CalculateTotalForce(tObserver, &tList, tBestIs, t * fm_c, tNPoints - 1));
                std::vector<Vec3<Acceleration> > tSaved(tNPoints);
                tSaved.swap(tCoulomb.m_accelerations[1]); // now zero accelerations
                TVector3 tVel = tInGeVPerFm(tCoulomb.CalculateTotalForce(tObserver, &tList, tBestIs, t * fm_c, tNPoints - 1));
                tSaved.swap(tCoulomb.m_accelerations[1]);
                TVector3 tAcc = tFull - tVel;
                // Exact retarded time: (t - t') = |x - z(t')|, the left minus right side decreases with t'
                double lo = t - 200, hi = t;
                for (int k = 0; k < 200; ++k) {
                    double m = 0.5 * (lo + hi);
                    double f = (t - m) - (tX - TVector3(0, 0, tSourceZ(m))).Mag();
                    (f > 0 ? lo : hi) = m;
                }
                double tR = 0.5 * (lo + hi);
                std::pair<TVector3, TVector3> tExact = ExactLienardWiechert(
                    kAlphaHbarC, tX - TVector3(0, 0, tSourceZ(tR)), TVector3(0, 0, tSourceV(tR)), TVector3(0, 0, tSourceA(tR)));
                tTObs.push_back(t);
                tTRet.push_back(tR);
                tAccX[0].push_back(tAcc.X()); tAccX[1].push_back(tExact.second.X());
                tAccZ[0].push_back(tAcc.Z()); tAccZ[1].push_back(tExact.second.Z());
                tVelX[0].push_back(tVel.X()); tVelX[1].push_back(tExact.first.X());
                tVelZ[0].push_back(tVel.Z()); tVelZ[1].push_back(tExact.first.Z());
            }
            for (int i = 0; i < 2; ++i) delete[] tBestIs[i];
            delete[] tBestIs;
            delete[] tCoulomb.m_coordinates;
            delete[] tCoulomb.m_velocities;
            delete[] tCoulomb.m_accelerations;
            delete[] tCoulomb.m_FirstStepExists;
            delete[] tCoulomb.m_LastStepExists;
            delete[] tCoulomb.m_CreationTimes;
            delete[] tCoulomb.m_DecayTimes;
            tCoulomb.Finalize();

            auto tMaxDiff = [](const std::vector<double> &a, const std::vector<double> &b) {
                double tMax = 0;
                for (size_t i = 0; i < a.size(); ++i) tMax = std::max(tMax, std::fabs(a[i] - b[i]));
                return tMax;
            };
            double tScaleAcc = std::max(MaxAbs(tAccX[1]), MaxAbs(tAccZ[1]));
            double tScaleVel = std::max(MaxAbs(tVelX[1]), MaxAbs(tVelZ[1]));
            snprintf(tBuf, sizeof(tBuf), "E: retarded times %.1f to %.1f fm/c (source speed up to %.3f), max |F_acc| / max |F_vel| = %.1f",
                     tTRet.front(), tTRet.back(), std::max(std::fabs(tSourceV(tTRet.front())), std::fabs(tSourceV(tTRet.back()))),
                     tScaleAcc / tScaleVel);
            Info(tBuf);
            Check("E: acceleration term dominates (max |F_acc| / max |F_vel| > 3)", tScaleAcc / tScaleVel > 3, 1, 0, false);
            Check("E: acceleration part, max |code - exact| / max |exact|",
                  std::max(tMaxDiff(tAccX[0], tAccX[1]), tMaxDiff(tAccZ[0], tAccZ[1])) / tScaleAcc, 0, 3e-3, false);
            Check("E: velocity part, max |code - exact| / max |exact|",
                  std::max(tMaxDiff(tVelX[0], tVelX[1]), tMaxDiff(tVelZ[0], tVelZ[1])) / tScaleVel, 0, 1e-3, false);
            AddCodeExact(tDir, "Fx_acc", "acceleration part of F_{x} on the observer", "F_{x} (GeV/fm)", tTObs, tAccX[0], tAccX[1]);
            AddCodeExact(tDir, "Fz_acc", "acceleration part of F_{z} on the observer", "F_{z} (GeV/fm)", tTObs, tAccZ[0], tAccZ[1]);
            AddCodeExact(tDir, "Fx_vel", "velocity part of F_{x} on the observer", "F_{x} (GeV/fm)", tTObs, tVelX[0], tVelX[1]);
            AddCodeExact(tDir, "Fz_vel", "velocity part of F_{z} on the observer", "F_{z} (GeV/fm)", tTObs, tVelZ[0], tVelZ[1]);
            AddGraph(tDir, "t_retarded", "Exact retarded time of the source;t observer (fm/c);t retarded (fm/c)", tTObs, tTRet);
        }
        if (chdir("..") != 0) { perror(".."); exit(2); }
    }

    // ----- F. Decay of the static source during the Rutherford scattering
    // The Z = 79 source decays into a daughter with the same charge and mass at rest at the same point, and a neutral
    // particle, close to the proton's closest approach. The proton must not notice: exact orbit of B.
    Info("F. Rutherford scattering of B, the source decays into an identical charged + a neutral daughter at closest approach");
    {
        const char *tDir = "F_decay_static";
        const double tDecayTimes[2] = {70.0, 70.05}; // on and between time steps
        for (int iDecay = 0; iDecay < 2; ++iDecay) {
            double tTD = tDecayTimes[iDecay];
            std::vector<PartInit> tInit = {{tProton, TVector3(tRuthB, 0, -tRuthZ0), TVector3(0, 0, tRuthP0)},
                                           {tHeavy, TVector3(0, 0, 0), TVector3(0, 0, 0), 0., -1, true},
                                           {tHeavyDaughter, TVector3(0, 0, 0), TVector3(0, 0, 0), tTD, 1},
                                           {tNeutron, TVector3(0, 0, 0), TVector3(0, 0, 0), tTD, 1}};
            std::string tSub = std::string(tDir) + "_" + std::to_string(iDecay);
            RunResult tRun = RunCoulomb(tSub.c_str(), tRuthSteps, tRuthStep, tInit, &tDB);
            std::string tWhat = Form("F (decay at %.2f fm/c): ", tTD);
            Check((tWhat + "code deflection vs exact orbit [rad]").c_str(), TVector3(0, 0, 1).Angle(tRun.finalMomenta[0]), tThetaExact, 3e-5);
            Check((tWhat + "|p| vs initial |p| [GeV]").c_str(), tRun.finalMomenta[0].Mag(), tRuthP0, 3e-5);
            Check((tWhat + "|p - p(B, no decay)| / |p|").c_str(), (tRun.finalMomenta[0] - tRunB.finalMomenta[0]).Mag() / tRuthP0, 0, 1e-5, false);
            Check((tWhat + "positions and times restored (1 = yes)").c_str(), tRun.restored, 1, 0, false);
            std::vector<double> tDrT, tDr, tT, tPx, tPxB;
            tDistanceFromOrbit(tRun.trajectory[0], tDrT, tDr);
            Check((tWhat + "max |r(code) - r(exact orbit)| [fm]").c_str(), MaxAbs(tDr), 0, 1e-3, false);
            for (size_t n = 0; n < tRun.fourMomenta[0].size(); ++n) {
                tT.push_back(tRun.trajectory[0][n].T());
                tPx.push_back(tRun.fourMomenta[0][n].Px());
                tPxB.push_back(tRunB.fourMomenta[0][n].Px());
            }
            std::string tLabel = Form("td%.2f", tTD);
            AddGraph(tDir, tLabel + "_dr_orbit", Form("decay at %.2f fm/c: r(code) - r(exact orbit) at the same polar angle;t (fm/c);#Deltar (fm)", tTD), tDrT, tDr);
            AddCodeExact(tDir, tLabel + "_px", Form("proton p_{x}, source decays at %.2f fm/c, exact = run B without decay", tTD),
                         "p_{x} (GeV)", tT, tPx, tPxB);
        }
    }

    // ----- G. Decay of a moving source passing a probe
    // A Z = 10 source moves with v = 0.5 along z at impact parameter b = 2 fm from a practically fixed unit-charge
    // probe, from 60 fm/c before to 60 fm/c after the closest approach; it decays at the closest approach into a
    // daughter with the same charge, mass and velocity, and a neutral particle. Exact transverse impulse on the probe
    // from the field of a uniformly moving charge: dp_x(t) = kappa/(b v) [g v s/sqrt(b^2 + g^2 v^2 s^2)]_{s=-60}^{t-60}.
    Info("G. Z = 10 source, v = 0.5, passes a probe at b = 2 fm and decays at closest approach: exact impulse");
    {
        const char *tDir = "G_decay_moving";
        const double tV = 0.5, tGamma = 1 / std::sqrt(1 - tV * tV), tB = 2.0, tTC = 60.0, tStep = 0.1;
        const int tSteps = 1200;
        const double tKappa = 10 * kAlphaHbarC;
        const double tPz = tGamma * tMoving->GetMass() * tV;
        auto tImpulse = [&](double t) {
            auto tF = [&](double s) { return tGamma * tV * s / std::sqrt(tB * tB + tGamma * tGamma * tV * tV * s * s); };
            return tKappa / (tB * tV) * (tF(t - tTC) - tF(-tTC));
        };
        std::vector<PartInit> tNoDecay = {{tProbe, TVector3(tB, 0, 0), TVector3()},
                                          {tMoving, TVector3(0, 0, -tV * tTC), TVector3(0, 0, tPz)}};
        RunResult tRunNoDecay = RunCoulomb((std::string(tDir) + "_nodecay").c_str(), tSteps, tStep, tNoDecay, &tDB);
        const double tTD = 60.05;
        std::vector<PartInit> tDecay = {{tProbe, TVector3(tB, 0, 0), TVector3()},
                                        {tMoving, TVector3(0, 0, -tV * tTC), TVector3(0, 0, tPz), 0., -1, true},
                                        {tMovingDaughter, TVector3(0, 0, tV * (tTD - tTC)), TVector3(0, 0, tPz), tTD, 1},
                                        {tNeutron, TVector3(0, 0, tV * (tTD - tTC)), TVector3(), tTD, 1}};
        RunResult tRun = RunCoulomb(tDir, tSteps, tStep, tDecay, &tDB);
        double tExact = tImpulse(tSteps * tStep);
        Check("G: no decay, probe p_x vs exact impulse [GeV]", tRunNoDecay.finalMomenta[0].X(), tExact, 3e-4);
        Check("G: decay, probe p_x vs exact impulse [GeV]", tRun.finalMomenta[0].X(), tExact, 3e-4);
        Check("G: decay, probe |p_z| / p_x (exact 0 by symmetry)", std::fabs(tRun.finalMomenta[0].Z()) / tExact, 0, 1e-4, false);
        // Tolerance: moving from z = -30 fm, the Float_t position of the father drifts ~1e-4 fm ahead of the exact line,
        // while the daughter is created at the exact position; the source jumping back changes the impulse by ~2e-5.
        Check("G: |p(decay) - p(no decay)| / p_x", (tRun.finalMomenta[0] - tRunNoDecay.finalMomenta[0]).Mag() / tExact, 0, 5e-5, false);
        std::vector<double> tT, tPx, tPxExact, tPxNoDecay;
        for (size_t n = 0; n < tRun.fourMomenta[0].size(); ++n) {
            double t = tRun.trajectory[0][n].T();
            tT.push_back(t);
            tPx.push_back(tRun.fourMomenta[0][n].Px());
            tPxExact.push_back(tImpulse(t));
            tPxNoDecay.push_back(tRunNoDecay.fourMomenta[0][n].Px());
        }
        AddCodeExact(tDir, "probe_px", "probe p_{x}, source decays at closest approach, exact = uniformly moving charge",
                     "p_{x} (GeV)", tT, tPx, tPxExact);
        AddCodeExact(tDir, "probe_px_vs_nodecay", "probe p_{x}, source decays at closest approach, exact = run without decay",
                     "p_{x} (GeV)", tT, tPx, tPxNoDecay);
    }

    // ----- H. Platonic solids of like charges released from rest
    // All vertices are equivalent, so the expansion is radial and self-similar: r'' = S kappa / (m r^2), with S from
    // the geometry. Heavy charges (100 GeV) keep the motion non-relativistic.
    Info("H. like charges (m = 100 GeV) at the vertices of Platonic solids, circumradius 2 fm, released from rest, 100 fm/c");
    {
        const char *tDir = "H_platonic";
        const double tPhi = (1 + std::sqrt(5.)) / 2;
        std::vector<std::pair<std::string, std::vector<TVector3> > > tSolids;
        tSolids.push_back({"tetrahedron", {{1, 1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, 1}}});
        tSolids.push_back({"octahedron", {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}}});
        std::vector<TVector3> tCube, tIcosahedron;
        for (int i = 0; i < 8; ++i) tCube.push_back(TVector3(i & 1 ? 1 : -1, i & 2 ? 1 : -1, i & 4 ? 1 : -1));
        for (int s1 : {-1, 1}) for (int s2 : {-1, 1}) {
            tIcosahedron.push_back(TVector3(0, s1, s2 * tPhi));
            tIcosahedron.push_back(TVector3(s1, s2 * tPhi, 0));
            tIcosahedron.push_back(TVector3(s2 * tPhi, 0, s1));
        }
        tSolids.push_back({"cube", tCube});
        tSolids.push_back({"icosahedron", tIcosahedron});
        const double tR0 = 2.0, tStep = 0.1, tM = tPlatonic->GetMass();
        const int tSteps = 1000;
        for (auto &tSolid : tSolids) {
            std::vector<TVector3> &tVertices = tSolid.second;
            for (TVector3 &tV : tVertices) tV = tV.Unit();
            double tSinAngle;
            double tS = ConfigurationS(tVertices, tSinAngle);
            std::vector<PartInit> tInit;
            for (const TVector3 &tV : tVertices) tInit.push_back({tPlatonic, tR0 * tV, TVector3()});
            std::string tName = tSolid.first;
            RunResult tRun = RunCoulomb((std::string(tDir) + "_" + tName).c_str(), tSteps, tStep, tInit, &tDB);
            double tPExact = NonRelRepulsionMomentum(tM, tS * kAlphaHbarC, tR0, tSteps * tStep);
            double tPMean = 0, tPMin = 1e30, tPMax = 0, tMaxAngle = 0;
            TVector3 tPSum;
            for (size_t i = 0; i < tVertices.size(); ++i) {
                const TVector3 &tP = tRun.finalMomenta[i];
                tPMean += tP.Mag() / tVertices.size();
                tPMin = std::min(tPMin, tP.Mag());
                tPMax = std::max(tPMax, tP.Mag());
                tMaxAngle = std::max(tMaxAngle, tP.Angle(tVertices[i]));
                tPSum += tP;
            }
            snprintf(tBuf, sizeof(tBuf), "H: %s, %zu charges, S = %.5f (symmetry: sin of force angle %.1e), final v = %.4f",
                     tName.c_str(), tVertices.size(), tS, tSinAngle, tPExact / tM);
            Info(tBuf);
            Check(("H: " + tName + ", mean |p| vs non-relativistic radial solution [GeV]").c_str(), tPMean, tPExact, 1e-3);
            Check(("H: " + tName + ", spread (max - min) / mean of |p|").c_str(), (tPMax - tPMin) / tPMean, 0, 1e-5, false);
            Check(("H: " + tName + ", max angle between p and radial direction [rad]").c_str(), tMaxAngle, 0, 1e-5, false);
            Check(("H: " + tName + ", |sum of p| / mean |p|").c_str(), tPSum.Mag() / tPMean, 0, 1e-5, false);
            std::vector<double> tT, tRCode, tRExact, tPrCode, tPrExact;
            for (size_t n = 0; n < tRun.trajectory[0].size(); ++n) {
                double t = tRun.trajectory[0][n].T();
                tT.push_back(t);
                tRCode.push_back(tRun.trajectory[0][n].Vect().Mag());
                tRExact.push_back(NonRelRepulsionDistance(tM, tS * kAlphaHbarC, tR0, t));
                tPrCode.push_back(tRun.fourMomenta[0][n].Vect().Dot(tVertices[0]));
                tPrExact.push_back(NonRelRepulsionMomentum(tM, tS * kAlphaHbarC, tR0, t));
            }
            AddCodeExact(tDir, tName + "_r", tName + ": distance of charge 0 from the centre, exact = non-relativistic", "r (fm)", tT, tRCode, tRExact);
            AddCodeExact(tDir, tName + "_pr", tName + ": radial momentum of charge 0, exact = non-relativistic", "p_{r} (GeV)", tT, tPrCode, tPrExact);
        }
    }

    // ----- I. Hexagon of charges, at rest and boosted to beta = 0.9 along its axis (Lorentz covariance)
    Info("I. hexagon of charges (m = 10 GeV), radius 1 fm, at rest and boosted to beta = 0.9 along z, 600 steps");
    {
        const char *tDir = "I_hexagon";
        const int tN = 6, tSteps = 600;
        const double tR0 = 1.0, tStep = 0.1, tBeta = 0.9, tGamma = 1 / std::sqrt(1 - tBeta * tBeta);
        std::vector<TVector3> tVertices;
        for (int i = 0; i < tN; ++i) tVertices.push_back(TVector3(std::cos(2 * M_PI * i / tN), std::sin(2 * M_PI * i / tN), 0));
        double tSinAngle;
        double tS = ConfigurationS(tVertices, tSinAngle);
        std::vector<PartInit> tRest, tLab;
        for (const TVector3 &tV : tVertices) {
            tRest.push_back({tSlow, tR0 * tV, TVector3()});
            tLab.push_back({tSlow, tR0 * tV, TVector3(0, 0, tGamma * tMSlow * tBeta)});
        }
        RunResult tRunRestHex = RunCoulomb((std::string(tDir) + "_rest").c_str(), tSteps, tStep, tRest, &tDB);
        RunResult tRunLab = RunCoulomb((std::string(tDir) + "_lab").c_str(), tSteps, tGamma * tStep, tLab, &tDB);
        double tPExact = NonRelRepulsionMomentum(tMSlow, tS * kAlphaHbarC, tR0, tSteps * tStep);
        double tMaxRel = 0, tMaxCov = 0;
        for (int i = 0; i < tN; ++i) {
            TVector3 tPRest = tRunRestHex.finalMomenta[i], tPLab = tRunLab.finalMomenta[i];
            tMaxRel = std::max(tMaxRel, std::fabs(tPRest.Mag() - tPExact) / tPExact);
            tMaxCov = std::max(tMaxCov, (TVector3(tPLab.X(), tPLab.Y(), 0) - TVector3(tPRest.X(), tPRest.Y(), 0)).Mag() / tPRest.Mag());
        }
        snprintf(tBuf, sizeof(tBuf), "I: S = %.5f, final v = %.4f", tS, tPExact / tMSlow);
        Info(tBuf);
        Check("I: rest frame, max | |p| - exact | / exact over the charges", tMaxRel, 0, 2e-3, false);
        Check("I: max |p_T(lab) - p_T(rest)| / |p| over the charges", tMaxCov, 0, 1e-5, false);
        std::vector<double> tT, tPrLab, tPrRest, tPrExact;
        for (size_t n = 0; n < tRunLab.fourMomenta[0].size(); ++n) {
            double t = tRunLab.trajectory[0][n].T();
            tT.push_back(t);
            tPrLab.push_back(tRunLab.fourMomenta[0][n].Vect().Dot(tVertices[0]));
            tPrRest.push_back(tRunRestHex.fourMomenta[0][n].Vect().Dot(tVertices[0]));
            tPrExact.push_back(NonRelRepulsionMomentum(tMSlow, tS * kAlphaHbarC, tR0, t / tGamma));
        }
        AddCodeExact(tDir, "pr_covariance", "radial momentum of charge 0 (lab), exact = rest frame run at t/#gamma", "p_{r} (GeV)", tT, tPrLab, tPrRest);
        AddCodeExact(tDir, "pr", "radial momentum of charge 0 (lab), exact = non-relativistic at t/#gamma", "p_{r} (GeV)", tT, tPrLab, tPrExact);
    }

    // ----- J. Random cluster of like charges vs instantaneous Coulomb N-body integration
    // No symmetry; the retarded interaction differs from the instantaneous one by O(v^2) effects.
    Info("J. 20 like charges (m = 20 GeV) at random in a sphere of 4 fm, from rest, 60 fm/c: instantaneous Coulomb N-body RK4");
    {
        const char *tDir = "J_cluster";
        const int tN = 20, tSteps = 600;
        const double tStep = 0.1, tM = tClusterType->GetMass();
        TRandom3 tRandom(4357);
        std::vector<TVector3> tPositions;
        while ((int)tPositions.size() < tN) {
            TVector3 tX(tRandom.Uniform(-4, 4), tRandom.Uniform(-4, 4), tRandom.Uniform(-4, 4));
            if (tX.Mag() > 4) continue;
            bool tTooClose = false;
            for (const TVector3 &tY : tPositions) if ((tX - tY).Mag() < 1.0) tTooClose = true;
            if (!tTooClose) tPositions.push_back(tX);
        }
        std::vector<PartInit> tInit;
        NBody tRef;
        for (const TVector3 &tX : tPositions) {
            tInit.push_back({tClusterType, tX, TVector3()});
            tRef.mass.push_back(tM);
            tRef.charge.push_back(1);
            tRef.x.push_back(tX);
            tRef.p.push_back(TVector3());
        }
        RunResult tRun = RunCoulomb(tDir, tSteps, tStep, tInit, &tDB);
        const int tSub = 20;
        std::vector<double> tT, tX0Code, tX0Ref, tPx0Code, tPx0Ref;
        for (int n = 0; n < tSteps; ++n) {
            tT.push_back(tRun.trajectory[0][n].T());
            tX0Code.push_back(tRun.trajectory[0][n].X());
            tX0Ref.push_back(tRef.x[0].X());
            tPx0Code.push_back(tRun.fourMomenta[0][n].Px());
            tPx0Ref.push_back(tRef.p[0].X());
            for (int k = 0; k < tSub; ++k) tRef.Step(tStep / tSub);
        }
        double tPMax = 0, tVMax = 0, tErrMax = 0;
        std::vector<double> tIndex, tRelErr;
        for (int i = 0; i < tN; ++i) {
            tPMax = std::max(tPMax, tRef.p[i].Mag());
            tVMax = std::max(tVMax, tRef.p[i].Mag() / std::sqrt(tM * tM + tRef.p[i].Mag2()));
        }
        for (int i = 0; i < tN; ++i) {
            double tErr = (tRun.finalMomenta[i] - tRef.p[i]).Mag() / tPMax;
            tErrMax = std::max(tErrMax, tErr);
            tIndex.push_back(i);
            tRelErr.push_back(tErr);
        }
        snprintf(tBuf, sizeof(tBuf), "J: max final speed %.4f, v^2 = %.1e (size of retardation effects)", tVMax, tVMax * tVMax);
        Info(tBuf);
        Check("J: max |p(code) - p(instantaneous)| / max |p| over the charges", tErrMax, 0, 3 * tVMax * tVMax, false);
        AddGraph(tDir, "final_p_diff", "|p(code) - p(instantaneous Coulomb)| / max|p| at the end;charge;relative difference", tIndex, tRelErr);
        AddCodeExact(tDir, "p0_x", "x of charge 0, exact = instantaneous Coulomb N-body", "x (fm)", tT, tX0Code, tX0Ref);
        AddCodeExact(tDir, "p0_px", "p_{x} of charge 0, exact = instantaneous Coulomb N-body", "p_{x} (GeV)", tT, tPx0Code, tPx0Ref);
    }

    // ----- L. A particle with a non-finite momentum must not affect the others
    // The Rutherford scattering of B with an additional charged particle with NaN momentum: it must be left out of the
    // propagation, so the proton must come out exactly as in B, and the bad particle must be left unchanged.
    Info("L. Rutherford scattering of B with an additional charged particle with NaN momentum");
    {
        const double tNaN = std::numeric_limits<double>::quiet_NaN();
        std::vector<PartInit> tInit = {{tProton, TVector3(tRuthB, 0, -tRuthZ0), TVector3(0, 0, tRuthP0)},
                                       {tHeavy, TVector3(0, 0, 0), TVector3(0, 0, 0)},
                                       {tProton, TVector3(-3, 1, 2), TVector3(tNaN, 0.1, 0.2)}};
        RunResult tRun = RunCoulomb("L_nan_particle", tRuthSteps, tRuthStep, tInit, &tDB);
        Check("L: proton |p - p(B)| / |p|", (tRun.finalMomenta[0] - tRunB.finalMomenta[0]).Mag() / tRuthP0, 0, 0, false);
        Check("L: bad particle left unchanged (p_y, p_z) [GeV]",
              std::fabs(tRun.finalMomenta[2].Y() - 0.1) + std::fabs(tRun.finalMomenta[2].Z() - 0.2), 0, 1e-7, false);
        Check("L: positions and times restored (1 = yes)", tRun.restored, 1, 0, false);
    }

    // ----- M. Close encounters, integrated with sub-steps (CoulombSubsteps = 1)
    // A proton passes a unit charge at impact parameters down to 0.01 fm, where the encounter lasts a fraction of the
    // time step of 0.1 fm/c. Reference: the exact relativistic Coulomb orbit in the rest frame of the source.
    // Without sub-steps (shown as INFO) the time step cannot follow the encounter.
    Info("M. close encounters of a proton with a unit charge, with sub-steps: exact relativistic orbits");
    {
        const double tM = tProton->GetMass(), tBeta = 0.6;
        const double tP0 = tM * tBeta / std::sqrt(1 - tBeta * tBeta);
        const double tZ0 = 5.0;
        // M1: fixed source; proton from (b, 0, -z0) along z, until it is back at the distance r0 (as in B)
        auto tFixedSource = [&](const char *aDir, double aB, ParticleType *aSource, double aStep, bool aSubsteps) {
            std::vector<PartInit> tInit = {{tProton, TVector3(aB, 0, -tZ0), TVector3(0, 0, tP0)},
                                           {aSource, TVector3(0, 0, 0), TVector3(0, 0, 0)}};
            int tSteps = (int)std::lround(2 * tZ0 / tBeta / aStep);
            sRunConfig = aSubsteps ? &sSubstepConfig : &sChecksConfig;
            RunResult tRun = RunCoulomb(aDir, tSteps, aStep, tInit, &tDB);
            sRunConfig = &sChecksConfig;
            return tRun;
        };
        struct FixedCase { const char *dir; double b; ParticleType *source; double kappa; const char *what; };
        for (const FixedCase &tCase : {FixedCase{"M1_fixed_b0.3", 0.3, tProbe, kAlphaHbarC, "repulsive, b = 0.3 fm"},
                                       FixedCase{"M1_fixed_b0.01", 0.01, tProbe, kAlphaHbarC, "repulsive, b = 0.01 fm"},
                                       FixedCase{"M1_fixed_b0.01_attractive", 0.01, tProbeNegative, -kAlphaHbarC, "attractive, b = 0.01 fm"}}) {
            const double tR0 = std::hypot(tCase.b, tZ0);
            const ExactOrbit tOrbit(tM, tCase.kappa, tR0, tP0, tCase.b);
            const double tTheta = tOrbit.DeflectionAtDistance(tR0);
            RunResult tOff = tFixedSource((std::string(tCase.dir) + "_nosubsteps").c_str(), tCase.b, tCase.source, 0.1, false);
            RunResult tRun = tFixedSource(tCase.dir, tCase.b, tCase.source, 0.1, true);
            double tThetaOff = TVector3(0, 0, 1).Angle(tOff.finalMomenta[0]);
            double tThetaCode = TVector3(0, 0, 1).Angle(tRun.finalMomenta[0]);
            snprintf(tBuf, sizeof(tBuf), "M1 %s: exact deflection %.5f rad; without sub-steps %.5f rad, |p| = %.5f GeV (initial %.5f)",
                     tCase.what, tTheta, tThetaOff, tOff.finalMomenta[0].Mag(), tP0);
            Info(tBuf);
            std::string tWhat = std::string("M1: ") + tCase.what + ", ";
            Check((tWhat + "deflection vs exact orbit [rad]").c_str(), tThetaCode, tTheta, 5e-4);
            Check((tWhat + "|p| vs initial |p| [GeV]").c_str(), tRun.finalMomenta[0].Mag(), tP0, 3e-5);
            Check((tWhat + "positions and times restored (1 = yes)").c_str(), tRun.restored, 1, 0, false);
            // Distance from the source vs the exact orbit at the same polar angle (as in B)
            const double tPsiPeri = std::atan2(tCase.b, -tZ0) - tOrbit.Phi(tR0);
            std::vector<double> tT, tRCode, tRExact;
            for (const TLorentzVector &tX : tRun.trajectory[0]) {
                tT.push_back(tX.T());
                tRCode.push_back(std::hypot(tX.X(), tX.Z()));
                tRExact.push_back(tOrbit.R(tPsiPeri - std::atan2(tX.X(), tX.Z())));
            }
            AddCodeExact(tCase.dir, "r", "distance from the source, exact = orbit at the same polar angle", "r (fm)", tT, tRCode, tRExact);
        }

        // M2: like a spectator, a source moving with beta_s = 0.63 passes a proton at rest at b = 0.02 fm.
        // In the rest frame of the source, the proton comes in with -beta_s: exact orbit there, boosted to the lab.
        const double tBetaS = 0.63, tGammaS = 1 / std::sqrt(1 - tBetaS * tBetaS), tBS = 0.02, tZS = 10.0;
        const double tMSource = tProbe->GetMass();
        const double tERest = tGammaS * tM, tPRest = tGammaS * tBetaS * tM; // proton in the rest frame of the source
        const double tThetaS = AsymptoticDeflection(tM, kAlphaHbarC, tERest, tBS);
        // Rest frame: from -z, deflected to +x; lab: boost with +beta_s
        const TVector3 tPExact(tPRest * std::sin(tThetaS), 0, tGammaS * (-tPRest * std::cos(tThetaS) + tBetaS * tERest));
        auto tMovingSource = [&](const char *aDir, double aStep, bool aSubsteps) {
            std::vector<PartInit> tInit = {{tProton, TVector3(tBS, 0, 0), TVector3(0, 0, 0)},
                                           {tProbe, TVector3(0, 0, -tZS), TVector3(0, 0, tMSource * tGammaS * tBetaS)}};
            // until the source is at +z0: in the rest frame the proton is then as far away as at the start
            int tSteps = (int)std::lround(2 * tZS / tBetaS / aStep);
            sRunConfig = aSubsteps ? &sSubstepConfig : &sChecksConfig;
            RunResult tRun = RunCoulomb(aDir, tSteps, aStep, tInit, &tDB);
            sRunConfig = &sChecksConfig;
            return tRun;
        };
        // Deflection in the rest frame of the source, from the final momentum in the lab
        auto tRestFrameDeflection = [&](const TVector3 &aP) {
            double tE = std::sqrt(tM * tM + aP.Mag2());
            TVector3 tPRestFinal(aP.X(), aP.Y(), tGammaS * (aP.Z() - tBetaS * tE));
            return TVector3(0, 0, -1).Angle(tPRestFinal);
        };
        RunResult tOff = tMovingSource("M2_moving_nosubsteps", 0.1, false);
        RunResult tRun = tMovingSource("M2_moving", 0.1, true);
        snprintf(tBuf, sizeof(tBuf), "M2: exact p = (%.5f, 0, %.5f) GeV; without sub-steps (%.5f, %.5f, %.5f) GeV",
                 tPExact.X(), tPExact.Z(), tOff.finalMomenta[0].X(), tOff.finalMomenta[0].Y(), tOff.finalMomenta[0].Z());
        Info(tBuf);
        Check("M2: rest-frame deflection vs exact orbit [rad]", tRestFrameDeflection(tRun.finalMomenta[0]), tThetaS, 3e-4);
        Check("M2: lab |p - p(exact)| / |p(exact)|", (tRun.finalMomenta[0] - tPExact).Mag() / tPExact.Mag(), 0, 3e-4, false);
        Check("M2: positions and times restored (1 = yes)", tRun.restored, 1, 0, false);
        // M3: the result must not depend on the time step
        for (double tStep : {0.05, 0.2}) {
            snprintf(tBuf, sizeof(tBuf), "M2_moving_step%.2f", tStep);
            RunResult tRunStep = tMovingSource(tBuf, tStep, true);
            snprintf(tBuf, sizeof(tBuf), "M3: step %.2f fm/c, |p - p(step 0.1)| / |p|", tStep);
            Check(tBuf, (tRunStep.finalMomenta[0] - tRun.finalMomenta[0]).Mag() / tRun.finalMomenta[0].Mag(), 0, 2e-4, false);
        }
        std::vector<double> tT, tXCode, tXOff;
        for (size_t n = 0; n < tRun.trajectory[0].size(); ++n) {
            tT.push_back(tRun.trajectory[0][n].T());
            tXCode.push_back(tRun.fourMomenta[0][n].Px());
            tXOff.push_back(tOff.fourMomenta[0][n].Px());
        }
        AddCodeExact("M2_moving", "px_substeps_vs_none", "p_{x} of the proton with sub-steps (code) and without (exact)", "p_{x} (GeV)",
                     tT, tXCode, tXOff);
    }

    // ----- K. Settings of the afterburner from the configuration
    // The defaults must be the values of the former constants, events.ini and events_template must contain all keys
    // with these values, all keys must be read, and missing keys must keep the defaults.
    Info("K. settings of the afterburner from the configuration: defaults, events.ini, events_template, all keys, no keys");
    {
        const std::vector<std::string> tKeys = {"CoulombUseSpectators", "CoulombExtrapolateSpectators", "CoulombExtrapolateFireball",
                                                "CoulombFindProblem", "CoulombLargeTerm", "CoulombDebugStep", "CoulombDebugEid",
                                                "CoulombDebugOtherEid", "CoulombSubsteps", "CoulombSubstepAccuracy",
                                                "CoulombSubstepMax", "CoulombSubstepDistance"};
        // The former constants, and the defaults of the later settings (sub-steps off)
        const std::vector<double> tFormerConstants = {1, 1, 1, 1, 1.0, -2, -2, -2, 0, 0.05, 1024, 1.0};
        auto tSettings = [](const CoulombAfterburner &aCoulomb) {
            return std::vector<double>{(double)aCoulomb.m_UseSpectators, (double)aCoulomb.m_ExtrapolateSpectators,
                                       (double)aCoulomb.m_ExtrapolateFireball, (double)aCoulomb.m_FindProblem,
                                       aCoulomb.m_LargeTerm.numerical_value_in(CoulombUnits::GeV / CoulombUnits::fm),
                                       (double)aCoulomb.m_DebugStep, (double)aCoulomb.m_DebugEid, (double)aCoulomb.m_DebugOtherEid,
                                       (double)aCoulomb.m_Substeps, aCoulomb.m_SubstepAccuracy, (double)aCoulomb.m_SubstepMax,
                                       aCoulomb.m_SubstepDistance.numerical_value_in(CoulombUnits::fm)};
        };
        // Number of settings different from the expected values, printing them
        auto tMismatches = [&](const char *aWhat, const std::vector<double> &aGot, const std::vector<double> &aExpected) {
            int tN = 0;
            for (size_t i = 0; i < tKeys.size(); ++i) {
                if (aGot[i] != aExpected[i]) {
                    printf("INFO  K: %s: %s = %g, expected %g\n", aWhat, tKeys[i].c_str(), aGot[i], aExpected[i]);
                    ++tN;
                }
            }
            return tN;
        };
        {
            CoulombAfterburner tCoulomb(1, 0.1, 0., nullptr, &tDB);
            Check("K: defaults different from the former constants", tMismatches("defaults", tSettings(tCoulomb), tFormerConstants), 0, 0, false);
        }
        if (tTopDir) {
            for (const char *tFile : {"events.ini", "events_template"}) {
                std::string tPath = std::string(tTopDir) + "/" + tFile;
                Configurator tConfig;
                Parser tParser(tPath.c_str());
                tParser.ReadINI(&tConfig);
                int tMissing = 0;
                for (const std::string &tKey : tKeys) {
                    if (!tConfig.HasParameter(tKey.c_str())) {
                        printf("INFO  K: %s: key %s missing\n", tFile, tKey.c_str());
                        ++tMissing;
                    }
                }
                Check((std::string("K: keys missing in ") + tFile).c_str(), tMissing, 0, 0, false);
                CoulombAfterburner tCoulomb(1, 0.1, 0., nullptr, &tDB);
                tCoulomb.ReadParameters(&tConfig);
                Check((std::string("K: settings from ") + tFile + " different from the former constants").c_str(),
                      tMismatches(tFile, tSettings(tCoulomb), tFormerConstants), 0, 0, false);
            }
        } else {
            Info("K: no repository directory given, events.ini and events_template not checked");
        }
        {
            const std::vector<std::string> tValues = {"0", "0", "0", "0", "2.5", "7", "11", "13", "1", "0.02", "64", "0.5"};
            Configurator tConfig;
            for (size_t i = 0; i < tKeys.size(); ++i) {
                Parameter tPar{tKeys[i].c_str(), tValues[i].c_str()};
                tConfig.AddParameter(&tPar);
            }
            CoulombAfterburner tCoulomb(1, 0.1, 0., nullptr, &tDB);
            tCoulomb.ReadParameters(&tConfig);
            Check("K: settings not taken from a configuration with all keys",
                  tMismatches("all keys", tSettings(tCoulomb), {0, 0, 0, 0, 2.5, 7, 11, 13, 1, 0.02, 64, 0.5}), 0, 0, false);
        }
        {
            Configurator tConfig;
            CoulombAfterburner tCoulomb(1, 0.1, 0., nullptr, &tDB);
            tCoulomb.ReadParameters(&tConfig);
            Check("K: defaults not kept with an empty configuration", tMismatches("no keys", tSettings(tCoulomb), tFormerConstants), 0, 0, false);
        }
    }

    WriteGraphs(tOutputFile);
    snprintf(tBuf, sizeof(tBuf), "Graphs written to %s", tOutputFile);
    Info(tBuf);

    printf("%s: %d failure(s)\n", sFailures ? "FAILED" : "ALL PASSED", sFailures);
    return sFailures != 0;
}
