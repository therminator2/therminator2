// Full-evolution checks of CoulombAfterburner::Apply on simple systems with known solutions.
// Build and run with checks/coulomb/run.sh (spectators are switched off in the copy that is tested).
//
// Units as in THERMINATOR: inside the afterburner lengths and times are in GeV^-1, momenta in GeV.
// Here positions and times are given in fm (fm/c) and converted with kHbarC.
//
// Besides the pass/fail checks, trajectories of the code and of the exact solutions are written as TGraphs
// to a ROOT file (first argument, default coulomb_checks.root), one directory per check:
//   <label>_code, <label>_exact, <label>_diff (= code - exact) as functions of t (fm/c),
//   and for B the paths in the scattering plane (path_code, path_rk4, path_exact).

#include "CoulombAfterburner.h"
#include "Event.h"
#include "Particle.h"
#include "ParticleDB.h"
#include "ParticleType.h"
#include "THGlobal.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMath.h"
#include "TVector3.h"

#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <utility>
#include <vector>

TString sModelINI;   // globals expected by the library code
TString sTimeStamp;

static const double kAlphaHbarC = 1.439964e-3; // GeV fm, alpha * hbar * c (CODATA)

static int sFailures = 0;

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

// ---------------------------------------------------------------------------------------------
// Running the afterburner

struct PartInit {
    ParticleType *type;
    TVector3 x; // fm
    TVector3 p; // GeV
};

struct RunResult {
    std::vector<TVector3> finalMomenta;                    // GeV
    std::vector<std::vector<TLorentzVector> > trajectory;  // (x, y, z, t) in fm, fm/c, at the time steps
    std::vector<std::vector<TLorentzVector> > fourMomenta; // (px, py, pz, E) in GeV, at the time steps
    bool restored;                                         // positions and times restored after Apply()
};

// Runs CoulombAfterburner::Apply on one event made of the given particles, all created at t = 0
static RunResult RunCoulomb(const char *aDir, int aSteps, double aStepFm, const std::vector<PartInit> &aInit, ParticleDB *aDB)
{
    mkdir(aDir, 0755);
    if (chdir(aDir) != 0) { perror(aDir); exit(2); }
    RunResult tResult;
    {
        CoulombAfterburner tCoulomb(aSteps, aStepFm, 0., nullptr, aDB);
        tCoulomb.SetKeepTrajectories(true);
        tCoulomb.Init();
        Event tEvent; // resets the eid counter
        for (const PartInit &tInit : aInit) {
            Particle tParticle(tInit.type);
            double tE = std::sqrt(tParticle.mass * tParticle.mass + tInit.p.Mag2());
            tParticle.SetParticlePX(tE, tInit.p.X(), tInit.p.Y(), tInit.p.Z(),
                                    0., tInit.x.X() / kHbarC, tInit.x.Y() / kHbarC, tInit.x.Z() / kHbarC);
            tEvent.GetParticleList()->push_back(tParticle);
        }
        tCoulomb.Apply(&tEvent);
        tResult.restored = true;
        size_t i = 0;
        for (Particle &tParticle : *tEvent.GetParticleList()) {
            tResult.finalMomenta.push_back(tParticle.VecP());
            TVector3 tX0 = (1 / kHbarC) * aInit[i].x;
            if (tParticle.T() != 0. || (tParticle.VecX() - tX0).Mag() > 1e-6 * (1 + tX0.Mag())) tResult.restored = false;
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

// Two equal charges released from rest at distance r0 (non-relativistic): time to reach distance x*r0 is
// t(x) = sqrt(mu r0^3 / (2 kappa)) [sqrt(x(x-1)) + ln(sqrt(x) + sqrt(x-1))].
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
    return std::sqrt(2 * aMu * aKappa / aR0 * (1 - 1 / x)); // momentum of each particle = mu * v_rel
}

// ---------------------------------------------------------------------------------------------

int main(int argc, char **argv)
{
    const char *tOutputFile = argc > 1 ? argv[1] : "coulomb_checks.root";

    ParticleDB tDB;
    ParticleType tProtonType = MakeType("proton", 2212, 0.938272, 1);   // also needed by the constructor
    ParticleType tNeutronType = MakeType("neutron", 2112, 0.939565, 0);
    ParticleType tHeavyType = MakeType("heavyZ79", 9079, 1.0e5, 79);   // practically fixed source
    ParticleType tSlowType = MakeType("slowQ1", 9001, 10.0, 1);        // heavy unit charge, slow motion
    tDB.AddParticleType(&tProtonType);
    tDB.AddParticleType(&tNeutronType);
    tDB.AddParticleType(&tHeavyType);
    tDB.AddParticleType(&tSlowType);
    ParticleType *tProton = tDB.FindByPID(2212); // pointers are valid only after all types are added
    ParticleType *tNeutron = tDB.FindByPID(2112);
    ParticleType *tHeavy = tDB.FindByPID(9079);
    ParticleType *tSlow = tDB.FindByPID(9001);
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

    // ----- B. Relativistic Rutherford scattering on a fixed Z = 79 charge
    Info("B. proton, T = 200 MeV, b = 4 fm, on a fixed Z = 79 charge: exact relativistic orbit");
    {
        const char *tDir = "B_rutherford";
        double tM = tProton->GetMass(); // as used by the code (float)
        double tKappa = 79 * kAlphaHbarC;
        double tB = 4.0, tZ0 = 40.0;
        double tP0 = std::sqrt(std::pow(tM + 0.200, 2) - tM * tM);
        double tR0 = std::hypot(tB, tZ0);
        ExactOrbit tOrbit(tM, tKappa, tR0, tP0, tB);
        // Time at which the proton is back at distance r0, from a fine RK4 integration
        State tS{TVector3(tB, 0, -tZ0), TVector3(0, 0, tP0)};
        double tH = 1e-3, t = 0, tPrevR = tR0;
        bool tPastPeriapsis = false;
        double tTReturn = 0;
        while (true) {
            State tNext = RK4Step(tS, tH, tM, tKappa);
            double tR = tNext.x.Mag();
            if (tR > tPrevR) tPastPeriapsis = true;
            if (tPastPeriapsis && tR >= tR0) { tTReturn = t + tH * (tR0 - tPrevR) / (tR - tPrevR); break; }
            tS = tNext; tPrevR = tR; t += tH;
        }
        double tStep = 0.1;
        int tSteps = (int)std::lround(tTReturn / tStep);
        double tT = tSteps * tStep;
        snprintf(tBuf, sizeof(tBuf), "B: back at r0 = %.2f fm after %.2f fm/c, run %d steps of %.1f fm/c", tR0, tTReturn, tSteps, tStep);
        Info(tBuf);

        std::vector<PartInit> tInit = {{tProton, TVector3(tB, 0, -tZ0), TVector3(0, 0, tP0)},
                                       {tHeavy, TVector3(0, 0, 0), TVector3(0, 0, 0)}};
        RunResult tRun = RunCoulomb(tDir, tSteps, tStep, tInit, &tDB);

        // RK4 reference at the time steps of the code (100 substeps per step) and at the final time
        int tSub = 100;
        State tRef{TVector3(tB, 0, -tZ0), TVector3(0, 0, tP0)};
        std::vector<State> tRefAtSteps;
        for (int n = 0; n < tSteps; ++n) {
            tRefAtSteps.push_back(tRef);
            for (int k = 0; k < tSub; ++k) tRef = RK4Step(tRef, tStep / tSub, tM, tKappa);
        }

        double tThetaExact = tOrbit.DeflectionAtDistance(tR0);
        double tThetaRef = TVector3(0, 0, 1).Angle(tRef.p);
        double tThetaCode = TVector3(0, 0, 1).Angle(tRun.finalMomenta[0]);
        Check("B: reference RK4 deflection vs exact orbit [rad]", tThetaRef, tThetaExact, 1e-4);
        Check("B: code deflection vs exact orbit [rad]", tThetaCode, tThetaExact, 3e-5);
        Check("B: code |p| vs initial |p| (energy conservation) [GeV]", tRun.finalMomenta[0].Mag(), tP0, 3e-5);
        Check("B: code |p - p_RK4| / |p|", (tRun.finalMomenta[0] - tRef.p).Mag() / tP0, 0, 1e-5, false);
        Check("B: positions and times restored (1 = yes)", tRun.restored, 1, 0, false);

        // Graphs. Scattering plane: horizontal z, vertical x; polar angle psi = atan2(x, z).
        // The orbit is traversed clockwise (L_y < 0); the start is at phi = -phi0, so psi = psi_peri - phi.
        double tPhi0 = tOrbit.Phi(tR0);
        double tPsiPeri = std::atan2(tB, -tZ0) - tPhi0;
        std::vector<double> tPathZ, tPathX;
        for (int k = 0; k <= 1000; ++k) {
            double tPhi = -tPhi0 + 2 * tPhi0 * k / 1000;
            double tR = tOrbit.R(tPhi), tPsi = tPsiPeri - tPhi;
            tPathZ.push_back(tR * std::cos(tPsi));
            tPathX.push_back(tR * std::sin(tPsi));
        }
        AddGraph(tDir, "path_exact", "Exact relativistic Coulomb orbit;z (fm);x (fm)", tPathZ, tPathX);
        std::vector<double> tTimes, tCodeZ, tCodeX, tRefZ, tRefX, tCodePx, tRefPx, tCodePz, tRefPz, tDrOrbit;
        const std::vector<TLorentzVector> &tTraj = tRun.trajectory[0];
        const std::vector<TLorentzVector> &tMom = tRun.fourMomenta[0];
        for (size_t n = 0; n < tTraj.size(); ++n) {
            tTimes.push_back(tTraj[n].T());
            tCodeZ.push_back(tTraj[n].Z());
            tCodeX.push_back(tTraj[n].X());
            tRefZ.push_back(tRefAtSteps[n].x.Z());
            tRefX.push_back(tRefAtSteps[n].x.X());
            tCodePx.push_back(tMom[n].Px());
            tRefPx.push_back(tRefAtSteps[n].p.X());
            tCodePz.push_back(tMom[n].Pz());
            tRefPz.push_back(tRefAtSteps[n].p.Z());
            // distance from the exact orbit at the same polar angle (independent of timing)
            double tPsi = std::atan2(tTraj[n].X(), tTraj[n].Z());
            tDrOrbit.push_back(std::hypot(tTraj[n].X(), tTraj[n].Z()) - tOrbit.R(tPsiPeri - tPsi));
        }
        AddGraph(tDir, "path_code", "Trajectory (code);z (fm);x (fm)", tCodeZ, tCodeX);
        AddGraph(tDir, "path_rk4", "Trajectory (RK4 of the exact equation of motion);z (fm);x (fm)", tRefZ, tRefX);
        AddGraph(tDir, "dr_orbit", "r(code) - r(exact orbit) at the same polar angle;t (fm/c);#Deltar (fm)", tTimes, tDrOrbit);
        AddCodeExact(tDir, "x", "proton x, exact = RK4", "x (fm)", tTimes, tCodeX, tRefX);
        AddCodeExact(tDir, "z", "proton z, exact = RK4", "z (fm)", tTimes, tCodeZ, tRefZ);
        AddCodeExact(tDir, "px", "proton p_{x}, exact = RK4", "p_{x} (GeV)", tTimes, tCodePx, tRefPx);
        AddCodeExact(tDir, "pz", "proton p_{z}, exact = RK4", "p_{z} (GeV)", tTimes, tCodePz, tRefPz);
        double tMaxDr = 0;
        for (double tDr : tDrOrbit) tMaxDr = std::max(tMaxDr, std::fabs(tDr));
        Check("B: max |r(code) - r(exact orbit)| along the path [fm]", tMaxDr, 0, 1e-3, false);
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
        // p_z(lab) = gamma beta E'(rest): the kinetic energy gained in the rest frame shows up in p_z.
        // Loose tolerance: the change per step (~1e-6 GeV) is comparable to the Float_t spacing at p_z = 20.6 GeV.
        double tDPzExpected = tGamma * tBeta * (std::sqrt(tMSlow * tMSlow + tPRest[1].Mag2()) - tMSlow);
        Check("D: change of p_z(lab) vs gamma beta x kinetic energy(rest) [GeV]", tP[1].Z() - tPz, tDPzExpected, 2e-2);
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

    WriteGraphs(tOutputFile);
    snprintf(tBuf, sizeof(tBuf), "Graphs written to %s", tOutputFile);
    Info(tBuf);

    printf("%s: %d failure(s)\n", sFailures ? "FAILED" : "ALL PASSED", sFailures);
    return sFailures != 0;
}
