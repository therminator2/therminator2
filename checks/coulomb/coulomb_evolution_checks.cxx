// Full-evolution checks of CoulombAfterburner::Apply on simple systems with known solutions.
// Build and run with checks/coulomb/run.sh (spectators are switched off in the copy that is tested).
//
// Units as in THERMINATOR: inside the afterburner lengths and times are in GeV^-1, momenta in GeV.
// Here positions and times are given in fm (fm/c) and converted with kHbarC.

#include "CoulombAfterburner.h"
#include "Event.h"
#include "Particle.h"
#include "ParticleDB.h"
#include "ParticleType.h"
#include "THGlobal.h"
#include "TMath.h"
#include "TVector3.h"

#include <cmath>
#include <cstdio>
#include <functional>
#include <sys/stat.h>
#include <unistd.h>
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

struct PartInit {
    ParticleType *type;
    TVector3 x; // fm
    TVector3 p; // GeV
};

// Runs CoulombAfterburner::Apply on one event made of the given particles, all created at t = 0,
// and returns their final momenta. aRestored tells if positions and times were restored afterwards.
static std::vector<TVector3> RunCoulomb(const char *aDir, int aSteps, double aStepFm,
                                        const std::vector<PartInit> &aInit, ParticleDB *aDB, bool &aRestored)
{
    mkdir(aDir, 0755);
    if (chdir(aDir) != 0) { perror(aDir); exit(2); }
    std::vector<TVector3> tMomenta;
    {
        CoulombAfterburner tCoulomb(aSteps, aStepFm, 0., nullptr, aDB);
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
        aRestored = true;
        size_t i = 0;
        for (Particle &tParticle : *tEvent.GetParticleList()) {
            tMomenta.push_back(tParticle.VecP());
            TVector3 tX0 = (1 / kHbarC) * aInit[i].x;
            if (tParticle.T() != 0. || (tParticle.VecX() - tX0).Mag() > 1e-6 * (1 + tX0.Mag())) aRestored = false;
            ++i;
        }
        tCoulomb.Finalize();
    }
    if (chdir("..") != 0) { perror(".."); exit(2); }
    return tMomenta;
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
// so u = A cos(g phi) - C. For a start at distance r0 with momentum p0 (impact parameter b), returns the angle
// between the momenta at the incoming and the outgoing point at the same distance r0 (mirror points of the orbit).
static double ExactDeflectionAtDistance(double aMass, double aKappa, double aR0, double aP0, double aB)
{
    double tE = std::sqrt(aMass * aMass + aP0 * aP0) + aKappa / aR0; // total energy
    double tL = aB * aP0;                                              // angular momentum
    double tG2 = 1 - aKappa * aKappa / (tL * tL);
    double tG = std::sqrt(tG2);
    double tC = tE * aKappa / (tL * tL * tG2);
    double tA = std::sqrt(tC * tC + (tE * tE - aMass * aMass) / (tL * tL * tG2));
    double tU0 = 1 / aR0;
    double tPhi0 = std::acos((tU0 + tC) / tA) / tG;              // outgoing point, periapsis at phi = 0
    double tDuDphi = -tA * tG * std::sin(tG * tPhi0);
    // p_r = -L u', p_phi = L u
    TVector3 tOut = (-tDuDphi) * TVector3(std::cos(tPhi0), std::sin(tPhi0), 0)
                    + tU0 * TVector3(-std::sin(tPhi0), std::cos(tPhi0), 0);
    TVector3 tIn(-tOut.X(), tOut.Y(), 0); // mirror image, time reversed
    return tIn.Angle(tOut);
}

// Two equal charges released from rest at distance r0 (non-relativistic): time to reach distance x*r0 is
// t(x) = sqrt(mu r0^3 / (2 kappa)) [sqrt(x(x-1)) + ln(sqrt(x) + sqrt(x-1))].
static double NonRelRepulsionMomentum(double aMu, double aKappa, double aR0, double aT)
{
    auto tOfX = [&](double x) {
        return std::sqrt(aMu * aR0 * aR0 * aR0 / (2 * aKappa)) * (std::sqrt(x * (x - 1)) + std::log(std::sqrt(x) + std::sqrt(x - 1)));
    };
    double lo = 1, hi = 2;
    while (tOfX(hi) < aT) hi *= 2;
    for (int i = 0; i < 200; ++i) { double m = 0.5 * (lo + hi); (tOfX(m) < aT ? lo : hi) = m; }
    double x = 0.5 * (lo + hi);
    return std::sqrt(2 * aMu * aKappa / aR0 * (1 - 1 / x)); // momentum of each particle = mu * v_rel
}

// ---------------------------------------------------------------------------------------------

int main()
{
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
    bool tRestored;

    // ----- A. No interaction: one charged and one neutral particle
    Info("A. one charged + one neutral particle: no force");
    {
        std::vector<PartInit> tInit = {{tProton, TVector3(1, 2, 3), TVector3(0.1, 0.2, 0.3)},
                                       {tNeutron, TVector3(-1, 0, 0), TVector3(0, 0, 0.5)}};
        std::vector<TVector3> tP = RunCoulomb("A_free", 50, 0.1, tInit, &tDB, tRestored);
        for (int i = 0; i < 2; ++i) {
            snprintf(tBuf, sizeof(tBuf), "A: |p_final - p_initial| of particle %d [GeV]", i);
            Check(tBuf, (tP[i] - tInit[i].p).Mag(), 0, 1e-6, false);
        }
        Check("A: positions and times restored (1 = yes)", tRestored, 1, 0, false);
    }

    // ----- B. Relativistic Rutherford scattering on a fixed Z = 79 charge
    Info("B. proton, T = 200 MeV, b = 4 fm, on a fixed Z = 79 charge: exact relativistic orbit");
    {
        double tM = tProton->GetMass(); // as used by the code (float)
        double tKappa = 79 * kAlphaHbarC;
        double tB = 4.0, tZ0 = 40.0;
        double tP0 = std::sqrt(std::pow(tM + 0.200, 2) - tM * tM);
        double tR0 = std::hypot(tB, tZ0);
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
        // RK4 reference state at exactly the final time of the code
        State tRef{TVector3(tB, 0, -tZ0), TVector3(0, 0, tP0)};
        int tRefSteps = (int)std::lround(tT / tH);
        for (int i = 0; i < tRefSteps; ++i) tRef = RK4Step(tRef, tT / tRefSteps, tM, tKappa);
        snprintf(tBuf, sizeof(tBuf), "B: back at r0 = %.2f fm after %.2f fm/c, run %d steps of %.1f fm/c", tR0, tTReturn, tSteps, tStep);
        Info(tBuf);

        std::vector<PartInit> tInit = {{tProton, TVector3(tB, 0, -tZ0), TVector3(0, 0, tP0)},
                                       {tHeavy, TVector3(0, 0, 0), TVector3(0, 0, 0)}};
        std::vector<TVector3> tP = RunCoulomb("B_rutherford", tSteps, tStep, tInit, &tDB, tRestored);
        double tThetaExact = ExactDeflectionAtDistance(tM, tKappa, tR0, tP0, tB);
        double tThetaRef = TVector3(0, 0, 1).Angle(tRef.p);
        double tThetaCode = TVector3(0, 0, 1).Angle(tP[0]);
        Check("B: reference RK4 deflection vs exact orbit [rad]", tThetaRef, tThetaExact, 1e-4);
        Check("B: code deflection vs exact orbit [rad]", tThetaCode, tThetaExact, 3e-5);
        Check("B: code |p| vs initial |p| (energy conservation) [GeV]", tP[0].Mag(), tP0, 3e-5);
        Check("B: code |p - p_RK4| / |p|", (tP[0] - tRef.p).Mag() / tP0, 0, 1e-5, false);
        Check("B: positions and times restored (1 = yes)", tRestored, 1, 0, false);
    }

    // ----- C. Two equal charges released from rest (slow, heavy)
    Info("C. two charges q = 1, m = 10 GeV, released from rest 1 fm apart, 60 fm/c: non-relativistic solution");
    std::vector<TVector3> tPRest;
    double tMSlow = tSlow->GetMass();
    double tKappa1 = kAlphaHbarC;
    int tStepsC = 600;
    double tStepC = 0.1;
    {
        std::vector<PartInit> tInit = {{tSlow, TVector3(-0.5, 0, 0), TVector3()},
                                       {tSlow, TVector3(0.5, 0, 0), TVector3()}};
        tPRest = RunCoulomb("C_rest", tStepsC, tStepC, tInit, &tDB, tRestored);
        double tPExact = NonRelRepulsionMomentum(tMSlow / 2, tKappa1, 1.0, tStepsC * tStepC);
        snprintf(tBuf, sizeof(tBuf), "C: final speed of each particle v = %.4f (relativistic corrections ~ v^2 = %.1e)",
                 tPExact / tMSlow, std::pow(tPExact / tMSlow, 2));
        Info(tBuf);
        Check("C: p_x of particle 2 vs non-relativistic solution [GeV]", tPRest[1].X(), tPExact, 5e-4);
        Check("C: p_x of particle 1 = -p_x of particle 2 [GeV]", tPRest[0].X(), -tPRest[1].X(), 1e-5);
        Check("C: transverse momenta |p_y| + |p_z| [GeV]",
              std::fabs(tPRest[0].Y()) + std::fabs(tPRest[0].Z()) + std::fabs(tPRest[1].Y()) + std::fabs(tPRest[1].Z()), 0, 1e-9, false);
    }

    // ----- D. The same pair moving with beta = 0.9 along z (Lorentz covariance)
    Info("D. pair of C boosted to beta = 0.9 along z, lab step = gamma x 0.1 fm/c: p_T(lab, t) = p_T(rest, t/gamma)");
    {
        double tBeta = 0.9, tGamma = 1 / std::sqrt(1 - tBeta * tBeta);
        double tPz = tGamma * tMSlow * tBeta;
        std::vector<PartInit> tInit = {{tSlow, TVector3(-0.5, 0, 0), TVector3(0, 0, tPz)},
                                       {tSlow, TVector3(0.5, 0, 0), TVector3(0, 0, tPz)}};
        std::vector<TVector3> tP = RunCoulomb("D_boosted", tStepsC, tGamma * tStepC, tInit, &tDB, tRestored);
        Check("D: p_x(lab) vs p_x(rest frame run C) [GeV]", tP[1].X(), tPRest[1].X(), 1e-5);
        Check("D: p_x(lab) vs non-relativistic solution in rest frame [GeV]", tP[1].X(),
              NonRelRepulsionMomentum(tMSlow / 2, tKappa1, 1.0, tStepsC * tStepC), 5e-4);
        // p_z(lab) = gamma beta E'(rest): the kinetic energy gained in the rest frame shows up in p_z.
        // Loose tolerance: the change per step (~1e-6 GeV) is comparable to the Float_t spacing at p_z = 20.6 GeV.
        double tDPzExpected = tGamma * tBeta * (std::sqrt(tMSlow * tMSlow + tPRest[1].Mag2()) - tMSlow);
        Check("D: change of p_z(lab) vs gamma beta x kinetic energy(rest) [GeV]", tP[1].Z() - tPz, tDPzExpected, 2e-2);
        Check("D: p_x of particle 1 = -p_x of particle 2 [GeV]", tP[0].X(), -tP[1].X(), 1e-5);
    }

    printf("%s: %d failure(s)\n", sFailures ? "FAILED" : "ALL PASSED", sFailures);
    return sFailures != 0;
}
