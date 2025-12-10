#include "CoulombAfterburner.h"
#include "Configurator.h"
#include "Parser.h"
#include "TVector3.h"
#include "TMath.h"
#include <limits>
#include <iostream>
#include "THGlobal.h"

using namespace std;
extern TString sModelINI;

static const bool USE_SPECTATORS = true;
static const bool EXTRAPOLATE_SPECTATORS = true;
static const bool EXTRAPOLATE_FIREBALL = false;

double CoulombAfterburner::Interval(TVector3 tX1, double tT1, TVector3 tX2, double tT2, bool debug) {
    if (debug) {
        cout << tT1 << " " << tT2 << " " << TMath::Power(tT1 - tT2, 2) << " " << (tX1 - tX2).Mag2() << " ";
        (tX1 - tX2).Print();
    }
    return TMath::Power(tT1 - tT2, 2) - (tX1 - tX2).Mag2();
}
double CoulombAfterburner::Interval(TVector3 tX1, double tT1, TLorentzVector tV, bool debug) {
    return this->Interval(tX1, tT1, tV.Vect(), tV.T(), debug);
}
double CoulombAfterburner::Interval(TLorentzVector tV1, TLorentzVector tV2, bool debug) {
    return this->Interval(tV1.Vect(), tV1.T(), tV2.Vect(), tV2.T(), debug);
}

TVector3 CoulombAfterburner::CalculateTotalForce(
    std::list<Particle>::iterator &tPartIter, std::list<Particle> *tParticles, int **tBestIs, double tCurrentTime,
    std::vector<TLorentzVector> *tCoordinates, std::vector<TVector3> *tVelocities, std::vector<TVector3> *tAccelerations)
{
    Int_t tQ = tPartIter->GetParticleType()->GetCharge();
    // TVector3 tForceTotal;
    TVector3 tTermCoulTotal;
    TVector3 tTermVeloTotal;
    TVector3 tTermAcceTotal;
    TVector3 tTermMagnTotal;
    TVector3 tPartPos(tPartIter->x, tPartIter->y, tPartIter->z); // 1/GeV

    TVector3 tPartMom(tPartIter->px, tPartIter->py, tPartIter->pz);
    double tP2 = tPartMom.Mag2();
    double tM2 = TMath::Power(tPartIter->GetParticleType()->GetMass(), 2);
    TVector3 tV(0, 0, 0);
    if (tP2 > 0 && TMath::Finite(tP2))
    {
        tV = 1. / TMath::Sqrt(tM2 / tP2 + 1) * tPartMom.Unit();
    }
    std::list<Particle>::iterator tPartOtherIter = tParticles->begin();
    for (; tPartOtherIter != tParticles->end(); ++tPartOtherIter)
    {

        if (tPartOtherIter->eid == tPartIter->eid)
            continue;

        Int_t tOtherQ = tPartOtherIter->GetParticleType()->GetCharge();
        if (tOtherQ == 0)
            continue;

        double tBestInterval;
        unsigned long tLastBestI = tBestIs[tPartIter->eid][tPartOtherIter->eid] > -1 ? tBestIs[tPartIter->eid][tPartOtherIter->eid] : 0;
        unsigned long tBestI = NearestInterval(tPartPos, tCurrentTime, tCoordinates[tPartOtherIter->eid], tLastBestI, tBestInterval, false);
        tBestIs[tPartIter->eid][tPartOtherIter->eid] = tBestI;

        // The best interval is at the end of trajectory of other particle and it's space-like:
        // the particles are not causally connected
        if ((tBestI + 1) == tCoordinates[tPartOtherIter->eid].size() && tBestInterval < 0)
        {
            continue;
        }

        if (tBestI >= 0 && tBestI < tCoordinates[tPartOtherIter->eid].size())
        {
            if (false && (tBestI - 2) >= 0 && (tBestI + 2) < tCoordinates[tPartOtherIter->eid].size())
            {
                double tIntervalDown2 = this->Interval(tPartPos, tCurrentTime, tCoordinates[tPartOtherIter->eid][tBestI - 2]);
                double tIntervalDown  = this->Interval(tPartPos, tCurrentTime, tCoordinates[tPartOtherIter->eid][tBestI - 1]);
                double tInterval      = this->Interval(tPartPos, tCurrentTime, tCoordinates[tPartOtherIter->eid][tBestI]);
                double tIntervalUp    = this->Interval(tPartPos, tCurrentTime, tCoordinates[tPartOtherIter->eid][tBestI + 1]);
                double tIntervalUp2   = this->Interval(tPartPos, tCurrentTime, tCoordinates[tPartOtherIter->eid][tBestI + 2]);
/*
                double tInterval = TMath::Power(tCurrentTime - tCoordinates[tPartOtherIter->eid][tBestI].T(), 2) - (tPartPos - tCoordinates[tPartOtherIter->eid][tBestI].Vect()).Mag2();
                double tIntervalDown2 = TMath::Power(tCurrentTime - tCoordinates[tPartOtherIter->eid][tBestI - 2].T(), 2) - (tPartPos - tCoordinates[tPartOtherIter->eid][tBestI - 2].Vect()).Mag2();
                double tIntervalDown = TMath::Power(tCurrentTime - tCoordinates[tPartOtherIter->eid][tBestI - 1].T(), 2) - (tPartPos - tCoordinates[tPartOtherIter->eid][tBestI - 1].Vect()).Mag2();
                double tIntervalUp = TMath::Power(tCurrentTime - tCoordinates[tPartOtherIter->eid][tBestI + 1].T(), 2) - (tPartPos - tCoordinates[tPartOtherIter->eid][tBestI + 1].Vect()).Mag2();
                double tIntervalUp2 = TMath::Power(tCurrentTime - tCoordinates[tPartOtherIter->eid][tBestI + 2].T(), 2) - (tPartPos - tCoordinates[tPartOtherIter->eid][tBestI + 2].Vect()).Mag2();
*/
                cout << " tBestI = " << tBestI << " " << tBestInterval << " "
                     << "Intervals " << tIntervalUp2 << " " << tIntervalUp << " " << tInterval << " " << tIntervalDown << " " << tIntervalDown2 << " " << endl
                     << tCurrentTime << " " << tCoordinates[tPartOtherIter->eid][tBestI].T() << " "
                     << (tPartPos - tCoordinates[tPartOtherIter->eid][tBestI].Vect()).Mag() << endl;

                cout << "X, V, A " << tCoordinates[tPartOtherIter->eid][tBestI].Vect().Mag() << " "
                     << tVelocities[tPartOtherIter->eid][tBestI].Mag() << "  " << endl;
//                     << tAccelerations[tPartOtherIter->eid][tBestI].Mag() << endl;
            }
            TVector3 tPartOtherPos = tCoordinates[tPartOtherIter->eid][tBestI].Vect();
            TVector3 tPartOtherVel = tVelocities[tPartOtherIter->eid][tBestI];
            TVector3 tPartOtherAcc = tAccelerations[tPartOtherIter->eid][tBestI];

            // Initial idea: True, when the null interval is NOT withing the trajectory of the particle from its creation
            // to its decay or to the end of simulation
            // Seems to not work, because at late times intervals jump from one point to another by more than the time step size,
            // therefore skipping, maybe will be fully removed
            /// if (TMath::Abs(tBestInterval) >= m_StepSize) {
            // Particle is primary AND null interval would be before freeze-out
            if (tPartOtherIter->pid == tPartOtherIter->fatherpid && tBestI == 0)
            {
                bool tOtherIsSpectator = TString(tPartOtherIter->GetParticleType()->GetName()).Contains("Spectator");
                if (tOtherIsSpectator && !EXTRAPOLATE_SPECTATORS) {
                    continue;
                } else if (!tOtherIsSpectator && !EXTRAPOLATE_FIREBALL) {
                    continue;
                }
                double tZeroTime = tCoordinates[tPartOtherIter->eid][0].T();
                TVector3 tX0 = tCoordinates[tPartOtherIter->eid][0].Vect();
                TVector3 tIniV = tVelocities[tPartOtherIter->eid][0];
                // double NearestPreFreezeoutIntervalTime(TVector3 &tX0, TVector3 &tV,TVector3 &tY,  double tT0, double tU, bool debug = false);
                double tTimePreFreeze = NearestPreFreezeoutIntervalTime(tX0, tIniV, tPartPos, tZeroTime, tCurrentTime/*, tPartIter->eid == 10 && (TString(tPartOtherIter->GetParticleType()->GetName()).Contains("Spectator"))*/);
                TVector3 tPosPreFreeze = tCoordinates[tPartOtherIter->eid][0].Vect() + (tTimePreFreeze-tZeroTime) * tIniV;
                tPartOtherPos = tPosPreFreeze;
                tPartOtherVel = tIniV;
                tPartOtherAcc = TVector3(0, 0, 0);
/*
                if (tPartIter->eid == 10 || !USE_SPECTATORS) {
                    if ((TString(tPartOtherIter->GetParticleType()->GetName()).Contains("Spectator") && tPartOtherVel.Z() > 0) || !USE_SPECTATORS) {
                        double tIntervalBestI = this->Interval(tPartPos, tCurrentTime, tCoordinates[tPartOtherIter->eid][tBestI]);
                        double tIntervalExtrap = this->Interval(tPartPos, tCurrentTime, tPartOtherPos, tTimePreFreeze);

                        cout << " extrapolating before freezeout " << tTimePreFreeze*kHbarC << " fm, current time " << tCurrentTime*kHbarC << " fm" << endl;
                        cout << " interval " << tIntervalExtrap << endl;
                        cout << " at best i = " << tBestI << " would be " << endl;
                        cout << " interval " << tIntervalBestI << " = ";
                        cout << (tPartPos - tCoordinates[tPartOtherIter->eid][tBestI].Vect()).Mag2() << " - ";
                        cout << TMath::Power(tCurrentTime - tCoordinates[tPartOtherIter->eid][tBestI].T(),  2) << endl;
                        tPartPos.Print();
                        tCoordinates[tPartOtherIter->eid][tBestI].Print();
                    } else if (tPartOtherIter->eid == 11) {
                       // cout << " extrapolating someone else before freezeout " << endl;
                       // cout << " " << tPartOtherPos.X() << " " << tPartOtherPos.Y() << " " << tPartOtherPos.Z() << endl;
                       // cout << " " << tPartOtherVel.X() << " " << tPartOtherVel.Y() << " " << tPartOtherVel.Z() << endl;
                    }

                }
                ///     } else {
                ///         continue;
                ///     }
*/
            } else if (tPartOtherIter->pid == tPartOtherIter->fatherpid && false) {
                if (tPartIter->eid == 10 || !USE_SPECTATORS) {
                    if ((TString(tPartOtherIter->GetParticleType()->GetName()).Contains("Spectator") && tPartOtherVel.Z() > 0) || !USE_SPECTATORS) {
                        cout << " not extrapolating before freezeout " << tCurrentTime << " " << tBestI << endl;
                        double tIntervalBestIup = this->Interval(tPartPos, tCurrentTime, tCoordinates[tPartOtherIter->eid][tBestI+1],true);
                        double tIntervalBestI = this->Interval(tPartPos, tCurrentTime, tCoordinates[tPartOtherIter->eid][tBestI], true);
                        double tIntervalBestIdown = this->Interval(tPartPos, tCurrentTime, tCoordinates[tPartOtherIter->eid][tBestI-1], true);
                        cout << " intervals " << tIntervalBestIup << " " << tIntervalBestI << " " << tIntervalBestIdown << endl;
                        cout << " " << tPartOtherPos.X() << " " << tPartOtherPos.Y() << " " << tPartOtherPos.Z() << endl;
                        cout << " " << tPartOtherVel.X() << " " << tPartOtherVel.Y() << " " << tPartOtherVel.Z() << endl;
                    }
                }
            }

            TVector3 tPosDiff = tPartPos - tPartOtherPos;
            TVector3 tDirection = tPosDiff.Unit();

            TVector3 v = tPartOtherVel; // ()
            TVector3 a = tPartOtherAcc; // GeV
            TVector3 vecR = tPosDiff;   // 1/GeV
            Double_t R = vecR.Mag();
            ;                             // 1/GeV
            TVector3 hatR = tDirection;   // ()
            TVector3 u = hatR - v;        // ()
            Double_t Rdotu = vecR.Dot(u); // 1/GeV
            Double_t Rdota = vecR.Dot(a); // ()
            Double_t v2 = v.Mag2();       // ()

            if (false && !TMath::Finite(Rdotu))
            {
                cout << endl
                     << ", Condition: " << (tBestInterval < m_StepSize) << " " << tBestInterval << " " << tBestI << " " << tCoordinates[tPartOtherIter->eid].size() << " " << tPartOtherIter->pid << " " << tPartOtherIter->eid << " " << tParticles->size() << " " << endl;
                cout << Rdotu << " "
                     << " " << v.X() << " " << v.Y() << " " << v.Z() << endl;
                cout << " " << vecR.X() << " " << vecR.Y() << " " << vecR.Z() << " " << u.X() << " " << u.Y() << " " << u.Z() << endl;
                cout << " " << tPartPos.X() << " " << tPartPos.Y() << " " << tPartPos.Z() << " " << tPartOtherPos.X() << " " << tPartOtherPos.Y() << " " << tPartOtherPos.Z() << endl;
                return TVector3(0, 0, 0);
            }
            TVector3 tTermCoul = (tOtherQ * qe * R / TMath::Power(Rdotu, 3) * (1 - v2)) * u; // sqrt(GeV*fm) / GeV * GeV**3 = sqrt(GeV*fm) * GeV**2
            TVector3 tTermVelo = (tOtherQ * qe * R / TMath::Power(Rdotu, 3) * Rdota) * u;    // sqrt(GeV*fm) / GeV * GeV**3 = sqrt(GeV*fm) * GeV**2
            TVector3 tTermAcce = (tOtherQ * qe * R / TMath::Power(Rdotu, 3) * (-Rdotu)) * a; // sqrt(GeV*fm) * GeV**2
            if (false && (tBestI - 2) >= 0 && (tBestI + 2) < tCoordinates[tPartOtherIter->eid].size())
            {
                cout << tTermCoul.Mag() << " "
                     << tTermVelo.Mag() << " "
                     << tTermAcce.Mag() << endl;
            }

            TVector3 tE = tTermCoul + tTermVelo + tTermAcce; // sqrt(GeV*fm) * GeV**2
            // Using BAC-CAB rule for tV x tB having c * tB = tDirection x tE
            TVector3 tVxB = hatR * tV.Dot(tE) - hatR.Dot(tV) * tE; // sqrt(GeV*fm) * GeV**2
            if (tPartIter->eid == 10 && TString(tPartOtherIter->GetParticleType()->GetName()).Contains("Spectator")) {
            //    cout << "COUL TOTAL " << tTermCoulTotal.Mag() << endl; 
            }
            tTermCoulTotal += tTermCoul;
            tTermVeloTotal += tTermVelo;
            tTermAcceTotal += tTermAcce;
            tTermMagnTotal += tVxB;
            if (tPartIter->eid == 10 && TString(tPartOtherIter->GetParticleType()->GetName()).Contains("Spectator")) {
            //    cout << "SPECTATOR CHARGE " << tPartOtherIter->GetParticleType()->GetCharge() << endl;
            //    cout << "COUL TOTAL " << tTermCoulTotal.Mag() << " CONTRIB " << tTermCoul.Mag() << endl; 
            }

        }
    }
    TVector3 tForceTotal = (tQ * qe / TMath::Power(kHbarC, 2)) * (tTermCoulTotal + tTermVeloTotal + tTermAcceTotal + tTermMagnTotal); // GeV*fm / ((GeV*fm)**2) * GeV**2 = GeV / fm

    return tForceTotal;
}

CoulombAfterburner::CoulombAfterburner(int aSteps, double aStepSize, Double_t aR, AbstractEventSaver *aEventSaver) 
    : m_nSteps(aSteps), m_StepSize(aStepSize / kHbarC), mR(aR), mEventSaver(aEventSaver)
{
}
CoulombAfterburner::~CoulombAfterburner()
{
}

TGraph *initImpBGraph() {
    Float_t xs[] { 185.870 + 185.642, 181.902 + 181.760, 184.599 + 184.364, 180.046 + 179.368, 183.434 + 182.991, 183.980 + 183.701, 184.294 + 183.154, 185.464 + 184.580, 180.894 + 180.582, 182.516 + 181.986, 178.901 + 178.906, 185.637 + 185.712, 177.916 + 177.302, 178.391 + 178.344, 177.536 + 177.389, 180.008 + 180.720, 172.162 + 172.490, 173.693 + 173.214, 175.419 + 175.353, 171.275 + 171.576, 168.162 + 167.389, 170.919 + 171.987, 169.970 + 168.991, 172.380 + 172.458, 167.680 + 167.973, 165.792 + 166.663, 159.967 + 159.881, 160.395 + 160.109, 159.053 + 161.153, 162.047 + 159.706, 155.590 + 155.196, 157.467 + 158.524, 153.694 + 150.953, 150.958 + 150.781, 151.297 + 150.622, 148.507 + 148.534, 150.576 + 151.296, 149.398 + 148.339, 140.584 + 139.497, 140.006 + 141.817, 141.453 + 141.444, 140.423 + 139.828, 132.688 + 130.706, 129.880 + 130.140, 132.745 + 132.244, 130.558 + 133.015, 128.724 + 127.144, 123.967 + 124.291, 122.225 + 122.776, 120.866 + 121.504, 121.731 + 122.749, 117.899 + 119.436, 115.054 + 116.298, 110.464 + 110.776, 109.516 + 110.150, 106.308 + 106.569, 108.114 + 107.820, 106.558 + 105.327, 105.752 + 107.747, 100.821 +  99.785, 97.137 +  98.315, 95.343 +  93.013, 94.033 +  95.004, 90.977 +  90.945, 92.298 +  92.903, 87.081 +  86.855, 87.430 +  86.606, 85.404 +  86.900, 87.092 +  86.733, 79.703 +  78.480, 77.819 +  75.356, 78.144 +  78.572, 76.003 +  76.257, 73.348 +  72.836, 73.849 +  73.083, 69.608 +  69.788, 67.071 +  67.497, 68.992 +  69.597, 65.500 +  64.320, 59.964 +  60.350, 60.009 +  59.817, 57.177 +  57.197, 55.319 +  55.347, 57.712 +  56.912, 51.972 +  52.519, 52.459 +  53.849, 49.385 +  49.521, 46.880 +  47.787, 48.240 +  46.627, 47.702 +  47.086, 42.209 +  42.466, 42.925 +  42.240, 40.645 +  39.476, 38.318 +  39.389, 38.285 +  38.486, 36.290 +  36.220, 33.978 +  34.830, 30.991 +  31.838, 31.793 +  31.406, 31.227 +  30.930, 29.945 +  29.247, 26.719 +  26.933, 27.071 +  26.197, 23.788 +  23.629, 23.166 +  23.339, 23.378 +  22.931, 21.188 +  21.000, 20.186 +  20.473, 18.152 +  17.882, 18.901 +  18.587, 16.549 +  16.860, 15.772 +  15.781, 14.534 +  15.124, 13.647 +  13.353, 12.323 +  12.683, 11.912 +  12.562, 11.230 +  10.695, 9.650 +   9.998, 9.744 +   9.642, 8.714 +   8.857, 8.015 +   7.976, 6.854 +   7.172, 7.120 +   7.135, 6.599 +   6.713, 5.862 +   6.019, 5.284 +   5.415, 4.996 +   5.038, 4.308 +   4.397, 4.053 +   3.915, 3.678 +   3.704, 3.718 +   3.672, 3.078 +   2.983, 2.691 +   2.589, 2.658 +   2.629, 2.202 +   2.177, 1.894 +   1.866, 1.885 +   1.812, 1.544 +   1.540, 1.477 +   1.411, 1.222 +   1.221, 1.073 +   1.047, 0.955 +   0.941, 0.896 +   0.842, 0.778 +   0.759, 0.683 +   0.666, 0.591 +   0.570, 0.543 +   0.498, 0.442 +   0.427, 0.363 +   0.353, 0.343 +   0.333, 0.281 +   0.268, 0.256 +   0.250, 0.213 +   0.203, 0.199 +   0.190, 0.158 +   0.155, 0.141 +   0.130, 0.118 +   0.110, 0.106 +   0.100, 0.092 +   0.087, 0.079 +   0.071, 0.063 +   0.056, 0.058 +   0.052, 0.049 +   0.044, 0.042 +   0.037, 0.035 +   0.032, 0.031 +   0.027, 0.026 +   0.023, 0.022 +   0.019, 0.018 +   0.016, 0.015 +   0.013, 0.012 +   0.011, 0.012 +   0.010, 0.010 +   0.009, 0.008 +   0.007, 0.007 +   0.006, 0.005 +   0.005, 0.005 +   0.004, 0.004 +   0.003, 0.003 +   0.003, 0.003 +   0.002, 0.002 +   0.002, 0.002 +   0.002, 0.002 +   0.001, 0.001 +   0.001, 0.001 +   0.001, 0.001 +   0.001, 0.001 +   0.001, 0.001 +   0.001, 0.001 +   0.000, 0.000 +   0.000, 0.000 +   0.000, 0.000 +   0.000, 0.000 +   0.000, 0.000 +   0.000, 0.000 +   0.000, 0.000 +   0.000, 0.000 +   0.000, 0.000 +   0.000, 0.000 +   0.000, 0.000 +   0.000 };
    Float_t ys[] { 0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05, 2.15, 2.25, 2.35, 2.45, 2.55, 2.65, 2.75, 2.85, 2.95, 3.05, 3.15, 3.25, 3.35, 3.45, 3.55, 3.65, 3.75, 3.85, 3.95, 4.05, 4.15, 4.25, 4.35, 4.45, 4.55, 4.65, 4.75, 4.85, 4.95, 5.05, 5.15, 5.25, 5.35, 5.45, 5.55, 5.65, 5.75, 5.85, 5.95, 6.05, 6.15, 6.25, 6.35, 6.45, 6.55, 6.65, 6.75, 6.85, 6.95, 7.05, 7.15, 7.25, 7.35, 7.45, 7.55, 7.65, 7.75, 7.85, 7.95, 8.05, 8.15, 8.25, 8.35, 8.45, 8.55, 8.65, 8.75, 8.85, 8.95, 9.05, 9.15, 9.25, 9.35, 9.45, 9.55, 9.65, 9.75, 9.85, 9.95, 10.05, 10.15, 10.25, 10.35, 10.45, 10.55, 10.65, 10.75, 10.85, 10.95, 11.05, 11.15, 11.25, 11.35, 11.45, 11.55, 11.65, 11.75, 11.85, 11.95, 12.05, 12.15, 12.25, 12.35, 12.45, 12.55, 12.65, 12.75, 12.85, 12.95, 13.05, 13.15, 13.25, 13.35, 13.45, 13.55, 13.65, 13.75, 13.85, 13.95, 14.05, 14.15, 14.25, 14.35, 14.45, 14.55, 14.65, 14.75, 14.85, 14.95, 15.05, 15.15, 15.25, 15.35, 15.45, 15.55, 15.65, 15.75, 15.85, 15.95, 16.05, 16.15, 16.25, 16.35, 16.45, 16.55, 16.65, 16.75, 16.85, 16.95, 17.05, 17.15, 17.25, 17.35, 17.45, 17.55, 17.65, 17.75, 17.85, 17.95, 18.05, 18.15, 18.25, 18.35, 18.45, 18.55, 18.65, 18.75, 18.85, 18.95, 19.05, 19.15, 19.25, 19.35, 19.45, 19.55, 19.65, 19.75, 19.85, 19.95 };
    return new TGraph(200, xs, ys);
}

TGraph *CoulombAfterburner::m_gImpBNPart = initImpBGraph();

bool CoulombAfterburner::ParticleExists(std::list<Particle>::iterator &tPartIter, double tDecayTime, double tTime)
{
    if (tPartIter->decayed)
    {
        return tPartIter->t <= tTime && tDecayTime > tTime;
    }
    else
    {
        return tPartIter->t <= tTime;
    }
}
unsigned long CoulombAfterburner::NearestInterval(TVector3 &tPos, double tTime, std::vector<TLorentzVector> &tTrajectory, unsigned long tLastBestI, double &tBestInterval, bool debug)
{
    int tBestI = 0;
    double tPrevInterval = 1000000;
    bool goUp = true;
    if (tLastBestI > 0 && tLastBestI < (tTrajectory.size() - 1))
    {
        double tIntervalDown = this->Interval(tPos, tTime, tTrajectory[tLastBestI - 1]);
        double tInterval     = this->Interval(tPos, tTime, tTrajectory[tLastBestI    ]);
        double tIntervalUp   = this->Interval(tPos, tTime, tTrajectory[tLastBestI + 1]);
/*
        double tInterval =     TMath::Power(tTime - tTrajectory[tLastBestI    ].T(), 2) - (tPos - tTrajectory[tLastBestI    ].Vect()).Mag2();
        double tIntervalDown = TMath::Power(tTime - tTrajectory[tLastBestI - 1].T(), 2) - (tPos - tTrajectory[tLastBestI - 1].Vect()).Mag2();
        double tIntervalUp =   TMath::Power(tTime - tTrajectory[tLastBestI + 1].T(), 2) - (tPos - tTrajectory[tLastBestI + 1].Vect()).Mag2();
*/
        if (TMath::Abs(tIntervalUp) > TMath::Abs(tInterval) && TMath::Abs(tIntervalDown) > TMath::Abs(tInterval))
        {
            tBestInterval = TMath::Sqrt(TMath::Abs(tInterval)); // Zostalo poprawione na sqrt(abs(...))
            return tLastBestI;
        }
        else
        {
            goUp = TMath::Abs(tIntervalUp) < TMath::Abs(tInterval);
        }
    }
    int tSign = 0;
    if (debug)
        cout << "NearestInterval: " << tLastBestI << " " << tTrajectory.size() << " " << goUp << endl;
    for (unsigned long i = tLastBestI; i < tTrajectory.size() && i >= 0; i = goUp ? i + 1 : i - 1)
    {
        double tInterval = this->Interval(tPos, tTime, tTrajectory[i ]);
/*
        double tInterval = TMath::Power(tTime - tTrajectory[i].T(), 2) - (tPos - tTrajectory[i].Vect()).Mag2();
*/
        int tNewSign = TMath::Sign(1.0, tInterval);
        if (debug)
        {
            cout << "NearestInterval: i = " << i << " " << tInterval << " " << tTime << " " << tTrajectory[i].T() << " " << (tPos - tTrajectory[i].Vect()).Mag2() << " " << (tPos - tTrajectory[i].Vect()).Mag() << " " << tPos.Mag() << " " << tTrajectory[i].Vect().Mag() << " " << tSign << endl;
        }
        if (i != tLastBestI)
        {
            if ((tSign != tNewSign && tSign != 0) || TMath::Abs(tInterval) > TMath::Abs(tPrevInterval))
            {
                tBestInterval = TMath::Sqrt(TMath::Abs(tPrevInterval)); // Zostalo poprawione na  zwracane poprzednie
                tBestI = goUp ? i - 1 : i + 1;                          // i indeks tak samo
                break;
            }
        }
        tPrevInterval = tInterval;
        tSign = tNewSign;
    }
    return tBestI;
}
double CoulombAfterburner::NearestPreFreezeoutIntervalTime(TVector3 &tX0, TVector3 &tV, TVector3 &tY, double tT0, double tU, bool debug)
{
    double tX0dotV = tX0.Dot(tV);
    double tX0dotY = tX0.Dot(tY);
    double tVdotY = tV.Dot(tY);
    double tX02 = tX0.Mag2();
    double tV2 = tV.Mag2();
    double tY2 = tY.Mag2();

    // A, B, C should be correct now, checked with CAS
    double tA = 1 - tV2;
    double tB = 2*(-tU + tV2*tT0 - tX0dotV + tVdotY);
    double tC = tU*tU - tX02 - tV2*tT0*tT0 - tY2 + 2*tX0dotV*tT0 + 2*tX0dotY - 2*tVdotY*tT0;
    double tDelta = tB * tB - 4 * tA * tC;
    double t1 = (-tB - TMath::Sqrt(tDelta)) / (2 * tA);
    double t2 = (-tB + TMath::Sqrt(tDelta)) / (2 * tA);
    if (debug) {
        cout << "PARAMS tT0 " << tT0 << " tU " << tU << endl;
        cout << "tV:  "; tV.Print();
        cout << "tX0: "; tX0.Print();
        cout << "tY:  "; tY.Print();
        cout << "tdiff " << t1-tT0 << " " << t2-tT0 << endl;
        cout << "dt " << t1-tU << " " << t2-tU << endl;
        cout << "tX1: "; (tX0 + tV*(t1-tT0)).Print();
        cout << "tX2: "; (tX0 + tV*(t2-tT0)).Print();
        cout << "dX1: "; (tY - (tX0 + tV*(t1-tT0))).Print();
        cout << "dX2: "; (tY - (tX0 + tV*(t2-tT0))).Print();
        cout << "|dx1|^2: " << (tY - (tX0 + tV*(t1-tT0))).Mag2() << endl;
        cout << "|dx2|^2: " << (tY - (tX0 + tV*(t2-tT0))).Mag2() << endl;
        cout << "t1 " << t1 << " t2 " << t2 << endl; 
        cout << "sol1 " << t1*t1*tA + t1*tB + tC << " sol2 " << t2*t2*tA + t2*tB + tC << endl;
        double int1 = this->Interval(tX0 + tV*(t1-tT0), t1, tY, tU, true);
        double int2 = this->Interval(tX0 + tV*(t2-tT0), t2, tY, tU, true);
        double int1_a = this->Interval(tX0 + tV*(t1+1/kHbarC-tT0), t1+1/kHbarC, tY, tU);
        double int1_b = this->Interval(tX0 + tV*(t1-1/kHbarC-tT0), t1-1/kHbarC, tY, tU);
        cout << "int1 " << int1 << " int2 " << int2 << endl;
        cout << "int1_a " << int1_a << " int1_b " << int1_b << endl;
    }
    // Absolute time coordinate !
    return t1 <= t2 ? t1 : t2;
}

Int_t CoulombAfterburner::TotalCharge(std::list<Particle> *tParticles)
{
    Int_t tResult = 0;
    for (const auto &tParticle : *tParticles)
    {
        // Fixed isospin means:
        // Don't add deuteron charge, because it's either caseB or we use deuterons from dedicated caseC list
        // Total charge is the charge of primary particles
        if (tParticle.pid != 80000 && tParticle.pid == tParticle.fatherpid && tParticle.GetParticleType() != nullptr)
        {
            tResult += tParticle.GetParticleType()->GetCharge();
        }
    }
    return tResult;
}

ParticleType *CoulombAfterburner::SpectatorType(std::list<Particle> *tParticles, Float_t &tSpectatorX) {

    ParticleType *tSpectatorType = new ParticleType;
    Float_t tA = 197.;
    Float_t tZ = 79.; // because they are used in fractions
    Int_t tTotalCharge = TotalCharge(tParticles);
    Int_t tSpectatorCharge = (2 * tZ - tTotalCharge)/2;
    cout << "TOTAL CHARGE " << tTotalCharge << " SPEC. CHARGE " << tSpectatorCharge << endl; 
    Int_t tNp = tSpectatorCharge;
    Int_t tNn = tNp * (tA - tZ) / tZ;
    Int_t tNpart = tTotalCharge * tA / tZ;
    Float_t tImpB = m_gImpBNPart->Eval(tNpart) / kHbarC; // Distances are in 1/MeV during calculations
    Float_t tSpectatorMass = tNp * 0.9382720 + tNn * 0.9395653;
    tSpectatorType->SetNumber(777777);
    tSpectatorType->SetName("Spectator");
    tSpectatorType->SetMass(tSpectatorMass);
    tSpectatorType->SetBarionN(tNp + tNn);
    tSpectatorType->SetI3((tNp - tNn)/2);

    // At Npart = 0, nuclei just grazing each other, spectator X would be the radius of the nucleus.
    // At Npart -> A, a head on collision, the center of mass of the spectators would also approach
    // the radius of the nucleus. so we assume that the center of mass will always be at x = tR
    Float_t tSpectatorXtoR = 1.;//tNpart / (2 * tA) + tImpB / (2 * tA) * (2 * tA - tNpart);
    Float_t tR = 1.2 * TMath::Power(tA, 1./3) / kHbarC;
    tSpectatorX = tSpectatorXtoR * tR;
    cout << "DEFINING SPECTATOR X " << tR << " x " << tSpectatorXtoR << " -> " << tSpectatorX << endl;

    return tSpectatorType;
}

ParticleType *CoulombAfterburner::AddSpectators(std::list<Particle> *tParticles, Float_t tEarliestTime)
{
    Float_t tSpectatorX;
    ParticleType *tSpectatorType = SpectatorType(tParticles, tSpectatorX);
    Float_t tSpectatorMass = tSpectatorType->GetMass();

    Float_t tSpectatorBeta = 0.63; // TODO: make an external parameter
    Float_t tSpectatorGamma = 1/TMath::Sqrt(1 - tSpectatorBeta*tSpectatorBeta);
    Float_t tSpectatorPz = tSpectatorGamma * tSpectatorMass * tSpectatorBeta;
    Float_t tSpectatorE = tSpectatorGamma * tSpectatorMass;
    Float_t tSpectatorZ = tSpectatorBeta * tEarliestTime;
    // Linear interpolation between:
    // tNpart = 2*tA -> tSpectatorX = tR
    // tNpart = 0 -> tSpectatorX = tImpB
    //Float_t tSpectatorX = tImpB * tA / (2 * tA - tNpart);

    if (USE_SPECTATORS) {
        Particle tSpectator1(tSpectatorType);
        tSpectator1.SetParticlePX(tSpectatorE, 0, 0, tSpectatorPz, tEarliestTime, tSpectatorX, 0, tSpectatorZ);
       // tSpectator1.SetParticlePX(tSpectatorE, 0, 0, tSpectatorPz, 0, tSpectatorX, 0, 0);
        tParticles->push_back(tSpectator1);

        Particle tSpectator2(tSpectatorType);
        tSpectator1.SetParticlePX(tSpectatorE, 0, 0, -tSpectatorPz, tEarliestTime, -tSpectatorX, 0, -tSpectatorZ);
       // tSpectator2.SetParticlePX(tSpectatorE, 0, 0, -tSpectatorPz, 0, -tSpectatorX, 0, 0);
        tParticles->push_back(tSpectator2);
    }
    return tSpectatorType;
}

void CoulombAfterburner::Apply(Event *tEvent)
{
    if (m_nSteps <= 0) return;

    bool debug = false;

    std::list<Particle> *tParticles = tEvent->GetParticleList();
    
    unsigned long N = tParticles->size() + 2;
    std::vector<TLorentzVector> *tCoordinates = new std::vector<TLorentzVector>[N];
    std::vector<TVector3> *tVelocities = new std::vector<TVector3>[N];
    std::vector<TVector3> *tAccelerations = new std::vector<TVector3>[N];
    double *tDecayTimes = new double[N];
    unsigned int *tFirstStepExists = new unsigned int[N];
    unsigned int *tLastStepExists = new unsigned int[N];
    double tEarliestTime = std::numeric_limits<double>::quiet_NaN();
    auto tPartIter = tParticles->begin();
    for (; tPartIter != tParticles->end(); ++tPartIter)
    {
        tCoordinates[tPartIter->eid] = std::vector<TLorentzVector>();
        if (std::isnan(tEarliestTime) || tPartIter->t < tEarliestTime)
        {
            tEarliestTime = tPartIter->t;
        }
        if (tPartIter->decayed)
        {
            std::list<Particle>::iterator tChildPartIter = tParticles->begin();
            for (; tChildPartIter != tParticles->end(); ++tChildPartIter)
            {
                if (tChildPartIter->fathereid == tPartIter->eid)
                {
                    tDecayTimes[tPartIter->eid] = tChildPartIter->t;
                    break;
                }
            }
        }
        unsigned int tFirstStep = (unsigned int)TMath::Floor((tPartIter->t - tEarliestTime) / m_StepSize);
        if (tFirstStep >= m_nSteps)
        {
            tFirstStep = m_nSteps - 1;
        }
        tFirstStepExists[tPartIter->eid] = tFirstStep;
        if (tPartIter->decayed)
        {
            unsigned int tLastStep = (unsigned int)TMath::Floor((tDecayTimes[tPartIter->eid] - tEarliestTime) / m_StepSize) - 1;
            if (tLastStep >= m_nSteps)
            {
                tLastStep = m_nSteps - 1;
            }
            tLastStepExists[tPartIter->eid] = tLastStep;
        }
        else
        {
            tLastStepExists[tPartIter->eid] = m_nSteps - 1;
        }
    }

    ParticleType *tSpectatorType = AddSpectators(tParticles, tEarliestTime);

    int **tBestIs = new int *[N];
    for (unsigned int ai = 0; ai < N; ++ai)
    {
        tBestIs[ai] = new int[N];
        for (unsigned int aj = 0; aj < N; ++aj)
        {
            tBestIs[ai][aj] = -1;
        }
    }

    TVector3 *tTotalForcesLast = NULL;
    for (unsigned int iStep = 0; iStep < m_nSteps; ++iStep)
    {

        // cout << "Step " << iStep << " of " << m_nSteps << " particles: " << tParticles->size() << endl;

        double tCurrentTime = tEarliestTime + iStep * m_StepSize;
        TVector3 *tTotalForces = new TVector3[N];
        std::list<Particle>::iterator tPartIter = tParticles->begin();
        for (; tPartIter != tParticles->end(); ++tPartIter)
        {
/*
            if (TString(tPartIter->GetParticleType()->GetName()).Contains("Spectator")) {
                cout << " at step " << iStep << " spectator is ";
                cout << "x(" << tPartIter->x << ", " << tPartIter->y << ", " << tPartIter->z << "), ";
                cout << "p(" << tPartIter->px << ", " << tPartIter->py << ", " << tPartIter->pz << ")" << endl;
            }
*/
            Int_t tQ = tPartIter->GetParticleType()->GetCharge();
            TVector3 tPartPos(tPartIter->x, tPartIter->y, tPartIter->z); // 1/GeV
            if (tQ != 0 && iStep >= tFirstStepExists[tPartIter->eid] && iStep <= tLastStepExists[tPartIter->eid])
            {
                tCoordinates[tPartIter->eid].push_back(TLorentzVector(tPartPos, tCurrentTime)); // 1/GeV
                TVector3 tPartMom(tPartIter->px, tPartIter->py, tPartIter->pz);
                double tP2 = tPartMom.Mag2();
                double tM2 = TMath::Power(tPartIter->GetParticleType()->GetMass(), 2);
                TVector3 tV(0, 0, 0);
                if (tP2 > 0 && TMath::Finite(tP2))
                {
                    tV = 1. / TMath::Sqrt(tM2 / tP2 + 1) * tPartMom.Unit();
                }
                tVelocities[tPartIter->eid].push_back(tV); // []
                if (debug && !TMath::Finite(tV.Mag()))
                {
                    cout << "In step " << iStep << " ERROR, infinite velocity " << tPartMom.X() << " " << tPartMom.Y() << " " << tPartMom.Z() << endl;
                    cout << "In step " << iStep << " ERROR, infinite velocity " << tV.X() << " " << tV.Y() << " " << tV.Z() << endl;
                }

                if (iStep == tFirstStepExists[tPartIter->eid])
                {
                    tAccelerations[tPartIter->eid].push_back(TVector3(0, 0, 0));
                }
                else
                {
                    if (tFirstStepExists[tPartIter->eid] < 0)
                    {
                        cerr << tParticles->size() << " " << N << " " << tPartIter->eid << " " << sizeof(tFirstStepExists) / sizeof(int *) << " " << tFirstStepExists[tPartIter->eid] << " " << iStep - tFirstStepExists[tPartIter->eid] << " " << iStep - tFirstStepExists[tPartIter->eid] - 1 << " " << iStep << " " << tVelocities[tPartIter->eid].size() << endl;
                    }
                    TVector3 &tV1 = tVelocities[tPartIter->eid][iStep - tFirstStepExists[tPartIter->eid] - 1];
                    TVector3 &tV2 = tVelocities[tPartIter->eid][iStep - tFirstStepExists[tPartIter->eid]];
                    TVector3 tAcce = 1. / (m_StepSize) * (tV2 - tV1);
                    if (debug && !TMath::Finite(tAcce.Z()))
                    {
                        cout << "Error, infinite acceleration " << iStep << " " << tVelocities[tPartIter->eid].size() << endl;
                        cout << "Error, infinite acceleration " << tAcce.X() << " " << tAcce.Y() << " " << tAcce.Z() << endl;
                        cout << "Error, infinite acceleration " << tV1.X() << " " << tV1.Y() << " " << tV1.Z() << endl;
                        cout << "Error, infinite acceleration " << tV2.X() << " " << tV2.Y() << " " << tV2.Z() << endl;
                        cout << "Error, infinite acceleration " << tFirstStepExists[tPartIter->eid] << " " << tLastStepExists[tPartIter->eid] << endl;
                    }
                    tAccelerations[tPartIter->eid].push_back(tAcce); // GeV
                }
            }
        }
        // Calculate force at old position
        tPartIter = tParticles->begin();
        for (; tPartIter != tParticles->end(); ++tPartIter)
        {

            TVector3 tForceTotal;
            // Possibly reuse calculation from the previous step.
            if (tTotalForcesLast != NULL) {
                tForceTotal = tTotalForcesLast[tPartIter->eid];
            } else {
                Int_t tQ = tPartIter->GetParticleType()->GetCharge();
                if (tQ != 0 && iStep >= tFirstStepExists[tPartIter->eid] && iStep <= tLastStepExists[tPartIter->eid])
                {
                    tForceTotal = CalculateTotalForce(tPartIter, tParticles, tBestIs, tCurrentTime,
                                                  tCoordinates, tVelocities, tAccelerations);
                }
            }
            tTotalForces[tPartIter->eid] = tForceTotal;
        }
        // Update positions (x_{n+1} = ...) using velocity and force (for acceleration) calculated at old position
        tPartIter = tParticles->begin();
        for (; tPartIter != tParticles->end(); ++tPartIter)
        {
            TVector3 tPartPos(tPartIter->x, tPartIter->y, tPartIter->z);
            TVector3 tPartMom(tPartIter->px, tPartIter->py, tPartIter->pz);

            TVector3 tForceTotal = tTotalForces[tPartIter->eid];
            TVector3 tImpulse = (m_StepSize * kHbarC) * tForceTotal; // GeV/fm * (1/GeV) * GeV*fm = GeV

            double tP2 = tPartMom.Mag2();
            double fM = tPartIter->GetParticleType()->GetMass();
            double tM2 = TMath::Power(fM, 2);

            if (tP2 > 0 && TMath::Finite(tP2))
            {                                                                           // guardian
                tPartPos += (m_StepSize / TMath::Sqrt(tM2 / tP2 + 1)) * tPartMom.Unit(); // 1/GeV
                            + (0.5 / fM / kHbarC) * TMath::Power(m_StepSize * kHbarC, 2) * tForceTotal; // 1/GeV (first 1/kHbarC should give this)
                if (debug && !TMath::Finite(tPartPos.Mag()))
                {
                    cout << "updating position: " << tM2 << " " << tP2 << " " << tPartMom.Mag() << endl;
                }
            }

            tPartIter->t  = tCurrentTime;
            tPartIter->x  = tPartPos.X();                                                                                                   
            tPartIter->y  = tPartPos.Y();
            tPartIter->z  = tPartPos.Z();

        }
        if (tTotalForcesLast != NULL) {                
            delete [] tTotalForcesLast;
        }
        tTotalForcesLast = tTotalForces;

        // Calculate force at new position
        tPartIter = tParticles->begin();
        for (; tPartIter != tParticles->end(); ++tPartIter)
        {

            TVector3 tForceTotal;
            Int_t tQ = tPartIter->GetParticleType()->GetCharge();
            if (tQ != 0 && iStep >= tFirstStepExists[tPartIter->eid] && iStep <= tLastStepExists[tPartIter->eid])
            {
                tForceTotal = CalculateTotalForce(tPartIter, tParticles, tBestIs, tCurrentTime,
                                                  tCoordinates, tVelocities, tAccelerations);
            }
            tTotalForces[tPartIter->eid] = tForceTotal;
        }
        // Update momenmta (p_{n+1} = ...) using average force from old and new position
        tPartIter = tParticles->begin();
        for (; tPartIter != tParticles->end(); ++tPartIter)
        {
            TVector3 tPartPos(tPartIter->x, tPartIter->y, tPartIter->z);
            TVector3 tPartMom(tPartIter->px, tPartIter->py, tPartIter->pz);

            TVector3 tForceTotal = 0.5 * (tTotalForces[tPartIter->eid] + tTotalForcesLast[tPartIter->eid]);
           // TVector3 tForceTotal = TVector3(0,0,0);
            TVector3 tImpulse = (m_StepSize * kHbarC) * tForceTotal; // GeV/fm * (1/GeV) * GeV*fm = GeV

            double tP2 = tPartMom.Mag2();
            double fM = tPartIter->GetParticleType()->GetMass();
            double tM2 = TMath::Power(fM, 2);

            tPartMom += tImpulse; // GeV

            double tE = TMath::Sqrt(tM2 + tP2);

            tPartIter->e = tE;
            tPartIter->px = tPartMom.X();
            tPartIter->py = tPartMom.Y();
            tPartIter->pz = tPartMom.Z();
        }
        if (mEventSaver != nullptr) {
            mEventSaver->Save(tEvent,NULL,tEvent->GetEventID());
        }
    }
    // Writing back to particles' objects original position, to not change HBT radii
    tPartIter = tParticles->begin();
    for (; tPartIter != tParticles->end(); ++tPartIter)
    {
        if (tCoordinates[tPartIter->eid].size() > 0)
        {
            Float_t tPartTime =  tCoordinates[tPartIter->eid].front().T();
            TVector3 tPartPos = tCoordinates[tPartIter->eid].front().Vect();
            TVector3 tPartPosFinal(tPartIter->x, tPartIter->y, tPartIter->z);
            tPartIter->t = tPartTime;
            tPartIter->x = tPartPos.X();
            tPartIter->y = tPartPos.Y();
            tPartIter->z = tPartPos.Z();
        }
    }
    delete[] tBestIs;

    delete[] tCoordinates;
    delete[] tVelocities;
    delete[] tAccelerations;
    delete[] tDecayTimes;
    delete[] tFirstStepExists;
    delete[] tLastStepExists;

    //tParticles->remove_if([](Particle &p) { return TString(p.GetParticleType()->GetName()).Contains("Spectator"); });
    //delete tSpectatorType;
}
