#include "CoulombAfterburner.h"
#include "Configurator.h"
#include "Parser.h"
#include "TVector3.h"
#include "TMath.h"
#include <limits>
#include <iostream>
#include "THGlobal.h"

using namespace std;
extern TString  sModelINI;


TVector3 CoulombAfterburner::CalculateTotalForce(
        std::list<Particle>::iterator &tPartIter, std::list<Particle> *tParticles, int ** tBestIs, double tCurrentTime,
        std::vector<TLorentzVector> *tCoordinates, std::vector<TVector3> *tVelocities, std::vector<TVector3> *tAccelerations
) {
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
    TVector3 tV(0,0,0);
    if (tP2 > 0 && TMath::Finite(tP2)) { 
        tV = 1. / TMath::Sqrt(tM2/tP2 + 1) * tPartMom.Unit();
    }
    std::list<Particle>::iterator tPartOtherIter = tParticles->begin();
    for ( ; tPartOtherIter != tParticles->end(); ++tPartOtherIter) {

        if (tPartOtherIter->eid == tPartIter->eid) continue;

        Int_t tOtherQ = tPartOtherIter->GetParticleType()->GetCharge();
        if (tOtherQ == 0) continue;

        double tBestInterval;
        int tLastBestI = tBestIs[tPartIter->eid][tPartOtherIter->eid] > -1 ? tBestIs[tPartIter->eid][tPartOtherIter->eid] : 0;
        int tBestI = NearestInterval(tPartPos, tCurrentTime, tCoordinates[tPartOtherIter->eid], tLastBestI, tBestInterval, false);
        tBestIs[tPartIter->eid][tPartOtherIter->eid] = tBestI;

        // The best interval is at the end of trajectory of other particle and it's space-like: 
        // the particles are not causally connected
        if ((tBestI + 1) == tCoordinates[tPartOtherIter->eid].size() && tBestInterval < 0) {
            continue;
        }

        if (tBestI >= 0 && tBestI < tCoordinates[tPartOtherIter->eid].size()) {
            if (false && (tBestI - 2) >= 0 && (tBestI + 2) < tCoordinates[tPartOtherIter->eid].size()) {
                double tInterval = TMath::Power(tCurrentTime - tCoordinates[tPartOtherIter->eid][tBestI].T(), 2) - (tPartPos  - tCoordinates[tPartOtherIter->eid][tBestI].Vect()).Mag2();
                double tIntervalDown2 = TMath::Power(tCurrentTime - tCoordinates[tPartOtherIter->eid][tBestI-2].T(), 2) - (tPartPos  - tCoordinates[tPartOtherIter->eid][tBestI-2].Vect()).Mag2();
                double tIntervalDown = TMath::Power(tCurrentTime - tCoordinates[tPartOtherIter->eid][tBestI-1].T(), 2) - (tPartPos  - tCoordinates[tPartOtherIter->eid][tBestI-1].Vect()).Mag2();
                double tIntervalUp = TMath::Power(tCurrentTime - tCoordinates[tPartOtherIter->eid][tBestI+1].T(), 2) - (tPartPos  - tCoordinates[tPartOtherIter->eid][tBestI+1].Vect()).Mag2();
                double tIntervalUp2 = TMath::Power(tCurrentTime - tCoordinates[tPartOtherIter->eid][tBestI+2].T(), 2) - (tPartPos  - tCoordinates[tPartOtherIter->eid][tBestI+2].Vect()).Mag2();
                cout << " tBestI = " << tBestI << " " << tBestInterval << " "
                    << "Intervals " << tIntervalUp2 << " " << tIntervalUp << " " << tInterval << " " << tIntervalDown << " " << tIntervalDown2 << " " << endl
                    << tCurrentTime << " " << tCoordinates[tPartOtherIter->eid][tBestI].T() << " "
                    << (tPartPos  - tCoordinates[tPartOtherIter->eid][tBestI].Vect()).Mag() << endl;

                cout << "X, V, A " << tCoordinates[tPartOtherIter->eid][tBestI].Vect().Mag() << " " 
                    << tVelocities[tPartOtherIter->eid][tBestI].Mag() << " " 
                    << tAccelerations[tPartOtherIter->eid][tBestI].Mag() << endl;
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
            if (tPartOtherIter->pid == tPartOtherIter->fatherpid && tBestI == 0) {
                TVector3 tX0 = tCoordinates[tPartOtherIter->eid][0].Vect() - tPartPos;
                TVector3 tIniV = tVelocities[tPartOtherIter->eid][0];
                double tTimePreFreeze = NearestPreFreezeoutIntervalTime(tX0, tIniV, tCurrentTime);
                TVector3 tPosPreFreeze = tCoordinates[tPartOtherIter->eid][0].Vect() + tTimePreFreeze * tIniV;
                tPartOtherPos = tPosPreFreeze;
                tPartOtherVel = tIniV;
                tPartOtherAcc = TVector3(0,0,0);
                ///     } else {
                ///         continue;
                ///     }
        }

        TVector3 tPosDiff = tPartPos - tPartOtherPos;
        TVector3 tDirection = tPosDiff.Unit();

        TVector3 v = tPartOtherVel; // ()
        TVector3 a = tPartOtherAcc; // GeV
        TVector3 vecR = tPosDiff;   // 1/GeV
        Double_t R = vecR.Mag();;   // 1/GeV
        TVector3 hatR = tDirection; // ()
        TVector3 u = hatR - v;      // ()
        Double_t Rdotu = vecR.Dot(u); // 1/GeV
        Double_t Rdota = vecR.Dot(a); // ()
        Double_t v2 = v.Mag2(); // ()

        if (false && !TMath::Finite(Rdotu)) {
            cout << endl << ", Condition: " << (tBestInterval < m_StepSize) << " " << tBestInterval << " " << tBestI << " " << tCoordinates[tPartOtherIter->eid].size() << " " << tPartOtherIter->pid << " " << tPartOtherIter->eid << " " << tParticles->size() << " " << endl;
            cout << Rdotu << " " << " " << v.X() << " " << v.Y() << " " << v.Z() << endl;
            cout << " " << vecR.X() << " " << vecR.Y() << " " << vecR.Z() << " " << u.X() << " " << u.Y() << " " << u.Z() << endl;
            cout << " " << tPartPos.X() << " " << tPartPos.Y() << " " << tPartPos.Z() << " " << tPartOtherPos.X() << " " << tPartOtherPos.Y() << " " << tPartOtherPos.Z() << endl;
            return TVector3(0, 0, 0);
        }
        TVector3 tTermCoul = (tOtherQ * qe * R / TMath::Power(Rdotu, 3) * (1-v2)) * u; // sqrt(GeV*fm) / GeV * GeV**3 = sqrt(GeV*fm) * GeV**2
        TVector3 tTermVelo = (tOtherQ * qe * R / TMath::Power(Rdotu, 3) * Rdota) * u;  // sqrt(GeV*fm) / GeV * GeV**3 = sqrt(GeV*fm) * GeV**2
        TVector3 tTermAcce = (tOtherQ * qe * R / TMath::Power(Rdotu, 3) * (-Rdotu)) * a; // sqrt(GeV*fm) * GeV**2
        if (false && (tBestI - 2) >= 0 && (tBestI + 2) < tCoordinates[tPartOtherIter->eid].size()) {
            cout << tTermCoul.Mag() << " " 
                << tTermVelo.Mag() << " " 
                << tTermAcce.Mag() << endl;
        }

        TVector3 tE = tTermCoul + tTermVelo + tTermAcce;  // sqrt(GeV*fm) * GeV**2
        // Using BAC-CAB rule for tV x tB having c * tB = tDirection x tE
        TVector3 tVxB = hatR * tV.Dot(tE) - hatR.Dot(tV) * tE; // sqrt(GeV*fm) * GeV**2
        tTermCoulTotal += tTermCoul;
        tTermVeloTotal += tTermVelo;
        tTermAcceTotal += tTermAcce;
        tTermMagnTotal += tVxB;
        }
    }
    TVector3 tForceTotal = (tQ * qe / TMath::Power(kHbarC, 2)) * (tTermCoulTotal + tTermVeloTotal + tTermAcceTotal + tTermMagnTotal); // GeV*fm / ((GeV*fm)**2) * GeV**2 = GeV / fm

    return tForceTotal;

}

CoulombAfterburner::CoulombAfterburner(int aSteps, double aStepSize) : m_nSteps(aSteps), m_StepSize(aStepSize / kHbarC) {
}
CoulombAfterburner::~CoulombAfterburner() {
}
bool CoulombAfterburner::ParticleExists(std::list<Particle>::iterator &tPartIter, double tDecayTime, double tTime) {
    if (tPartIter->decayed) {
        return tPartIter->t <= tTime && tDecayTime > tTime;
    } else {
        return tPartIter->t <= tTime;
    }
}
int CoulombAfterburner::NearestInterval(TVector3 &tPos, double tTime, std::vector<TLorentzVector> &tTrajectory, int tLastBestI, double &tBestInterval, bool debug) {
    int tBestI = 0;
    double tPrevInterval = 1000000;
    bool goUp = true;
    if (tLastBestI > 0 && tLastBestI < (tTrajectory.size() - 1)) {
        double tInterval = TMath::Power(tTime - tTrajectory[tLastBestI].T(), 2) - (tPos  - tTrajectory[tLastBestI].Vect()).Mag2();
        double tIntervalDown = TMath::Power(tTime - tTrajectory[tLastBestI-1].T(), 2) - (tPos  - tTrajectory[tLastBestI-1].Vect()).Mag2();
        double tIntervalUp = TMath::Power(tTime - tTrajectory[tLastBestI+1].T(), 2) - (tPos  - tTrajectory[tLastBestI+1].Vect()).Mag2();
        if (TMath::Abs(tIntervalUp) > TMath::Abs(tInterval) && TMath::Abs(tIntervalDown) > TMath::Abs(tInterval)) {
            tBestInterval = TMath::Sqrt(TMath::Abs(tInterval)); // Zostalo poprawione na sqrt(abs(...))
            return tLastBestI;
        } else {
            goUp = TMath::Abs(tIntervalUp) < TMath::Abs(tInterval);
        }
    }
    int tSign = 0;
    if (debug) cout << "NearestInterval: " << tLastBestI << " " << tTrajectory.size() << " " << goUp << endl;
    for (int i = tLastBestI; i < tTrajectory.size() && i >= 0 ; i = goUp ? i + 1 : i - 1) {
        double tInterval = TMath::Power(tTime - tTrajectory[i].T(), 2) - (tPos  - tTrajectory[i].Vect()).Mag2();
        int tNewSign = TMath::Sign(1.0, tInterval);
        if (debug) {
            cout << "NearestInterval: i = " << i << " " << tInterval << " " << tTime << " " << tTrajectory[i].T() << " " << (tPos  - tTrajectory[i].Vect()).Mag2() << " " << (tPos  - tTrajectory[i].Vect()).Mag() << " " << tPos.Mag() << " " << tTrajectory[i].Vect().Mag() << " " << tSign << endl;
        }
        if (i == tLastBestI) {
            tPrevInterval = tInterval;
        } else {
            // if ((tSign != TMath::Sign(1.0, tPrevInterval - tInterval)) && tSign != 0) {
            if (tSign != tNewSign && tSign != 0) {
                tBestInterval = TMath::Sqrt(TMath::Abs(tPrevInterval)); // Zostalo poprawione na  zwracane poprzednie
                tBestI = goUp ? i - 1 : i + 1; // i indeks tak samo
                break;
            }
            // tSign = TMath::Sign(1.0, tPrevInterval - tInterval);
            tPrevInterval = tInterval;
        }
        tSign = tNewSign;
        }
        return tBestI;
    }
    double CoulombAfterburner::NearestPreFreezeoutIntervalTime(TVector3 &tX0, TVector3 &tV, double tCurrentTime) {
        double tX0dotV = tX0.Dot(tV);
        double tV2 = tV.Mag2();
        // A, B, C should be correct
        double tA = 1 - tV2;
        double tB = -2*(tX0dotV + tCurrentTime);
        double tC = tCurrentTime * tCurrentTime - tX0.Mag2();
        double tDelta = tB*tB - 4*tA*tC;
        double t1 = (-tB - TMath::Sqrt(tDelta)) / (2*tA);
        double t2 = (-tB + TMath::Sqrt(tDelta)) / (2*tA);
        return t1 <= 0 ? t1 : t2;
    }
    void CoulombAfterburner::Apply(Event *tEvent)
    {
        bool debug = false;

        Int_t tPrimaryCharge = 0.;
        std::list<Particle>::iterator tPartIter = tEvent->GetParticleList()->begin();
        for ( ; tPartIter != tEvent->GetParticleList()->end(); ++tPartIter) {
            if (tPartIter->pid == tPartIter->fatherpid) {
                tPrimaryCharge += tPartIter->GetParticleType()->GetCharge();
            }
        }
        std::list<Particle>* tParticles = tEvent->GetParticleList();
        Int_t N = tParticles->size();
        TVector3 *tInitialVelocities = new TVector3[N];;
        std::vector<TLorentzVector> *tCoordinates = new std::vector<TLorentzVector>[N];
        std::vector<TVector3> *tVelocities = new std::vector<TVector3>[N];
        std::vector<TVector3> *tAccelerations = new std::vector<TVector3>[N];
        double *tDecayTimes = new double[N];
        unsigned int *tFirstStepExists = new unsigned int[N];
        unsigned int *tLastStepExists  = new unsigned int[N];
        double tEarliestTime = std::numeric_limits<double>::quiet_NaN();
        int favoriteEID = -1;
        double favoriteE0 = -1.;
        tPartIter = tParticles->begin();
        for ( ; tPartIter != tParticles->end(); ++tPartIter) {
            tCoordinates[tPartIter->eid] = std::vector<TLorentzVector>();
            if (std::isnan(tEarliestTime) || tPartIter->t < tEarliestTime) {
                tEarliestTime = tPartIter->t;
            }
            if (tPartIter->decayed) {
                std::list<Particle>::iterator tChildPartIter = tParticles->begin();
                for ( ; tChildPartIter != tParticles->end(); ++tChildPartIter) {
                    if (tChildPartIter->fathereid == tPartIter->eid) {
                        tDecayTimes[tPartIter->eid] = tChildPartIter->t;
                        break;
                    }
                }
            }
            unsigned int tFirstStep = (unsigned int)TMath::Floor((tPartIter->t - tEarliestTime)/m_StepSize);
            if (tFirstStep >= m_nSteps) {
                tFirstStep = m_nSteps - 1;
            }
            tFirstStepExists[tPartIter->eid] = tFirstStep;
            if (tPartIter->decayed) {
                unsigned int tLastStep = (unsigned int)TMath::Floor((tDecayTimes[tPartIter->eid] - tEarliestTime)/m_StepSize) - 1;
                if (tLastStep >= m_nSteps) {
                    tLastStep = m_nSteps - 1;
                }
                tLastStepExists[tPartIter->eid] = tLastStep;
            } else {
                tLastStepExists[tPartIter->eid] = m_nSteps - 1;
            }
            if (tPartIter->pid == 211) {
                if (favoriteEID == -1) {
                    favoriteEID = tPartIter->eid;
                    favoriteE0 = tPartIter->e;
                }
            }
        }

        int ** tBestIs = new int*[N];
        for (unsigned int ai = 0; ai < N; ++ai) {
            tBestIs[ai] = new int[N];
            for (unsigned int aj = 0; aj < N; ++aj) {
                tBestIs[ai][aj] = -1;
            }
        }

        TVector3 * tTotalForcesLast = NULL;
        for (unsigned int iStep = 0; iStep < m_nSteps; ++iStep) {

           // cout << "Step " << iStep << " of " << m_nSteps << " particles: " << tParticles->size() << endl;

            double tCurrentTime = tEarliestTime + iStep * m_StepSize;

            TVector3 * tTotalForces = new TVector3[N];
            std::list<Particle>::iterator tPartIter = tParticles->begin();
            for ( ; tPartIter != tParticles->end(); ++tPartIter) {
                Int_t tQ = tPartIter->GetParticleType()->GetCharge();
                TVector3 tPartPos(tPartIter->x, tPartIter->y, tPartIter->z); // 1/GeV
                if (tQ != 0 && iStep >= tFirstStepExists[tPartIter->eid] && iStep <= tLastStepExists[tPartIter->eid]) {   
                    tCoordinates[tPartIter->eid].push_back(TLorentzVector(tPartPos, tCurrentTime)); // 1/GeV
                    TVector3 tPartMom(tPartIter->px, tPartIter->py, tPartIter->pz);
                    double tP2 = tPartMom.Mag2();
                    double tM2 = TMath::Power(tPartIter->GetParticleType()->GetMass(), 2);
                    TVector3 tV(0,0,0);
                    if (tP2 > 0 && TMath::Finite(tP2)) { 
                        tV = 1. / TMath::Sqrt(tM2/tP2 + 1) * tPartMom.Unit();
                    }
                    tVelocities[tPartIter->eid].push_back(tV); // []
                    if (debug && !TMath::Finite(tV.Mag())) {
                        cout << "In step " << iStep << " ERROR, infinite velocity " << tPartMom.X() << " " << tPartMom.Y() << " " << tPartMom.Z() << endl;
                        cout << "In step " << iStep << " ERROR, infinite velocity " << tV.X() << " " << tV.Y() << " " << tV.Z() << endl;
                    }

                    if (iStep == tFirstStepExists[tPartIter->eid]) {
                        tAccelerations[tPartIter->eid].push_back(TVector3(0,0,0));
                    } else {
                        if (tFirstStepExists[tPartIter->eid] < 0) {
                            cerr << tParticles->size() << " " << N << " " << tPartIter->eid << " " << sizeof(tFirstStepExists)/sizeof(int*) << " " << tFirstStepExists[tPartIter->eid] << " " << iStep-tFirstStepExists[tPartIter->eid] << " " << iStep-tFirstStepExists[tPartIter->eid]-1 << " " << iStep << " " << tVelocities[tPartIter->eid].size() << endl;
                        }
                        TVector3 &tV1 = tVelocities[tPartIter->eid][iStep-tFirstStepExists[tPartIter->eid]-1];
                        TVector3 &tV2 = tVelocities[tPartIter->eid][iStep-tFirstStepExists[tPartIter->eid]];
                        TVector3 tAcce = 1./(m_StepSize) * (tV2 - tV1);
                        if (debug && !TMath::Finite(tAcce.Z())) {
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
            tPartIter = tParticles->begin();
            for ( ; tPartIter != tParticles->end(); ++tPartIter) {

                TVector3 tForceTotal;
                Int_t tQ = tPartIter->GetParticleType()->GetCharge();
                if (tQ != 0 && iStep >= tFirstStepExists[tPartIter->eid] && iStep <= tLastStepExists[tPartIter->eid]) {
                    tForceTotal = CalculateTotalForce(tPartIter, tParticles, tBestIs, tCurrentTime,
                            tCoordinates, tVelocities, tAccelerations
                            );
                }
                tTotalForces[tPartIter->eid] = tForceTotal;

            }
            tPartIter = tParticles->begin();
            for ( ; tPartIter != tParticles->end(); ++tPartIter) {
                TVector3 tPartPos(tPartIter->x, tPartIter->y, tPartIter->z);
                TVector3 tPartMom(tPartIter->px, tPartIter->py, tPartIter->pz);

                TVector3 tForceTotal = tTotalForces[tPartIter->eid];
                TVector3 tImpulse = (m_StepSize * kHbarC) * tForceTotal; // GeV/fm * (1/GeV) * GeV*fm = GeV
                double zeff = tForceTotal.Mag() * tPartPos.Mag() * kHbarC * tPartPos.Mag() * kHbarC / (qe*qe);

                double tP2 = tPartMom.Mag2();
                double fM = tPartIter->GetParticleType()->GetMass();
                double tM2 = TMath::Power(fM, 2);

                if (tP2 > 0 && TMath::Finite(tP2)) { // guardian
                    tPartPos += (m_StepSize / TMath::Sqrt(tM2/tP2 + 1)) * tPartMom.Unit() // 1/GeV
                        + (0.5/fM) * TMath::Power(m_StepSize * kHbarC, 2) * tForceTotal;
                    if (debug && !TMath::Finite(tPartPos.Mag())) {
                        cout << "updating position: " << tM2 << " " << tP2 << " " << tPartMom.Mag() << endl;
                    }
                }
                tPartMom += tImpulse; // GeV
                //// cout << tPartIter->pid << " " << tImpulse.X() << " " << tImpulse.Y() << " " << tImpulse.Z() << " " << tPartMom.X() << " " << tPartMom.Y() << " " << tPartMom.Z() << endl;

                double tE = TMath::Sqrt( tM2 + tP2 );
                tPartIter->x  = tPartPos.X();
                tPartIter->y  = tPartPos.Y();
                tPartIter->z  = tPartPos.Z();
                tPartIter->e  = tE;
                tPartIter->px = tPartMom.X();
                tPartIter->py = tPartMom.Y();
                tPartIter->pz = tPartMom.Z();
                }
                delete [] tTotalForcesLast;
                tTotalForcesLast = tTotalForces;
            }
            // Writing back to particles' objects original position, to not change HBT radii
            tPartIter = tParticles->begin();
            for ( ; tPartIter != tParticles->end(); ++tPartIter) {
                if (tCoordinates[tPartIter->eid].size() > 0) {
                    TVector3 tPartPos = tCoordinates[tPartIter->eid].front().Vect();
                    TVector3 tPartPosFinal(tPartIter->x, tPartIter->y, tPartIter->z);
                    tPartIter->x  = tPartPos.X();
                    tPartIter->y  = tPartPos.Y();
                    tPartIter->z  = tPartPos.Z();
                    if (tPartIter->eid == favoriteEID) {
                        //// cout << "+> " << tPartIter->e << " - " << favoriteE0 << " = " << tPartIter->e - favoriteE0 << endl;
                        //// cout << "-> " << tPartPos.Mag() << " " << tPartPosFinal.Mag() << endl;
                    }
                }
            }
            delete [] tBestIs;

            delete [] tCoordinates;
            delete [] tVelocities;
            delete [] tAccelerations;
            delete [] tDecayTimes;
            delete [] tFirstStepExists;
            delete [] tLastStepExists;
        }
