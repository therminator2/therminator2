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

CoulombAfterburner::CoulombAfterburner(int aSteps, double aStepSize) : m_nSteps(aSteps), m_StepSize(aStepSize / kHbarC) {
/*
    m_fileHistOut = new TFile("controlHists.root","recreate");
    m_hDistPresentRetarded = new TProfile("hDistPresentRetarded","hDistPresentRetarded",400,0,400);
    m_hImpulseTimeStep = new TProfile("hImpulseTimeStep","hImpulseTimeStep",4000,0,4000);
    m_hWorkTimeStep = new TProfile("hWorkTimeStep","hWorkTimeStep",4000,0,4000);
    m_hProtonsTimeStep = new TProfile("hProtonsTimeStep","hProtonsTimeStep",4000,0,4000);
    m_hProtonsInsideTimeStep = new TProfile("hProtonsInsideTimeStep","hProtonsInsideTimeStep",4000,0,4000);
    m_hInteractionsTimeStep = new TProfile("hInteractionsTimeStep","hInteractionsTimeStep",4000,0,4000);
    m_hElseTimeStep = new TProfile("hElseTimeStep","hElseTimeStep",4000,0,4000);
    m_hElseBestInegTimeStep = new TProfile("hElseBestInegTimeStep","hElseBestInegTimeStep",4000,0,4000);
    m_hElseBestIltSizeTimeStep = new TProfile("hElseBestIltSizeTimeStep","hElseBestIltSizeTimeStep",4000,0,4000);
    m_hContinueTimeStep = new TProfile("hContinueTimeStep","hContinueTimeStep",4000,0,4000);
    m_hNetChargeTimeStep = new TProfile("hNetChargeTimeStep","hNetChargeTimeStep",4000,0,4000);
    m_hZeffTimeStep = new TProfile("hZeffTimeStep","hZeffTimeStep",4000,0,4000);
    m_hZeffR = new TProfile("hZeffR","hZeffR",400,0,400);
    m_gPositionTimeStep_1 = new TGraph(aSteps);
    m_gPositionTimeStep_1->SetName("gPositionTimeStep_1");
    m_gPositionTimeStep_2 = new TGraph(aSteps);
    m_gPositionTimeStep_2->SetName("gPositionTimeStep_2");
*/
}
CoulombAfterburner::~CoulombAfterburner() {
/*
    m_fileHistOut->cd();
    m_hDistPresentRetarded->Write();
    m_hImpulseTimeStep->Write();
    m_hWorkTimeStep->Write();
    m_hProtonsTimeStep->Write();
    m_hProtonsInsideTimeStep->Write();
    m_hInteractionsTimeStep->Write();
    m_hElseTimeStep->Write();
    m_hElseBestInegTimeStep->Write();
    m_hElseBestIltSizeTimeStep->Write();
    m_hContinueTimeStep->Write();
    m_hNetChargeTimeStep->Write();
    m_hZeffTimeStep->Write();
    m_hZeffR->Write();
    m_gPositionTimeStep_1->Write();
    m_gPositionTimeStep_2->Write();
    m_fileHistOut->Close();
*/
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
/*
        int tMissingCharge = 2*79 - tPrimaryCharge;
        int tMissingBarionN = int(1.0 * tMissingCharge * 197 / 79);
        double tMissingI3 = tMissingCharge - tMissingBarionN / 2;
        double tMp = 9.382e-01;
        double tMn = 9.395e-01;
        double tMissingM = tMp * tMissingCharge + tMn * (tMissingBarionN - tMissingCharge);

        ParticleType *tSpecType = new ParticleType;
        double tSpecM = tMissingM / 2;
        tSpecType->SetMass(tSpecM);
        tSpecType->SetBarionN(tMissingBarionN / 2);
        tSpecType->SetI3(tMissingI3 / 2);
        double tSpecV = 0.63;
        double tSpecGamma = 1/TMath::Sqrt(1-tSpecV*tSpecV);

        Configurator* tModelParam;
        Parser*       tParser;

        tModelParam = new Configurator;
        tParser     = new Parser(sModelINI.Data());                                                                                                                                  
        tParser->ReadINI(tModelParam);
        delete tParser;

        double tSpecT;
        try {
            tSpecT          = tModelParam->GetParameter("T0").Atof() / kHbarC;
        }
        catch (TString &tError) {
            PRINT_MESSAGE("<Model_SR::ReadParameters>\tCaught exception " << tError);
            PRINT_MESSAGE("\tDid not find one of the necessary model parameters.");
            exit(_ERROR_CONFIG_PARAMETER_NOT_FOUND_);
        }
        delete tModelParam;

        double tSpecE = tSpecGamma * tSpecM;
        double tSpecP = TMath::Sqrt(tSpecE*tSpecE - tSpecM*tSpecM);
        Particle *tSpec1 = new Particle(tSpecType);
        Particle *tSpec2 = new Particle(tSpecType);
        tSpec1->SetParticlePX(tSpecE, 0, 0, +tSpecP, tSpecT, 0, 0, +tSpecT*tSpecV);
        tSpec2->SetParticlePX(tSpecE, 0, 0, -tSpecP, tSpecT, 0, 0, -tSpecT*tSpecV);
        tEvent->GetParticleList()->push_back(*tSpec1);
        tEvent->GetParticleList()->push_back(*tSpec2);
*/
        static const double qe = TMath::Sqrt(1.4399764e-3); // sqrt(GeV * fm)

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
/*
            TVector3 tPartPos(tPartIter->x, tPartIter->y, tPartIter->z); // 1/GeV
            TVector3 tPartMom(tPartIter->px, tPartIter->py, tPartIter->pz);
            double tP2 = tPartMom.Mag2();
            double tM2 = TMath::Power(tPartIter->GetParticleType()->GetMass(), 2);
            double tDeltaE = 0.016;
            double tDeltaP = TMath::Sqrt(TMath::Power(tPartIter->e + tDeltaE, 2) + tM2);
            TVector3 tImpulse = tDeltaP*tPartPos.Unit()*tPartIter->GetParticleType()->GetCharge();
            tPartMom += tImpulse; // GeV

            double tE = TMath::Sqrt( tM2 + tP2 );
            tPartIter->e  = tE;
            tPartIter->px = tPartMom.X();
            tPartIter->py = tPartMom.Y();
            tPartIter->pz = tPartMom.Z();
*/
        }

        int ** tBestIs = new int*[N];
        for (unsigned int ai = 0; ai < N; ++ai) {
            tBestIs[ai] = new int[N];
            for (unsigned int aj = 0; aj < N; ++aj) {
                tBestIs[ai][aj] = -1;
            }
        }

        for (unsigned int iStep = 0; iStep < m_nSteps; ++iStep) {

           // cout << "Step " << iStep << " of " << m_nSteps << " particles: " << tParticles->size() << endl;

            double tCurrentTime = tEarliestTime + iStep * m_StepSize;

            TVector3 * tTotalForces = new TVector3[N];
            TVector3 * tTotalTermsCoul = new TVector3[N];
            TVector3 * tTotalTermsVelo = new TVector3[N];
            TVector3 * tTotalTermsAcce = new TVector3[N];
            TVector3 * tTotalTermsMagn = new TVector3[N];
            int * nProtons = new int[N];
            int * nProtonsInside = new int[N];
            int * nInteractions = new int[N];
            int * nElses = new int[N];
            int * nContinues = new int[N];
            double * netCharges = new double[N];
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

                Int_t tQ = tPartIter->GetParticleType()->GetCharge();
                TVector3 tForceTotal;
                TVector3 tTermCoulTotal;
                TVector3 tTermVeloTotal;
                TVector3 tTermAcceTotal;
                TVector3 tTermMagnTotal;
                TVector3 tPartPos(tPartIter->x, tPartIter->y, tPartIter->z); // 1/GeV

                int nProton = 0;
                int nProtonInside = 0;
                int n = 0;
                int nElse = 0;
                int nElseBestIneg = 0;
                int nElseBestIltSize = 0;
                int nContinue = 0;
                double netCharge = 0.;
                if (tQ != 0 && iStep >= tFirstStepExists[tPartIter->eid] && iStep <= tLastStepExists[tPartIter->eid]) {   

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
                        if (false) {
                            if (tFirstStepExists[tPartOtherIter->eid] <= iStep && tLastStepExists[tPartOtherIter->eid] > iStep) {
                                TVector3 tPartOtherPos(tPartOtherIter->x, tPartOtherIter->y, tPartOtherIter->z); // 1/GeV
                                TVector3 tPosDiff = tPartPos - tPartOtherPos;
                                TVector3 vecR = tPosDiff;   // 1/GeV
                                Double_t R = vecR.Mag();;   // 1/GeV
                                TVector3 tTermCoul = (tOtherQ * qe / TMath::Power(R, 3)) * vecR; // sqrt(GeV*fm) / GeV * GeV**3 = sqrt(GeV*fm) * GeV**2
                                tTermCoulTotal += tTermCoul;
                                ++n;
                                netCharge += tOtherQ;
                                if (tPartOtherIter->pid == 2212) {
                                    ++nProton;
                                    if (tPartOtherPos.Mag() < tPartPos.Mag()) {
                                        ++nProtonInside;
                                    }
                                }
                            }
                        } else {

                            double tBestInterval;
                            int tLastBestI = tBestIs[tPartIter->eid][tPartOtherIter->eid] > -1 ? tBestIs[tPartIter->eid][tPartOtherIter->eid] : 0;
                            if (debug && tPartIter->eid == favoriteEID) {
                                cout << "Before calling NearestInterval, tPartPos.Mag() = " << tPartPos.Mag() << endl;
                            }
                            int tBestI = NearestInterval(tPartPos, tCurrentTime, tCoordinates[tPartOtherIter->eid], tLastBestI, tBestInterval, false && debug && tPartIter->eid == favoriteEID);
                            tBestIs[tPartIter->eid][tPartOtherIter->eid] = tBestI;

                            // The best interval is at the end of trajectory of other particle and it's space-like: 
                            // the particles are not causally connected
                            if ((tBestI + 1) == tCoordinates[tPartOtherIter->eid].size() && tBestInterval < 0) {
                                ++nContinue;
                                continue;
                            }

                            if (tBestI >= 0 && tBestI < tCoordinates[tPartOtherIter->eid].size()) {
                                if (debug && (tBestI - 2) >= 0 && (tBestI + 2) < tCoordinates[tPartOtherIter->eid].size()) {
                                    double tInterval = TMath::Power(tCurrentTime - tCoordinates[tPartOtherIter->eid][tBestI].T(), 2) - (tPartPos  - tCoordinates[tPartOtherIter->eid][tBestI].Vect()).Mag2();
                                    double tIntervalDown2 = TMath::Power(tCurrentTime - tCoordinates[tPartOtherIter->eid][tBestI-2].T(), 2) - (tPartPos  - tCoordinates[tPartOtherIter->eid][tBestI-2].Vect()).Mag2();
                                    double tIntervalDown = TMath::Power(tCurrentTime - tCoordinates[tPartOtherIter->eid][tBestI-1].T(), 2) - (tPartPos  - tCoordinates[tPartOtherIter->eid][tBestI-1].Vect()).Mag2();
                                    double tIntervalUp = TMath::Power(tCurrentTime - tCoordinates[tPartOtherIter->eid][tBestI+1].T(), 2) - (tPartPos  - tCoordinates[tPartOtherIter->eid][tBestI+1].Vect()).Mag2();
                                    double tIntervalUp2 = TMath::Power(tCurrentTime - tCoordinates[tPartOtherIter->eid][tBestI+2].T(), 2) - (tPartPos  - tCoordinates[tPartOtherIter->eid][tBestI+2].Vect()).Mag2();
                                    cout << "In step " << iStep << " tBestI = " << tBestI << " " << tBestInterval << " "
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

                            if (debug && !TMath::Finite(Rdotu)) {
                                cout << endl << "Step " << iStep << ", Condition: " << (tBestInterval < m_StepSize) << " " << tBestInterval << " " << tBestI << " " << tCoordinates[tPartOtherIter->eid].size() << " " << tPartOtherIter->pid << " " << tPartOtherIter->eid << " " << tParticles->size() << " " << tEvent->GetEventID() << endl;
                                cout << Rdotu << " " << " " << v.X() << " " << v.Y() << " " << v.Z() << endl;
                                cout << " " << vecR.X() << " " << vecR.Y() << " " << vecR.Z() << " " << u.X() << " " << u.Y() << " " << u.Z() << endl;
                                cout << " " << tPartPos.X() << " " << tPartPos.Y() << " " << tPartPos.Z() << " " << tPartOtherPos.X() << " " << tPartOtherPos.Y() << " " << tPartOtherPos.Z() << endl;
                                return;
                            }
                            TVector3 tTermCoul = (tOtherQ * qe * R / TMath::Power(Rdotu, 3) * (1-v2)) * u; // sqrt(GeV*fm) / GeV * GeV**3 = sqrt(GeV*fm) * GeV**2
                            TVector3 tTermVelo = (tOtherQ * qe * R / TMath::Power(Rdotu, 3) * Rdota) * u;  // sqrt(GeV*fm) / GeV * GeV**3 = sqrt(GeV*fm) * GeV**2
                            TVector3 tTermAcce = (tOtherQ * qe * R / TMath::Power(Rdotu, 3) * (-Rdotu)) * a; // sqrt(GeV*fm) * GeV**2
                            if (debug && (tBestI - 2) >= 0 && (tBestI + 2) < tCoordinates[tPartOtherIter->eid].size()) {
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
                            ++n;
                            netCharge += tOtherQ;
                            if (tPartOtherIter->pid == 2212) {
                                ++nProton;
                                if (tPartOtherPos.Mag() < tPartPos.Mag()) {
                                    ++nProtonInside;
                                }
                            }
                            } else {
                                ++nElse;
                                if (!(tBestI >= 0 )) {
                                    ++nElseBestIneg;
                                }
                                if (tPartIter->eid == favoriteEID && !(tBestI < tCoordinates[tPartOtherIter->eid].size())) {
                                    // cout << iStep << " " << favoriteEID << " " << tPartOtherIter->eid << " " << tPartOtherIter->pid << " " << tPartOtherIter->fatherpid << " " << tBestI << " " << tCoordinates[tPartOtherIter->eid].size() << endl;
                                    ++nElseBestIltSize;
                                }
                            }
                        }
                    }
                }
                nProtons[tPartIter->eid] = nProton;
                nProtonsInside[tPartIter->eid] = nProtonInside;
                nInteractions[tPartIter->eid] = n;
                nElses[tPartIter->eid] = nElse;
                nContinues[tPartIter->eid] = nContinue;
                netCharges[tPartIter->eid] = netCharge;
                tTotalForces[tPartIter->eid] = tForceTotal;
                tTotalTermsCoul[tPartIter->eid] = tTermCoulTotal;
                tTotalTermsVelo[tPartIter->eid] = tTermVeloTotal;
                tTotalTermsAcce[tPartIter->eid] = tTermAcceTotal;
                tTotalTermsMagn[tPartIter->eid] = tTermMagnTotal;
                if (tPartIter->pid == 211 && tPartIter->fatherpid == 211) {
                    /**
                    m_hProtonsTimeStep->Fill(iStep,nProton);
                    m_hProtonsInsideTimeStep->Fill(iStep,nProtonInside);
                    m_hInteractionsTimeStep->Fill(iStep,n);
                    m_hElseTimeStep->Fill(iStep,nElse);
                    m_hElseBestInegTimeStep->Fill(iStep,nElseBestIneg);
                    m_hElseBestIltSizeTimeStep->Fill(iStep,nElseBestIltSize);
                    m_hContinueTimeStep->Fill(iStep,nContinue);
                    m_hNetChargeTimeStep->Fill(iStep,netCharge);
                    **/
                }
                /**
                if (tPartIter->GetParticleType()->GetCharge() >= 0) {
                    m_gPositionTimeStep_1->SetPoint(iStep,tPartPos.Z()*kHbarC,iStep*m_StepSize*kHbarC);
                } else {
                    m_gPositionTimeStep_2->SetPoint(iStep,tPartPos.Z()*kHbarC,iStep*m_StepSize*kHbarC);
                }
                **/
            }
            tPartIter = tParticles->begin();
            for ( ; tPartIter != tParticles->end(); ++tPartIter) {
                TVector3 tPartPos(tPartIter->x, tPartIter->y, tPartIter->z);
                TVector3 tPartMom(tPartIter->px, tPartIter->py, tPartIter->pz);
                TVector3 tTermCoulTotal = tTotalTermsCoul[tPartIter->eid];
                TVector3 tTermVeloTotal = tTotalTermsVelo[tPartIter->eid];
                TVector3 tTermAcceTotal = tTotalTermsAcce[tPartIter->eid];
                TVector3 tTermMagnTotal = tTotalTermsMagn[tPartIter->eid];
                Int_t tQ = tPartIter->GetParticleType()->GetCharge();
                TVector3 tForceTotal = (tQ * qe / TMath::Power(kHbarC, 2)) * (tTermCoulTotal + tTermVeloTotal + tTermAcceTotal + tTermMagnTotal); // GeV*fm / ((GeV*fm)**2) * GeV**2 = GeV / fm
                /**
                if (tPartIter->pid == 211 && tPartIter->fatherpid == 211) {
                    m_hImpulseTimeStep->Fill(iStep,(m_StepSize * kHbarC) * tForceTotal.Mag());
                }
                **/
                // TVector3 tForceTotal = tTotalForces[tPartIter->eid];
                TVector3 tImpulse = (m_StepSize * kHbarC) * tForceTotal; // GeV/fm * (1/GeV) * GeV*fm = GeV
                double zeff = tForceTotal.Mag() * tPartPos.Mag() * kHbarC * tPartPos.Mag() * kHbarC / (qe*qe);
                /**
                if (tPartIter->pid == 211 && tPartIter->fatherpid==211) {
                    m_hZeffTimeStep->Fill(iStep,zeff);
                    m_hZeffR->Fill(tPartPos.Mag()*kHbarC,zeff);
                }
                **/

                if (debug && tPartIter->eid == favoriteEID) {
                    // if (tPartIter->eid == favoriteEID) {
                    cout << tPartPos.Mag() * kHbarC << " " << tPartMom.Mag() << " " << tQ * qe * tTermCoulTotal.Mag() << " " << tQ * qe * tTermVeloTotal.Mag() << " " << tQ * qe * tTermAcceTotal.Mag() << " " << tQ * qe * tTermMagnTotal.Mag() << " " << tForceTotal.Mag() << " " << nInteractions[tPartIter->eid] << " " << zeff << " " << netCharges[tPartIter->eid] << endl;
                }
                double tP2 = tPartMom.Mag2();
                double tM2 = TMath::Power(tPartIter->GetParticleType()->GetMass(), 2);
                if (tP2 > 0 && TMath::Finite(tP2)) { // guardian
                    tPartPos += (m_StepSize / TMath::Sqrt(tM2/tP2 + 1)) * tPartMom.Unit(); // 1/GeV
                    if (debug && !TMath::Finite(tPartPos.Mag())) {
                        cout << "updating position: " << tM2 << " " << tP2 << " " << tPartMom.Mag() << endl;
                    }
                }
                tPartMom += tImpulse; // GeV
                //// cout << tPartIter->pid << " " << tImpulse.X() << " " << tImpulse.Y() << " " << tImpulse.Z() << " " << tPartMom.X() << " " << tPartMom.Y() << " " << tPartMom.Z() << endl;

                double tE = TMath::Sqrt( tM2 + tP2 );
                /**
                if (tPartIter->pid == 211 && tPartIter->fatherpid == 211) {
                    m_hWorkTimeStep->Fill(iStep,tE-tPartIter->e);
                }
                **/
                tPartIter->x  = tPartPos.X();
                tPartIter->y  = tPartPos.Y();
                tPartIter->z  = tPartPos.Z();
                tPartIter->e  = tE;
                tPartIter->px = tPartMom.X();
                tPartIter->py = tPartMom.Y();
                tPartIter->pz = tPartMom.Z();
                }
                delete [] tTotalTermsCoul;
                delete [] tTotalTermsVelo;
                delete [] tTotalTermsAcce;
                delete [] tTotalTermsMagn;
                delete [] tTotalForces;
                delete [] nElses;
                delete [] nProtonsInside;
                delete [] nProtons;
                delete [] nInteractions;
                delete [] netCharges;
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
