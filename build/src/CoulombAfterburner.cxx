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

TString TVector3toString(TVector3 &v) {
    return Form("[%f, %f, %f]", v.X(), v.Y(), v.Z());
}

TVector3 CoulombAfterburner::CalculateTotalForce(
        std::list<Particle>::iterator &tPartIter, std::list<Particle> *tParticles, int **tBestIs, double tCurrentTime,
        Int_t iStep)
{
    if (IsSpectator(tPartIter)) {
        return TVector3();
    }
    //quantity tQe = 2 * hep::meter;
    Int_t tQ = tPartIter->GetParticleType()->GetCharge();
    // TVector3 tForceTotal;
    TVector3 tTermCoulTotal;
    TVector3 tTermVeloTotal;
    TVector3 tTermAcceTotal;
    TVector3 tTermMagnTotal;
    TVector3 tPartPos = tPartIter->VecX(); // 1/GeV

    TVector3 tPartMom = tPartIter->VecP();
    TVector3 tV(0, 0, 0);
    double tP = tPartIter->GetP();
    double tE = tPartIter->E();
    if (tP > 0 && TMath::Finite(tP))
    {
        tV = tP / tE * tPartMom.Unit();
    }
    // double tP2 = tPartMom.Mag2();
    // double tM2 = TMath::Power(tPartIter->mass, 2);
    // if (tP2 > 0 && TMath::Finite(tP2))
    // {
    //     tV = 1. / TMath::Sqrt(tM2 / tP2 + 1) * tPartMom.Unit();
    // }
    std::list<Particle>::iterator tPartOtherIter = tParticles->begin();
    for (; tPartOtherIter != tParticles->end(); ++tPartOtherIter)
    {

        if (tPartOtherIter->eid == tPartIter->eid)
            continue;

        Int_t tOtherQ = tPartOtherIter->GetParticleType()->GetCharge();
        if (tOtherQ == 0)
            continue;

        // Skip interaction between a father and a child
        if (tPartIter->eid == tPartOtherIter->fathereid || tPartOtherIter->eid == tPartIter->fathereid) {
            continue;
        }

        // Skip interaction between decay products at the step when they were created
        if (tPartIter->fathereid != -1 && tPartIter->fathereid == tPartOtherIter->fathereid && iStep == m_FirstStepExists[tPartOtherIter->eid]) {
            continue;
        }

        TVector3 tPartOtherPos, tPartOtherVel, tPartOtherAcc;
        if (!RetardedSource(tPartIter, tPartOtherIter, tPartPos, tCurrentTime, tBestIs, iStep,
                    tPartOtherPos, tPartOtherVel, tPartOtherAcc,
                    m_FindProblem && tPartIter->eid == m_DebugEid && iStep == m_DebugStep && tPartOtherIter->eid == m_DebugOtherEid))
        {
            continue;
        }
        if (m_FindProblem && tPartIter->eid == m_DebugEid && iStep == m_DebugStep && tPartOtherIter->eid == m_DebugOtherEid) {
            cout << "v_x here = " << tPartOtherVel.X() << endl;
        }

        TVector3 tPosDiff = tPartPos - tPartOtherPos;
        TVector3 tDirection = tPosDiff.Unit();
        if (m_FindProblem && tPartIter->eid == m_DebugEid && iStep == m_DebugStep && tPartOtherIter->eid == m_DebugOtherEid) {
            cout << "v_x here = " << tPartOtherVel.X() << endl;
        }

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

        if (m_FindProblem && tPartIter->eid == m_DebugEid && iStep == m_DebugStep/* && tPartOtherIter->eid == m_DebugOtherEid*/)
       // if (!TMath::Finite(Rdotu))
        {
            cout << endl << "Size of other coords " << m_coordinates[tPartOtherIter->eid].size()
                << " other pid " << tPartOtherIter->pid << " other eid " << tPartOtherIter->eid
                << " other father eid " << tPartOtherIter->fathereid << " other father pid " << tPartOtherIter->fatherpid
                << " n parts " << tParticles->size() << " " << endl;
            cout << "Other exists since step " << m_FirstStepExists[tPartOtherIter->eid] << endl;
            cout << "first times " << m_coordinates[tPartIter->eid][0].T() << " " << m_coordinates[tPartOtherIter->eid][0].T() << endl;
            cout << "1/(R.u)**3 = " << 1/TMath::Power(Rdotu, 3) << "other Q " << tOtherQ << endl;
            cout << "R.u = " << Rdotu << " |v| = " << v.Mag() << ", R.a = " << Rdota << endl;
            cout << "R = (" << vecR.X() << "," << vecR.Y() << "," << vecR.Z() << "), u = " 
                << u.X() << "," << u.Y() << "," << u.Z() << ")" << endl;
            cout << "pos = (" << tPartPos.X() << "," << tPartPos.Y() << "," << tPartPos.Z() << "), other = (" 
                << tPartOtherPos.X() << "," << tPartOtherPos.Y() << "," << tPartOtherPos.Z() << ")" << endl;
            cout << "R/(R.u)**3 = " << R / TMath::Power(Rdotu, 3) << endl;
            cout << "R/(R.u)**3 (1-v2) = " << R / TMath::Power(Rdotu, 3) * (1-v2) << endl;
            cout << "R/(R.u)**3 (1-v2) u_x = " << R / TMath::Power(Rdotu, 3) * (1-v2) * u.X() << endl;
            cout << "kQ R/(R.u)**3 (1-v2) u_x = " << k * tOtherQ * R / TMath::Power(Rdotu, 3) * (1-v2) * u.X() << endl;
            Int_t vel_size = m_velocities[tPartOtherIter->eid].size();
            for (Int_t vel_i = 0; vel_i < vel_size; ++vel_i) {
                auto vel = m_velocities[tPartOtherIter->eid][vel_i];
                if (vel.Mag() > 1) {
                    cout << "---- v(" << vel_i << " of " << vel_size << ") = (" << vel.X() << "," << vel.Y() << "," << vel.Z() << "), |v| = " << vel.Mag() << endl;
                }
            }
           // return TVector3(0, 0, 0);
        }
        // TVector3 tTermCoul = (tOtherQ * qe * R / TMath::Power(Rdotu, 3) * (1 - v2)) * u; // sqrt(GeV*fm) / GeV * GeV**3 = sqrt(GeV*fm) * GeV**2
        // TVector3 tTermVelo = (tOtherQ * qe * R / TMath::Power(Rdotu, 3) * Rdota) * u;    // sqrt(GeV*fm) / GeV * GeV**3 = sqrt(GeV*fm) * GeV**2
        // TVector3 tTermAcce = (tOtherQ * qe * R / TMath::Power(Rdotu, 3) * (-Rdotu)) * a; // sqrt(GeV*fm) * GeV**2
        TVector3 tTermCoul = k * (tQ / TMath::Power(kHbarC, 2)) * (tOtherQ * R / TMath::Power(Rdotu, 3) * (1 - v2)) * u; // GeV*fm / GeV * GeV**3 = GeV*fm * GeV**2
        TVector3 tTermVelo = k * (tQ / TMath::Power(kHbarC, 2)) * (tOtherQ * R / TMath::Power(Rdotu, 3) * Rdota) * u;    // GeV*fm / GeV * GeV**3 = GeV*fm * GeV**2
        TVector3 tTermAcce = k * (tQ / TMath::Power(kHbarC, 2)) * (tOtherQ * R / TMath::Power(Rdotu, 3) * (-Rdotu)) * a; // GeV*fm * GeV**2

        TVector3 tE = tTermCoul + tTermVelo + tTermAcce; // GeV*fm * GeV**2
                                                         // Using BAC-CAB rule for tV x tB having c * tB = tDirection x tE
        TVector3 tTermMagn = hatR * tV.Dot(tE) - hatR.Dot(tV) * tE; // GeV*fm * GeV**2
        if (tPartIter->eid == 10 && TString(tPartOtherIter->GetParticleType()->GetName()).Contains("Spectator")) {
            //    cout << "COUL TOTAL " << tTermCoulTotal.Mag() << endl; 
        }

        // Storing all pairs at all steps would take many GB, keep only the large contributions
        if (m_FindProblem && tTermCoul.Mag() > m_LargeTerm) {
            m_termsCoul->back()[std::make_pair(tPartIter->eid, tPartOtherIter->eid)] = tTermCoul;
            m_termsVelo->back()[std::make_pair(tPartIter->eid, tPartOtherIter->eid)] = tTermVelo;
            m_termsAcce->back()[std::make_pair(tPartIter->eid, tPartOtherIter->eid)] = tTermAcce;
            m_termsMagn->back()[std::make_pair(tPartIter->eid, tPartOtherIter->eid)] = tTermMagn;
            m_Rs->back()[std::make_pair(tPartIter->eid, tPartOtherIter->eid)] = vecR;
            m_us->back()[std::make_pair(tPartIter->eid, tPartOtherIter->eid)] = u;
            m_vs->back()[std::make_pair(tPartIter->eid, tPartOtherIter->eid)] = v;
        }

        if (m_FindProblem && tPartIter->eid == m_DebugEid && iStep == m_DebugStep) {
            cout << " -- other pid = " << tPartOtherIter->pid << " other eid = " << tPartOtherIter->eid << endl;;
            cout << " -- R = " << R << endl;
            cout << " --- TermCoul = (" << tTermCoul.X() << "," << tTermCoul.Y() << "," << tTermCoul.Z() << ")" << endl;
            cout << " --- TermVelo = (" << tTermVelo.X() << "," << tTermVelo.Y() << "," << tTermVelo.Z() << ")" << endl;
            cout << " --- TermAcce = (" << tTermAcce.X() << "," << tTermAcce.Y() << "," << tTermAcce.Z() << ")" << endl;
            cout << " --- TermMagn = (" << tTermMagn.X() << "," << tTermMagn.Y() << "," << tTermMagn.Z() << ")" << endl;
            cout << " ------- " << endl;
        }
        tTermCoulTotal += tTermCoul;
        tTermVeloTotal += tTermVelo;
        tTermAcceTotal += tTermAcce;
        tTermMagnTotal += tTermMagn;
        if (m_FindProblem && tPartIter->eid == m_DebugEid && iStep == m_DebugStep) {
            cout << " --- TermCoulTotal = (" << tTermCoulTotal.X() << "," << tTermCoulTotal.Y() << "," << tTermCoulTotal.Z() << ")" << endl;
            cout << " --- TermVeloTotal = (" << tTermVeloTotal.X() << "," << tTermVeloTotal.Y() << "," << tTermVeloTotal.Z() << ")" << endl;
            cout << " --- TermAcceTotal = (" << tTermAcceTotal.X() << "," << tTermAcceTotal.Y() << "," << tTermAcceTotal.Z() << ")" << endl;
            cout << " --- TermMagnTotal = (" << tTermMagnTotal.X() << "," << tTermMagnTotal.Y() << "," << tTermMagnTotal.Z() << ")" << endl;
            cout << " ------- " << endl;
        }
        if (tPartIter->eid == 10 && TString(tPartOtherIter->GetParticleType()->GetName()).Contains("Spectator")) {
            //    cout << "SPECTATOR CHARGE " << tPartOtherIter->GetParticleType()->GetCharge() << endl;
            //    cout << "COUL TOTAL " << tTermCoulTotal.Mag() << " CONTRIB " << tTermCoul.Mag() << endl; 
        }

    }
    TVector3 tForceTotal = tTermCoulTotal + tTermVeloTotal + tTermAcceTotal + tTermMagnTotal; // GeV**3 * fm / (GeV**2 * fm) = GeV / fm

    return tForceTotal;
}

std::pair<TVector3, Double_t> CoulombAfterburner::CalculatePotentials(
        std::list<Particle>::iterator &tPartIter, std::list<Particle> *tParticles, int **tBestIs, double tCurrentTime,
        Int_t iStep)
{
    if (IsSpectator(tPartIter)) {
        return std::make_pair(TVector3(), 0);
    }
    Double_t tTermScalarTotal = 0;
    TVector3 tTermVectorTotal;
    TVector3 tPartPos = tPartIter->VecX(); // 1/GeV
    std::list<Particle>::iterator tPartOtherIter = tParticles->begin();
    for (; tPartOtherIter != tParticles->end(); ++tPartOtherIter)
    {

        if (tPartOtherIter->eid == tPartIter->eid)
            continue;

        Int_t tOtherQ = tPartOtherIter->GetParticleType()->GetCharge();
        if (tOtherQ == 0)
            continue;

        TVector3 tPartOtherPos, tPartOtherVel, tPartOtherAcc;
        if (!RetardedSource(tPartIter, tPartOtherIter, tPartPos, tCurrentTime, tBestIs, iStep,
                    tPartOtherPos, tPartOtherVel, tPartOtherAcc))
        {
            continue;
        }

        TVector3 v = tPartOtherVel;                // ()
        TVector3 vecR = tPartPos - tPartOtherPos;  // 1/GeV
        Double_t R = vecR.Mag();                   // 1/GeV

        Double_t tTermScalar = k * tOtherQ / (R - vecR.Dot(v)) / kHbarC; // GeV * fm / (1/GeV) / (GeV * fm) = GeV
        TVector3 tTermVector = tTermScalar * v; // GeV
        tTermScalarTotal += tTermScalar;
        tTermVectorTotal += tTermVector;
    }

    return std::make_pair(tTermVectorTotal, tTermScalarTotal);
}

bool CoulombAfterburner::IsSpectator(std::list<Particle>::iterator &tPartIter)
{
    Int_t tPDG = TMath::Abs(tPartIter->GetParticleType()->GetPDGCode());
    return tPDG == 72212 || tPDG == 72112;
}

bool CoulombAfterburner::HasFiniteKinematics(const Particle &aParticle)
{
    return TMath::Finite(aParticle.E()) && TMath::Finite(aParticle.Px()) && TMath::Finite(aParticle.Py())
        && TMath::Finite(aParticle.Pz()) && TMath::Finite(aParticle.mass) && TMath::Finite(aParticle.T())
        && TMath::Finite(aParticle.X()) && TMath::Finite(aParticle.Y()) && TMath::Finite(aParticle.Z());
}

// Finds the position, velocity and acceleration of tPartOtherIter at the point where its trajectory crosses
// the past light cone of the space-time point (tPos, tTime). Returns false if there is no such point,
// i.e. the two particles are not causally connected.
// tTime may be up to one time step after the last recorded point of the trajectories (force at the new position).
bool CoulombAfterburner::RetardedSource(std::list<Particle>::iterator &tPartIter, std::list<Particle>::iterator &tPartOtherIter,
        TVector3 &tPos, double tTime, int **tBestIs, Int_t iStep,
        TVector3 &tOtherPos, TVector3 &tOtherVel, TVector3 &tOtherAcc, bool debug)
{
    Int_t tOtherEid = tPartOtherIter->eid;
    std::vector<TLorentzVector> &tTrajectory = m_coordinates[tOtherEid];
    std::vector<TVector3> &tVelocities = m_velocities[tOtherEid];
    std::vector<TVector3> &tAccelerations = m_accelerations[tOtherEid];
    unsigned long tSize = tTrajectory.size();
    if (tSize == 0) {
        return false;
    }

    double tBestInterval;
    unsigned long tLastBestI = tBestIs[tPartIter->eid][tOtherEid] > -1 ? tBestIs[tPartIter->eid][tOtherEid] : 0;
    unsigned long tBestI = NearestInterval(tPos, tTime, tTrajectory, tLastBestI, tBestInterval, debug);
    tBestIs[tPartIter->eid][tOtherEid] = tBestI;
    if (debug) {
        cout << "RetardedSource: other eid " << tOtherEid << " bestI " << tBestI << " of " << tSize
            << " best interval " << tBestInterval << endl;
    }

    // Accelerations are backward differences of the velocities, a_k = (v_k - v_{k-1}) / dt, so they belong to the
    // middle of the time step, t_k - dt/2: interpolate linearly between these midpoints. a_0 is not used, it is set
    // to 0 at the first step, when there is no previous velocity.
    auto tAccelerationAt = [&](double aT) {
        if (tSize < 2) return TVector3(0, 0, 0);
        double tK = (aT - tTrajectory[0].T()) / m_StepSize + 0.5;
        if (tK <= 1) return tAccelerations[1];
        if (tK >= tSize - 1) return tAccelerations[tSize - 1];
        long k = (long)TMath::Floor(tK);
        double w = tK - k;
        return (1 - w) * tAccelerations[k] + w * tAccelerations[k + 1];
    };

    // The whole trajectory is inside the past light cone, the crossing would be after its last point: extrapolate
    // forwards from the last point. Either the trajectory is recorded only up to the current step, or the particle
    // decayed between its last time step and its decay time.
    if (tBestI + 1 == tSize && tBestInterval > 0) {
        bool tOtherExistsNow = (Int_t)m_FirstStepExists[tOtherEid] <= iStep && iStep <= (Int_t)m_LastStepExists[tOtherEid];
        TVector3 tX0 = tTrajectory[tBestI].Vect();
        TVector3 tV0 = tVelocities[tBestI];
        double tT0 = tTrajectory[tBestI].T();
        double tTimeRet = NearestPreFreezeoutIntervalTime(tX0, tV0, tPos, tT0, tTime);
        if (!TMath::Finite(tTimeRet)) {
            return false;
        }
        if (!tOtherExistsNow && tTimeRet >= m_DecayTimes[tOtherEid]) {
            // The other particle decayed before the crossing, its field is taken over by its decay products
            return false;
        }
        tOtherPos = tX0 + (tTimeRet - tT0) * tV0;
        tOtherVel = tV0;
        tOtherAcc = tAccelerationAt(tTimeRet);
        return true;
    }

    bool tIsPrimary = tPartOtherIter->pid == tPartOtherIter->fatherpid;
    // The whole trajectory of a decay product is outside the past light cone, the crossing would be before its first
    // point: extrapolate backwards to its creation, between the creation time and the first time step
    if (tBestI == 0 && tBestInterval < 0 && !tIsPrimary) {
        TVector3 tX0 = tTrajectory[0].Vect();
        TVector3 tV0 = tVelocities[0];
        double tT0 = tTrajectory[0].T();
        double tTimeRet = NearestPreFreezeoutIntervalTime(tX0, tV0, tPos, tT0, tTime);
        if (!TMath::Finite(tTimeRet) || tTimeRet < m_CreationTimes[tOtherEid]) {
            // The decay product is not created yet, as seen from (tPos, tTime)
            return false;
        }
        tOtherPos = tX0 + (tTimeRet - tT0) * tV0;
        tOtherVel = tV0;
        tOtherAcc = tAccelerationAt(tTimeRet);
        return true;
    }

    tOtherPos = tTrajectory[tBestI].Vect();
    tOtherVel = tVelocities[tBestI];
    tOtherAcc = tAccelerations[tBestI];

    if (tIsPrimary && (tBestI == 0 || tBestI == 1))
    {
        // Primary particle and the crossing is close to or before freeze-out:
        // extrapolate the trajectory backwards from freeze-out with the initial velocity
        bool tOtherIsSpectator = IsSpectator(tPartOtherIter);
        if (tOtherIsSpectator && !m_ExtrapolateSpectators) {
            return false;
        } else if (!tOtherIsSpectator && !m_ExtrapolateFireball) {
            return false;
        }
        double tZeroTime = tTrajectory[0].T();
        TVector3 tX0 = tTrajectory[0].Vect();
        TVector3 tIniV = tVelocities[0];
        double tTimePreFreeze = NearestPreFreezeoutIntervalTime(tX0, tIniV, tPos, tZeroTime, tTime);
        // No solution for extrapolated past trajectory of the other particle crossing the past light cone
        // of the current particle. They are not causaly connected
        if (!TMath::Finite(tTimePreFreeze)) {
            return false;
        }
        tOtherPos = tX0 + (tTimePreFreeze - tZeroTime) * tIniV;
        tOtherVel = tIniV;
        tOtherAcc = TVector3(0, 0, 0);
    } else {
        // The crossing is between points tZeroI-1 and tZeroI of the trajectory. A non-negative interval at tBestI
        // means tBestI is before the crossing. tZeroI >= 1, because tBestI == 0 with a negative interval was handled above.
        unsigned long tZeroI = (tBestInterval >= 0 && tBestI + 1 < tSize) ? tBestI + 1 : tBestI;
        if (tZeroI == 0) {
            return true;
        }
        double tZeroTime = tTrajectory[tZeroI].T();
        TVector3 tX0 = tTrajectory[tZeroI].Vect();
        TVector3 tZeroV = tVelocities[tZeroI];
        double tTimeBest = NearestPreFreezeoutIntervalTime(tX0, tZeroV, tPos, tZeroTime, tTime);
        if (TMath::Finite(tTimeBest)) {
            tOtherPos = tX0 + (tTimeBest - tZeroTime) * tZeroV;
            // Weight of point tZeroI: 0 at its predecessor, 1 at tZeroI itself
            double tInterpolatingWeight = (tTimeBest - tZeroTime + m_StepSize) / m_StepSize;
            tInterpolatingWeight = TMath::Max(0., TMath::Min(1., tInterpolatingWeight));
            tOtherVel =
                tZeroV * tInterpolatingWeight +
                tVelocities[tZeroI-1] * (1-tInterpolatingWeight);
            tOtherAcc = tAccelerationAt(tTimeBest);
            if (debug) {
                cout << "RetardedSource: zeroI " << tZeroI << " weight " << tInterpolatingWeight
                    << " v_x " << tOtherVel.X() << endl;
            }
        }
    }
    return true;
}


CoulombAfterburner::CoulombAfterburner(int aSteps, double aStepSize, Double_t aR, AbstractEventSaver *aEventSaver, ParticleDB *tPartDB) 
    : m_nSteps(aSteps), m_StepSize(aStepSize / kHbarC), mR(aR), mEventSaver(aEventSaver), m_hNpartImpB(nullptr), mKeepTrajectories(false),
      m_UseSpectators(true), m_ExtrapolateSpectators(true), m_ExtrapolateFireball(true),
      m_FindProblem(true), m_LargeTerm(1.0), m_DebugStep(-2), m_DebugEid(-2), m_DebugOtherEid(-2)
{
    m_mcg=new TGlauberMC("Au","Au",23.8,0);
    ParticleType *tProtonType = tPartDB->FindByPID(2212);
    ParticleType *tNeutronType = tPartDB->FindByPID(2112);
    mProtonSpectatorAType = new ParticleType(*tProtonType);
    mNeutronSpectatorAType = new ParticleType(*tNeutronType);
    mProtonSpectatorAType->SetPDGCode(72212);
    mNeutronSpectatorAType->SetPDGCode(72112);
    mProtonSpectatorBType = new ParticleType(*tProtonType);
    mNeutronSpectatorBType = new ParticleType(*tNeutronType);
    mProtonSpectatorBType->SetPDGCode(-72212);
    mNeutronSpectatorBType->SetPDGCode(-72112);
}
CoulombAfterburner::~CoulombAfterburner()
{
    delete m_mcg;
    delete mProtonSpectatorAType;
    delete mNeutronSpectatorAType;
    delete mProtonSpectatorBType;
    delete mNeutronSpectatorBType;

}
void CoulombAfterburner::ReadParameters(Configurator *aMainConfig)
{
    if (aMainConfig->HasParameter("CoulombUseSpectators")) {
        m_UseSpectators = aMainConfig->GetParameter("CoulombUseSpectators").Atoi() != 0;
    }
    if (aMainConfig->HasParameter("CoulombExtrapolateSpectators")) {
        m_ExtrapolateSpectators = aMainConfig->GetParameter("CoulombExtrapolateSpectators").Atoi() != 0;
    }
    if (aMainConfig->HasParameter("CoulombExtrapolateFireball")) {
        m_ExtrapolateFireball = aMainConfig->GetParameter("CoulombExtrapolateFireball").Atoi() != 0;
    }
    if (aMainConfig->HasParameter("CoulombFindProblem")) {
        m_FindProblem = aMainConfig->GetParameter("CoulombFindProblem").Atoi() != 0;
    }
    if (aMainConfig->HasParameter("CoulombLargeTerm")) {
        m_LargeTerm = aMainConfig->GetParameter("CoulombLargeTerm").Atof();
    }
    if (aMainConfig->HasParameter("CoulombDebugStep")) {
        m_DebugStep = aMainConfig->GetParameter("CoulombDebugStep").Atoi();
    }
    if (aMainConfig->HasParameter("CoulombDebugEid")) {
        m_DebugEid = aMainConfig->GetParameter("CoulombDebugEid").Atoi();
    }
    if (aMainConfig->HasParameter("CoulombDebugOtherEid")) {
        m_DebugOtherEid = aMainConfig->GetParameter("CoulombDebugOtherEid").Atoi();
    }
}
void CoulombAfterburner::Init() {
   // cout << "BEGIN INIT " << endl;
    m_fileHistOut = new TFile("controlHists.root","recreate");
    m_hTotalHamil = new TH1F("hTotalHamil","hTotalHamil",m_nSteps,0,m_nSteps*m_StepSize);
    m_hTotalMomMS = new TH1F("hTotalMomMS","hTotalMomMS",m_nSteps,0,m_nSteps*m_StepSize);
    m_hTotalPx = new TH1F("hTotalPx","hTotalPx",m_nSteps,0,m_nSteps*m_StepSize);
    m_hTotalPy = new TH1F("hTotalPy","hTotalPy",m_nSteps,0,m_nSteps*m_StepSize);
    m_hTotalPz = new TH1F("hTotalPz","hTotalPz",m_nSteps,0,m_nSteps*m_StepSize);
    m_hTotalAx = new TH1F("hTotalAx","hTotalAx",m_nSteps,0,m_nSteps*m_StepSize);
    m_hTotalAy = new TH1F("hTotalAy","hTotalAy",m_nSteps,0,m_nSteps*m_StepSize);
    m_hTotalAz = new TH1F("hTotalAz","hTotalAz",m_nSteps,0,m_nSteps*m_StepSize);
    m_hTotalV = new TH1F("hTotalV","hTotalV",m_nSteps,0,m_nSteps*m_StepSize);
    m_hTotalRadPower = new TH1F("hTotalRadPower","hTotalRadPower",m_nSteps,0,m_nSteps*m_StepSize);
    m_hTotalPdotA = new TH1F("hTotalPdotA","hTotalPdotA",m_nSteps,0,m_nSteps*m_StepSize);
    m_pTotalHamil = new TProfile("pTotalHamil","pTotalHamil",m_nSteps,0,m_nSteps*m_StepSize);
    m_pTotalMomMS = new TProfile("pTotalMomMS","pTotalMomMS",m_nSteps,0,m_nSteps*m_StepSize);
    m_pTotalPx = new TProfile("pTotalPx","pTotalPx",m_nSteps,0,m_nSteps*m_StepSize);
    m_pTotalPy = new TProfile("pTotalPy","pTotalPy",m_nSteps,0,m_nSteps*m_StepSize);
    m_pTotalPz = new TProfile("pTotalPz","pTotalPz",m_nSteps,0,m_nSteps*m_StepSize);
    m_pTotalAx = new TProfile("pTotalAx","pTotalAx",m_nSteps,0,m_nSteps*m_StepSize);
    m_pTotalAy = new TProfile("pTotalAy","pTotalAy",m_nSteps,0,m_nSteps*m_StepSize);
    m_pTotalAz = new TProfile("pTotalAz","pTotalAz",m_nSteps,0,m_nSteps*m_StepSize);
    m_pTotalV = new TProfile("pTotalV","pTotalV",m_nSteps,0,m_nSteps*m_StepSize);
    m_pTotalRadPower = new TProfile("pTotalRadPower","pTotalRadPower",m_nSteps,0,m_nSteps*m_StepSize);
    m_pTotalPdotA = new TProfile("pTotalPdotA","pTotalPdotA",m_nSteps,0,m_nSteps*m_StepSize);
    m_h2TotalHamil = new TH2F("h2TotalHamil","h2TotalHamil",m_nSteps,0,m_nSteps*m_StepSize,400,-20,20);
    m_h2TotalMomMS = new TH2F("h2TotalMomMS","h2TotalMomMS",m_nSteps,0,m_nSteps*m_StepSize,400,-20,20);
    m_h2TotalPx = new TH2F("h2TotalPx","h2TotalPx",m_nSteps,0,m_nSteps*m_StepSize,400,-2,2);
    m_h2TotalPy = new TH2F("h2TotalPy","h2TotalPy",m_nSteps,0,m_nSteps*m_StepSize,400,-2,2);
    m_h2TotalPz = new TH2F("h2TotalPz","h2TotalPz",m_nSteps,0,m_nSteps*m_StepSize,400,-2,2);
    m_h2TotalAx = new TH2F("h2TotalAx","h2TotalAx",m_nSteps,0,m_nSteps*m_StepSize,400,-2,2);
    m_h2TotalAy = new TH2F("h2TotalAy","h2TotalAy",m_nSteps,0,m_nSteps*m_StepSize,400,-2,2);
    m_h2TotalAz = new TH2F("h2TotalAz","h2TotalAz",m_nSteps,0,m_nSteps*m_StepSize,400,-2,2);
    m_h2TotalV = new TH2F("h2TotalV","h2TotalV",m_nSteps,0,m_nSteps*m_StepSize,400,-2,2);
    m_h2TotalRadPower = new TH2F("h2TotalRadPower","h2TotalRadPower",m_nSteps,0,m_nSteps*m_StepSize,400,-2,2);
    m_h2TotalPdotA = new TH2F("h2TotalPdotA","h2TotalPdotA",m_nSteps,0,m_nSteps*m_StepSize,400,-2,2);
   // cout << "MIDDLE INIT " << endl;

    m_termsCoul = new std::vector< std::map< std::pair<int, int>, TVector3> >();
    m_termsVelo = new std::vector< std::map< std::pair<int, int>, TVector3> >();
    m_termsAcce = new std::vector< std::map< std::pair<int, int>, TVector3> >();
    m_termsMagn = new std::vector< std::map< std::pair<int, int>, TVector3> >();
    m_Rs = new std::vector< std::map< std::pair<int, int>, TVector3> >();
    m_us = new std::vector< std::map< std::pair<int, int>, TVector3> >();
    m_vs = new std::vector< std::map< std::pair<int, int>, TVector3> >();
   // cout << "END INIT " << endl;

}
void CoulombAfterburner::Finalize() {
    cout << "Before writing hists to file" << endl;
    m_fileHistOut->cd();
    m_hTotalHamil->Write();
    m_hTotalMomMS->Write();
    m_hTotalPx->Write();
    m_hTotalPy->Write();
    m_hTotalPz->Write();
    m_hTotalAx->Write();
    m_hTotalAy->Write();
    m_hTotalAz->Write();
    m_hTotalV->Write();
    m_hTotalRadPower->Write();
    m_hTotalPdotA->Write();
    m_pTotalHamil->Write();
    m_pTotalMomMS->Write();
    m_pTotalPx->Write();
    m_pTotalPy->Write();
    m_pTotalPz->Write();
    m_pTotalAx->Write();
    m_pTotalAy->Write();
    m_pTotalAz->Write();
    m_pTotalV->Write();
    m_pTotalRadPower->Write();
    m_pTotalPdotA->Write();
    m_h2TotalHamil->Write();
    m_h2TotalMomMS->Write();
    m_h2TotalPx->Write();
    m_h2TotalPy->Write();
    m_h2TotalPz->Write();
    m_h2TotalAx->Write();
    m_h2TotalAy->Write();
    m_h2TotalAz->Write();
    m_h2TotalV->Write();
    m_h2TotalRadPower->Write();
    m_h2TotalPdotA->Write();
    m_fileHistOut->Close();
    cout << "After writing hists to file" << endl;

    delete m_termsCoul;
    delete m_termsVelo;
    delete m_termsAcce;
    delete m_termsMagn;
    delete m_Rs;
    delete m_us;
    delete m_vs;
}

TGraph *initImpBGraph() {
    Double_t xs[] { 185.870 + 185.642, 181.902 + 181.760, 184.599 + 184.364, 180.046 + 179.368, 183.434 + 182.991, 183.980 + 183.701, 184.294 + 183.154, 185.464 + 184.580, 180.894 + 180.582, 182.516 + 181.986, 178.901 + 178.906, 185.637 + 185.712, 177.916 + 177.302, 178.391 + 178.344, 177.536 + 177.389, 180.008 + 180.720, 172.162 + 172.490, 173.693 + 173.214, 175.419 + 175.353, 171.275 + 171.576, 168.162 + 167.389, 170.919 + 171.987, 169.970 + 168.991, 172.380 + 172.458, 167.680 + 167.973, 165.792 + 166.663, 159.967 + 159.881, 160.395 + 160.109, 159.053 + 161.153, 162.047 + 159.706, 155.590 + 155.196, 157.467 + 158.524, 153.694 + 150.953, 150.958 + 150.781, 151.297 + 150.622, 148.507 + 148.534, 150.576 + 151.296, 149.398 + 148.339, 140.584 + 139.497, 140.006 + 141.817, 141.453 + 141.444, 140.423 + 139.828, 132.688 + 130.706, 129.880 + 130.140, 132.745 + 132.244, 130.558 + 133.015, 128.724 + 127.144, 123.967 + 124.291, 122.225 + 122.776, 120.866 + 121.504, 121.731 + 122.749, 117.899 + 119.436, 115.054 + 116.298, 110.464 + 110.776, 109.516 + 110.150, 106.308 + 106.569, 108.114 + 107.820, 106.558 + 105.327, 105.752 + 107.747, 100.821 +  99.785, 97.137 +  98.315, 95.343 +  93.013, 94.033 +  95.004, 90.977 +  90.945, 92.298 +  92.903, 87.081 +  86.855, 87.430 +  86.606, 85.404 +  86.900, 87.092 +  86.733, 79.703 +  78.480, 77.819 +  75.356, 78.144 +  78.572, 76.003 +  76.257, 73.348 +  72.836, 73.849 +  73.083, 69.608 +  69.788, 67.071 +  67.497, 68.992 +  69.597, 65.500 +  64.320, 59.964 +  60.350, 60.009 +  59.817, 57.177 +  57.197, 55.319 +  55.347, 57.712 +  56.912, 51.972 +  52.519, 52.459 +  53.849, 49.385 +  49.521, 46.880 +  47.787, 48.240 +  46.627, 47.702 +  47.086, 42.209 +  42.466, 42.925 +  42.240, 40.645 +  39.476, 38.318 +  39.389, 38.285 +  38.486, 36.290 +  36.220, 33.978 +  34.830, 30.991 +  31.838, 31.793 +  31.406, 31.227 +  30.930, 29.945 +  29.247, 26.719 +  26.933, 27.071 +  26.197, 23.788 +  23.629, 23.166 +  23.339, 23.378 +  22.931, 21.188 +  21.000, 20.186 +  20.473, 18.152 +  17.882, 18.901 +  18.587, 16.549 +  16.860, 15.772 +  15.781, 14.534 +  15.124, 13.647 +  13.353, 12.323 +  12.683, 11.912 +  12.562, 11.230 +  10.695, 9.650 +   9.998, 9.744 +   9.642, 8.714 +   8.857, 8.015 +   7.976, 6.854 +   7.172, 7.120 +   7.135, 6.599 +   6.713, 5.862 +   6.019, 5.284 +   5.415, 4.996 +   5.038, 4.308 +   4.397, 4.053 +   3.915, 3.678 +   3.704, 3.718 +   3.672, 3.078 +   2.983, 2.691 +   2.589, 2.658 +   2.629, 2.202 +   2.177, 1.894 +   1.866, 1.885 +   1.812, 1.544 +   1.540, 1.477 +   1.411, 1.222 +   1.221, 1.073 +   1.047, 0.955 +   0.941, 0.896 +   0.842, 0.778 +   0.759, 0.683 +   0.666, 0.591 +   0.570, 0.543 +   0.498, 0.442 +   0.427, 0.363 +   0.353, 0.343 +   0.333, 0.281 +   0.268, 0.256 +   0.250, 0.213 +   0.203, 0.199 +   0.190, 0.158 +   0.155, 0.141 +   0.130, 0.118 +   0.110, 0.106 +   0.100, 0.092 +   0.087, 0.079 +   0.071, 0.063 +   0.056, 0.058 +   0.052, 0.049 +   0.044, 0.042 +   0.037, 0.035 +   0.032, 0.031 +   0.027, 0.026 +   0.023, 0.022 +   0.019, 0.018 +   0.016, 0.015 +   0.013, 0.012 +   0.011, 0.012 +   0.010, 0.010 +   0.009, 0.008 +   0.007, 0.007 +   0.006, 0.005 +   0.005, 0.005 +   0.004, 0.004 +   0.003, 0.003 +   0.003, 0.003 +   0.002, 0.002 +   0.002, 0.002 +   0.002, 0.002 +   0.001, 0.001 +   0.001, 0.001 +   0.001, 0.001 +   0.001, 0.001 +   0.001, 0.001 +   0.001, 0.001 +   0.000, 0.000 +   0.000, 0.000 +   0.000, 0.000 +   0.000, 0.000 +   0.000, 0.000 +   0.000, 0.000 +   0.000, 0.000 +   0.000, 0.000 +   0.000, 0.000 +   0.000, 0.000 +   0.000, 0.000 +   0.000 };
    Double_t ys[] { 0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05, 2.15, 2.25, 2.35, 2.45, 2.55, 2.65, 2.75, 2.85, 2.95, 3.05, 3.15, 3.25, 3.35, 3.45, 3.55, 3.65, 3.75, 3.85, 3.95, 4.05, 4.15, 4.25, 4.35, 4.45, 4.55, 4.65, 4.75, 4.85, 4.95, 5.05, 5.15, 5.25, 5.35, 5.45, 5.55, 5.65, 5.75, 5.85, 5.95, 6.05, 6.15, 6.25, 6.35, 6.45, 6.55, 6.65, 6.75, 6.85, 6.95, 7.05, 7.15, 7.25, 7.35, 7.45, 7.55, 7.65, 7.75, 7.85, 7.95, 8.05, 8.15, 8.25, 8.35, 8.45, 8.55, 8.65, 8.75, 8.85, 8.95, 9.05, 9.15, 9.25, 9.35, 9.45, 9.55, 9.65, 9.75, 9.85, 9.95, 10.05, 10.15, 10.25, 10.35, 10.45, 10.55, 10.65, 10.75, 10.85, 10.95, 11.05, 11.15, 11.25, 11.35, 11.45, 11.55, 11.65, 11.75, 11.85, 11.95, 12.05, 12.15, 12.25, 12.35, 12.45, 12.55, 12.65, 12.75, 12.85, 12.95, 13.05, 13.15, 13.25, 13.35, 13.45, 13.55, 13.65, 13.75, 13.85, 13.95, 14.05, 14.15, 14.25, 14.35, 14.45, 14.55, 14.65, 14.75, 14.85, 14.95, 15.05, 15.15, 15.25, 15.35, 15.45, 15.55, 15.65, 15.75, 15.85, 15.95, 16.05, 16.15, 16.25, 16.35, 16.45, 16.55, 16.65, 16.75, 16.85, 16.95, 17.05, 17.15, 17.25, 17.35, 17.45, 17.55, 17.65, 17.75, 17.85, 17.95, 18.05, 18.15, 18.25, 18.35, 18.45, 18.55, 18.65, 18.75, 18.85, 18.95, 19.05, 19.15, 19.25, 19.35, 19.45, 19.55, 19.65, 19.75, 19.85, 19.95 };
    return new TGraph(200, xs, ys);
}

TGraph *CoulombAfterburner::m_gImpBNPart = initImpBGraph();

bool CoulombAfterburner::ParticleExists(std::list<Particle>::iterator &tPartIter, double tDecayTime, double tTime)
{
    if (tPartIter->decayed)
    {
        return tPartIter->T() <= tTime && tDecayTime > tTime;
    }
    else
    {
        return tPartIter->T() <= tTime;
    }
}
// Along a trajectory with |v| < 1 the interval to the point (tPos, tTime) decreases monotonically: it is positive
// (time-like) before the crossing with the past light cone and negative (space-like) after it.
// Points later than tTime are in the future and are treated as space-like.
// Starting from the previous best index, find the last time-like and the first space-like point and return
// the one with smaller |interval|, with tBestInterval = sign(interval) * sqrt(|interval|).
// If the whole trajectory is space-like, returns 0 with tBestInterval < 0;
// if the whole trajectory is time-like, returns size - 1 with tBestInterval > 0.
unsigned long CoulombAfterburner::NearestInterval(TVector3 &tPos, double tTime, std::vector<TLorentzVector> &tTrajectory, unsigned long tLastBestI, double &tBestInterval, bool debug)
{
    long tSize = tTrajectory.size();
    if (tSize == 0) {
        tBestInterval = std::numeric_limits<double>::quiet_NaN();
        return 0;
    }
    auto tIsTimeLikePast = [&](long j) {
        return tTrajectory[j].T() <= tTime && this->Interval(tPos, tTime, tTrajectory[j]) >= 0;
    };
    long i = (long)tLastBestI < tSize ? (long)tLastBestI : tSize - 1;
    long tLastTimeLike; // -1 if none
    if (tIsTimeLikePast(i)) {
        tLastTimeLike = i;
        while (tLastTimeLike + 1 < tSize && tIsTimeLikePast(tLastTimeLike + 1)) {
            ++tLastTimeLike;
        }
    } else {
        tLastTimeLike = i - 1;
        while (tLastTimeLike >= 0 && !tIsTimeLikePast(tLastTimeLike)) {
            --tLastTimeLike;
        }
    }
    long tFirstSpaceLike = tLastTimeLike + 1; // tSize if none

    long tBestI;
    if (tLastTimeLike < 0) {
        tBestI = 0;
    } else if (tFirstSpaceLike >= tSize) {
        tBestI = tSize - 1;
    } else {
        double tIntervalTimeLike  = this->Interval(tPos, tTime, tTrajectory[tLastTimeLike]);
        double tIntervalSpaceLike = this->Interval(tPos, tTime, tTrajectory[tFirstSpaceLike]);
        tBestI = TMath::Abs(tIntervalTimeLike) <= TMath::Abs(tIntervalSpaceLike) ? tLastTimeLike : tFirstSpaceLike;
    }
    double tInterval = this->Interval(tPos, tTime, tTrajectory[tBestI]);
    tBestInterval = (tInterval >= 0 ? +1 : -1) * TMath::Sqrt(TMath::Abs(tInterval));
    if (debug) {
        cout << "NearestInterval: start " << tLastBestI << " size " << tSize << " last time-like " << tLastTimeLike
            << " best " << tBestI << " interval " << tBestInterval << endl;
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
    // if (!TMath::Finite(t1) || !TMath::Finite(t2)) {
    if (debug) {
        cout << "PARAMS tT0 " << tT0 << " tU " << tU << endl;
        cout << "A: " << tA << " B: " << tB << " C: " << tC << " Delta: " << tDelta << endl;
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
        if (tParticle.pid != 80000 && tParticle.pid / 10000 != 7 && tParticle.pid == tParticle.fatherpid && tParticle.GetParticleType() != nullptr)
        {
            // cout << "partial charge " << tResult << " " << tParticle.pid << " " << tParticles->size() << endl;
            tResult += tParticle.GetParticleType()->GetCharge();
        }
    }
    return tResult;
}

Double_t CoulombAfterburner::SampleImpB(std::list<Particle> *tParticles) {
    Double_t tA = 197.;
    Double_t tZ = 79.; // because they are used in fractions
    Int_t tTotalCharge = TotalCharge(tParticles);
    Int_t tNpart = tTotalCharge * tA / tZ;
    cout << "total charge " << tTotalCharge << " npart " << tNpart << " all parts " << tParticles->size() << endl;    

    if (m_hImpBprojs.count(tNpart) <= 0) {
        TString tHistName = Form("%s_%i", m_hNpartImpB->GetName(), tNpart);
        Int_t tBin = m_hNpartImpB->GetYaxis()->FindBin(tNpart);
        TH1F *tProj = (TH1F*)m_hNpartImpB->ProjectionX(tHistName, tBin, tBin);
        m_hImpBprojs[tNpart] = tProj;
        if (tProj->Integral() <= 0) {
            cout << "ERROR : integral 0 for projection in bin " << tBin << " Npart " << tNpart << endl;
        }
    }
    Double_t tImpB = m_hImpBprojs.at(tNpart)->GetRandom() / kHbarC;;

    cout << "End sampling b" << endl;
    // Double_t tImpB = m_gImpBNPart->Eval(tNpart) / kHbarC; // Distances are in 1/MeV during calculations
    return tImpB;
}

ParticleType *CoulombAfterburner::SpectatorType(std::list<Particle> *tParticles, Double_t &tSpectatorX) {

    ParticleType *tSpectatorType = new ParticleType;
    Double_t tA = 197.;
    Double_t tZ = 79.; // because they are used in fractions
    Int_t tTotalCharge = TotalCharge(tParticles);
    Int_t tSpectatorCharge = (2 * tZ - tTotalCharge)/2;
    cout << "TOTAL CHARGE " << tTotalCharge << " SPEC. CHARGE " << tSpectatorCharge << endl; 
    Int_t tNp = tSpectatorCharge;
    Int_t tNn = tNp * (tA - tZ) / tZ;
    Int_t tNpart = tTotalCharge * tA / tZ;
    Double_t tImpB = m_gImpBNPart->Eval(tNpart) / kHbarC; // Distances are in 1/MeV during calculations
    Double_t tSpectatorMass = tNp * 0.9382720 + tNn * 0.9395653;
    tSpectatorType->SetNumber(777777);
    tSpectatorType->SetName("Spectator");
    tSpectatorType->SetMass(tSpectatorMass);
    tSpectatorType->SetBarionN(tNp + tNn);
    tSpectatorType->SetI3((tNp - tNn)/2);

    // At Npart = 0, nuclei just grazing each other, spectator X would be the radius of the nucleus.
    // At Npart -> A, a head on collision, the center of mass of the spectators would also approach
    // the radius of the nucleus. so we assume that the center of mass will always be at x = tR
    Double_t tSpectatorXtoR = 1.;//tNpart / (2 * tA) + tImpB / (2 * tA) * (2 * tA - tNpart);
    Double_t tR = 1.2 * TMath::Power(tA, 1./3) / kHbarC;
    tSpectatorX = tSpectatorXtoR * tR;
    cout << "DEFINING SPECTATOR X " << tR << " x " << tSpectatorXtoR << " -> " << tSpectatorX << endl;

    return tSpectatorType;
}

Int_t CoulombAfterburner::AddSpectators(std::list<Particle> *tParticles, Double_t tEarliestTime)
{
    Int_t tNspec = 0;

    if (m_hNpartImpB == nullptr) {
        cout << "POINTER " << m_mcg << endl;
        m_mcg->SetMinDistance(0.9);
        m_mcg->SetNodeDistance(0);
        m_mcg->SetCalcLength(0);
        m_mcg->SetCalcArea(0);
        m_mcg->SetCalcCore(0);
        m_mcg->SetDetail(99);
        // m_mcg->Run(100000);
        m_mcg->Run(10000);
        TNtuple  *nt=m_mcg->GetNtuple();
        m_hNpartImpB = new TH2F("hNpartImpB","hNpartImpB",100,0,20,400,0,400);
        nt->Project("hNpartImpB","Npart:B");
    }
    Double_t tImpB = SampleImpB(tParticles) * kHbarC;
    m_mcg->SetBmin(tImpB);
    m_mcg->SetBmax(tImpB);
    m_mcg->Run(1);
    TObjArray* tNucleons=m_mcg->GetNucleons();
    Int_t nNucls=tNucleons->GetEntries();
    for (Int_t iNucl=0; iNucl<nNucls; ++iNucl) {
        TGlauNucleon *tNucl=(TGlauNucleon *)tNucleons->At(iNucl);
        Int_t tSpectatorDir = +1;
        if (tNucl->IsSpectator()) {
            ParticleType *tSpectatorType;
            if (tNucl->IsInNucleusA()) {
                if (tNucl->IsProton()) {
                    tSpectatorType = mProtonSpectatorAType;
                } else {                              
                    tSpectatorType = mNeutronSpectatorAType;
                }
            } else {
                if (tNucl->IsProton()) {
                    tSpectatorType = mProtonSpectatorBType;
                } else {                               
                    tSpectatorType = mNeutronSpectatorBType;
                }
                tSpectatorDir = -1;
            }
            // Double_t tSpectatorMass = tSpectatorType->GetMass();
            Double_t tSpectatorMass = (79 * 0.9382720 + (197-79) * 0.9395653)/197;

            Double_t tSpectatorBeta = 0.63; // TODO: make an external parameter
            Double_t tSpectatorGamma = 1/TMath::Sqrt(1 - tSpectatorBeta*tSpectatorBeta);
            Double_t tSpectatorPz = tSpectatorDir * tSpectatorGamma * tSpectatorMass * tSpectatorBeta;
            Double_t tSpectatorE = tSpectatorGamma * tSpectatorMass;
            Double_t tSpectatorShift = tSpectatorDir * tSpectatorBeta * tEarliestTime;
            Double_t tSpectatorX = tNucl->GetX() / kHbarC;
            Double_t tSpectatorY = tNucl->GetY() / kHbarC;
            Double_t tSpectatorZ = tNucl->GetZ() / kHbarC / tSpectatorGamma + tSpectatorShift; // Lorentz contraction and shift

            Particle tSpectator(tSpectatorType);
            tSpectator.SetParticlePX(tSpectatorE, 0, 0, tSpectatorPz, tEarliestTime, tSpectatorX, tSpectatorY, tSpectatorZ);
            tParticles->push_back(tSpectator);
            ++tNspec;
        }
    }
    cout << "End adding spectators" << endl;
    return tNspec;
}

void CoulombAfterburner::Apply(Event *tEvent)
{
    if (m_nSteps <= 0) return;

    cout << " KUKU 1 " << endl;

    bool debug = false;

    std::list<Particle> *tParticles = tEvent->GetParticleList();
    double tEarliestTime = std::numeric_limits<double>::quiet_NaN();
    for (auto tPartIter = tParticles->begin(); tPartIter != tParticles->end(); ++tPartIter) {
        if (!HasFiniteKinematics(*tPartIter))
            continue;
        if (std::isnan(tEarliestTime) || tPartIter->T() < tEarliestTime)
        {
            tEarliestTime = tPartIter->T();
        }
    }

    if (m_UseSpectators) {
        AddSpectators(tParticles, tEarliestTime);
    }
    // All arrays below are indexed by eid
    unsigned long N = 0;
    for (auto tPartIter = tParticles->begin(); tPartIter != tParticles->end(); ++tPartIter) {
        if ((unsigned long)tPartIter->eid + 1 > N) {
            N = tPartIter->eid + 1;
        }
    }
    m_coordinates = new std::vector<TLorentzVector>[N];
    m_fourmomenta = new std::vector<TLorentzVector>[N];
    m_velocities = new std::vector<TVector3>[N];
    m_accelerations = new std::vector<TVector3>[N];
    m_pids = new std::map<int, int>();
    m_DecayTimes = new double[N];
    m_CreationTimes = new double[N];
    m_FirstStepExists = new unsigned int[N];
    m_LastStepExists = new unsigned int[N];
    TLorentzVector *tOriginalCoordinates = new TLorentzVector[N];

    for (auto tPartIter = tParticles->begin(); tPartIter != tParticles->end(); ++tPartIter)
    {
        m_coordinates[tPartIter->eid] = std::vector<TLorentzVector>();
        tOriginalCoordinates[tPartIter->eid] = TLorentzVector(tPartIter->VecX(), tPartIter->T());
        m_CreationTimes[tPartIter->eid] = tPartIter->T();
        m_DecayTimes[tPartIter->eid] = std::numeric_limits<double>::infinity();
        if (tPartIter->decayed)
        {
            std::list<Particle>::iterator tChildPartIter = tParticles->begin();
            for (; tChildPartIter != tParticles->end(); ++tChildPartIter)
            {
                if (tChildPartIter->fathereid == tPartIter->eid)
                {
                    m_DecayTimes[tPartIter->eid] = tChildPartIter->T();
                    break;
                }
            }
        }
        // A particle exists at the time steps t with creation time <= t < decay time, so that the first step
        // of a decay product directly follows the last step of its father
        double tFirstStep = TMath::Min((double)m_nSteps, TMath::Ceil((tPartIter->T() - tEarliestTime) / m_StepSize));
        double tLastStep = m_nSteps - 1;
        if (TMath::Finite(m_DecayTimes[tPartIter->eid]))
        {
            tLastStep = TMath::Min(tLastStep, TMath::Ceil((m_DecayTimes[tPartIter->eid] - tEarliestTime) / m_StepSize) - 1);
        }
        if (!HasFiniteKinematics(*tPartIter))
        {
            // A single particle with NaN momentum or position would turn the forces on all charged particles to NaN
            cout << "<CoulombAfterburner::Apply>\tParticle eid " << tPartIter->eid << ", pid " << tPartIter->pid
                << " has a non-finite momentum, mass or position: not propagated and not a source of the field" << endl;
            tFirstStep = m_nSteps;
            tLastStep = 0;
        }
        else if (tFirstStep > tLastStep)
        {
            // The particle doesn't exist at any time step
            tFirstStep = m_nSteps;
            tLastStep = 0;
        }
        m_FirstStepExists[tPartIter->eid] = (unsigned int)tFirstStep;
        m_LastStepExists[tPartIter->eid] = (unsigned int)tLastStep;
    }
    cout << "End setting initial coordinates" << endl;
    // if (m_UseSpectators) {
    //     Int_t tNspec = AddSpectators(tParticles, tEarliestTime);
    // }
    int **tBestIs = new int *[N];
    for (unsigned int ai = 0; ai < N; ++ai)
    {
        tBestIs[ai] = new int[N];
        for (unsigned int aj = 0; aj < N; ++aj)
        {
            tBestIs[ai][aj] = -1;
        }
    }
    m_termsCoul->clear();
    m_termsVelo->clear();
    m_termsAcce->clear();
    m_termsMagn->clear();
    m_Rs->clear();
    m_us->clear();
    m_vs->clear();

    // Reference values from the first step, to monitor conservation laws
    Double_t tTotalHamil_0 = 0;
    Double_t tTotalMomMS_0 = 0;
    Double_t tTotalPx_0 = 0;
    Double_t tTotalPy_0 = 0;
    Double_t tTotalPz_0 = 0;
    Double_t tTotalAx_0 = 0;
    Double_t tTotalAy_0 = 0;
    Double_t tTotalAz_0 = 0;
    Double_t tTotalV_0 = 0;
    Double_t tTotalPdotA_0 = 0;
    Double_t tRadiatedEnergy = 0;
    Double_t tRadiatedEnergy_0 = 0;

    TVector3 *tTotalForcesLast = NULL;
    for (unsigned int iStep = 0; iStep < m_nSteps; ++iStep)
    {
       // cout << "1. Step " << iStep << " of " << m_nSteps << " particles: " << tParticles->size() << endl;

        m_termsCoul->push_back(std::map< std::pair<int, int>, TVector3>());
        m_termsVelo->push_back(std::map< std::pair<int, int>, TVector3>());
        m_termsAcce->push_back(std::map< std::pair<int, int>, TVector3>());
        m_termsMagn->push_back(std::map< std::pair<int, int>, TVector3>());
        m_Rs->push_back(std::map< std::pair<int, int>, TVector3>());
        m_us->push_back(std::map< std::pair<int, int>, TVector3>());
        m_vs->push_back(std::map< std::pair<int, int>, TVector3>());

        double tCurrentTime = tEarliestTime + iStep * m_StepSize;
        TVector3 *tTotalForces = new TVector3[N];
        std::list<Particle>::iterator tPartIter = tParticles->begin();
        for (; tPartIter != tParticles->end(); ++tPartIter)
        {
            /*
               if (TString(tPartIter->GetParticleType()->GetName()).Contains("Spectator")) {
               cout << " at step " << iStep << " spectator is ";
               cout << "x(" << tPartIter->X() << ", " << tPartIter->Y() << ", " << tPartIter->Z() << "), ";
               cout << "p(" << tPartIter->Px() << ", " << tPartIter->Py() << ", " << tPartIter->Pz() << ")" << endl;
               }
               */
            if (iStep >= m_FirstStepExists[tPartIter->eid] && iStep <= m_LastStepExists[tPartIter->eid])
            {
                if (iStep == m_FirstStepExists[tPartIter->eid])
                {
                    // Propagate freely from the creation time to the first time step
                    TVector3 tPos = tPartIter->VecX();
                    double tP = tPartIter->GetP();
                    if (tP > 0 && TMath::Finite(tP))
                    {
                        tPos += ((tCurrentTime - tPartIter->T()) * tP / tPartIter->E()) * tPartIter->VecP().Unit();
                    }
                    tPartIter->SetTXYZ(tCurrentTime, tPos.X(), tPos.Y(), tPos.Z());
                }
                Int_t tQ = tPartIter->GetParticleType()->GetCharge();
                TVector3 tPartPos = tPartIter->VecX(); // 1/GeV
                TVector3 tPartMom = tPartIter->VecP();
                double tE = tPartIter->E();
                m_coordinates[tPartIter->eid].push_back(TLorentzVector(tPartPos, tCurrentTime)); // 1/GeV
                m_fourmomenta[tPartIter->eid].push_back(TLorentzVector(tPartMom, tE)); // 1/GeV
                (*m_pids)[tPartIter->eid] = tPartIter->pid;
                if (tQ != 0)
                {
                    // cout << m_coordinates->size() << " " << tPartIter->eid << " " << tPartIter->pid << " " << tParticles->size() << endl;
                    // cout << m_coordinates[tPartIter->eid].size() << endl;
                    // cout << "===" << endl;
                    // m_coordinates[tPartIter->eid].push_back(TLorentzVector(tPartPos, tCurrentTime)); // 1/GeV
                    // double tP2 = tPartMom.Mag2();
                    // double tM2 = TMath::Power(tPartIter->mass, 2);
                    TVector3 tV(0, 0, 0);
                    // if (tP2 > 0 && TMath::Finite(tP2))
                    // {
                    //     tV = 1. / TMath::Sqrt(tM2 / tP2 + 1) * tPartMom.Unit();
                    // }
                    double tP = tPartIter->GetP();
                    if (tP > 0 && TMath::Finite(tP))
                    {
                        tV = tP / tE * tPartMom.Unit();
                    }
                    m_velocities[tPartIter->eid].push_back(tV); // []

                    if (iStep == m_FirstStepExists[tPartIter->eid])
                    {
                        m_accelerations[tPartIter->eid].push_back(TVector3(0, 0, 0));
                    }
                    else
                    {
                        if (m_FirstStepExists[tPartIter->eid] < 0)
                        {
                            cerr << tParticles->size() << " " << N << " " << tPartIter->eid << " " << sizeof(m_FirstStepExists) / sizeof(int *) << " " << m_FirstStepExists[tPartIter->eid] << " " << iStep - m_FirstStepExists[tPartIter->eid] << " " << iStep - m_FirstStepExists[tPartIter->eid] - 1 << " " << iStep << " " << m_velocities[tPartIter->eid].size() << endl;
                        }
                        TVector3 &tV1 = m_velocities[tPartIter->eid][iStep - m_FirstStepExists[tPartIter->eid] - 1];
                        TVector3 &tV2 = m_velocities[tPartIter->eid][iStep - m_FirstStepExists[tPartIter->eid]];
                        TVector3 tAcce = 1. / (m_StepSize) * (tV2 - tV1);
                        m_accelerations[tPartIter->eid].push_back(tAcce); // GeV
                    }
                }
            }
        }
        // Calculate force at old position
        tPartIter = tParticles->begin();
        for (; tPartIter != tParticles->end(); ++tPartIter)
        {

            TVector3 tForceTotal;
            // For spectator nucleons, the force stays zero
            // Reuse the force at the new position from the previous step, unless the particle has just been created
            Int_t tQ = tPartIter->GetParticleType()->GetCharge();
            if (tQ != 0 && iStep >= m_FirstStepExists[tPartIter->eid] && iStep <= m_LastStepExists[tPartIter->eid])
            {
                if (tTotalForcesLast != NULL && iStep != m_FirstStepExists[tPartIter->eid]) {
                    tForceTotal = tTotalForcesLast[tPartIter->eid];
                } else {
                    tForceTotal = CalculateTotalForce(tPartIter, tParticles, tBestIs, tCurrentTime, iStep);
                }
            }
            if (m_FindProblem && tPartIter->eid == m_DebugEid) {
                cout << "total force at old pos, step " << iStep << ", eid = " << tPartIter->eid
                    << ", F = (" << tForceTotal.X() << "," << tForceTotal.Y() << "," << tForceTotal.Z() << ") "  << endl;
            }
            tTotalForces[tPartIter->eid] = tForceTotal;
        }
        // Update positions x_{n+1} = x_n + v_n dt + a_n dt^2 / 2, using the force calculated at old position
        tPartIter = tParticles->begin();
        for (; tPartIter != tParticles->end(); ++tPartIter)
        {
            if (iStep < m_FirstStepExists[tPartIter->eid] || iStep > m_LastStepExists[tPartIter->eid])
                continue;

            TVector3 tPartPos = tPartIter->VecX();
            TVector3 tPartMom = tPartIter->VecP();
            TVector3 tForceTotal = tTotalForces[tPartIter->eid];

            TVector3 tV(0, 0, 0);
            double tP = tPartIter->GetP();
            double tE = tPartIter->E();
            if (tP > 0 && TMath::Finite(tP))
            {                                                                           // guardian
                tV = tP / tE * tPartMom.Unit();
            }
            // Relativistic acceleration dv/dt = (F - v (v.F)) / E
            TVector3 tAcc = (kHbarC / tE) * (tForceTotal - tV.Dot(tForceTotal) * tV); // GeV*fm * GeV/fm / GeV = GeV
            tPartPos += m_StepSize * tV + (0.5 * m_StepSize * m_StepSize) * tAcc; // 1/GeV
            if (debug && !TMath::Finite(tPartPos.Mag()))
            {
                cout << "updating position: " << tE << " " << tP << " " << tPartMom.Mag() << endl;
            }

            tPartIter->SetTXYZ(tCurrentTime + m_StepSize, tPartPos.X(), tPartPos.Y(), tPartPos.Z());
        }

        // Calculate force at new position x_{n+1} and time t_{n+1}
        TVector3 *tTotalForcesNew = new TVector3[N];
        tPartIter = tParticles->begin();
        for (; tPartIter != tParticles->end(); ++tPartIter)
        {

            TVector3 tForceTotal;
            Int_t tQ = tPartIter->GetParticleType()->GetCharge();
            if (tQ != 0 && iStep >= m_FirstStepExists[tPartIter->eid] && iStep <= m_LastStepExists[tPartIter->eid])
            {
                tForceTotal = CalculateTotalForce(tPartIter, tParticles, tBestIs, tCurrentTime + m_StepSize, iStep);
            }
            if (m_FindProblem && tPartIter->eid == m_DebugEid) {
                cout << "total force at new pos, step " << iStep << ", eid = " << tPartIter->eid
                    << ", F = (" << tForceTotal.X() << "," << tForceTotal.Y() << "," << tForceTotal.Z() << ") "  << endl;
            }
            tTotalForcesNew[tPartIter->eid] = tForceTotal;
        }
        // Update momenta p_{n+1} = p_n + (F_n + F_{n+1}) dt / 2, using average force from old and new position
        tPartIter = tParticles->begin();
        for (; tPartIter != tParticles->end(); ++tPartIter)
        {
            if (iStep < m_FirstStepExists[tPartIter->eid] || iStep > m_LastStepExists[tPartIter->eid])
                continue;

            TVector3 tPartPos = tPartIter->VecX();
            TVector3 tPartMom = tPartIter->VecP();
            TVector3 tForceTotal = 0.5 * (tTotalForces[tPartIter->eid] + tTotalForcesNew[tPartIter->eid]);
            TVector3 tImpulse = (m_StepSize * kHbarC) * tForceTotal; // GeV/fm * (1/GeV) * GeV*fm = GeV
            if (m_FindProblem && tPartIter->eid == m_DebugEid) {
                cout << "total force " << iStep << ", eid = " << tPartIter->eid
                    << ", F = (" << tForceTotal.X() << "," << tForceTotal.Y() << "," << tForceTotal.Z() << ") "  << endl;
            }

            double fM = tPartIter->mass;
            double tM2 = fM*fM;

            tPartMom += tImpulse; // GeV

            double tP2 = tPartMom.Mag2();
            double tE = TMath::Sqrt(tM2 + tP2);
            if (m_FindProblem && tPartIter->eid == m_DebugEid) {
                cout << "before updating momenta " << iStep << ", eid = " << tPartIter->eid
                    << ", p = (" << tPartIter->Px() << "," << tPartIter->Py() << "," << tPartIter->Pz() << "), |p| = " 
                    << tPartIter->GetP() << ", E = " << tPartIter->E() << ", mass = " << tPartIter->mass << endl;
            }

            tPartIter->SetEPxPyPz(tE, tPartMom.X(), tPartMom.Y(), tPartMom.Z());
            if (m_FindProblem && tPartIter->eid == m_DebugEid) {
                cout << "after updating momenta " << iStep << ", eid = " << tPartIter->eid
                    << ", p = (" << tPartIter->Px() << "," << tPartIter->Py() << "," << tPartIter->Pz() << "), |p| = " 
                    << tPartIter->GetP() << ", E = " << tPartIter->E() << ", mass = " << tPartIter->mass << endl;
            }

            TVector3 tVnew = tPartMom.Mag() / tE * tPartMom.Unit();
            if (m_FindProblem && (!TMath::Finite(tVnew.Mag()) || tVnew.Mag() >= 1))
            {
                cout << "In step " << iStep << " ERROR, wrong velocity eid " << tPartIter->eid << ", pid " << tPartIter->pid << endl;
                cout << "In step " << iStep << " p = " << tPartMom.X() << " " << tPartMom.Y() << " " << tPartMom.Z() << endl;
                cout << "In step " << iStep << " v = " << tVnew.X() << " " << tVnew.Y() << " " << tVnew.Z() << endl;
                cout << "In step " << iStep << " |p| = " << tPartMom.Mag() << ", E = " << tE << endl;
                for (auto mapping : m_termsCoul->back()) {
                    if (mapping.first.first == tPartIter->eid && mapping.second.Mag() > m_LargeTerm) {
                        cout << "possible reason: " << mapping.second.X() << " " << mapping.second.Y() << " " << mapping.second.Z() << endl;
                        cout << "other eid: " << mapping.first.second << endl;
                    } 
                }
            }
        }
        // The force at the new position is the force at the old position in the next step
        if (tTotalForcesLast != NULL) {
            delete [] tTotalForcesLast;
        }
        tTotalForcesLast = tTotalForcesNew;

        // Calculating total energy, momentum and potentials at time t_{n+1}
        Double_t tTotalHamil = 0;
        Double_t tTotalMomMS = 0;
        Double_t tTotalPx = 0;
        Double_t tTotalPy = 0;
        Double_t tTotalPz = 0;
        Double_t tTotalAx = 0;
        Double_t tTotalAy = 0;
        Double_t tTotalAz = 0;
        Double_t tTotalV = 0;
        Double_t tTotalPdotA = 0;
        Double_t tTotalRadPower = 0;

        tPartIter = tParticles->begin();
        for (; tPartIter != tParticles->end(); ++tPartIter)
        {
            if (iStep < m_FirstStepExists[tPartIter->eid] || iStep > m_LastStepExists[tPartIter->eid])
                continue;

            Int_t tQ = tPartIter->GetParticleType()->GetCharge();
            TVector3 tPartMom = tPartIter->VecP(); // kinetic momentum

            auto tPotentials = CalculatePotentials(tPartIter, tParticles, tBestIs, tCurrentTime + m_StepSize, iStep);
            TVector3 tPotA = tPotentials.first;
            Double_t tPotV = tPotentials.second;
            double fM = tPartIter->mass;
            double tM2 = fM*fM;

            Double_t tEkin = TMath::Sqrt(tM2 + tPartMom.Mag2()); // GeV
            // Each pair contributes to the potentials of both particles, hence 1/2 in the potential energy
            Double_t tH = tEkin + 0.5*tQ*tPotV; // GeV
            tTotalHamil += tH;
            tTotalMomMS += tEkin;
            tTotalPx += tPartIter->Px();
            tTotalPy += tPartIter->Py();
            tTotalPz += tPartIter->Pz();
            tTotalAx += tQ*tPotA.X();
            tTotalAy += tQ*tPotA.Y();
            tTotalAz += tQ*tPotA.Z();
            tTotalV  += 0.5*tQ*tPotV;
            tTotalPdotA  += 2*tQ*tPartMom.Dot(tPotA);

            if (tQ != 0 && m_velocities[tPartIter->eid].size() > 0 && m_accelerations[tPartIter->eid].size() > 0) {
                TVector3 tV = m_velocities[tPartIter->eid].back();
                double tBeta2 = tV.Mag2();
                double tGamma2 = 1/(1-tBeta2);
                TVector3 tA = m_accelerations[tPartIter->eid].back();
                double tRadPower = (2./3) * k * tQ * tQ * TMath::Power(tGamma2, 3) * (tA.Mag2() - tV.Cross(tA).Mag2()); // GeV*fm * GeV**2
                if (tRadPower < 0) {
                    PrintProblem(tPartIter->eid, iStep, "");
                }
                tTotalRadPower += tRadPower;
            }
            if (m_FindProblem && tPartIter->eid == m_DebugEid) {
                cout << "at end of step " << iStep << ", eid = " << tPartIter->eid
                    << ", p = (" << tPartIter->Px() << "," << tPartIter->Py() << "," << tPartIter->Pz() << "), |p| = "
                    << tPartIter->GetP() << ", E = " << tPartIter->E() << ", mass = " << tPartIter->mass << endl;
            }
        }
        tRadiatedEnergy += tTotalRadPower / kHbarC * m_StepSize; // GeV*fm * GeV**2 / (GeV*fm) / GeV = GeV

        if (iStep == 0) {
            tTotalHamil_0 = tTotalHamil;
            tTotalMomMS_0 = tTotalMomMS;
            tTotalPx_0 = tTotalPx;
            tTotalPy_0 = tTotalPy;
            tTotalPz_0 = tTotalPz;
            tTotalAx_0 = tTotalAx;
            tTotalAy_0 = tTotalAy;
            tTotalAz_0 = tTotalAz;
            tTotalV_0 = tTotalV;
            tTotalPdotA_0 = tTotalPdotA;
            tRadiatedEnergy_0 = tRadiatedEnergy;
        } else {
            m_hTotalHamil->SetBinContent(iStep+1, tTotalHamil);
            m_hTotalMomMS->SetBinContent(iStep+1, tTotalMomMS);
            m_hTotalPx->SetBinContent(iStep+1, tTotalPx);
            m_hTotalPy->SetBinContent(iStep+1, tTotalPy);
            m_hTotalPz->SetBinContent(iStep+1, tTotalPz);
            m_hTotalAx->SetBinContent(iStep+1, tTotalAx);
            m_hTotalAy->SetBinContent(iStep+1, tTotalAy);
            m_hTotalAz->SetBinContent(iStep+1, tTotalAz);
            m_hTotalV->SetBinContent(iStep+1, tTotalV);
            // RadPower histograms contain the energy radiated since the first step
            m_hTotalRadPower->SetBinContent(iStep+1, tRadiatedEnergy - tRadiatedEnergy_0);
            m_hTotalPdotA->SetBinContent(iStep+1, tTotalPdotA);
            m_pTotalHamil->Fill((iStep+0.1)*m_StepSize, tTotalHamil-tTotalHamil_0);
            m_pTotalMomMS->Fill((iStep+0.1)*m_StepSize, tTotalMomMS-tTotalMomMS_0);
            m_pTotalPx->Fill((iStep+0.1)*m_StepSize, tTotalPx-tTotalPx_0);
            m_pTotalPy->Fill((iStep+0.1)*m_StepSize, tTotalPy-tTotalPy_0);
            m_pTotalPz->Fill((iStep+0.1)*m_StepSize, tTotalPz-tTotalPz_0);
            m_pTotalAx->Fill((iStep+0.1)*m_StepSize, tTotalAx-tTotalAx_0);
            m_pTotalAy->Fill((iStep+0.1)*m_StepSize, tTotalAy-tTotalAy_0);
            m_pTotalAz->Fill((iStep+0.1)*m_StepSize, tTotalAz-tTotalAz_0);
            m_pTotalV->Fill((iStep+0.1)*m_StepSize, tTotalV-tTotalV_0);
            m_pTotalRadPower->Fill((iStep+0.1)*m_StepSize, tRadiatedEnergy - tRadiatedEnergy_0);
            m_pTotalPdotA->Fill((iStep+0.1)*m_StepSize, tTotalPdotA-tTotalPdotA_0);
            m_h2TotalHamil->Fill((iStep+0.1)*m_StepSize, tTotalHamil-tTotalHamil_0);
            m_h2TotalMomMS->Fill((iStep+0.1)*m_StepSize, tTotalMomMS-tTotalMomMS_0);
            m_h2TotalPx->Fill((iStep+0.1)*m_StepSize, tTotalPx-tTotalPx_0);
            m_h2TotalPy->Fill((iStep+0.1)*m_StepSize, tTotalPy-tTotalPy_0);
            m_h2TotalPz->Fill((iStep+0.1)*m_StepSize, tTotalPz-tTotalPz_0);
            m_h2TotalAx->Fill((iStep+0.1)*m_StepSize, tTotalAx-tTotalAx_0);
            m_h2TotalAy->Fill((iStep+0.1)*m_StepSize, tTotalAy-tTotalAy_0);
            m_h2TotalAz->Fill((iStep+0.1)*m_StepSize, tTotalAz-tTotalAz_0);
            m_h2TotalV->Fill((iStep+0.1)*m_StepSize, tTotalV-tTotalV_0);
            m_h2TotalRadPower->Fill((iStep+0.1)*m_StepSize, tRadiatedEnergy - tRadiatedEnergy_0);
            m_h2TotalPdotA->Fill((iStep+0.1)*m_StepSize, tTotalPdotA-tTotalPdotA_0);
        }
        if (tTotalForces != NULL) {
            delete [] tTotalForces;
        }
    }
    if (tTotalForcesLast != NULL) {
        delete [] tTotalForcesLast;
    }
    // Writing back to particles' objects original position, to not change HBT radii
    for (auto tPartIter = tParticles->begin(); tPartIter != tParticles->end(); ++tPartIter)
    {
        TLorentzVector &tOrig = tOriginalCoordinates[tPartIter->eid];
        tPartIter->SetTXYZ(tOrig.T(), tOrig.X(), tOrig.Y(), tOrig.Z());
    }
    for (unsigned long ai = 0; ai < N; ++ai)
    {
        delete[] tBestIs[ai];
    }
    delete[] tBestIs;

    if (mKeepTrajectories) {
        mKeptCoordinates.assign(m_coordinates, m_coordinates + N);
        mKeptFourMomenta.assign(m_fourmomenta, m_fourmomenta + N);
    }
    delete[] m_coordinates;
    delete[] m_fourmomenta;
    delete[] m_velocities;
    delete[] m_accelerations;
    delete[] m_DecayTimes;
    delete[] m_CreationTimes;
    delete[] m_FirstStepExists;
    delete[] m_LastStepExists;
    delete[] tOriginalCoordinates;
    delete m_pids;

    //tParticles->remove_if([](Particle &p) { return TString(p.GetParticleType()->GetName()).Contains("Spectator"); });
    //delete tSpectatorType;
    cout << "Before saving event" << endl;
    if (mEventSaver != nullptr) {
        mEventSaver->Save(tEvent,NULL,tEvent->GetEventID());
    }
    cout << "After saving event" << endl;


}

void CoulombAfterburner::PrintProblem(Int_t eid, Int_t step, TString indent, Int_t depth) {
    // Avoid endless recursion when two particles point at each other
    if (depth > 3) {
        return;
    }
    Int_t firstStep = m_FirstStepExists[eid];
    Int_t index = step - firstStep;
    if (index < 0 || index >= (Int_t)m_velocities[eid].size() || index >= (Int_t)m_fourmomenta[eid].size()) {
        cout << indent.Data() << "problem in step " << step << ", eid = " << eid << ": particle does not exist at this step" << endl;
        return;
    }
    TVector3 tV = m_velocities[eid].at(index);
    double tBeta2 = tV.Mag2();
    double tGamma2 = 1/(1-tBeta2);
    TVector3 tA = m_accelerations[eid].at(index);
    TLorentzVector tPartFourMom = m_fourmomenta[eid].at(index);
    cout << indent.Data() << "problem in step " << step << ", pid = " << (*m_pids)[eid] << ", eid = " << eid
        << ", gamma^2 = " << tGamma2
        << ", |a| = " << tA.Mag() << ", |v| = (" << tV.Mag() << ", a^2 - (v x a)^2 = "
        << tA.Mag2() - tV.Cross(tA).Mag2() << endl;
    cout << indent.Data() << "|p| = " << tPartFourMom.Vect().Mag() << ", E = " << tPartFourMom.E() << endl;
    for (int is = 0; is < m_termsCoul->size(); ++is) {
        for (auto mapping : m_termsCoul->at(is)) {
            if (mapping.first.first == eid && mapping.second.Mag() > m_LargeTerm) {
                cout << indent.Data() << "possible reason: " << mapping.second.X() << " " << mapping.second.Y() << " " << mapping.second.Z() << endl;
                cout << indent.Data() << "step number " << is << ", other eid: " << mapping.first.second << endl;
                cout << indent.Data() << "R: " << m_Rs->at(is)[std::make_pair(mapping.first.first, mapping.first.second)].Mag() << endl;
                cout << indent.Data() << "u: " << m_us->at(is)[std::make_pair(mapping.first.first, mapping.first.second)].Mag() << endl;
                cout << indent.Data() << "v: " << m_vs->at(is)[std::make_pair(mapping.first.first, mapping.first.second)].Mag() << endl;
                PrintProblem(mapping.first.second, is, indent + "**** ", depth + 1);
            }
        }
    }

}
