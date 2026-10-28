#include "CoulombAfterburner.h"
#include "Configurator.h"
#include "Parser.h"
#include "TVector3.h"
#include "TMath.h"
#include <functional>
#include <limits>
#include <map>
#include <iostream>
#include "THGlobal.h"

using namespace std;
using namespace CoulombUnits;
// CoulombUnits::Energy is always written in full: THERMINATOR has a class ::Energy
extern TString sModelINI;

bool CoulombUnits::CheckInternalUnits()
{
    return TMath::Abs(Length(1. * inverse_gev_length).numerical_value_in(fm) - kHbarC) < 1e-15;
}

Area CoulombAfterburner::Interval(const Vec3<Length> &tX1, Time tT1, const Vec3<Length> &tX2, Time tT2, bool debug) {
    if (debug) {
        cout << tT1 << " " << tT2 << " " << pow<2>((tT1 - tT2) * c) << " " << Mag2(tX1 - tX2) << " " << (tX1 - tX2) << endl;
    }
    return pow<2>((tT1 - tT2) * c) - Mag2(tX1 - tX2);
}
Area CoulombAfterburner::Interval(const Vec3<Length> &tX1, Time tT1, const SpaceTimePoint &tPoint, bool debug) {
    return Interval(tX1, tT1, tPoint.x, tPoint.t, debug);
}

std::vector<TLorentzVector> CoulombAfterburner::GetTrajectory(Int_t aEid) const
{
    std::vector<TLorentzVector> tResult;
    for (const SpaceTimePoint &tPoint : mKeptCoordinates.at(aEid)) {
        tResult.push_back(TLorentzVector(tPoint.x.x.numerical_value_in(inverse_gev_length), tPoint.x.y.numerical_value_in(inverse_gev_length),
                                         tPoint.x.z.numerical_value_in(inverse_gev_length), tPoint.t.numerical_value_in(inverse_gev_time)));
    }
    return tResult;
}

std::vector<TLorentzVector> CoulombAfterburner::GetFourMomenta(Int_t aEid) const
{
    std::vector<TLorentzVector> tResult;
    for (const FourMomentum &tMomentum : mKeptFourMomenta.at(aEid)) {
        tResult.push_back(TLorentzVector(tMomentum.p.x.numerical_value_in(GeV / c), tMomentum.p.y.numerical_value_in(GeV / c),
                                         tMomentum.p.z.numerical_value_in(GeV / c), tMomentum.e.numerical_value_in(GeV)));
    }
    return tResult;
}

Vec3<Force> CoulombAfterburner::CalculateTotalForce(
        std::list<Particle>::iterator &tPartIter, std::list<Particle> *tParticles, int **tBestIs, Time tCurrentTime,
        Int_t iStep, Vec3<Force> *aMagneticForce)
{
    if (aMagneticForce != nullptr) {
        *aMagneticForce = Vec3<Force>();
    }
    if (IsSpectator(tPartIter)) {
        return Vec3<Force>();
    }
    Int_t tQ = tPartIter->GetParticleType()->GetCharge();
    Vec3<Force> tTermCoulTotal;
    Vec3<Force> tTermVeloTotal;
    Vec3<Force> tTermAcceTotal;
    Vec3<Force> tTermMagnTotal;
    Vec3<Length> tPartPos = PositionOf(*tPartIter);
    Vec3<Dimensionless> tBeta = VelocityOf(*tPartIter) / c;

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
        if (AreSiblings(*tPartIter, *tPartOtherIter) && iStep == (Int_t)m_FirstStepExists[tPartOtherIter->eid]) {
            continue;
        }

        // The force within a close pair is integrated with sub-steps, see IntegrateCloseGroup()
        if (IsClosePair(tPartIter->eid, tPartOtherIter->eid)) {
            continue;
        }

        Vec3<Length> tPartOtherPos;
        Vec3<Speed> tPartOtherVel;
        Vec3<Acceleration> tPartOtherAcc;
        if (!RetardedSource(tPartIter, tPartOtherIter, tPartPos, tCurrentTime, tBestIs, iStep,
                    tPartOtherPos, tPartOtherVel, tPartOtherAcc,
                    m_FindProblem && tPartIter->eid == m_DebugEid && iStep == m_DebugStep && tPartOtherIter->eid == m_DebugOtherEid))
        {
            continue;
        }
        if (m_FindProblem && tPartIter->eid == m_DebugEid && iStep == m_DebugStep && tPartOtherIter->eid == m_DebugOtherEid) {
            cout << "v_x here = " << tPartOtherVel.x << endl;
        }

        // Lienard-Wiechert field of the other particle at its retarded position, velocity and acceleration
        Vec3<Length> vecR = tPartPos - tPartOtherPos;
        Length R = Mag(vecR);
        Vec3<Dimensionless> hatR = UnitVector(vecR);
        Vec3<Dimensionless> v = tPartOtherVel / c;
        Vec3<Dimensionless> u = hatR - v;
        Length Rdotu = Dot(vecR, u);
        Dimensionless Rdota = Dot(vecR, tPartOtherAcc) / square(c);
        Dimensionless v2 = Mag2(v);

        // Short-distance cut-off between decay products of the same father. The distance is the one in the rest frame
        // of the source, gamma (R - R.v/c), the same for both of them.
        if (AreSiblings(*tPartIter, *tPartOtherIter) && Rdotu / sqrt(1 - v2) < m_SiblingMinDistance) {
            continue;
        }

        if (m_FindProblem && tPartIter->eid == m_DebugEid && iStep == m_DebugStep/* && tPartOtherIter->eid == m_DebugOtherEid*/)
        {
            cout << endl << "Size of other coords " << m_coordinates[tPartOtherIter->eid].size()
                << " other pid " << tPartOtherIter->pid << " other eid " << tPartOtherIter->eid
                << " other father eid " << tPartOtherIter->fathereid << " other father pid " << tPartOtherIter->fatherpid
                << " n parts " << tParticles->size() << " " << endl;
            cout << "Other exists since step " << m_FirstStepExists[tPartOtherIter->eid] << endl;
            cout << "first times " << m_coordinates[tPartIter->eid][0].t << " " << m_coordinates[tPartOtherIter->eid][0].t << endl;
            cout << "1/(R.u)**3 = " << 1. / pow<3>(Rdotu) << " other Q " << tOtherQ << endl;
            cout << "R.u = " << Rdotu << " |v| = " << Mag(v) << ", R.a/c^2 = " << Rdota << endl;
            cout << "R = " << vecR << ", u = " << u << endl;
            cout << "pos = " << tPartPos << ", other = " << tPartOtherPos << endl;
            cout << "R/(R.u)**3 = " << R / pow<3>(Rdotu) << endl;
            cout << "R/(R.u)**3 (1-v2) = " << R / pow<3>(Rdotu) * (1 - v2) << endl;
            cout << "R/(R.u)**3 (1-v2) u_x = " << R / pow<3>(Rdotu) * (1 - v2) * u.x << endl;
            cout << "kQ R/(R.u)**3 (1-v2) u_x = " << kCoulomb * tOtherQ * R / pow<3>(Rdotu) * (1 - v2) * u.x << endl;
            Int_t vel_size = m_velocities[tPartOtherIter->eid].size();
            for (Int_t vel_i = 0; vel_i < vel_size; ++vel_i) {
                const Vec3<Speed> &vel = m_velocities[tPartOtherIter->eid][vel_i];
                if (Mag(vel) > 1. * c) {
                    cout << "---- v(" << vel_i << " of " << vel_size << ") = " << vel << ", |v| = " << Mag(vel) << endl;
                }
            }
        }
        // E = k q' R / (R.u)^3 [(1 - v^2) u + (R.a) u / c^2 - (R.u) a / c^2], force q E + q v x B with c B = hatR x E
        Force tPrefactor = kCoulomb * tQ * tOtherQ * R / pow<3>(Rdotu);
        Vec3<Force> tTermCoul = tPrefactor * (1 - v2) * u;
        Vec3<Force> tTermVelo = tPrefactor * Rdota * u;
        Vec3<Force> tTermAcce = tPrefactor * (-Rdotu) * tPartOtherAcc / square(c);
        Vec3<Force> tE = tTermCoul + tTermVelo + tTermAcce;
        // Using BAC-CAB rule for beta x (hatR x E)
        Vec3<Force> tTermMagn = hatR * Dot(tBeta, tE) - Dot(hatR, tBeta) * tE;

        // Storing all pairs at all steps would take many GB, keep only the large contributions
        if (m_FindProblem && Mag(tTermCoul) > m_LargeTerm) {
            m_termsCoul->back()[std::make_pair(tPartIter->eid, tPartOtherIter->eid)] = tTermCoul;
            m_termsVelo->back()[std::make_pair(tPartIter->eid, tPartOtherIter->eid)] = tTermVelo;
            m_termsAcce->back()[std::make_pair(tPartIter->eid, tPartOtherIter->eid)] = tTermAcce;
            m_termsMagn->back()[std::make_pair(tPartIter->eid, tPartOtherIter->eid)] = tTermMagn;
            m_Rs->back()[std::make_pair(tPartIter->eid, tPartOtherIter->eid)] = vecR;
            m_us->back()[std::make_pair(tPartIter->eid, tPartOtherIter->eid)] = u;
            m_vs->back()[std::make_pair(tPartIter->eid, tPartOtherIter->eid)] = tPartOtherVel;
        }

        if (m_FindProblem && tPartIter->eid == m_DebugEid && iStep == m_DebugStep) {
            cout << " -- other pid = " << tPartOtherIter->pid << " other eid = " << tPartOtherIter->eid << endl;;
            cout << " -- R = " << R << endl;
            cout << " --- TermCoul = " << tTermCoul << endl;
            cout << " --- TermVelo = " << tTermVelo << endl;
            cout << " --- TermAcce = " << tTermAcce << endl;
            cout << " --- TermMagn = " << tTermMagn << endl;
            cout << " ------- " << endl;
        }
        tTermCoulTotal += tTermCoul;
        tTermVeloTotal += tTermVelo;
        tTermAcceTotal += tTermAcce;
        tTermMagnTotal += tTermMagn;
        if (aMagneticForce != nullptr) {
            *aMagneticForce += Cross(hatR, tE);
        }
        if (m_FindProblem && tPartIter->eid == m_DebugEid && iStep == m_DebugStep) {
            cout << " --- TermCoulTotal = " << tTermCoulTotal << endl;
            cout << " --- TermVeloTotal = " << tTermVeloTotal << endl;
            cout << " --- TermAcceTotal = " << tTermAcceTotal << endl;
            cout << " --- TermMagnTotal = " << tTermMagnTotal << endl;
            cout << " ------- " << endl;
        }
    }
    return tTermCoulTotal + tTermVeloTotal + tTermAcceTotal + tTermMagnTotal;
}

std::pair<Vec3<CoulombUnits::Energy>, CoulombUnits::Energy> CoulombAfterburner::CalculatePotentials(
        std::list<Particle>::iterator &tPartIter, std::list<Particle> *tParticles, int **tBestIs, Time tCurrentTime,
        Int_t iStep)
{
    if (IsSpectator(tPartIter)) {
        return std::make_pair(Vec3<CoulombUnits::Energy>(), CoulombUnits::Energy::zero());
    }
    CoulombUnits::Energy tTermScalarTotal = CoulombUnits::Energy::zero();
    Vec3<CoulombUnits::Energy> tTermVectorTotal;
    Vec3<Length> tPartPos = PositionOf(*tPartIter);
    std::list<Particle>::iterator tPartOtherIter = tParticles->begin();
    for (; tPartOtherIter != tParticles->end(); ++tPartOtherIter)
    {

        if (tPartOtherIter->eid == tPartIter->eid)
            continue;

        Int_t tOtherQ = tPartOtherIter->GetParticleType()->GetCharge();
        if (tOtherQ == 0)
            continue;

        Vec3<Length> tPartOtherPos;
        Vec3<Speed> tPartOtherVel;
        Vec3<Acceleration> tPartOtherAcc;
        if (!RetardedSource(tPartIter, tPartOtherIter, tPartPos, tCurrentTime, tBestIs, iStep,
                    tPartOtherPos, tPartOtherVel, tPartOtherAcc))
        {
            continue;
        }

        Vec3<Dimensionless> v = tPartOtherVel / c;
        Vec3<Length> vecR = tPartPos - tPartOtherPos;
        Length R = Mag(vecR);

        // Lienard-Wiechert potentials of a unit charge: phi = k q' / (R - R.v/c), A c = phi v / c
        CoulombUnits::Energy tTermScalar = kCoulomb * tOtherQ / (R - Dot(vecR, v));
        Vec3<CoulombUnits::Energy> tTermVector = tTermScalar * v;
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
        const Vec3<Length> &tPos, Time tTime, int **tBestIs, Int_t iStep,
        Vec3<Length> &tOtherPos, Vec3<Speed> &tOtherVel, Vec3<Acceleration> &tOtherAcc, bool debug)
{
    Int_t tOtherEid = tPartOtherIter->eid;
    std::vector<SpaceTimePoint> &tTrajectory = m_coordinates[tOtherEid];
    std::vector<Vec3<Speed> > &tVelocities = m_velocities[tOtherEid];
    std::vector<Vec3<Acceleration> > &tAccelerations = m_accelerations[tOtherEid];
    unsigned long tSize = tTrajectory.size();
    if (tSize == 0) {
        return false;
    }

    Length tBestInterval;
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
    auto tAccelerationAt = [&](Time aT) {
        if (tSize < 2) return Vec3<Acceleration>();
        double tK = ((aT - tTrajectory[0].t) / m_StepSize).numerical_value_in(one) + 0.5;
        if (tK <= 1) return tAccelerations[1];
        if (tK >= tSize - 1) return tAccelerations[tSize - 1];
        long k = (long)TMath::Floor(tK);
        double w = tK - k;
        return Vec3<Acceleration>((1 - w) * tAccelerations[k] + w * tAccelerations[k + 1]);
    };

    // The whole trajectory is inside the past light cone, the crossing would be after its last point: extrapolate
    // forwards from the last point. Either the trajectory is recorded only up to the current step, or the particle
    // decayed between its last time step and its decay time.
    if (tBestI + 1 == tSize && tBestInterval > Length::zero()) {
        bool tOtherExistsNow = (Int_t)m_FirstStepExists[tOtherEid] <= iStep && iStep <= (Int_t)m_LastStepExists[tOtherEid];
        const Vec3<Length> &tX0 = tTrajectory[tBestI].x;
        const Vec3<Speed> &tV0 = tVelocities[tBestI];
        Time tT0 = tTrajectory[tBestI].t;
        Time tTimeRet = NearestPreFreezeoutIntervalTime(tX0, tV0, tPos, tT0, tTime);
        if (!IsFinite(tTimeRet)) {
            return false;
        }
        if (!tOtherExistsNow && tTimeRet >= m_DecayTimes[tOtherEid]) {
            // The other particle decayed before the crossing, its field is taken over by its decay products
            return false;
        }
        tOtherPos = tX0 + tV0 * (tTimeRet - tT0);
        tOtherVel = tV0;
        tOtherAcc = tAccelerationAt(tTimeRet);
        return true;
    }

    bool tIsPrimary = tPartOtherIter->pid == tPartOtherIter->fatherpid;
    // The whole trajectory of a decay product is outside the past light cone, the crossing would be before its first
    // point: extrapolate backwards to its creation, between the creation time and the first time step
    if (tBestI == 0 && tBestInterval < Length::zero() && !tIsPrimary) {
        const Vec3<Length> &tX0 = tTrajectory[0].x;
        const Vec3<Speed> &tV0 = tVelocities[0];
        Time tT0 = tTrajectory[0].t;
        Time tTimeRet = NearestPreFreezeoutIntervalTime(tX0, tV0, tPos, tT0, tTime);
        if (!IsFinite(tTimeRet) || tTimeRet < m_CreationTimes[tOtherEid]) {
            // The decay product is not created yet, as seen from (tPos, tTime)
            return false;
        }
        tOtherPos = tX0 + tV0 * (tTimeRet - tT0);
        tOtherVel = tV0;
        tOtherAcc = tAccelerationAt(tTimeRet);
        return true;
    }

    tOtherPos = tTrajectory[tBestI].x;
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
        Time tZeroTime = tTrajectory[0].t;
        const Vec3<Length> &tX0 = tTrajectory[0].x;
        const Vec3<Speed> &tIniV = tVelocities[0];
        Time tTimePreFreeze = NearestPreFreezeoutIntervalTime(tX0, tIniV, tPos, tZeroTime, tTime);
        // No solution for extrapolated past trajectory of the other particle crossing the past light cone
        // of the current particle. They are not causaly connected
        if (!IsFinite(tTimePreFreeze)) {
            return false;
        }
        tOtherPos = tX0 + tIniV * (tTimePreFreeze - tZeroTime);
        tOtherVel = tIniV;
        tOtherAcc = Vec3<Acceleration>();
    } else {
        // The crossing is between points tZeroI-1 and tZeroI of the trajectory. A non-negative interval at tBestI
        // means tBestI is before the crossing. tZeroI >= 1, because tBestI == 0 with a negative interval was handled above.
        unsigned long tZeroI = (tBestInterval >= Length::zero() && tBestI + 1 < tSize) ? tBestI + 1 : tBestI;
        if (tZeroI == 0) {
            return true;
        }
        Time tZeroTime = tTrajectory[tZeroI].t;
        const Vec3<Length> &tX0 = tTrajectory[tZeroI].x;
        const Vec3<Speed> &tZeroV = tVelocities[tZeroI];
        Time tTimeBest = NearestPreFreezeoutIntervalTime(tX0, tZeroV, tPos, tZeroTime, tTime);
        if (IsFinite(tTimeBest)) {
            tOtherPos = tX0 + tZeroV * (tTimeBest - tZeroTime);
            // Weight of point tZeroI: 0 at its predecessor, 1 at tZeroI itself
            double tInterpolatingWeight = ((tTimeBest - tZeroTime + m_StepSize) / m_StepSize).numerical_value_in(one);
            tInterpolatingWeight = TMath::Max(0., TMath::Min(1., tInterpolatingWeight));
            tOtherVel = tZeroV * tInterpolatingWeight + tVelocities[tZeroI-1] * (1-tInterpolatingWeight);
            tOtherAcc = tAccelerationAt(tTimeBest);
            if (debug) {
                cout << "RetardedSource: zeroI " << tZeroI << " weight " << tInterpolatingWeight
                    << " v_x " << tOtherVel.x << endl;
            }
        }
    }
    return true;
}


CoulombAfterburner::CoulombAfterburner(int aSteps, double aStepSize, Double_t aR, AbstractEventSaver *aEventSaver, ParticleDB *tPartDB)
    : m_nSteps(aSteps), m_StepSize(aStepSize * fm_c), mR(aR), mEventSaver(aEventSaver), m_hNpartImpB(nullptr), mKeepTrajectories(false),
      m_UseSpectators(true), m_ExtrapolateSpectators(true), m_ExtrapolateFireball(true),
      m_FindProblem(true), m_LargeTerm(1.0 * GeV / fm), m_DebugStep(-2), m_DebugEid(-2), m_DebugOtherEid(-2),
      m_Substeps(false), m_SubstepAccuracy(0.05), m_SubstepMax(1024), m_SubstepDistance(1.0 * fm),
      m_SiblingMinDistance(1.0 * fm)
{
    if (!CheckInternalUnits()) {
        cout << "<CoulombAfterburner::CoulombAfterburner>\tCoulombUnits::inverse_gev_length is not kHbarC fm" << endl;
        exit(1);
    }
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
        m_LargeTerm = aMainConfig->GetParameter("CoulombLargeTerm").Atof() * GeV / fm;
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
    if (aMainConfig->HasParameter("CoulombSubsteps")) {
        m_Substeps = aMainConfig->GetParameter("CoulombSubsteps").Atoi() != 0;
    }
    if (aMainConfig->HasParameter("CoulombSubstepAccuracy")) {
        m_SubstepAccuracy = aMainConfig->GetParameter("CoulombSubstepAccuracy").Atof();
    }
    if (aMainConfig->HasParameter("CoulombSubstepMax")) {
        m_SubstepMax = TMath::Max(1, aMainConfig->GetParameter("CoulombSubstepMax").Atoi());
    }
    if (aMainConfig->HasParameter("CoulombSubstepDistance")) {
        m_SubstepDistance = aMainConfig->GetParameter("CoulombSubstepDistance").Atof() * fm;
    }
    if (aMainConfig->HasParameter("CoulombSiblingMinDistance")) {
        m_SiblingMinDistance = aMainConfig->GetParameter("CoulombSiblingMinDistance").Atof() * fm;
    }
}
void CoulombAfterburner::Init() {
    // The time axes of the control histograms are in GeV^-1
    const double tTimeRange = (m_nSteps * m_StepSize).numerical_value_in(inverse_gev_time);
    m_fileHistOut = new TFile("controlHists.root","recreate");
    m_hTotalHamil = new TH1F("hTotalHamil","hTotalHamil",m_nSteps,0,tTimeRange);
    m_hTotalMomMS = new TH1F("hTotalMomMS","hTotalMomMS",m_nSteps,0,tTimeRange);
    m_hTotalPx = new TH1F("hTotalPx","hTotalPx",m_nSteps,0,tTimeRange);
    m_hTotalPy = new TH1F("hTotalPy","hTotalPy",m_nSteps,0,tTimeRange);
    m_hTotalPz = new TH1F("hTotalPz","hTotalPz",m_nSteps,0,tTimeRange);
    m_hTotalAx = new TH1F("hTotalAx","hTotalAx",m_nSteps,0,tTimeRange);
    m_hTotalAy = new TH1F("hTotalAy","hTotalAy",m_nSteps,0,tTimeRange);
    m_hTotalAz = new TH1F("hTotalAz","hTotalAz",m_nSteps,0,tTimeRange);
    m_hTotalV = new TH1F("hTotalV","hTotalV",m_nSteps,0,tTimeRange);
    m_hTotalRadPower = new TH1F("hTotalRadPower","hTotalRadPower",m_nSteps,0,tTimeRange);
    m_hTotalPdotA = new TH1F("hTotalPdotA","hTotalPdotA",m_nSteps,0,tTimeRange);
    m_pTotalHamil = new TProfile("pTotalHamil","pTotalHamil",m_nSteps,0,tTimeRange);
    m_pTotalMomMS = new TProfile("pTotalMomMS","pTotalMomMS",m_nSteps,0,tTimeRange);
    m_pTotalPx = new TProfile("pTotalPx","pTotalPx",m_nSteps,0,tTimeRange);
    m_pTotalPy = new TProfile("pTotalPy","pTotalPy",m_nSteps,0,tTimeRange);
    m_pTotalPz = new TProfile("pTotalPz","pTotalPz",m_nSteps,0,tTimeRange);
    m_pTotalAx = new TProfile("pTotalAx","pTotalAx",m_nSteps,0,tTimeRange);
    m_pTotalAy = new TProfile("pTotalAy","pTotalAy",m_nSteps,0,tTimeRange);
    m_pTotalAz = new TProfile("pTotalAz","pTotalAz",m_nSteps,0,tTimeRange);
    m_pTotalV = new TProfile("pTotalV","pTotalV",m_nSteps,0,tTimeRange);
    m_pTotalRadPower = new TProfile("pTotalRadPower","pTotalRadPower",m_nSteps,0,tTimeRange);
    m_pTotalPdotA = new TProfile("pTotalPdotA","pTotalPdotA",m_nSteps,0,tTimeRange);
    m_h2TotalHamil = new TH2F("h2TotalHamil","h2TotalHamil",m_nSteps,0,tTimeRange,400,-20,20);
    m_h2TotalMomMS = new TH2F("h2TotalMomMS","h2TotalMomMS",m_nSteps,0,tTimeRange,400,-20,20);
    m_h2TotalPx = new TH2F("h2TotalPx","h2TotalPx",m_nSteps,0,tTimeRange,400,-2,2);
    m_h2TotalPy = new TH2F("h2TotalPy","h2TotalPy",m_nSteps,0,tTimeRange,400,-2,2);
    m_h2TotalPz = new TH2F("h2TotalPz","h2TotalPz",m_nSteps,0,tTimeRange,400,-2,2);
    m_h2TotalAx = new TH2F("h2TotalAx","h2TotalAx",m_nSteps,0,tTimeRange,400,-2,2);
    m_h2TotalAy = new TH2F("h2TotalAy","h2TotalAy",m_nSteps,0,tTimeRange,400,-2,2);
    m_h2TotalAz = new TH2F("h2TotalAz","h2TotalAz",m_nSteps,0,tTimeRange,400,-2,2);
    m_h2TotalV = new TH2F("h2TotalV","h2TotalV",m_nSteps,0,tTimeRange,400,-2,2);
    m_h2TotalRadPower = new TH2F("h2TotalRadPower","h2TotalRadPower",m_nSteps,0,tTimeRange,400,-2,2);
    m_h2TotalPdotA = new TH2F("h2TotalPdotA","h2TotalPdotA",m_nSteps,0,tTimeRange,400,-2,2);

    m_termsCoul = new PairTerms<Force>();
    m_termsVelo = new PairTerms<Force>();
    m_termsAcce = new PairTerms<Force>();
    m_termsMagn = new PairTerms<Force>();
    m_Rs = new PairTerms<Length>();
    m_us = new PairTerms<Dimensionless>();
    m_vs = new PairTerms<Speed>();
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

// Along a trajectory with |v| < c the interval to the point (tPos, tTime) decreases monotonically: it is positive
// (time-like) before the crossing with the past light cone and negative (space-like) after it.
// Points later than tTime are in the future and are treated as space-like.
// Starting from the previous best index, find the last time-like and the first space-like point and return
// the one with smaller |interval|, with tBestInterval = sign(interval) * sqrt(|interval|).
// If the whole trajectory is space-like, returns 0 with tBestInterval < 0;
// if the whole trajectory is time-like, returns size - 1 with tBestInterval > 0.
unsigned long CoulombAfterburner::NearestInterval(const Vec3<Length> &tPos, Time tTime, std::vector<SpaceTimePoint> &tTrajectory,
                                                  unsigned long tLastBestI, Length &tBestInterval, bool debug)
{
    long tSize = tTrajectory.size();
    if (tSize == 0) {
        tBestInterval = NaN<Length>();
        return 0;
    }
    auto tIsTimeLikePast = [&](long j) {
        return tTrajectory[j].t <= tTime && Interval(tPos, tTime, tTrajectory[j]) >= Area::zero();
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
        Area tIntervalTimeLike  = Interval(tPos, tTime, tTrajectory[tLastTimeLike]);
        Area tIntervalSpaceLike = Interval(tPos, tTime, tTrajectory[tFirstSpaceLike]);
        tBestI = abs(tIntervalTimeLike) <= abs(tIntervalSpaceLike) ? tLastTimeLike : tFirstSpaceLike;
    }
    Area tInterval = Interval(tPos, tTime, tTrajectory[tBestI]);
    tBestInterval = (tInterval >= Area::zero() ? 1. : -1.) * Length(sqrt(abs(tInterval)));
    if (debug) {
        cout << "NearestInterval: start " << tLastBestI << " size " << tSize << " last time-like " << tLastTimeLike
            << " best " << tBestI << " interval " << tBestInterval << endl;
    }
    return tBestI;
}

// Time t at which the straight line x(t) = X0 + V (t - T0) crosses the past light cone of the point (Y, U),
// c^2 (U - t)^2 = |Y - x(t)|^2. With D = Y - X0 + V T0 this is
// (c^2 - V^2) t^2 - 2 (c^2 U - V.D) t + c^2 U^2 - D^2 = 0, and the earlier of the two solutions is the retarded time.
// NaN if there is no solution.
Time CoulombAfterburner::NearestPreFreezeoutIntervalTime(const Vec3<Length> &tX0, const Vec3<Speed> &tV, const Vec3<Length> &tY,
                                                         Time tT0, Time tU, bool debug)
{
    const auto tC2 = pow<2>(1. * c);
    Vec3<Length> tD = tY - tX0 + tV * tT0;
    auto tA = tC2 - Mag2(tV);
    auto tB = -2. * (tC2 * tU - Dot(tV, tD));
    auto tC = tC2 * pow<2>(tU) - Mag2(tD);
    auto tDelta = pow<2>(tB) - 4. * tA * tC;
    if (tDelta < decltype(tDelta)::zero()) {
        return NaN<Time>();
    }
    Time t1 = (-tB - sqrt(tDelta)) / (2. * tA);
    Time t2 = (-tB + sqrt(tDelta)) / (2. * tA);
    if (debug) {
        cout << "PARAMS tT0 " << tT0 << " tU " << tU << endl;
        cout << "A: " << tA << " B: " << tB << " C: " << tC << " Delta: " << tDelta << endl;
        cout << "tV:  " << tV << endl;
        cout << "tX0: " << tX0 << endl;
        cout << "tY:  " << tY << endl;
        cout << "tdiff " << t1-tT0 << " " << t2-tT0 << endl;
        cout << "dt " << t1-tU << " " << t2-tU << endl;
        cout << "tX1: " << Vec3<Length>(tX0 + tV*(t1-tT0)) << endl;
        cout << "tX2: " << Vec3<Length>(tX0 + tV*(t2-tT0)) << endl;
        cout << "|dx1|^2: " << Mag2(tY - (tX0 + tV*(t1-tT0))) << endl;
        cout << "|dx2|^2: " << Mag2(tY - (tX0 + tV*(t2-tT0))) << endl;
        cout << "t1 " << t1 << " t2 " << t2 << endl;
        cout << "int1 " << Interval(tX0 + tV*(t1-tT0), t1, tY, tU, true) << " int2 " << Interval(tX0 + tV*(t2-tT0), t2, tY, tU, true) << endl;
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
            tResult += tParticle.GetParticleType()->GetCharge();
        }
    }
    return tResult;
}

Length CoulombAfterburner::SampleImpB(std::list<Particle> *tParticles) {
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
    // The impact parameters of the Glauber MC are in fm
    Length tImpB = m_hImpBprojs.at(tNpart)->GetRandom() * fm;

    cout << "End sampling b" << endl;
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
    [[maybe_unused]] Length tImpB = m_gImpBNPart->Eval(tNpart) * fm;
    Mass tSpectatorMass = (tNp * 0.9382720 + tNn * 0.9395653) * GeV / square(c);
    tSpectatorType->SetNumber(777777);
    tSpectatorType->SetName("Spectator");
    tSpectatorType->SetMass(tSpectatorMass.numerical_value_in(GeV / square(c)));
    tSpectatorType->SetBarionN(tNp + tNn);
    tSpectatorType->SetI3((tNp - tNn)/2);

    // At Npart = 0, nuclei just grazing each other, spectator X would be the radius of the nucleus.
    // At Npart -> A, a head on collision, the center of mass of the spectators would also approach
    // the radius of the nucleus. so we assume that the center of mass will always be at x = tR
    Double_t tSpectatorXtoR = 1.;//tNpart / (2 * tA) + tImpB / (2 * tA) * (2 * tA - tNpart);
    Length tR = 1.2 * TMath::Power(tA, 1./3) * fm;
    tSpectatorX = (tSpectatorXtoR * tR).numerical_value_in(inverse_gev_length); // GeV^-1, like the positions of particles
    cout << "DEFINING SPECTATOR X " << tR << " x " << tSpectatorXtoR << " -> " << tSpectatorX << " GeV^-1" << endl;

    return tSpectatorType;
}

Int_t CoulombAfterburner::AddSpectators(std::list<Particle> *tParticles, Time tEarliestTime)
{
    Int_t tNspec = 0;
    // The spectators get the eids after the largest one in the event: particles read from a file keep their eids,
    // while the counter of Particle starts again from 0 with every event
    Int_t tNextEid = 0;
    for (const Particle &tParticle : *tParticles) {
        tNextEid = TMath::Max(tNextEid, tParticle.eid + 1);
    }

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
    // TGlauberMC works in fm
    Length tImpB = SampleImpB(tParticles);
    m_mcg->SetBmin(tImpB.numerical_value_in(fm));
    m_mcg->SetBmax(tImpB.numerical_value_in(fm));
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
            Mass tSpectatorMass = (79 * 0.9382720 + (197-79) * 0.9395653)/197 * GeV / square(c);

            Double_t tSpectatorBeta = 0.63; // TODO: make an external parameter
            Double_t tSpectatorGamma = 1/TMath::Sqrt(1 - tSpectatorBeta*tSpectatorBeta);
            Momentum tSpectatorPz = tSpectatorDir * tSpectatorGamma * tSpectatorBeta * tSpectatorMass * c;
            CoulombUnits::Energy tSpectatorE = tSpectatorGamma * tSpectatorMass * square(c);
            Length tSpectatorShift = tSpectatorDir * tSpectatorBeta * tEarliestTime * c;
            // Lorentz contraction and shift
            Vec3<Length> tSpectatorPos(tNucl->GetX() * fm, tNucl->GetY() * fm, tNucl->GetZ() * fm / tSpectatorGamma + tSpectatorShift);

            Particle tSpectator(tSpectatorType);
            tSpectator.SetParticlePX(tSpectatorE.numerical_value_in(GeV), 0, 0, tSpectatorPz.numerical_value_in(GeV / c),
                                     tEarliestTime.numerical_value_in(inverse_gev_time), tSpectatorPos.x.numerical_value_in(inverse_gev_length),
                                     tSpectatorPos.y.numerical_value_in(inverse_gev_length), tSpectatorPos.z.numerical_value_in(inverse_gev_length));
            tSpectator.eid = tNextEid++;
            tParticles->push_back(tSpectator);
            ++tNspec;
        }
    }
    cout << "End adding spectators" << endl;
    return tNspec;
}

bool CoulombAfterburner::AreSiblings(const Particle &aParticle1, const Particle &aParticle2)
{
    // Primary particles and spectators have fathereid -1
    return aParticle1.fathereid != -1 && aParticle1.fathereid == aParticle2.fathereid;
}

bool CoulombAfterburner::IsClosePair(Int_t aEid1, Int_t aEid2) const
{
    if (m_ClosePairs.empty()) return false;
    return m_ClosePairs.count(std::make_pair(TMath::Min(aEid1, aEid2), TMath::Max(aEid1, aEid2))) > 0;
}

Vec3<Force> CoulombAfterburner::UniformMotionForce(Int_t aCharges, const Vec3<Length> &aR, const Vec3<Speed> &aSourceVelocity,
        const Vec3<Speed> &aTestVelocity)
{
    // Field of a uniformly moving charge at its present position (the Lienard-Wiechert field without acceleration):
    // E = k q' (1 - beta^2) R / (R^2 - |beta x R|^2)^(3/2), c B = beta x E, force q E + q v x B
    Vec3<Dimensionless> tBetaSource = aSourceVelocity / c;
    Vec3<Dimensionless> tBetaTest = aTestVelocity / c;
    Length tD = Length(sqrt(Mag2(aR) - Mag2(Cross(tBetaSource, aR))));
    Vec3<Force> tE = kCoulomb * aCharges * (1 - Mag2(tBetaSource)) * aR / pow<3>(tD);
    // Using BAC-CAB rule for beta_test x (beta_source x E)
    return tE + tBetaSource * Dot(tBetaTest, tE) - tE * Dot(tBetaTest, tBetaSource);
}

// Finds the pairs of charged particles whose encounter during the time step [t_n, t_n + dt] is too short to be
// followed with the time step, and joins pairs with common particles into groups. All particles are at t_n.
// The encounter time is the shorter of the passage time d/|w| (closest distance d of the straight-line trajectories
// during the step, relative velocity w) and the time sqrt(d^3 mu / k|q q'|) in which the Coulomb force alone changes
// the distance by about d (reduced mass mu, with E/c^2 as the mass; the mass of the particle for a spectator, which
// does not recoil). A pair is close if the step is longer than the encounter time by more than 1/m_SubstepAccuracy;
// the number of sub-steps is then dt / (m_SubstepAccuracy * encounter time), rounded up to a power of 2.
std::vector<CoulombAfterburner::CloseGroup> CoulombAfterburner::FindCloseGroups(std::list<Particle> *tParticles, Int_t iStep)
{
    m_ClosePairs.clear();
    std::vector<CloseGroup> tGroups;
    if (!m_Substeps) {
        return tGroups;
    }
    std::vector<std::list<Particle>::iterator> tCharged;
    for (auto tPartIter = tParticles->begin(); tPartIter != tParticles->end(); ++tPartIter) {
        if (tPartIter->GetParticleType()->GetCharge() != 0
            && iStep >= (Int_t)m_FirstStepExists[tPartIter->eid] && iStep <= (Int_t)m_LastStepExists[tPartIter->eid]) {
            tCharged.push_back(tPartIter);
        }
    }
    // Groups by union-find over the indices in tCharged
    std::vector<size_t> tParent(tCharged.size());
    for (size_t i = 0; i < tParent.size(); ++i) tParent[i] = i;
    std::function<size_t(size_t)> tRoot = [&](size_t i) { return tParent[i] == i ? i : (tParent[i] = tRoot(tParent[i])); };
    struct ClosePair { size_t a, b; int substeps; };
    std::vector<ClosePair> tPairs;

    for (size_t a = 0; a < tCharged.size(); ++a) {
        for (size_t b = a + 1; b < tCharged.size(); ++b) {
            std::list<Particle>::iterator tA = tCharged[a], tB = tCharged[b];
            bool tSpectatorA = IsSpectator(tA), tSpectatorB = IsSpectator(tB);
            if (tSpectatorA && tSpectatorB) continue;
            // Pairs without interaction in CalculateTotalForce stay without it
            if (tA->eid == tB->fathereid || tB->eid == tA->fathereid) continue;
            if (AreSiblings(*tA, *tB)
                && (iStep == (Int_t)m_FirstStepExists[tA->eid] || iStep == (Int_t)m_FirstStepExists[tB->eid])) continue;

            Vec3<Length> tD = PositionOf(*tA) - PositionOf(*tB);
            Vec3<Speed> tW = VelocityOf(*tA) - VelocityOf(*tB);
            auto tW2 = Mag2(tW);
            Time tS = Time::zero();
            if (tW2 > decltype(tW2)::zero()) {
                tS = Time(-Dot(tD, tW) / tW2);
                if (tS < Time::zero()) tS = Time::zero();
                if (tS > m_StepSize) tS = m_StepSize;
            }
            Length tDMin = Mag(tD + tW * tS);
            if (AreSiblings(*tA, *tB)) {
                // Decay products interact only beyond m_SiblingMinDistance: not close if they stay within it during
                // the whole step (the distance is largest at the ends of the step), else not closer than it
                if (Mag(tD) < m_SiblingMinDistance && Mag(tD + tW * m_StepSize) < m_SiblingMinDistance) continue;
                if (tDMin < m_SiblingMinDistance) tDMin = m_SiblingMinDistance;
            }
            if (tDMin > m_SubstepDistance) continue;

            CoulombUnits::Energy tEA = EnergyOf(*tA), tEB = EnergyOf(*tB);
            Mass tMu = tSpectatorA ? Mass(tEB / square(c)) : tSpectatorB ? Mass(tEA / square(c)) : Mass(tEA * tEB / ((tEA + tEB) * square(c)));
            Int_t tCharges = TMath::Abs(tA->GetParticleType()->GetCharge() * tB->GetParticleType()->GetCharge());
            Time tEncounter = Time(sqrt(pow<3>(tDMin) * tMu / (kCoulomb * tCharges)));
            if (tW2 > decltype(tW2)::zero()) {
                Time tPassage = Time(tDMin / Mag(tW));
                if (tPassage < tEncounter) tEncounter = tPassage;
            }
            double tSubsteps = (m_StepSize / (m_SubstepAccuracy * tEncounter)).numerical_value_in(one);
            if (!(tSubsteps > 1)) continue;
            int tN = 1;
            while (tN < tSubsteps && tN < m_SubstepMax) tN *= 2;
            tN = TMath::Min(tN, m_SubstepMax);

            tPairs.push_back(ClosePair{a, b, tN});
            m_ClosePairs.insert(std::make_pair(TMath::Min(tA->eid, tB->eid), TMath::Max(tA->eid, tB->eid)));
            tParent[tRoot(a)] = tRoot(b);
        }
    }

    std::map<size_t, size_t> tGroupOfRoot;
    std::vector<bool> tAdded(tCharged.size(), false);
    for (const ClosePair &tPair : tPairs) {
        size_t tR = tRoot(tPair.a);
        if (tGroupOfRoot.find(tR) == tGroupOfRoot.end()) {
            tGroupOfRoot[tR] = tGroups.size();
            tGroups.push_back(CloseGroup{{}, 1});
        }
        CloseGroup &tGroup = tGroups[tGroupOfRoot[tR]];
        tGroup.substeps = TMath::Max(tGroup.substeps, tPair.substeps);
        for (size_t i : {tPair.a, tPair.b}) {
            if (!tAdded[i]) {
                tGroup.members.push_back(tCharged[i]);
                tAdded[i] = true;
            }
        }
    }
    return tGroups;
}

// Velocity Verlet with the force split into the far part (all particles outside the close pairs, retarded, computed
// with the time step) and the close part (within the close pairs of the group, computed at every sub-step):
// p += F_far(t_n) dt/2; sub-steps of velocity Verlet with F_close; and p += F_far(t_n + dt) dt/2 in Apply().
// The close part uses the field of uniform motion at the simultaneous position of the partner: the light travel
// time over the distances of close pairs is short, and the acceleration term is negligible there.
// Spectators move uniformly, they are only sources of the field.
void CoulombAfterburner::IntegrateCloseGroup(const CloseGroup &aGroup, const Vec3<Force> *tFarForces, Time tCurrentTime)
{
    const size_t tN = aGroup.members.size();
    std::vector<Vec3<Length> > tX(tN);
    std::vector<Vec3<Momentum> > tP(tN);
    std::vector<Vec3<Speed> > tV(tN);
    std::vector<Vec3<Force> > tF(tN);
    std::vector<Mass> tM(tN);
    std::vector<Int_t> tQ(tN), tEid(tN);
    std::vector<std::vector<bool> > tSiblings(tN, std::vector<bool>(tN));
    std::vector<bool> tSpectator(tN);
    for (size_t i = 0; i < tN; ++i) {
        std::list<Particle>::iterator tPartIter = aGroup.members[i];
        tX[i] = PositionOf(*tPartIter);
        tP[i] = MomentumOf(*tPartIter);
        tM[i] = MassOf(*tPartIter);
        tQ[i] = tPartIter->GetParticleType()->GetCharge();
        tEid[i] = tPartIter->eid;
        tSpectator[i] = IsSpectator(tPartIter);
        tV[i] = VelocityOf(*tPartIter);
        if (!tSpectator[i]) {
            tP[i] += 0.5 * tFarForces[tEid[i]] * m_StepSize;
        }
        for (size_t j = 0; j < i; ++j) {
            tSiblings[i][j] = tSiblings[j][i] = AreSiblings(*tPartIter, *aGroup.members[j]);
        }
    }
    auto tUpdateVelocities = [&]() {
        for (size_t i = 0; i < tN; ++i) {
            if (!tSpectator[i]) tV[i] = VelocityOf(tP[i], EnergyOf(tM[i], tP[i]));
        }
    };
    auto tCloseForces = [&]() {
        for (size_t i = 0; i < tN; ++i) {
            tF[i] = Vec3<Force>();
            if (tSpectator[i]) continue;
            for (size_t j = 0; j < tN; ++j) {
                if (j == i || !IsClosePair(tEid[i], tEid[j])) continue;
                // Short-distance cut-off between decay products, with the distance in the rest frame of the source,
                // gamma sqrt(R^2 - |v/c x R|^2) for the simultaneous R
                if (tSiblings[i][j]) {
                    Vec3<Length> tR = tX[i] - tX[j];
                    Vec3<Dimensionless> tBeta = tV[j] / c;
                    if (Mag2(tR) - Mag2(Cross(tBeta, tR)) < (1 - Mag2(tBeta)) * pow<2>(m_SiblingMinDistance)) continue;
                }
                tF[i] += UniformMotionForce(tQ[i] * tQ[j], tX[i] - tX[j], tV[j], tV[i]);
            }
        }
    };

    const Time tH = m_StepSize / aGroup.substeps;
    tUpdateVelocities();
    tCloseForces();
    for (int k = 0; k < aGroup.substeps; ++k) {
        for (size_t i = 0; i < tN; ++i) {
            if (!tSpectator[i]) tP[i] += 0.5 * tF[i] * tH;
        }
        tUpdateVelocities();
        for (size_t i = 0; i < tN; ++i) {
            tX[i] += tV[i] * tH;
        }
        // The magnetic part depends on the velocities at the end of the sub-step: evaluate the force once with the
        // velocities after the first half kick, and again with the velocities predicted from it (second order)
        tCloseForces();
        std::vector<Vec3<Momentum> > tPHalf = tP;
        for (size_t i = 0; i < tN; ++i) {
            if (!tSpectator[i]) tP[i] = tPHalf[i] + 0.5 * tF[i] * tH;
        }
        tUpdateVelocities();
        tCloseForces();
        for (size_t i = 0; i < tN; ++i) {
            if (!tSpectator[i]) tP[i] = tPHalf[i] + 0.5 * tF[i] * tH;
        }
        tUpdateVelocities();
    }

    for (size_t i = 0; i < tN; ++i) {
        std::list<Particle>::iterator tPartIter = aGroup.members[i];
        SetSpaceTime(*tPartIter, tCurrentTime + m_StepSize, tX[i]);
        if (!tSpectator[i]) {
            SetFourMomentum(*tPartIter, EnergyOf(tM[i], tP[i]), tP[i]);
        }
    }
}

void CoulombAfterburner::Apply(Event *tEvent)
{
    if (m_nSteps <= 0) return;

    cout << " KUKU 1 " << endl;

    bool debug = false;

    std::list<Particle> *tParticles = tEvent->GetParticleList();
    Time tEarliestTime = NaN<Time>();
    for (auto tPartIter = tParticles->begin(); tPartIter != tParticles->end(); ++tPartIter) {
        if (!HasFiniteKinematics(*tPartIter))
            continue;
        if (!IsFinite(tEarliestTime) || TimeOf(*tPartIter) < tEarliestTime)
        {
            tEarliestTime = TimeOf(*tPartIter);
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
    m_coordinates = new std::vector<SpaceTimePoint>[N];
    m_fourmomenta = new std::vector<FourMomentum>[N];
    m_velocities = new std::vector<Vec3<Speed> >[N];
    m_accelerations = new std::vector<Vec3<Acceleration> >[N];
    m_pids = new std::map<int, int>();
    m_DecayTimes = new Time[N];
    m_CreationTimes = new Time[N];
    m_FirstStepExists = new unsigned int[N];
    m_LastStepExists = new unsigned int[N];
    // The original coordinates as stored in the particles (GeV^-1), restored unchanged at the end
    TLorentzVector *tOriginalCoordinates = new TLorentzVector[N];

    for (auto tPartIter = tParticles->begin(); tPartIter != tParticles->end(); ++tPartIter)
    {
        m_coordinates[tPartIter->eid] = std::vector<SpaceTimePoint>();
        tOriginalCoordinates[tPartIter->eid] = TLorentzVector(tPartIter->VecX(), tPartIter->T());
        m_CreationTimes[tPartIter->eid] = TimeOf(*tPartIter);
        m_DecayTimes[tPartIter->eid] = std::numeric_limits<double>::infinity() * fm_c;
        if (tPartIter->decayed)
        {
            std::list<Particle>::iterator tChildPartIter = tParticles->begin();
            for (; tChildPartIter != tParticles->end(); ++tChildPartIter)
            {
                if (tChildPartIter->fathereid == tPartIter->eid)
                {
                    m_DecayTimes[tPartIter->eid] = TimeOf(*tChildPartIter);
                    break;
                }
            }
        }
        // A particle exists at the time steps t with creation time <= t < decay time, so that the first step
        // of a decay product directly follows the last step of its father
        double tFirstStep = TMath::Min((double)m_nSteps,
                                       TMath::Ceil(((TimeOf(*tPartIter) - tEarliestTime) / m_StepSize).numerical_value_in(one)));
        double tLastStep = m_nSteps - 1;
        if (IsFinite(m_DecayTimes[tPartIter->eid]))
        {
            tLastStep = TMath::Min(tLastStep, TMath::Ceil(((m_DecayTimes[tPartIter->eid] - tEarliestTime) / m_StepSize).numerical_value_in(one)) - 1);
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

    // Sums over the particles, to monitor conservation laws; p.A in GeV/c * GeV
    using MomentumTimesEnergy = quantity<GeV / c * GeV, double>;
    // Reference values from the first step
    CoulombUnits::Energy tTotalHamil_0 = CoulombUnits::Energy::zero();
    CoulombUnits::Energy tTotalMomMS_0 = CoulombUnits::Energy::zero();
    Vec3<Momentum> tTotalP_0;
    Vec3<CoulombUnits::Energy> tTotalA_0;
    CoulombUnits::Energy tTotalV_0 = CoulombUnits::Energy::zero();
    MomentumTimesEnergy tTotalPdotA_0 = MomentumTimesEnergy::zero();
    CoulombUnits::Energy tRadiatedEnergy = CoulombUnits::Energy::zero();
    CoulombUnits::Energy tRadiatedEnergy_0 = CoulombUnits::Energy::zero();

    Vec3<Force> *tTotalForcesLast = NULL;
    // Particles in close pairs in the previous step, whose force from that step lacks the close pairs
    std::vector<char> tCloseLast(N, 0);
    long tNClosePairSteps = 0;
    int tMaxSubsteps = 0;
    for (unsigned int iStep = 0; iStep < m_nSteps; ++iStep)
    {
        m_termsCoul->emplace_back();
        m_termsVelo->emplace_back();
        m_termsAcce->emplace_back();
        m_termsMagn->emplace_back();
        m_Rs->emplace_back();
        m_us->emplace_back();
        m_vs->emplace_back();

        Time tCurrentTime = tEarliestTime + iStep * m_StepSize;
        Vec3<Force> *tTotalForces = new Vec3<Force>[N];
        std::list<Particle>::iterator tPartIter = tParticles->begin();
        for (; tPartIter != tParticles->end(); ++tPartIter)
        {
            if (iStep >= m_FirstStepExists[tPartIter->eid] && iStep <= m_LastStepExists[tPartIter->eid])
            {
                if (iStep == m_FirstStepExists[tPartIter->eid])
                {
                    // Propagate freely from the creation time to the first time step
                    SetSpaceTime(*tPartIter, tCurrentTime,
                                 PositionOf(*tPartIter) + VelocityOf(*tPartIter) * (tCurrentTime - TimeOf(*tPartIter)));
                }
                Int_t tQ = tPartIter->GetParticleType()->GetCharge();
                m_coordinates[tPartIter->eid].push_back(SpaceTimePoint{tCurrentTime, PositionOf(*tPartIter)});
                m_fourmomenta[tPartIter->eid].push_back(FourMomentum{EnergyOf(*tPartIter), MomentumOf(*tPartIter)});
                (*m_pids)[tPartIter->eid] = tPartIter->pid;
                if (tQ != 0)
                {
                    m_velocities[tPartIter->eid].push_back(VelocityOf(*tPartIter));

                    if (iStep == m_FirstStepExists[tPartIter->eid])
                    {
                        m_accelerations[tPartIter->eid].push_back(Vec3<Acceleration>());
                    }
                    else
                    {
                        const Vec3<Speed> &tV1 = m_velocities[tPartIter->eid][iStep - m_FirstStepExists[tPartIter->eid] - 1];
                        const Vec3<Speed> &tV2 = m_velocities[tPartIter->eid][iStep - m_FirstStepExists[tPartIter->eid]];
                        m_accelerations[tPartIter->eid].push_back((tV2 - tV1) / m_StepSize);
                    }
                }
            }
        }
        // Close pairs of this step: their mutual force is left out of the force below and integrated with sub-steps
        std::vector<CloseGroup> tCloseGroups = FindCloseGroups(tParticles, iStep);
        std::vector<char> tClose(N, 0);
        for (const CloseGroup &tGroup : tCloseGroups) {
            for (const std::list<Particle>::iterator &tMember : tGroup.members) {
                tClose[tMember->eid] = 1;
            }
            tMaxSubsteps = TMath::Max(tMaxSubsteps, tGroup.substeps);
        }
        tNClosePairSteps += m_ClosePairs.size();

        // Calculate force at old position
        tPartIter = tParticles->begin();
        for (; tPartIter != tParticles->end(); ++tPartIter)
        {

            Vec3<Force> tForceTotal;
            // For spectator nucleons, the force stays zero
            // Reuse the force at the new position from the previous step, unless the particle has just been created,
            // or it is in a close pair in this or the previous step (the close pairs left out of the force differ)
            Int_t tQ = tPartIter->GetParticleType()->GetCharge();
            if (tQ != 0 && iStep >= m_FirstStepExists[tPartIter->eid] && iStep <= m_LastStepExists[tPartIter->eid])
            {
                if (tTotalForcesLast != NULL && iStep != m_FirstStepExists[tPartIter->eid]
                    && !tClose[tPartIter->eid] && !tCloseLast[tPartIter->eid]) {
                    tForceTotal = tTotalForcesLast[tPartIter->eid];
                } else {
                    tForceTotal = CalculateTotalForce(tPartIter, tParticles, tBestIs, tCurrentTime, iStep);
                }
            }
            if (m_FindProblem && tPartIter->eid == m_DebugEid) {
                cout << "total force at old pos, step " << iStep << ", eid = " << tPartIter->eid << ", F = " << tForceTotal << endl;
            }
            tTotalForces[tPartIter->eid] = tForceTotal;
        }
        // Update positions x_{n+1} = x_n + v_n dt + a_n dt^2 / 2, using the force calculated at old position
        tPartIter = tParticles->begin();
        for (; tPartIter != tParticles->end(); ++tPartIter)
        {
            if (iStep < m_FirstStepExists[tPartIter->eid] || iStep > m_LastStepExists[tPartIter->eid])
                continue;
            // Particles in close groups are advanced below
            if (tClose[tPartIter->eid])
                continue;

            const Vec3<Force> &tForceTotal = tTotalForces[tPartIter->eid];
            Vec3<Speed> tV = VelocityOf(*tPartIter);
            CoulombUnits::Energy tE = EnergyOf(*tPartIter);
            // Relativistic acceleration dv/dt = (F - v (v.F) / c^2) c^2 / E
            Vec3<Acceleration> tAcc = (tForceTotal - tV * Dot(tV, tForceTotal) / square(c)) * square(c) / tE;
            Vec3<Length> tPartPos = PositionOf(*tPartIter) + tV * m_StepSize + 0.5 * tAcc * pow<2>(m_StepSize);
            if (debug && !IsFinite(tPartPos))
            {
                cout << "updating position: " << tE << " " << Mag(MomentumOf(*tPartIter)) << endl;
            }

            SetSpaceTime(*tPartIter, tCurrentTime + m_StepSize, tPartPos);
        }
        for (const CloseGroup &tGroup : tCloseGroups) {
            IntegrateCloseGroup(tGroup, tTotalForces, tCurrentTime);
        }

        // Calculate force at new position x_{n+1} and time t_{n+1}, still with the old velocity; the magnetic force
        // with the new velocity is calculated below, from q c B
        Vec3<Force> *tTotalForcesNew = new Vec3<Force>[N];
        std::vector<Vec3<Force> > tMagneticForcesNew(N);
        tPartIter = tParticles->begin();
        for (; tPartIter != tParticles->end(); ++tPartIter)
        {

            Vec3<Force> tForceTotal;
            Int_t tQ = tPartIter->GetParticleType()->GetCharge();
            if (tQ != 0 && iStep >= m_FirstStepExists[tPartIter->eid] && iStep <= m_LastStepExists[tPartIter->eid])
            {
                tForceTotal = CalculateTotalForce(tPartIter, tParticles, tBestIs, tCurrentTime + m_StepSize, iStep,
                                                  &tMagneticForcesNew[tPartIter->eid]);
            }
            if (m_FindProblem && tPartIter->eid == m_DebugEid) {
                cout << "total force at new pos, step " << iStep << ", eid = " << tPartIter->eid << ", F = " << tForceTotal << endl;
            }
            tTotalForcesNew[tPartIter->eid] = tForceTotal;
        }
        // Update momenta p_{n+1} = p_n + (F_n + F_{n+1}) dt / 2, using average force from old and new position.
        // F_{n+1} contains the magnetic force with v_{n+1}, which depends on p_{n+1}: it is calculated with the
        // velocity predicted with v_n in F_{n+1}, which makes the error of the magnetic force second order.
        tPartIter = tParticles->begin();
        for (; tPartIter != tParticles->end(); ++tPartIter)
        {
            if (iStep < m_FirstStepExists[tPartIter->eid] || iStep > m_LastStepExists[tPartIter->eid])
                continue;

            // Momentum after the half kick with the force at the old position. For close groups, this kick and the
            // forces within the group are already in the momentum (IntegrateCloseGroup)
            Vec3<Momentum> tPartMomHalf = MomentumOf(*tPartIter);
            if (!tClose[tPartIter->eid]) {
                tPartMomHalf += 0.5 * tTotalForces[tPartIter->eid] * m_StepSize;
            }
            Vec3<Force> &tForceNew = tTotalForcesNew[tPartIter->eid];
            Vec3<Momentum> tPartMomPred = tPartMomHalf + 0.5 * tForceNew * m_StepSize;
            Vec3<Speed> tVPred = VelocityOf(tPartMomPred, EnergyOf(MassOf(*tPartIter), tPartMomPred));
            tForceNew += Cross((tVPred - VelocityOf(*tPartIter)) / c, tMagneticForcesNew[tPartIter->eid]);
            if (m_FindProblem && tPartIter->eid == m_DebugEid) {
                cout << "total force " << iStep << ", eid = " << tPartIter->eid << ", F = "
                    << (tClose[tPartIter->eid] ? 0.5 * tForceNew : 0.5 * (tTotalForces[tPartIter->eid] + tForceNew)) << endl;
            }

            Vec3<Momentum> tPartMom = tPartMomHalf + 0.5 * tForceNew * m_StepSize;
            CoulombUnits::Energy tE = EnergyOf(MassOf(*tPartIter), tPartMom);
            if (m_FindProblem && tPartIter->eid == m_DebugEid) {
                cout << "before updating momenta " << iStep << ", eid = " << tPartIter->eid
                    << ", p = " << MomentumOf(*tPartIter) << ", |p| = " << Mag(MomentumOf(*tPartIter))
                    << ", E = " << EnergyOf(*tPartIter) << ", mass = " << MassOf(*tPartIter) << endl;
            }

            SetFourMomentum(*tPartIter, tE, tPartMom);
            if (m_FindProblem && tPartIter->eid == m_DebugEid) {
                cout << "after updating momenta " << iStep << ", eid = " << tPartIter->eid
                    << ", p = " << MomentumOf(*tPartIter) << ", |p| = " << Mag(MomentumOf(*tPartIter))
                    << ", E = " << EnergyOf(*tPartIter) << ", mass = " << MassOf(*tPartIter) << endl;
            }

            Vec3<Speed> tVnew = VelocityOf(tPartMom, tE);
            if (m_FindProblem && (!IsFinite(Mag(tVnew)) || Mag(tVnew) >= 1. * c))
            {
                cout << "In step " << iStep << " ERROR, wrong velocity eid " << tPartIter->eid << ", pid " << tPartIter->pid << endl;
                cout << "In step " << iStep << " p = " << tPartMom << endl;
                cout << "In step " << iStep << " v = " << tVnew << endl;
                cout << "In step " << iStep << " |p| = " << Mag(tPartMom) << ", E = " << tE << endl;
                for (auto mapping : m_termsCoul->back()) {
                    if (mapping.first.first == tPartIter->eid && Mag(mapping.second) > m_LargeTerm) {
                        cout << "possible reason: " << mapping.second << endl;
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
        tCloseLast.swap(tClose);

        // Calculating total energy, momentum and potentials at time t_{n+1}
        CoulombUnits::Energy tTotalHamil = CoulombUnits::Energy::zero();
        CoulombUnits::Energy tTotalMomMS = CoulombUnits::Energy::zero();
        Vec3<Momentum> tTotalP;
        Vec3<CoulombUnits::Energy> tTotalA;
        CoulombUnits::Energy tTotalV = CoulombUnits::Energy::zero();
        MomentumTimesEnergy tTotalPdotA = MomentumTimesEnergy::zero();
        Power tTotalRadPower = Power::zero();

        tPartIter = tParticles->begin();
        for (; tPartIter != tParticles->end(); ++tPartIter)
        {
            if (iStep < m_FirstStepExists[tPartIter->eid] || iStep > m_LastStepExists[tPartIter->eid])
                continue;

            Int_t tQ = tPartIter->GetParticleType()->GetCharge();
            Vec3<Momentum> tPartMom = MomentumOf(*tPartIter); // kinetic momentum

            auto tPotentials = CalculatePotentials(tPartIter, tParticles, tBestIs, tCurrentTime + m_StepSize, iStep);
            const Vec3<CoulombUnits::Energy> &tPotA = tPotentials.first;
            CoulombUnits::Energy tPotV = tPotentials.second;

            CoulombUnits::Energy tEkin = EnergyOf(MassOf(*tPartIter), tPartMom);
            // Each pair contributes to the potentials of both particles, hence 1/2 in the potential energy
            CoulombUnits::Energy tH = tEkin + 0.5*tQ*tPotV;
            tTotalHamil += tH;
            tTotalMomMS += tEkin;
            tTotalP += tPartMom;
            tTotalA += tQ * tPotA;
            tTotalV += 0.5*tQ*tPotV;
            tTotalPdotA += 2*tQ*Dot(tPartMom, tPotA);

            if (tQ != 0 && m_velocities[tPartIter->eid].size() > 0 && m_accelerations[tPartIter->eid].size() > 0) {
                Vec3<Dimensionless> tBeta = m_velocities[tPartIter->eid].back() / c;
                double tGamma2 = 1 / (1 - Mag2(tBeta)).numerical_value_in(one);
                const Vec3<Acceleration> &tA = m_accelerations[tPartIter->eid].back();
                // Lienard's generalisation of the Larmor formula
                Power tRadPower = (2./3) * kCoulomb * tQ * tQ * TMath::Power(tGamma2, 3) * (Mag2(tA) - Mag2(Cross(tBeta, tA))) / pow<3>(1. * c);
                if (tRadPower < Power::zero()) {
                    PrintProblem(tPartIter->eid, iStep, "");
                }
                tTotalRadPower += tRadPower;
            }
            if (m_FindProblem && tPartIter->eid == m_DebugEid) {
                cout << "at end of step " << iStep << ", eid = " << tPartIter->eid
                    << ", p = " << tPartMom << ", |p| = " << Mag(tPartMom) << ", E = " << EnergyOf(*tPartIter)
                    << ", mass = " << MassOf(*tPartIter) << endl;
            }
        }
        tRadiatedEnergy += tTotalRadPower * m_StepSize;

        if (iStep == 0) {
            tTotalHamil_0 = tTotalHamil;
            tTotalMomMS_0 = tTotalMomMS;
            tTotalP_0 = tTotalP;
            tTotalA_0 = tTotalA;
            tTotalV_0 = tTotalV;
            tTotalPdotA_0 = tTotalPdotA;
            tRadiatedEnergy_0 = tRadiatedEnergy;
        } else {
            // Histograms: energies and potentials in GeV, momenta in GeV/c, p.A in GeV^2/c, time axis in GeV^-1
            const double tT = ((iStep + 0.1) * m_StepSize).numerical_value_in(inverse_gev_time);
            const double tHamil = tTotalHamil.numerical_value_in(GeV), tHamil_0 = tTotalHamil_0.numerical_value_in(GeV);
            const double tMomMS = tTotalMomMS.numerical_value_in(GeV), tMomMS_0 = tTotalMomMS_0.numerical_value_in(GeV);
            const Vec3<Momentum> tDP = tTotalP - tTotalP_0;
            const Vec3<CoulombUnits::Energy> tDA = tTotalA - tTotalA_0;
            const double tDV = (tTotalV - tTotalV_0).numerical_value_in(GeV);
            const double tDPdotA = (tTotalPdotA - tTotalPdotA_0).numerical_value_in(GeV / c * GeV);
            const double tRadiated = (tRadiatedEnergy - tRadiatedEnergy_0).numerical_value_in(GeV);
            m_hTotalHamil->SetBinContent(iStep+1, tHamil);
            m_hTotalMomMS->SetBinContent(iStep+1, tMomMS);
            m_hTotalPx->SetBinContent(iStep+1, tTotalP.x.numerical_value_in(GeV / c));
            m_hTotalPy->SetBinContent(iStep+1, tTotalP.y.numerical_value_in(GeV / c));
            m_hTotalPz->SetBinContent(iStep+1, tTotalP.z.numerical_value_in(GeV / c));
            m_hTotalAx->SetBinContent(iStep+1, tTotalA.x.numerical_value_in(GeV));
            m_hTotalAy->SetBinContent(iStep+1, tTotalA.y.numerical_value_in(GeV));
            m_hTotalAz->SetBinContent(iStep+1, tTotalA.z.numerical_value_in(GeV));
            m_hTotalV->SetBinContent(iStep+1, tTotalV.numerical_value_in(GeV));
            // RadPower histograms contain the energy radiated since the first step
            m_hTotalRadPower->SetBinContent(iStep+1, tRadiated);
            m_hTotalPdotA->SetBinContent(iStep+1, tTotalPdotA.numerical_value_in(GeV / c * GeV));
            m_pTotalHamil->Fill(tT, tHamil-tHamil_0);
            m_pTotalMomMS->Fill(tT, tMomMS-tMomMS_0);
            m_pTotalPx->Fill(tT, tDP.x.numerical_value_in(GeV / c));
            m_pTotalPy->Fill(tT, tDP.y.numerical_value_in(GeV / c));
            m_pTotalPz->Fill(tT, tDP.z.numerical_value_in(GeV / c));
            m_pTotalAx->Fill(tT, tDA.x.numerical_value_in(GeV));
            m_pTotalAy->Fill(tT, tDA.y.numerical_value_in(GeV));
            m_pTotalAz->Fill(tT, tDA.z.numerical_value_in(GeV));
            m_pTotalV->Fill(tT, tDV);
            m_pTotalRadPower->Fill(tT, tRadiated);
            m_pTotalPdotA->Fill(tT, tDPdotA);
            m_h2TotalHamil->Fill(tT, tHamil-tHamil_0);
            m_h2TotalMomMS->Fill(tT, tMomMS-tMomMS_0);
            m_h2TotalPx->Fill(tT, tDP.x.numerical_value_in(GeV / c));
            m_h2TotalPy->Fill(tT, tDP.y.numerical_value_in(GeV / c));
            m_h2TotalPz->Fill(tT, tDP.z.numerical_value_in(GeV / c));
            m_h2TotalAx->Fill(tT, tDA.x.numerical_value_in(GeV));
            m_h2TotalAy->Fill(tT, tDA.y.numerical_value_in(GeV));
            m_h2TotalAz->Fill(tT, tDA.z.numerical_value_in(GeV));
            m_h2TotalV->Fill(tT, tDV);
            m_h2TotalRadPower->Fill(tT, tRadiated);
            m_h2TotalPdotA->Fill(tT, tDPdotA);
        }
        if (tTotalForces != NULL) {
            delete [] tTotalForces;
        }
    }
    if (tTotalForcesLast != NULL) {
        delete [] tTotalForcesLast;
    }
    m_ClosePairs.clear();
    if (m_Substeps) {
        cout << "<CoulombAfterburner::Apply>\tClose pairs integrated with sub-steps: " << tNClosePairSteps
            << " pair-steps, at most " << tMaxSubsteps << " sub-steps in a step" << endl;
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
    Vec3<Dimensionless> tBeta = m_velocities[eid].at(index) / c;
    double tGamma2 = 1 / (1 - Mag2(tBeta)).numerical_value_in(one);
    const Vec3<Acceleration> &tA = m_accelerations[eid].at(index);
    const FourMomentum &tPartFourMom = m_fourmomenta[eid].at(index);
    cout << indent.Data() << "problem in step " << step << ", pid = " << (*m_pids)[eid] << ", eid = " << eid
        << ", gamma^2 = " << tGamma2
        << ", |a| = " << Mag(tA) << ", |v|/c = " << Mag(tBeta) << ", a^2 - (v/c x a)^2 = "
        << Mag2(tA) - Mag2(Cross(tBeta, tA)) << endl;
    cout << indent.Data() << "|p| = " << Mag(tPartFourMom.p) << ", E = " << tPartFourMom.e << endl;
    for (int is = 0; is < (int)m_termsCoul->size(); ++is) {
        for (auto mapping : m_termsCoul->at(is)) {
            if (mapping.first.first == eid && Mag(mapping.second) > m_LargeTerm) {
                cout << indent.Data() << "possible reason: " << mapping.second << endl;
                cout << indent.Data() << "step number " << is << ", other eid: " << mapping.first.second << endl;
                cout << indent.Data() << "R: " << Mag(m_Rs->at(is)[std::make_pair(mapping.first.first, mapping.first.second)]) << endl;
                cout << indent.Data() << "u: " << Mag(m_us->at(is)[std::make_pair(mapping.first.first, mapping.first.second)]) << endl;
                cout << indent.Data() << "v: " << Mag(m_vs->at(is)[std::make_pair(mapping.first.first, mapping.first.second)]) << endl;
                PrintProblem(mapping.first.second, is, indent + "**** ", depth + 1);
            }
        }
    }

}
