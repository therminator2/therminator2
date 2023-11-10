#ifndef COULOMB_AFTERBURNER
#define COULOMB_AFTERBURNER

#include "AbstractAfterburner.h"
#include "TLorentzVector.h"

#include "TFile.h"
#include "TGraph.h"
#include "TProfile.h"

class CoulombAfterburner : public AbstractAfterburner {
  public:
    CoulombAfterburner(int, double);
    virtual ~CoulombAfterburner();
    virtual void Apply(Event *);

  protected:
    void ReadParameters();
    unsigned int m_nSteps;
    Double_t m_StepSize;

  private:
    bool ParticleExists(std::list<Particle>::iterator &tPartIter, double tDecayTime, double tTime);
    int NearestInterval(TVector3 &tPos, double tTime, std::vector<TLorentzVector> &tTrajectory, int tLastBestI, double &tInterval, bool debug = false);
    double NearestPreFreezeoutIntervalTime(TVector3 &tX0, TVector3 &tV, double tTime0);
    TVector3 CalculateTotalForce(std::list<Particle>::iterator &, std::list<Particle> *, int ** tBestIs, double tCurrentTime,
            std::vector<TLorentzVector> *tCoordinates, std::vector<TVector3> *tVelocities, std::vector<TVector3> *tAccelerations
            );

    TFile *m_fileHistOut;
    TProfile *m_hDistPresentRetarded;
    TProfile *m_hImpulseTimeStep;
    TProfile *m_hWorkTimeStep;
    TProfile *m_hProtonsTimeStep;
    TProfile *m_hProtonsInsideTimeStep;
    TProfile *m_hInteractionsTimeStep;
    TProfile *m_hElseTimeStep;
    TProfile *m_hElseBestInegTimeStep;
    TProfile *m_hElseBestIltSizeTimeStep;
    TProfile *m_hContinueTimeStep;
    TProfile *m_hNetChargeTimeStep;
    TProfile *m_hZeffTimeStep;
    TProfile *m_hZeffR;
    TGraph *m_gPositionTimeStep_1;
    TGraph *m_gPositionTimeStep_2;

    static constexpr double qe = sqrt(1.4399764e-3); // sqrt(GeV * fm)                                                                  

};

#endif
