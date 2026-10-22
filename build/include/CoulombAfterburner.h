#ifndef COULOMB_AFTERBURNER
#define COULOMB_AFTERBURNER

#include "AbstractAfterburner.h"
#include "AbstractEventSaver.h"
#include "ParticleDB.h"
#include "TGlauberMC.h"
#include "TLorentzVector.h"

#include "TFile.h"
#include "TGraph.h"
#include "TProfile.h"

class CoulombAfterburner : public AbstractAfterburner {
  public:
    CoulombAfterburner(int, double, double, AbstractEventSaver *, ParticleDB *);
    virtual ~CoulombAfterburner();
    virtual void Init();
    virtual void Apply(Event *);
    virtual void Finalize();
    static ParticleType *SpectatorType(std::list<Particle> *, Double_t &tSpectatorX);
    // For checks: keep the recorded trajectory (t, x, y, z in GeV^-1) and four-momenta (GeV) of every particle
    // at all time steps of the last event, indexed by eid; off by default
    void SetKeepTrajectories(bool aKeep) { mKeepTrajectories = aKeep; }
    const std::vector<TLorentzVector> &GetTrajectory(Int_t aEid) const { return mKeptCoordinates.at(aEid); }
    const std::vector<TLorentzVector> &GetFourMomenta(Int_t aEid) const { return mKeptFourMomenta.at(aEid); }

  protected:
    void ReadParameters();
    unsigned int m_nSteps;
    Double_t m_StepSize;

  private:
    static Int_t TotalCharge(std::list<Particle> *);
    static TGraph *m_gImpBNPart;

    bool ParticleExists(std::list<Particle>::iterator &tPartIter, double tDecayTime, double tTime);
    unsigned long NearestInterval(TVector3 &tPos, double tTime, std::vector<TLorentzVector> &tTrajectory, unsigned long tLastBestI, double &tInterval, bool debug = false);
    double NearestPreFreezeoutIntervalTime(TVector3 &tX0, TVector3 &tV,TVector3 &tY,  double tT0, double tU, bool debug = false);
    TVector3 CalculateTotalForce(std::list<Particle>::iterator &, std::list<Particle> *, int ** tBestIs, double tCurrentTime,
            Int_t iStep
            );
    std::pair<TVector3, Double_t>  CalculatePotentials(std::list<Particle>::iterator &, std::list<Particle> *, int ** tBestIs,
            double tCurrentTime, Int_t iStep
            );
    bool RetardedSource(std::list<Particle>::iterator &tPartIter, std::list<Particle>::iterator &tPartOtherIter,
            TVector3 &tPos, double tTime, int **tBestIs, Int_t iStep,
            TVector3 &tOtherPos, TVector3 &tOtherVel, TVector3 &tOtherAcc, bool debug = false);
    static bool IsSpectator(std::list<Particle>::iterator &);
    double Interval(TVector3, double, TVector3, double, bool debug = false);
    double Interval(TVector3, double, TLorentzVector, bool debug = false);
    double Interval(TLorentzVector, TLorentzVector, bool debug = false);
    Int_t AddSpectators(std::list<Particle> *, Double_t tEarliestTime);
    Double_t SampleImpB(std::list<Particle> *);

    Double_t mR;
    AbstractEventSaver *mEventSaver;

    TFile *m_fileHistOut;

    TH1F *m_hTotalHamil;
    TH1F *m_hTotalMomMS;
    TH1F *m_hTotalPx;
    TH1F *m_hTotalPy;                                         
    TH1F *m_hTotalPz;
    TH1F *m_hTotalAx;
    TH1F *m_hTotalAy;                                         
    TH1F *m_hTotalAz;
    TH1F *m_hTotalV;
    TH1F *m_hTotalRadPower;
    TH1F *m_hTotalPdotA;
    TProfile *m_pTotalHamil;
    TProfile *m_pTotalMomMS;
    TProfile *m_pTotalPx;
    TProfile *m_pTotalPy;
    TProfile *m_pTotalPz;
    TProfile *m_pTotalAx;
    TProfile *m_pTotalAy;
    TProfile *m_pTotalAz;
    TProfile *m_pTotalV;
    TProfile *m_pTotalRadPower;
    TProfile *m_pTotalPdotA;
    TH2F *m_h2TotalHamil;
    TH2F *m_h2TotalMomMS;
    TH2F *m_h2TotalPx;
    TH2F *m_h2TotalPy;
    TH2F *m_h2TotalPz;
    TH2F *m_h2TotalAx;
    TH2F *m_h2TotalAy;
    TH2F *m_h2TotalAz;
    TH2F *m_h2TotalV;
    TH2F *m_h2TotalRadPower;
    TH2F *m_h2TotalPdotA;

/*
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
*/
    TGlauberMC *m_mcg;
    ParticleType *mProtonSpectatorAType;
    ParticleType *mNeutronSpectatorAType;
    ParticleType *mProtonSpectatorBType;
    ParticleType *mNeutronSpectatorBType;

    TH2F *m_hNpartImpB;
    std::map<Int_t, TH1F *> m_hImpBprojs;

//    static csonstexpr double qe = sqrt(1.4399764e-3); // sqrt(GeV * fm)                                                                  
    static constexpr double epsilon_0 = 55.26349406; // e^2 / (GeV * fm) = 1 / (GeV * fm) - we assume the unit charge = 1;
    static constexpr double k = 1/(4*TMath::Pi()*epsilon_0); // GeV * fm

  unsigned int *m_FirstStepExists;
  unsigned int *m_LastStepExists;

  std::vector< std::map< std::pair< int, int>, TVector3> > *m_termsCoul;
  std::vector< std::map< std::pair< int, int>, TVector3> > *m_termsVelo;
  std::vector< std::map< std::pair< int, int>, TVector3> > *m_termsAcce;
  std::vector< std::map< std::pair< int, int>, TVector3> > *m_termsMagn;
  std::vector< std::map< std::pair< int, int>, TVector3> > *m_Rs;
  std::vector< std::map< std::pair< int, int>, TVector3> > *m_us;
  std::vector< std::map< std::pair< int, int>, TVector3> > *m_vs;
  std::map<int, int> *m_pids;

  std::vector<TLorentzVector> *m_coordinates;
  std::vector<TLorentzVector> *m_fourmomenta;
  std::vector<TVector3> *m_velocities;
  std::vector<TVector3> *m_accelerations;

   void PrintProblem(Int_t eid, Int_t step, TString indent, Int_t depth = 0);

  bool mKeepTrajectories;
  std::vector<std::vector<TLorentzVector> > mKeptCoordinates;
  std::vector<std::vector<TLorentzVector> > mKeptFourMomenta;

};

#endif
