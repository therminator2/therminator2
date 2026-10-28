#ifndef COULOMB_AFTERBURNER
#define COULOMB_AFTERBURNER

#include "AbstractAfterburner.h"
#include "AbstractEventSaver.h"
#include "Configurator.h"
#include "CoulombUnits.h"
#include "ParticleDB.h"
#include "TGlauberMC.h"
#include "TLorentzVector.h"

#include "TFile.h"
#include "TGraph.h"
#include "TProfile.h"

#include <set>

// Quantities with units (CoulombUnits.h) are used throughout; numbers without units only where the class exchanges
// data with the rest of THERMINATOR, which uses GeV^-1 for positions and times, and GeV for energies and momenta.
class CoulombAfterburner : public AbstractAfterburner {
  public:
    // Number of time steps, length of a time step in fm/c
    CoulombAfterburner(int, double, double, AbstractEventSaver *, ParticleDB *);
    virtual ~CoulombAfterburner();
    virtual void Init();
    virtual void Apply(Event *);
    virtual void Finalize();
    static ParticleType *SpectatorType(std::list<Particle> *, Double_t &tSpectatorX);
    // For checks: keep the recorded trajectory (t, x, y, z in GeV^-1) and four-momenta (GeV) of every particle
    // at all time steps of the last event, indexed by eid; off by default
    void SetKeepTrajectories(bool aKeep) { mKeepTrajectories = aKeep; }
    std::vector<TLorentzVector> GetTrajectory(Int_t aEid) const;
    std::vector<TLorentzVector> GetFourMomenta(Int_t aEid) const;
    // Optional settings from the main configuration (events.ini); keys that are not present keep the defaults
    void ReadParameters(Configurator *aMainConfig);

    // A point of a recorded trajectory and the four-momentum there
    struct SpaceTimePoint {
        CoulombUnits::Time t;
        CoulombUnits::Vec3<CoulombUnits::Length> x;
    };
    struct FourMomentum {
        CoulombUnits::Energy e;
        CoulombUnits::Vec3<CoulombUnits::Momentum> p;
    };

  protected:
    unsigned int m_nSteps;
    CoulombUnits::Time m_StepSize;

  private:
    static Int_t TotalCharge(std::list<Particle> *);
    static TGraph *m_gImpBNPart;

    unsigned long NearestInterval(const CoulombUnits::Vec3<CoulombUnits::Length> &tPos, CoulombUnits::Time tTime,
            std::vector<SpaceTimePoint> &tTrajectory, unsigned long tLastBestI, CoulombUnits::Length &tInterval, bool debug = false);
    CoulombUnits::Time NearestPreFreezeoutIntervalTime(const CoulombUnits::Vec3<CoulombUnits::Length> &tX0,
            const CoulombUnits::Vec3<CoulombUnits::Speed> &tV, const CoulombUnits::Vec3<CoulombUnits::Length> &tY,
            CoulombUnits::Time tT0, CoulombUnits::Time tU, bool debug = false);
    // Force q E + q v x B on the particle, with its current velocity v. If aMagneticForce is given, it is set to
    // q c B (a force), so that the force for another velocity v' is the result + (v' - v) / c x q c B.
    CoulombUnits::Vec3<CoulombUnits::Force> CalculateTotalForce(std::list<Particle>::iterator &, std::list<Particle> *,
            int ** tBestIs, CoulombUnits::Time tCurrentTime, Int_t iStep,
            CoulombUnits::Vec3<CoulombUnits::Force> *aMagneticForce = nullptr);
    // Vector potential (as phi v / c, an energy, like phi) and scalar potential of unit charges, at the particle
    std::pair<CoulombUnits::Vec3<CoulombUnits::Energy>, CoulombUnits::Energy> CalculatePotentials(std::list<Particle>::iterator &,
            std::list<Particle> *, int ** tBestIs, CoulombUnits::Time tCurrentTime, Int_t iStep);
    bool RetardedSource(std::list<Particle>::iterator &tPartIter, std::list<Particle>::iterator &tPartOtherIter,
            const CoulombUnits::Vec3<CoulombUnits::Length> &tPos, CoulombUnits::Time tTime, int **tBestIs, Int_t iStep,
            CoulombUnits::Vec3<CoulombUnits::Length> &tOtherPos, CoulombUnits::Vec3<CoulombUnits::Speed> &tOtherVel,
            CoulombUnits::Vec3<CoulombUnits::Acceleration> &tOtherAcc, bool debug = false);
    // Close pairs: pairs of charged particles whose encounter during the current time step is too short for the step
    // (see FindCloseGroups()). Their mutual force is left out of CalculateTotalForce() and integrated with sub-steps.
    struct CloseGroup {
        std::vector<std::list<Particle>::iterator> members;
        int substeps;
    };
    std::vector<CloseGroup> FindCloseGroups(std::list<Particle> *, Int_t iStep);
    bool IsClosePair(Int_t aEid1, Int_t aEid2) const;
    // Advances a group over one time step: half kick with the far force at the old positions, sub-steps with the
    // forces within the group; the half kick with the far force at the new positions follows in Apply()
    void IntegrateCloseGroup(const CloseGroup &, const CoulombUnits::Vec3<CoulombUnits::Force> *tFarForces,
            CoulombUnits::Time tCurrentTime);
    // Force on a charge moving with aTestVelocity from a source moving uniformly with aSourceVelocity, at the
    // simultaneous relative position aR = x(test) - x(source); aCharges = q(test) q(source)
    static CoulombUnits::Vec3<CoulombUnits::Force> UniformMotionForce(Int_t aCharges, const CoulombUnits::Vec3<CoulombUnits::Length> &aR,
            const CoulombUnits::Vec3<CoulombUnits::Speed> &aSourceVelocity, const CoulombUnits::Vec3<CoulombUnits::Speed> &aTestVelocity);
    static bool IsSpectator(std::list<Particle>::iterator &);
    static bool HasFiniteKinematics(const Particle &);
    // Minkowski interval c^2 (t1 - t2)^2 - |x1 - x2|^2
    static CoulombUnits::Area Interval(const CoulombUnits::Vec3<CoulombUnits::Length> &tX1, CoulombUnits::Time tT1,
            const CoulombUnits::Vec3<CoulombUnits::Length> &tX2, CoulombUnits::Time tT2, bool debug = false);
    static CoulombUnits::Area Interval(const CoulombUnits::Vec3<CoulombUnits::Length> &tX1, CoulombUnits::Time tT1,
            const SpaceTimePoint &tPoint, bool debug = false);
    Int_t AddSpectators(std::list<Particle> *, CoulombUnits::Time tEarliestTime);
    CoulombUnits::Length SampleImpB(std::list<Particle> *);

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

  unsigned int *m_FirstStepExists;
  // Creation and decay times of the particles (infinity for particles that do not decay), indexed by eid
  CoulombUnits::Time *m_CreationTimes;
  CoulombUnits::Time *m_DecayTimes;
  unsigned int *m_LastStepExists;

  // Large pair contributions for every step, (eid, eid of the other particle) -> vector
  template<typename Q>
  using PairTerms = std::vector< std::map< std::pair< int, int>, CoulombUnits::Vec3<Q> > >;
  PairTerms<CoulombUnits::Force> *m_termsCoul;
  PairTerms<CoulombUnits::Force> *m_termsVelo;
  PairTerms<CoulombUnits::Force> *m_termsAcce;
  PairTerms<CoulombUnits::Force> *m_termsMagn;
  PairTerms<CoulombUnits::Length> *m_Rs;
  PairTerms<CoulombUnits::Dimensionless> *m_us;
  PairTerms<CoulombUnits::Speed> *m_vs;
  std::map<int, int> *m_pids;

  // Recorded trajectories, four-momenta, velocities and accelerations at the time steps, indexed by eid
  std::vector<SpaceTimePoint> *m_coordinates;
  std::vector<FourMomentum> *m_fourmomenta;
  std::vector<CoulombUnits::Vec3<CoulombUnits::Speed> > *m_velocities;
  std::vector<CoulombUnits::Vec3<CoulombUnits::Acceleration> > *m_accelerations;

   void PrintProblem(Int_t eid, Int_t step, TString indent, Int_t depth = 0);

  bool mKeepTrajectories;
  std::vector<std::vector<SpaceTimePoint> > mKeptCoordinates;
  std::vector<std::vector<FourMomentum> > mKeptFourMomenta;

  // Settings, see ReadParameters()
  bool m_UseSpectators;         // add the spectator nucleons of a Glauber MC event
  bool m_ExtrapolateSpectators; // extrapolate spectators back in time before their creation when needed
  bool m_ExtrapolateFireball;   // extrapolate primary fireball particles back in time before freeze-out when needed
  bool m_FindProblem;           // debugging printouts, and storing of large pair contributions for every step
  CoulombUnits::Force m_LargeTerm; // with m_FindProblem, pair contributions with a Coulomb term larger than this
                                   // are stored, so that a particle getting faster than light can be traced back
  int m_DebugStep;              // time step, eid and eid of the other particle for detailed debugging printouts
  int m_DebugEid;
  int m_DebugOtherEid;
  bool m_Substeps;              // integrate the forces within close pairs with sub-steps of the time step
  double m_SubstepAccuracy;     // a sub-step is this fraction of the encounter time of a close pair
  int m_SubstepMax;             // maximum number of sub-steps in a time step
  CoulombUnits::Length m_SubstepDistance; // pairs further apart than this during the time step are not close pairs

  // Close pairs of the current time step, (smaller eid, larger eid)
  std::set<std::pair<Int_t, Int_t> > m_ClosePairs;

};

#endif
