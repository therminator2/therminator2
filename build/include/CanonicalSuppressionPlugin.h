#ifndef CANONICAL_SUPPRESSION_PLUGIN
#define CANONICAL_SUPPRESSION_PLUGIN

#include "AbstractPlugin.h"
#include "TLorentzVector.h"
#include "utils/Grid2D.h"
#include "utils/Option.h"

#include "TFile.h"
#include "TGraph.h"
#include "TProfile.h"

class CanonicalSuppressionPlugin : public AbstractPlugin {
  public:
    CanonicalSuppressionPlugin();
    virtual ~CanonicalSuppressionPlugin();
    virtual void Apply(Event *);

  protected:
    void ReadParameters();
    Double_t m_Rc;

  private:
    Utils::Grid2D<Float_t> grid;
    std::map<Int_t, std::vector<UInt_t>> eidToIs;
    std::map<Int_t, std::vector<UInt_t>> eidToJs;

    const Utils::Option<Particle> &&findParticleByEid(Event *tEvent, Int_t tEid);
    bool brokenMultistrange(Event *tEvent, Int_t tEid);
    bool brokenMultistrange(Particle &);
};

/*
// Maybe still will be used if the lambda-in-lambda approach will not work, but then eidToIs and eidToJs
// still have to be added to the predicate class:
class ParticleBreaksConservationLaw {
  private:
    Utils::Grid2D<Float_t> grid;
  public:
    ParticleBreaksConservationLaw(Utils::Grid2D<Float_t> grid) : grid(grid) {}
    bool operator()(Particle &part);
};
*/
#endif
