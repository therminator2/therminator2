#ifndef CANONICAL_SUPPRESSION_AFTERBURNER
#define CANONICAL_SUPPRESSION_AFTERBURNER

#include "AbstractAfterburner.h"
#include "TLorentzVector.h"
#include "utils/Grid2D.h"

#include "TFile.h"
#include "TGraph.h"
#include "TProfile.h"

class CanonicalSuppressionAfterburner : public AbstractAfterburner {
  public:
    CanonicalSuppressionAfterburner(Double_t);
    virtual ~CanonicalSuppressionAfterburner();
    virtual void Apply(Event *);

  protected:
    void ReadParameters();
    Double_t m_Rc;

  private:
    Utils::Grid2D<Float_t> grid;
    std::map<Int_t, std::vector<UInt_t>> eidToIs;
    std::map<Int_t, std::vector<UInt_t>> eidToJs;
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
