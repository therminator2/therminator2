#ifndef RADIUS_ADJUSTER
#define RADIUS_ADJUSTER

#include "Integrator.h"
#include "Model_SR.h"

class RadiusAdjuster {
  public:
    virtual ~RadiusAdjuster();
    virtual void Apply(Integrator *, Model_SR *, double);

  protected:
    void ReadParameters();
    unsigned int m_nSteps;
    Double_t m_StepSize;

  private:
    Integrator *mInteg;
    Model_SR *mModel;
    static const Double_t mEps;

};

#endif
