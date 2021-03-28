/********************************************************************************
 *                                                                              *
 *             THERMINATOR 2: THERMal heavy-IoN generATOR 2                     *
 *                                                                              *
 * Version:                                                                     *
 *      Release, 2.0.3, 1 February 2011                                         *
 *                                                                              *
 * Authors:                                                                     *
 *      Radoslaw Ryblewski   (Radoslaw.Ryblewski@ifj.edu.pl)                    *
 *                                                                              *
 *                                                                              *
 *                                                                              *
 *                                                                              *
 * Project homepage:                                                            *
 *      http://therminator2.ifj.edu.pl/                                         *
 *                                                                              *
 * For the detailed description of the program and further references           *
 * to the description of the model please refer to                              *
 * http://arxiv.org/abs/1102.0273                                               *
 *                                                                              *
 * This code can be freely used and redistributed. However if you decide to     *
 * make modifications to the code, please, inform the authors.                  *
 * Any publication of results obtained using this code must include the         *
 * reference to arXiv:1102.0273 and the published version of it, when           *
 * available.                                                                   *
 *                                                                              *
 ********************************************************************************/

#ifndef _TH2_MODEL_SR_H_
  #define _TH2_MODEL_SR_H_

#include "Model.h"
#include "Thermodynamics.h"

class Model_SR
: public Model
{
  public:
    Model_SR();
    Model_SR(TRandom2* aRandom);
    ~Model_SR();
    
    double GetIntegrand(ParticleType* aPartType, bool);
    void   AddParameterBranch(TTree* aTree);
    
  protected:
      /*
      double mRMax;
    double mVR;  
    double mh;
    double mA0;
    double mT0;
    */
    double mR;
    double mA;  
    double mT0;
    double mH;
    double mDel;
    double mEps;
    double mGammaS; 

    Thermodynamics* mThermo;
    
  private:
    void Description();
    void ReadParameters();
    bool PointInGrid(Vector3D *aV, double aX, double aY, double aZ, unsigned int &aI, unsigned int &aJ, unsigned int &aK);
    void GetGridAxisPoint(double aMin, double aMax, double aPts, double aVal, unsigned int &aRoundN, double &aRoundError);
};

#endif

/*! @file Model_SR.h
 * @brief Definition of Model_SR class. Simmens-Rasmussen freeze-out Model.
 */
/*! @class Model_SR
 *  
 *  NEEDS THE UPDATE (RR)
 * 
 */
