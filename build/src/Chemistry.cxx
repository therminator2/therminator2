/********************************************************************************
 *                                                                              *
 *             THERMINATOR 2: THERMal heavy-IoN generATOR 2                     *
 *                                                                              *
 * Version:                                                                     *
 *      Release, 2.0.3, 1 February 2011                                         *
 *                                                                              *
 * Authors:                                                                     *
 *      Mikolaj Chojnacki   (Mikolaj.Chojnacki@ifj.edu.pl)                      *
 *      Adam Kisiel         (kisiel@if.pw.edu.pl)                               *
 *      Wojciech Broniowski (Wojciech.Broniowski@ifj.edu.pl)                    *
 *      Wojciech Florkowski (Wojciech.Florkowski@ifj.edu.pl)                    *
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

#include <TMath.h>
#include "Chemistry.h"
#include "THGlobal.h"

using namespace TMath;

Chemistry::Chemistry()
: mMuBConst(0.0),	mMuIConst(0.0),		mMuQConst(0.0),		mMuSConst(0.0),		mMuCConst(0.0),
  mMuBVar(0),		mMuIVar(0),		mMuQVar(0),		mMuSVar(0),		mMuCVar(0),
  mLambdaQConst(1.0),	mLambdaIConst(1.0),	mLambdaSConst(1.0),	mLambdaCConst(1.0),
  mGammaQConst(1.0),				mGammaSConst(1.0),	mGammaCConst(1.0),
  mLambdaQVar(0),	mLambdaIVar(0),		mLambdaSVar(0),		mLambdaCVar(0),
  mGammaQVar(0),				mGammaSVar(0),		mGammaCVar(0),
  mChemistryType(-1),   mUseMuQ(false)
{
}

Chemistry::~Chemistry()
{
  delete mMuBVar;
  delete mMuIVar;
  delete mMuQVar;
  delete mMuSVar;
  delete mMuCVar;
  delete mLambdaQVar;
  delete mLambdaIVar;
  delete mLambdaSVar;
  delete mLambdaCVar;
  delete mGammaQVar;
  delete mGammaSVar;
  delete mGammaCVar;
}

double Chemistry::GetChemicalPotential(ParticleType* aPartType) const
{
  if (mUseMuQ) {
    return (aPartType->GetBarionN()  * mMuBConst +
	    aPartType->GetCharge()   * mMuQConst +
	    aPartType->GetStrangeN() * mMuSConst +
	    aPartType->GetCharmN()   * mMuCConst );
  } else {
    return (aPartType->GetBarionN()  * mMuBConst +
	    aPartType->GetI3()       * mMuIConst +
	    aPartType->GetStrangeN() * mMuSConst +
	    aPartType->GetCharmN()   * mMuCConst );
  }
}

double Chemistry::GetChemicalPotential(ParticleType* aPartType, double aX, double aY, double aZ) const
{
  if (mUseMuQ) {
    return (aPartType->GetBarionN()  * mMuBVar->Interpolate(aX, aY, aZ) +
	    aPartType->GetCharge()   * mMuQVar->Interpolate(aX, aY, aZ) +
	    aPartType->GetStrangeN() * mMuSVar->Interpolate(aX, aY, aZ) +
	    aPartType->GetCharmN()   * mMuCVar->Interpolate(aX, aY, aZ) );
  } else {
    return (aPartType->GetBarionN()  * mMuBVar->Interpolate(aX, aY, aZ) +
	    aPartType->GetI3()       * mMuIVar->Interpolate(aX, aY, aZ) +
	    aPartType->GetStrangeN() * mMuSVar->Interpolate(aX, aY, aZ) +
	    aPartType->GetCharmN()   * mMuCVar->Interpolate(aX, aY, aZ) );
  }
}

double Chemistry::GetFugacity(ParticleType* aPartType) const
{
  return mLambdaIConst
	  * Power(mGammaQConst,  aPartType->GetNumberQ() + aPartType->GetNumberAQ())
	  * Power(mLambdaQConst, aPartType->GetNumberQ() - aPartType->GetNumberAQ())
	  * Power(mGammaSConst,  aPartType->GetNumberS() + aPartType->GetNumberAS())
	  * Power(mLambdaSConst, aPartType->GetNumberS() - aPartType->GetNumberAS())
	  * Power(mGammaCConst,  aPartType->GetNumberC() + aPartType->GetNumberAC())
	  * Power(mLambdaCConst, aPartType->GetNumberC() - aPartType->GetNumberAC());
}

double Chemistry::GetFugacity(ParticleType* aPartType, double aX, double aY, double aZ) const
{
  return mLambdaIVar->Interpolate(aX, aY, aZ)
	  * Power(mGammaQVar ->Interpolate(aX, aY, aZ), aPartType->GetNumberQ() + aPartType->GetNumberAQ())
	  * Power(mLambdaQVar->Interpolate(aX, aY, aZ), aPartType->GetNumberQ() - aPartType->GetNumberAQ())
	  * Power(mGammaSVar ->Interpolate(aX, aY, aZ), aPartType->GetNumberS() + aPartType->GetNumberAS())
	  * Power(mLambdaSVar->Interpolate(aX, aY, aZ), aPartType->GetNumberS() - aPartType->GetNumberAS())
	  * Power(mGammaCVar ->Interpolate(aX, aY, aZ), aPartType->GetNumberC() + aPartType->GetNumberAC())
	  * Power(mLambdaCVar->Interpolate(aX, aY, aZ), aPartType->GetNumberC() - aPartType->GetNumberAC());
}

int    Chemistry::GetChemistryType() const	{ return mChemistryType; }
bool   Chemistry::GetUseMuQ() const	{ return mUseMuQ; }
double Chemistry::GetMuB() const		{ return mMuBConst; }
double Chemistry::GetMuI() const		{ return mMuIConst; }
double Chemistry::GetMuQ() const		{ return mMuQConst; }
double Chemistry::GetMuS() const		{ return mMuSConst; }
double Chemistry::GetMuC() const		{ return mMuCConst; }
double Chemistry::GetLambdaQ() const		{ return mLambdaQConst; }
double Chemistry::GetLambdaI() const		{ return mLambdaIConst; }
double Chemistry::GetLambdaS() const		{ return mLambdaSConst; }
double Chemistry::GetLambdaC() const		{ return mLambdaCConst; }
double Chemistry::GetGammaQ() const		{ return mGammaQConst; }
double Chemistry::GetGammaS() const		{ return mGammaSConst; }
double Chemistry::GetGammaC() const		{ return mGammaCConst; }

void Chemistry::SetChemistry(double aMuB, double aMuI, double aMuS, double aMuC)
{
  mChemistryType = 0;
  mUseMuQ = false;
  mMuBConst = aMuB; mMuIConst = aMuI; mMuSConst = aMuS; mMuCConst = aMuC;
}
void Chemistry::SetChemistry(Vector3D* aMuB, Vector3D* aMuI, Vector3D* aMuS, Vector3D* aMuC)
{
  mChemistryType = 1;
  mUseMuQ = false;
  mMuBVar = aMuB; mMuIVar = aMuI; mMuSVar = aMuS; ; mMuCVar = aMuC;
}
void Chemistry::SetChemistryQ(double aMuB, double aMuQ, double aMuS, double aMuC)
{
  mChemistryType = 0;
  mUseMuQ = true;
  mMuBConst = aMuB; mMuQConst = aMuQ; mMuSConst = aMuS; mMuCConst = aMuC;
  PRINT_MESSAGE("MuQ set and it will be used")
}
void Chemistry::SetChemistryQ(Vector3D* aMuB, Vector3D* aMuQ, Vector3D* aMuS, Vector3D* aMuC)
{
  mChemistryType = 1;
  mUseMuQ = true;
  mMuBVar = aMuB; mMuQVar = aMuQ; mMuSVar = aMuS; ; mMuCVar = aMuC;
  PRINT_MESSAGE("MuQ set and it will be used")
}
void Chemistry::SetChemistry(double aLambdaQ, double aLambdaI, double aLambdaS, double aLambdaC, double aGammaQ, double aGammaS,  double aGammaC)
{
  mChemistryType = 2;
  mLambdaQConst = aLambdaQ; mLambdaIConst = aLambdaI; mLambdaSConst = aLambdaS; mLambdaCConst = aLambdaC; 
  mGammaQConst  = aGammaQ;                            mGammaSConst  = aGammaS;  mGammaCConst  = aGammaC;

}
void Chemistry::SetGammas(double aGammaQ, double aGammaS,  double aGammaC)
{
  mGammaQConst  = aGammaQ;                            mGammaSConst  = aGammaS;  mGammaCConst  = aGammaC;
}
void Chemistry::SetChemistry(Vector3D* aLambdaQ, Vector3D* aLambdaI, Vector3D* aLambdaS, Vector3D* aLambdaC, Vector3D* aGammaQ, Vector3D* aGammaS, Vector3D* aGammaC)
{
  mChemistryType = 3;
  mLambdaQVar = aLambdaQ; mLambdaIVar = aLambdaI; mLambdaSVar = aLambdaS; mLambdaCVar = aLambdaC;
  mGammaQVar  = aGammaQ;                          mGammaCVar  = aGammaC;  mGammaCVar  = aGammaC;
  
}
