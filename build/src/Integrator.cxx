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

#include <stdio.h>
#include <fstream>
#include <TDatime.h>
#include <TString.h>
#include "Integrator.h"
#include <TMath.h>
#include "THGlobal.h"

// User-End Models
#include "Model_KrakowSFO.h"
#include "Model_BlastWave.h"
#include "Model_BWA.h"
#include "Model_Lhyquid3D.h"
#include "Model_Lhyquid2DBI.h"
#include "Model_SR.h"
/**********************************************************************************************
 * // [2] Include your Model_* class header file
 * #include "Model_Example.h"
 **********************************************************************************************/

using namespace std;
using namespace TMath;

extern void     AddLogEntry(const char* aEntry);
extern void     CopyINIFile();
extern TString  sEventDIR;
extern int      sModel;
extern int      sSeed;

Integrator::Integrator(ParticleDB *mDB)
: mNSamples(0), mRandom(0), mFOModel(0), mDB(mDB)
{
}

Integrator::Integrator(int aNSamples, ParticleDB *mDB)
: mNSamples(aNSamples), mDB(mDB)
{ 
  mRandom = new TRandom2();
#ifdef _ROOT_4_
  mRandom->SetSeed2(41321 + sSeed, 8457);
#else
  mRandom->SetSeed(41321 + sSeed);
#endif

  switch (sModel) {
// HYDRO INSPIRED
    case 0:
      mFOModel = new Model_KrakowSFO(mRandom);
      break;
    case 1:
      mFOModel = new Model_BlastWave(mRandom);
      break;
    case 2:
    case 3:
    case 4:
    case 5:
    case 6:
      mFOModel = new Model_BWA(mRandom);
      break;
// HYDRO BASED
    case 10:
      mFOModel = new Model_Lhyquid3D(mRandom);
      break;
    case 11:
      mFOModel = new Model_Lhyquid2DBI(mRandom);
      break;
 // USER DEFINED
    case 99:
      Model_SR *tModel_SR = new Model_SR(mRandom);
      mFOModel = tModel_SR;
      double tProtonsReq = tModel_SR->GetProtonsReq();
      if (tProtonsReq > 0) {
          AdjustRadius(tModel_SR, tProtonsReq);
      }
      break;
/**********************************************************************************************
 * // [3] add a case for your Model_*
 * case 99:
 *   mFOModel = new Model_Example(mRandom);
 *   break;
 **********************************************************************************************/
  };
  CopyINIFile();
}

Integrator::~Integrator()
{
  delete mRandom;
  delete mFOModel;
};

Model* Integrator::GetModel()
{
  return mFOModel;
}

void Integrator::GenerateParticles(ParticleType* aPartType, int aPartCount, list<Particle>* aParticles)
{
  int	   tIter = 0;
  double   tFMax;
  double   tVal;
  double   tValTest;
  Particle* tParticle; 
  
  PRINT_DEBUG_3("Integrator::GenerateParticles\t "<< aPartType->GetName() <<" B:"<< aPartType->GetBarionN() <<" I3:"<< aPartType->GetI3() <<" S:"<< aPartType->GetStrangeN()<<" C:"<< aPartType->GetCharmN());
  tFMax = aPartType->GetMaxIntegrand();
  while (tIter < aPartCount) {
    tVal      = mFOModel->GetIntegrand(aPartType,true,false).first;
    tValTest  = mRandom->Rndm() * tFMax;

 //   cout<<pdg<<"  tFMax : "<<tFMax<<"  tVal : "<<tVal<<"  tValTest : "<<tValTest<<" tIter :  "<<tIter<<"  aPartCount :  "<<aPartCount<<endl;
    if (tValTest<tVal) {
    //    cout<<"przeszedl test !!!"<<endl;
	tParticle = new Particle(aPartType);
	mFOModel->SetParticlePX(tParticle);

//	Float_t masa = sqrt(pow(tParticle->e,2) - pow(tParticle->px,2) - pow(tParticle->py,2) - pow(tParticle->pz,2));
//	if (Abs(pdg) == 2224 || Abs(pdg) == 2214 || Abs(pdg) == 2114 || Abs(pdg) == 1114)  cout<<"masa delty : "<<masa<<"  ";
	aParticles->push_back(*tParticle);
	tIter++;
	delete tParticle;
    }
//    cout<<"przeszedl test"<<endl;
  }
//    cout<<"przeszedl test 2 "<<endl;

}

void Integrator::SetMultiplicities()
{
  // Make or read table with probabilities
  ifstream tFileIn;
  ofstream tFileOut;
  TDatime  tDate;
  char     tBuff[2*kFileNameMaxChar];
  char     tMultiName[kFileNameMaxChar];
  char     tPart[100];
  double   tMaxInt;
  double   tMulti;

  sprintf(tMultiName,"%sfmultiplicity_%s.txt",sEventDIR.Data(),mFOModel->GetHash());

  tFileIn.open(tMultiName);
  if ((tFileIn) && (tFileIn.is_open())) {
    PRINT_DEBUG_1("<Integrator::SetMultiplicities>\tReading Max Integrand and Multiplicity values from " << tMultiName);
    while (!tFileIn.eof()) {
      *tPart = 0;
      tFileIn >> tPart;
      if(!(*tPart) || (*tPart == '#')) {
	tFileIn.getline(tPart,100);
	continue;
      }
      tFileIn >> tMaxInt >> tMulti;
      mDB->GetParticleType(tPart)->SetMaxIntegrand(tMaxInt);
      mDB->GetParticleType(tPart)->SetMultiplicity(tMulti);
      PRINT_DEBUG_2("\t"<<tPart << " " << tMaxInt << " " << tMulti);
    }
    tFileIn.close();
    sprintf(tBuff,"[input]\t%s",tMultiName);
    AddLogEntry(tBuff);
  } else {
    char tTempName[kFileNameMaxChar];
    sprintf(tTempName,"%sfmultiplicity_%s.tmp",sEventDIR.Data(),mFOModel->GetHash());
    tFileOut.open(tTempName);
    if ((tFileOut) && (tFileOut.is_open())) {
      PRINT_DEBUG_1("<Integrator::SetMultiplicities>\tMax Integrand and Multiplicity file " << tMultiName << " not found.");
      tDate.Set();
      PRINT_MESSAGE("["<<tDate.AsSQLString() << "]\tCalculating Max Integrand and Multiplicity");
      tFileOut << mFOModel->GetDescription();
      tFileOut << endl;
      tFileOut << "# Particle name\tMax integrand\tMultiplicity" << endl;
      for(int tIter = 0; tIter < mDB->GetParticleTypeCount(); tIter++) {
        tMulti  = Integrate(mDB->GetParticleType(tIter));
        tMaxInt = mDB->GetParticleType(tIter)->GetMaxIntegrand();
        tFileOut << mDB->GetParticleType(tIter)->GetName() << '\t' << tMaxInt << '\t' << tMulti << endl;
	cout << "\r\tparticle ("<< (tIter + 1) << "/" << mDB->GetParticleTypeCount() << "): " << mDB->GetParticleType(tIter)->GetName();
	cout.flush();
      }
      cout << endl;
      tFileOut.close();
      if(rename(tTempName,tMultiName) != 0) {
        PRINT_MESSAGE("<Integrator::SetMultiplicities>\tUnable to rename temp file to " << tMultiName);
        exit(_ERROR_GENERAL_FILE_NOT_FOUND_);
      }
    } else {
      PRINT_MESSAGE("<Integrator::SetMultiplicities>\tUnable to create file " << tMultiName);
      exit(_ERROR_GENERAL_FILE_NOT_FOUND_);
    }
    sprintf(tBuff,"[output]\t%s\tfmultiplicity_%s.txt",sEventDIR.Data(),mFOModel->GetHash());
    AddLogEntry(tBuff);
  }
}

void Integrator::Randomize()
{
  TDatime tDate;
#ifdef _ROOT_4_
  mRandom->SetSeed2(tDate.Get(), (tDate.Get() % 11) * 7 + (tDate.Get() / 7));
#else
  mRandom->SetSeed(tDate.Get());
#endif
}

double Integrator::Integrate(ParticleType* aPartType)
{
  double tMaxInt;
  double tMulti;
  double tVolume;
  double tVal;
  int    tIter;
  
  tMaxInt = 0.0;
  tMulti  = 0.0;

//  mRandom->SetSeed(41321);

  if (strcmp(aPartType->GetName(),"Dl1232mnb")==0 ||
      strcmp(aPartType->GetName(),"Dl1232zrb")==0 ||
      strcmp(aPartType->GetName(),"Dl1232plb")==0 ||
      strcmp(aPartType->GetName(),"Dl1232ppb")==0 ||
      strcmp(aPartType->GetName(),"Dl1232min")==0 ||
      strcmp(aPartType->GetName(),"Dl1232zer")==0 ||
      strcmp(aPartType->GetName(),"Dl1232plu")==0 ||
      strcmp(aPartType->GetName(),"Dl1232plp")==0)
  {

    int oversamp = 10;

    for (tIter = 0; tIter < oversamp*mNSamples; tIter++) 
    {
      std::pair<double, double> tIntegrands = mFOModel->GetIntegrand(aPartType,true);
      tVal = tIntegrands.first;
      if (tVal>tMaxInt)
        tMaxInt = tVal;

      tMulti += tVal;
      
    }
    tMulti *= mFOModel->GetHyperCubeVolume() / (1.0 * oversamp * mNSamples);
    aPartType->SetMaxIntegrand(tMaxInt);
    aPartType->SetMultiplicity(tMulti);
  }
  else if(strcmp(aPartType->GetName(),"Lm2350zer")==0 ||
          strcmp(aPartType->GetName(),"Lm2350zrb")==0 ||
          strcmp(aPartType->GetName(),"Lm2110zer")==0 ||
          strcmp(aPartType->GetName(),"Lm2110zrb")==0 ||
          strcmp(aPartType->GetName(),"Lm2100zer")==0 ||
          strcmp(aPartType->GetName(),"Lm2100zrb")==0 ||
          strcmp(aPartType->GetName(),"Lm1890zer")==0 ||
          strcmp(aPartType->GetName(),"Lm1890zrb")==0 ||
          strcmp(aPartType->GetName(),"Lm1830zer")==0 ||
          strcmp(aPartType->GetName(),"Lm1830zrb")==0 ||
          strcmp(aPartType->GetName(),"Lm1820zer")==0 ||
          strcmp(aPartType->GetName(),"Lm1820zrb")==0 ||
          strcmp(aPartType->GetName(),"Lm1810zer")==0 ||
          strcmp(aPartType->GetName(),"Lm1810zrb")==0 ||
          strcmp(aPartType->GetName(),"Lm1800zer")==0 ||
          strcmp(aPartType->GetName(),"Lm1800zrb")==0 ||
          strcmp(aPartType->GetName(),"Lm1690zer")==0 ||
          strcmp(aPartType->GetName(),"Lm1690zrb")==0 ||
          strcmp(aPartType->GetName(),"Lm1670zer")==0 ||
          strcmp(aPartType->GetName(),"Lm1670zrb")==0 ||
          strcmp(aPartType->GetName(),"Lm1600zer")==0 ||
          strcmp(aPartType->GetName(),"Lm1600zrb")==0 ||
          strcmp(aPartType->GetName(),"Lm1520zer")==0 ||
          strcmp(aPartType->GetName(),"Lm1520zrb")==0 ||
          strcmp(aPartType->GetName(),"Lm1405zer")==0 ||
          strcmp(aPartType->GetName(),"Lm1405zrb")==0 ||
          strcmp(aPartType->GetName(),"Lm1115zer")==0 ||
          strcmp(aPartType->GetName(),"Lm1115zrb")==0 ) //new hack - JK
  {
    int oversamp = 10;
    
    for (tIter = 0; tIter < oversamp*mNSamples; tIter++) 
    {
      tVal = mFOModel->GetIntegrand(aPartType,true).first;

      if (tVal>tMaxInt) 
        tMaxInt = tVal;

      tMulti += tVal;
    }
    tMulti *= mFOModel->GetHyperCubeVolume() / (1.0 * oversamp * mNSamples);
    aPartType->SetMaxIntegrand(tMaxInt);
    aPartType->SetMultiplicity(tMulti);
  }
  else
  {
    /// END_HACK_RR
    for (tIter = 0; tIter < mNSamples; tIter++) 
    {
      std::pair<double, double> tIntegrands = mFOModel->GetIntegrand(aPartType,true);
      tVal = tIntegrands.first;

      if (tVal>tMaxInt) 
        tMaxInt = tVal;

      tMulti += tVal;

      tVal = tIntegrands.second;
      tVolume += tVal;
    }
    tMulti *= mFOModel->GetHyperCubeVolume() / (1.0 * mNSamples);
    tVolume *= 2 * TMath::Pi() * TMath::Pi() * mFOModel->GetHyperCubeVolume() / (1.0 * mNSamples) * TMath::Power(kHbarC, 3)/3;
    aPartType->SetMaxIntegrand(tMaxInt);
    aPartType->SetMultiplicity(tMulti);
  }

  return tMulti;
}

void Integrator::AdjustRadius(Model_SR * aModel, double aRequiredProtons) {
  static const Double_t tEps = 0.05;
  ParticleType *tProtonType = this->mDB->GetParticleType("pr0938plu");
  Double_t tMult = this->Integrate(tProtonType);
  Int_t tIteration = 1;
  while (TMath::Abs(tMult - aRequiredProtons) > tEps) {
    Double_t tOldR = aModel->GetR();
    Double_t tF = TMath::Power(aRequiredProtons / tMult, 1./3);
    tF = (tF - 1)/tIteration + 1;
    Double_t tNewR = tOldR * tF;
    cout << "Iteration " << tIteration << ", required mult. = " << aRequiredProtons << ", actual mult = " << tMult;
    cout << ", changing radius from " << tOldR << " to " << tNewR << endl;
    aModel->SetR(tNewR);
    tMult = this->Integrate(tProtonType);
    tIteration++;
  }
}
