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

#include <sstream>
#include <TMath.h>
#include "THGlobal.h"
#include "Configurator.h"
#include "Parser.h"
#include "Hypersurface_Library.h"
#include "Model_SR.h"
#include "Vector3D.h"

using namespace TMath;
using namespace std;

extern Configurator* sMainConfig;
extern TString	sModelINI;
extern TString	sEventDIR;
extern TString	sTimeStamp;
extern int	sModel;
extern int	sRandomize;
extern int	sIntegrateSample;

Model_SR::Model_SR()
: Model(),
  mThermo(0),
  mR(0.0),  mH(0.0), mA(0.0), mT0(0.0), mGammaS(0.0)
{
}

Model_SR::Model_SR(TRandom2* aRandom)
: Model(aRandom) 
{
  mName = "SR";
  mThermo = new Thermodynamics();
  ReadParameters();
  Description();
  mHyperCube = mR * 2*TMath::Pi()*Pi() * 2*TMath::Pi()*Pi();
}

Model_SR::~Model_SR()
{
  delete mThermo;
}

double Model_SR::GetIntegrand(ParticleType* aPartType, bool finiteWidth)
{
  double dSIGMAdotP, Integrand;
  double Spin, Gs, Statistics;
  double Temp, Mu;
  double R, Phi, Theta;
  double P, PhiP, ThetaP;
  double dPdZet, Ep, kappa;
    double kappau, kappax;
  double Mt, Pt, RapP;

// type of statistics: Bose-Einstein or Fermi-Dirac
  Spin = aPartType->GetSpin();
  Gs	= 2 * Spin + 1;
  Statistics = ( (Spin - static_cast<int>(Spin)) < 0.01 ? -1.0 : +1.0 );
  Temp	= mThermo->GetTemperature();
  Mu = mThermo->GetChemicalPotential(aPartType);
//    mRandom->SetSeed(41321);
  
// generate random spatial position 
  R	    = mR * mRandom->Rndm();
  Phi	= 2.0 * Pi() * mRandom->Rndm();  
  Theta	= Pi() * mRandom->Rndm();  
    
// generate random momentum position 
  {
    double Zet = mRandom->Rndm();
    P	= Zet / (1.0 - Zet);	// 0 <= p <= Infinity
    dPdZet	= 1.0 / ( (1.0 - Zet) * (1.0 - Zet) );
  }
  PhiP	= 2.0 * Pi() * mRandom->Rndm(); 
  ThetaP	= Pi() * mRandom->Rndm();  
  
// other variables

  double spectralFunctionWeight;
  double M;
//  finiteWidth = false;
  GetParticleMass(aPartType, finiteWidth,M,spectralFunctionWeight);

  Ep	= Hypot(M,P);
    // for spherical
//   kappa = Cos(Theta) * Cos(ThetaP) + Sin(Theta) * Sin(ThetaP) * Cos(Phi - PhiP);
    // for spheroidal
     kappau = Cos(Theta) * Cos(ThetaP) * sqrt(1+mDel) + Sin(Theta) * Sin(ThetaP) * Cos(Phi - PhiP)* sqrt(1-mDel);
     kappax = Cos(Theta) * Cos(ThetaP) + Sin(Theta) * Sin(ThetaP) * Cos(Phi - PhiP) * (sqrt(1+mEps) / sqrt(1-mEps));
    


//   Float_t tVR = 0.5;  // TO JEST TYLKO TO TESTU NA STALY PRZYEPLYW !!!!
  Float_t tVR = TanH(mH * R);

  // invariants
  double Lgamma     =  1.0/sqrt(1-(1+mDel*Cos(2*Theta))*tVR*tVR);
  double UdotP      = (Ep - P * tVR * kappau) * Lgamma;
  dSIGMAdotP = (1-mEps) * R*R * Sin(Theta) * (Ep * sqrt(1+mEps) - P * mA * kappax);
  
// disable particle emission back to the hydro region
  if(dSIGMAdotP < 0.0) 
    dSIGMAdotP = 0.0;
// particle X and P coordinates - required to be initiated
  Xt = mT0 + mA * R;
//  Xt = mT0 + tVR * R; /// ta linijka byla chyba zanim dodalismy Hubbla
  Xx = Sqrt(1-mEps) * R * Cos(Phi) * Sin(Theta);
  Xy = Sqrt(1-mEps) * R * Sin(Phi) * Sin(Theta);
  Xz = Sqrt(1+mEps) * R * Cos(Theta);
  
  Pe = Ep;
  Px = P * Cos(PhiP) * Sin(ThetaP);
  Py = P * Sin(PhiP) * Sin(ThetaP);
  Pz = P * Cos(ThetaP);

  Pt  = Hypot(Px,Py);
  Mt  = Hypot(M,Pt);
  RapP = 1./2*TMath::Log((Pe+Pz)/(Pe-Pz));

  Float_t Ekin  = Mt * CosH(RapP);
  Float_t Px0   = Pt * Cos(PhiP);
  Float_t Py0   = Pt * Sin(PhiP);
  Float_t Pz0   = Mt * SinH(RapP);

  Px = Px0;
  Py = Py0;
  Pz = Pz0;
  Pe = Ekin;


  double fugacity 
      = Power(mThermo->GetGammaQ(),  aPartType->GetNumberQ() + aPartType->GetNumberAQ())
      * Power(mThermo->GetGammaS(),  aPartType->GetNumberS() + aPartType->GetNumberAS())
      * Exp(mThermo->GetChemicalPotential(aPartType) / Temp);
  double invFugacity = 1. / fugacity;

  int pdg = aPartType->GetPDGCode();

  double T          = mThermo->GetTemperature();
//  double Mu         = mThermo->GetChemicalPotential(aPartType);
  double Upsilon    = pow(mGammaS, aPartType->GetNumberS()+ aPartType->GetNumberAS()) * Exp(Mu/T);
  double DP = P * P * Sin(ThetaP) * dPdZet / Ep;
  double F  = (Gs / kTwoPi3) / (Exp(UdotP/T)/Upsilon  + Statistics);

  Integrand = F * DP * dSIGMAdotP;

  if (false) {
 // if (TString(aPartType->GetName()).Contains("Ka0492plu")) {
 // if (TString(aPartType->GetName()).Contains("pi0139plu")) {
    cout << fugacity << "\t" << Integrand << "\t" << mThermo->GetChemicalPotential(aPartType) << "\t"
         << UdotP << "\t" << Exp(UdotP/Temp) << "\t"
         << Power(mThermo->GetGammaS(),  aPartType->GetNumberS() + aPartType->GetNumberAS()) << "\t"
         << mThermo->GetGammaS() << "\t" << aPartType->GetNumberS() << "\t" << aPartType->GetNumberAS() << endl;
  }
// integrand
  return  Integrand;
  }

void Model_SR::Description()
{
  ostringstream oss;
  oss << "##################################################"<< endl;
  oss << MODEL_NAME(mName);
  oss << "# - freeze-out Cart. time  : " <<MODEL_PAR_DESC(mT0  * kHbarC,	"[fm]");
  oss << "# - radial size            : " <<MODEL_PAR_DESC(mR   * kHbarC,	"[fm]");
  oss << "# - Hubble velocity        : " <<MODEL_PAR_DESC(mH,		"[c]");
  oss << "# - hypersurface slope (A) : " <<MODEL_PAR_DESC(mA,		"[1]");
    
  oss << "# - delta                  : " <<MODEL_PAR_DESC(mDel,         "[1]");
  oss << "# - epsilon                : " <<MODEL_PAR_DESC(mEps,         "[1]");
  oss << "# - gamma_S                : " <<MODEL_PAR_DESC(mGammaS,		"[1]");
  /*
  oss << "# - freeze-out Cart. time  : " <<MODEL_PAR_DESC(mT0   * kHbarC,	"[fm]");
  oss << "# - radial size            : " <<MODEL_PAR_DESC(mRMax   * kHbarC,	"[fm]");
  oss << "# - radial velocity        : " <<MODEL_PAR_DESC(mVR,		"[c]");
  oss << "# - h const for Hubble expansion : " <<MODEL_PAR_DESC(mh, "");
  oss << "# - A const freeze-out geom : " <<MODEL_PAR_DESC(mA0, "");
  oss << "# - radial velocity        : " <<MODEL_PAR_DESC(mVR,		"[c]");
  */
  oss << "# - freeze-out temperature : " <<MODEL_PAR_DESC(mThermo->GetTemperature() * 1000.0,	"[MeV]");
  if ((mThermo->GetChemistryType() == 0) || (mThermo->GetChemistryType() == 1)) 
  {
    oss << "# - chem. potential Mu_B   : " <<MODEL_PAR_DESC(mThermo->GetMuB() * 1000.0,	"[MeV]");
    oss << "# - chem. potential Mu_I3  : " <<MODEL_PAR_DESC(mThermo->GetMuI() * 1000.0,	"[MeV]");
    oss << "# - chem. potential Mu_S   : " <<MODEL_PAR_DESC(mThermo->GetMuS() * 1000.0,	"[MeV]");
    oss << "# - chem. potential Mu_C   : " <<MODEL_PAR_DESC(mThermo->GetMuC() * 1000.0,	"[MeV]");
  }
  else
  {
    oss << "# - fugacity Lambda_I3     : " <<MODEL_PAR_DESC(mThermo->GetLambdaI(), "[1]");
    oss << "# - fugacity Lambda_Q      : " <<MODEL_PAR_DESC(mThermo->GetLambdaQ(), "[1]");
    oss << "# - fugacity Lambda_S      : " <<MODEL_PAR_DESC(mThermo->GetLambdaS(), "[1]");
    oss << "# - fugacity Lambda_C      : " <<MODEL_PAR_DESC(mThermo->GetLambdaC(), "[1]");   
    oss << "# - fugacity Gamma_Q       : " <<MODEL_PAR_DESC(mThermo->GetGammaQ(),  "[1]");
    oss << "# - fugacity Gamma_S       : " <<MODEL_PAR_DESC(mThermo->GetGammaS(),  "[1]");
    oss << "# - fugacity Gamma_C       : " <<MODEL_PAR_DESC(mThermo->GetGammaC(),  "[1]");   
  }
  oss << "# Parameters hash (CRC32)  : " <<MODEL_PAR_DESC(mHash,		"");
  oss << "# Integration samples      : " <<MODEL_PAR_DESC(sIntegrateSample,	"");
  oss << "# Random seed              : " <<MODEL_PAR_DESC((sRandomize ? "yes" : "no"),"");
  oss << "# Generation date          : " <<sTimeStamp<<" #"<<endl;
  oss << "##################################################"<< endl;
  mDescription = oss.str();
}

void Model_SR::AddParameterBranch(TTree* aTree)  
{
  Model_t_SR tPar;
  /*
  tPar.T0        = mT0  * kHbarC;
  tPar.RMax      = mRMax  * kHbarC;
  tPar.VR        = mVR;
  tPar.h         = mh;
  tPar.A0        = mA0;
  */
  tPar.T0        = mT0  * kHbarC;
  tPar.R         = mR   * kHbarC;
  tPar.H         = mH;
  tPar.A         = mA;

  tPar.Temp      = mThermo->GetTemperature() * 1000.0;
  tPar.Chem      = mThermo->GetChemistryType();
  tPar.MuB       = mThermo->GetMuB() * 1000.0;
  tPar.MuI       = mThermo->GetMuI() * 1000.0;
  tPar.MuS       = mThermo->GetMuS() * 1000.0;
  tPar.MuC       = mThermo->GetMuC() * 1000.0; 
  
  tPar.LambdaQ       = mThermo->GetLambdaQ(); 
  tPar.LambdaI       = mThermo->GetLambdaI(); 
  tPar.LambdaS       = mThermo->GetLambdaS(); 
  tPar.LambdaC       = mThermo->GetLambdaC(); 
  
  tPar.GammaQ       = mThermo->GetGammaQ();  
  tPar.GammaS       = mThermo->GetGammaS(); 
  tPar.GammaC       = mThermo->GetGammaC(); 
 
  cout << "PRINTOUT HERE : " << _MODEL_T_FORMAT_SR_ << endl;
  aTree->Branch(_MODEL_T_BRANCH_, &tPar, _MODEL_T_FORMAT_SR_)->Fill();
}

void Model_SR::ReadParameters()
{
  Hypersurface_Library *tLib;
  Configurator*	tModelParam;
  Parser*	tParser;

  tModelParam = new Configurator;
  tParser     = new Parser(sModelINI.Data());
  tParser->ReadINI(tModelParam);
  delete tParser;

  try {
    mT0  	 = tModelParam->GetParameter("T0").Atof() / kHbarC;		// [GeV^-1]
    /*
    mRMax	 = tModelParam->GetParameter("RMax").Atof() / kHbarC;		// [GeV^-1]
    mVR		 = tModelParam->GetParameter("VR").Atof();			// [c]
    mh		 = tModelParam->GetParameter("h").Atof();			//
    mA0		 = tModelParam->GetParameter("A0").Atof();
    mR    	 = tModelParam->GetParameter("R").Atof() / kHbarC;		// [GeV^-1]
    */
    mDel         = tModelParam->GetParameter("del").Atof();            // [1]
    mEps         = tModelParam->GetParameter("eps").Atof();            // [1]
    mH		 = tModelParam->GetParameter("H").Atof();			// [1]
    mA		 = tModelParam->GetParameter("A").Atof();			// [1]
    mGammaS	 = tModelParam->GetParameter("GammaS").Atof();			// [1]

    //
    mThermo->SetTemperature(tModelParam->GetParameter("Temperature").Atof() * 0.001);	// [GeV]
    if (tModelParam->GetParameter("Chemistry") == "chemical_potential")
    {
      if (tModelParam->HasParameter("MuQ")) {
        mThermo->SetChemistryQ  (tModelParam->GetParameter("MuB").Atof() * 0.001,
			        tModelParam->GetParameter("MuQ").Atof() * 0.001,
			        tModelParam->GetParameter("MuS").Atof() * 0.001,
			        tModelParam->GetParameter("MuC").Atof() * 0.001);	// [GeV]
      } else {
        mThermo->SetChemistry  (tModelParam->GetParameter("MuB").Atof() * 0.001,
			        tModelParam->GetParameter("MuI").Atof() * 0.001,
			        tModelParam->GetParameter("MuS").Atof() * 0.001,
			        tModelParam->GetParameter("MuC").Atof() * 0.001);	// [GeV]
      }
      mThermo->SetGammas(
                            tModelParam->GetParameter("GammaQ").Atof(),
			    tModelParam->GetParameter("GammaS").Atof(),
			    tModelParam->GetParameter("GammaC").Atof()); // [1]
    }
    else if (tModelParam->GetParameter("Chemistry") == "gamma_lambda")
    {
      mThermo->SetChemistry(
				tModelParam->GetParameter("LambdaQ").Atof(),
			    tModelParam->GetParameter("LambdaI").Atof(),
			    tModelParam->GetParameter("LambdaS").Atof(),
			    tModelParam->GetParameter("LambdaC").Atof(),
                            tModelParam->GetParameter("GammaQ").Atof(),
			    tModelParam->GetParameter("GammaS").Atof(),
			    tModelParam->GetParameter("GammaC").Atof()); // [1]

      cout << "HERE PARAMETERS ARE : " << endl;
      cout << tModelParam->GetParameter("LambdaQ").Atof() << endl;
      cout << tModelParam->GetParameter("LambdaI").Atof() << endl;
      cout << tModelParam->GetParameter("LambdaS").Atof() << endl;
      cout << tModelParam->GetParameter("LambdaC").Atof() << endl;
      cout << tModelParam->GetParameter("GammaQ").Atof() << endl;
      cout << tModelParam->GetParameter("GammaS").Atof() << endl;
      cout << tModelParam->GetParameter("GammaC").Atof() << endl;
      cout << mThermo->GetGammaS() << endl;

    }
    else
    {
      PRINT_MESSAGE("\tDid not find one of the necessary model parameters.");
      exit(_ERROR_CONFIG_PARAMETER_NOT_FOUND_);
    }

    TString tRadiusXML = tModelParam->GetParameter("RadiusFile");
    tLib    = new Hypersurface_Library;
    tParser = new Parser(tRadiusXML.Data());
    tParser->ReadXML(tLib);
    delete tParser;
    Vector3D *tVectorR = tLib->GetXMLTag("VECTOR3D", "name","R")->GetXMLVector3D();
    if (
            mH < tVectorR->GetXMin() || mH > tVectorR->GetXMax() ||
            mEps < tVectorR->GetYMin() || mEps > tVectorR->GetYMax() ||
            mDel < tVectorR->GetZMin() || mDel > tVectorR->GetZMax()
       ) {
        cout << "The values of (H, epsilon, delta) = (" << mH << "," << mEps << "," << mDel << ") outside of the grid (" 
            << tVectorR->GetXMin() << ":" << tVectorR->GetXMax() << "," 
            << tVectorR->GetYMin() << ":" << tVectorR->GetYMax() << "," 
            << tVectorR->GetZMin() << ":" << tVectorR->GetZMax() << ")"
            << endl;
        exit(_ERROR_CONFIG_PARAMETER_NOT_FOUND_);
    }
    unsigned int tI,tJ,tK;
    if (PointInGrid(tVectorR, mH, mEps, mDel, tI, tJ, tK)) {
      cout << "(i,j,k) = (" << tI << "," << tJ << "," << tK << ")" << endl;
      mR = tVectorR->operator()(tI, tJ, tK) / kHbarC;;
      cout << "Found in grid R = " << mR * kHbarC << " fm for (H, epsilon, delta) = (" << mH << "," << mEps << "," << mDel << ")" << endl;
    } else {
      mR = tVectorR->Interpolate(mH, mEps, mDel) / kHbarC;;
      cout << "Interpolated in grid R = " << mR * kHbarC << " fm for (H, epsilon, delta) = (" << mH << "," << mEps << "," << mDel << ")" << endl;
    }
  } catch (TString tError) {
    PRINT_MESSAGE("<Model_SR::ReadParameters>\tCaught exception " << tError);
    PRINT_MESSAGE("\tDid not find one of the necessary model parameters.");
    exit(_ERROR_CONFIG_PARAMETER_NOT_FOUND_);
  }
  
// calculate parameter hash 
  ostringstream oss;
  oss << sModel;
  oss << mT0;
  oss << mR;
  oss << mH;
  oss << mA;
  oss << mGammaS;
  oss << mThermo->GetTemperature() ;
  if ((mThermo->GetChemistryType() == 0) || (mThermo->GetChemistryType() == 1))
  {
    oss << mThermo->GetMuB() << mThermo->GetMuI() << mThermo->GetMuS() << mThermo->GetMuC();
  }
  else
  {
    oss << mThermo->GetLambdaQ() << mThermo->GetLambdaI() << mThermo->GetLambdaS() << mThermo->GetLambdaC();
    oss << mThermo->GetGammaQ()                           << mThermo->GetGammaS()  << mThermo->GetGammaC();
  }
  CalculateHash(TString(oss.str()));

// create event subdirectory if needed
  try {
    sEventDIR += tModelParam->GetParameter("EventSubDir");
    CreateEventSubDir();
  } catch (TString tError) {
  }
  
  delete tModelParam;
}

bool Model_SR::PointInGrid(Vector3D *aV, double aX, double aY, double aZ, unsigned int &aI, unsigned int &aJ, unsigned int &aK) {
  double tSmall = 0.001;

  unsigned int tRoundNX, tRoundNY, tRoundNZ;
  double tRoundErrorX, tRoundErrorY, tRoundErrorZ;

  GetGridAxisPoint(aV->GetXMin(), aV->GetXMax(), aV->GetXPts(), aX, tRoundNX, tRoundErrorX);
  GetGridAxisPoint(aV->GetYMin(), aV->GetYMax(), aV->GetYPts(), aY, tRoundNY, tRoundErrorY);
  GetGridAxisPoint(aV->GetZMin(), aV->GetZMax(), aV->GetZPts(), aZ, tRoundNZ, tRoundErrorZ);

  cout << "tSmall = " << tSmall << endl;
  if (tRoundErrorX < tSmall && tRoundErrorY < tSmall && tRoundErrorZ < tSmall) {
    aI = tRoundNX;
    aJ = tRoundNY;
    aK = tRoundNZ;
    return true;
  } else {
    return false;
  }
}

void Model_SR::GetGridAxisPoint(double aMin, double aMax, double aPts, double aVal, unsigned int &aRoundN, double &aRoundError) {
  double tStep = (aMax - aMin)/(aPts-1);
  double tN = (aVal - aMin) / tStep;
  aRoundN = TMath::Nint(tN);
  aRoundError = TMath::Abs(tN - aRoundN);
  cout << tN << " " << aRoundN << " " << aRoundError << endl;
}
