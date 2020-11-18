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

#include <sstream>
#include <TMath.h>
#include "THGlobal.h"
#include "Configurator.h"
#include "Parser.h"
#include "Model_Example.h"

using namespace TMath;
using namespace std;

extern Configurator* sMainConfig;
extern TString	sModelINI;
extern TString	sEventDIR;
extern TString	sTimeStamp;
extern int	sModel;
extern int	sRandomize;
extern int	sIntegrateSample;

Model_Example::Model_Example()
: Model(),
  mThermo(0),
  mParameter1(0.0), mParameter2(0.0), mParameter3(0.0)
{
}

Model_Example::Model_Example(TRandom2* aRandom)
: Model(aRandom)
{
  mName = "Example";
  mThermo = new Thermodynamics();
  ReadParameters();
  Description();
  mHyperCube = 1.0;
}

Model_Example::~Model_Example()
{
  delete mThermo;
}

<<<<<<< HEAD
double Model_Example::GetIntegrand(ParticleType* aPartType, bool)
=======
double Model_Example::GetIntegrand(ParticleType* aPartType)
>>>>>>> 6eef3d1fc1f3ac17aedbca11e4a34ae9c4017af9
{
  double Integrand, dSigmaP, PdotU, tStatistics;
  double Spin;
  double Tau, Rho, PhiS, RapS;	// space-time coordinates 
  double Mt,  Pt,  PhiP, RapP;	// four-momentum
  double dPt;
  
  Spin	= aPartType->GetSpin();
  
// Generate spacial components, transform form [fm] to [GeV^-1]
  Tau   = 1.0;
  Rho	= 1.0        * mRandom->Rndm();
  PhiS	= 2.0 * Pi() * mRandom->Rndm();
  RapS	= 1.0        * mRandom->Rndm();

// Generate momentum components
  {
    double Zet = mRandom->Rndm();
    Pt  = Zet / (1.0 - Zet);
    dPt = 1.0 / ( (1.0 - Zet) * (1.0 - Zet) );
  }
  PhiP	= 2.0 * Pi() * mRandom->Rndm();
  RapP	= 1.0        * mRandom->Rndm();
  Mt	= Hypot(aPartType->GetMass(), Pt);

// Invariants
  PdotU	= 1.0;
  
  dSigmaP = 1.0;

// Type of statistics Bose-Einstein or Fermi-Dirac
  tStatistics = ( (Spin - static_cast<int>(Spin)) < 0.01 ? -1.0 : +1.0 );
  
// Value of the integrand
  Integrand = (2.0 * Spin + 1.0) * 1.0 / kTwoPi3 * Pt * dPt * dSigmaP * PdotU;

// Return values
  Xt = Tau * CosH(RapS);
  Xx = Rho * Cos(PhiS);
  Xy = Rho * Sin(PhiS);
  Xz = Tau * SinH(RapS);
  Pe = Mt  * CosH(RapP);
  Px = Pt  * Cos(PhiP);
  Py = Pt  * Sin(PhiP);
  Pz = Mt  * SinH(RapP);
  return Integrand;
}

void Model_Example::Description()
{
  ostringstream oss;
  oss << "##################################################"<< endl;
  oss << MODEL_NAME(mName);
  oss << "# - parameter #1           : " <<MODEL_PAR_DESC(mParameter1,"[unit]");
  oss << "# - parameter #2           : " <<MODEL_PAR_DESC(mParameter2,"[unit]");
  oss << "# - parameter #3           : " <<MODEL_PAR_DESC(mParameter3,"[unit]");
  oss << "# Parameters hash (CRC32)  : " <<MODEL_PAR_DESC(mHash,	"");
  oss << "# Integration samples      : " <<MODEL_PAR_DESC(sIntegrateSample,"");
  oss << "# Random seed              : " <<MODEL_PAR_DESC((sRandomize ? "yes" : "no"),"");
  oss << "# Generation date          : " <<sTimeStamp<<" #"<<endl;
  oss << "##################################################"<< endl;
  mDescription = oss.str();
}

void Model_Example::AddParameterBranch(TTree* aTree)
{
  Model_t tPar;
  
  tPar.dummy = 0.0;
  aTree->Branch(_MODEL_T_BRANCH_, &tPar, _MODEL_T_FORMAT_)->Fill();
}

void Model_Example::ReadParameters()
{
  Configurator* tModelParam;
  Parser*      tParser;
  
  tModelParam = new Configurator;
  tParser     = new Parser(sModelINI.Data());
  tParser->ReadINI(tModelParam);
  delete tParser;
  
  try {
    mParameter1	= (tModelParam->GetParameter("Parameter1")).Atof();
    mParameter2	= (tModelParam->GetParameter("Parameter2")).Atof();
    mParameter3	= (tModelParam->GetParameter("Parameter3")).Atof();
  } catch (TString tError) {
    PRINT_MESSAGE("<Model_Example::ReadParameters>\tCaught exception " << tError);
    PRINT_MESSAGE("\tDid not find one of the necessary model parameters.");
    exit(_ERROR_CONFIG_PARAMETER_NOT_FOUND_);
  }
  
// calculate parameter hash
  ostringstream oss;
  oss << sModel;
  oss << mParameter1 << mParameter2 << mParameter3;
  CalculateHash(TString(oss.str()));

// create event subdirectory if needed
  try {
    sEventDIR += tModelParam.GetParameter("EventSubDir");
    CreateEventSubDir();
  } catch (TString tError) {
  }
  
  delete tModelParam;
}
