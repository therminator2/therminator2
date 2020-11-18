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
#include "Model_Lhyquid2DBI.h"

using namespace TMath;
using namespace std;

extern TString	sTimeStamp;
extern int	sModel;
extern int	sRandomize;
extern int	sIntegrateSample;

Model_Lhyquid2DBI::Model_Lhyquid2DBI()
: Model(), Hypersurface_Lhyquid2D()
{
}

Model_Lhyquid2DBI::Model_Lhyquid2DBI(TRandom2* aRandom)
: Model(aRandom), Hypersurface_Lhyquid2D()
{ 
  mName = "Lhyquid 2+1 BInv";
  ReadParameters();
  Description();
  mHyperCube = mRapPRange * 2*TMath::Pi() * 1.0 * GetHyperCubeSpatialVolume();
}

Model_Lhyquid2DBI::~Model_Lhyquid2DBI()
{
}

<<<<<<< HEAD
double Model_Lhyquid2DBI::GetIntegrand(ParticleType* aPartType, bool)
=======
double Model_Lhyquid2DBI::GetIntegrand(ParticleType* aPartType)
>>>>>>> 6eef3d1fc1f3ac17aedbca11e4a34ae9c4017af9
{
  double dSigmaP;
  double Spin, Statistics;
  double Mt,  Pt,  PhiP, RapP;
  double dPt;

// Type of statistics Bose-Einstein or Fermi-Dirac
  Spin	= aPartType->GetSpin();
  Statistics = ( (Spin - static_cast<int>(Spin)) < 0.01 ? -1.0 : +1.0 );
// Generate random position on the hypersurface
  SetPositionOnHypersurface(mRandom, &Xt, &Xx, &Xy, &Xz);
// Generate momentum components
  {
    double Zet = mRandom->Rndm();
    Pt	= Zet / (1.0 - Zet);	// 0 <= pT <= Infinity
    dPt	= 1.0 / ( (1.0 - Zet) * (1.0 - Zet) );
  }
  PhiP	= 2.0 * Pi() * mRandom->Rndm();
  RapP	= mRapPRange * mRandom->Rndm() - 0.5 * mRapPRange;
  Mt	= Hypot(aPartType->GetMass(),Pt);
// d Sigma_\mu p^\mu
  dSigmaP = GetDSigmaP(Mt, Pt, PhiP, RapP);   
#ifdef _MODEL_LHYQUID_ONLY_BACK_FLOW_
// only particle emission directed back to the hydro-region
  if(dSigmaP >= 0.0) 
    dSigmaP = 0.0;
  else
    dSigmaP = -dSigmaP;
#else
// disable particle emission back to the hydro-region
  if(dSigmaP < 0.0) 
    dSigmaP = 0.0;
#endif
// particle P coordinates
  Pe = Mt  * CosH(RapP);
  Px = Pt  * Cos(PhiP);
  Py = Pt  * Sin(PhiP);
  Pz = Mt  * SinH(RapP);
// integrand  
  return (2.0 * Spin + 1.0) * 1.0 / kTwoPi3 * Pt * dPt * dSigmaP * 1.0 / (Exp( ( GetPdotU(Mt, Pt, PhiP, RapP) - mThermo->GetChemicalPotential(aPartType) ) / mThermo->GetTemperature() ) + Statistics );
}

void Model_Lhyquid2DBI::Description()
{
  ostringstream oss;
  char tCentrality[10];
  sprintf(tCentrality,"%g-%g",mCentralityMin,mCentralityMax);
  
  oss << "##################################################"<< endl;
  oss << MODEL_NAME(mName);
  oss << "# - rapidity range         : " <<MODEL_PAR_DESC(mRapPRange,		"[units]");
  oss << "# - spatial rapidity range : " <<MODEL_PAR_DESC(mRapSRange,		"[units]");
  oss << "# - initial proper time    : " <<MODEL_PAR_DESC(mTauI * kHbarC,	"[fm]");
  oss << "# - freeze-out temperature : " <<MODEL_PAR_DESC(mThermo->GetTemperature() * 1000.0,	"[MeV]");
  oss << "# - chem. potential Mu_B   : " <<MODEL_PAR_DESC(mThermo->GetMuB() * 1000.0,	"[MeV]");
  oss << "# - chem. potential Mu_I3  : " <<MODEL_PAR_DESC(mThermo->GetMuI() * 1000.0,	"[MeV]");
  oss << "# - chem. potential Mu_S   : " <<MODEL_PAR_DESC(mThermo->GetMuS() * 1000.0,	"[MeV]");
  oss << "# - chem. potential Mu_C   : " <<MODEL_PAR_DESC(mThermo->GetMuC() * 1000.0,	"[MeV]");
  oss << "# - device name            : " <<MODEL_PAR_DESC(mDeviceName,	"");
  oss << "# - colliding particles    : " <<MODEL_PAR_DESC(mCollidingSystem,	"");
  oss << "# - sqrt(s_NN)             : " <<MODEL_PAR_DESC(mCollidingEnergy,	"[GeV]");
  oss << "# - centrality             : " <<MODEL_PAR_DESC(tCentrality,	"[%]");
  oss << "# - impact parameter       : " <<MODEL_PAR_DESC(mImpactParameter,	"[fm]");
  oss << "# - init. central temp.    : " <<MODEL_PAR_DESC(mTempI,		"[MeV]");
  oss << "# Parameters hash (CRC32)  : " <<MODEL_PAR_DESC(mHash,		"");
  oss << "# Integration samples      : " <<MODEL_PAR_DESC(sIntegrateSample,	"");
  oss << "# Random seed              : " <<MODEL_PAR_DESC((sRandomize ? "yes" : "no"),"");
  oss << "# Generation date          : " <<sTimeStamp<<" #"<<endl;
  oss << "##################################################"<< endl;
  mDescription = oss.str();
}

void Model_Lhyquid2DBI::AddParameterBranch(TTree* aTree)
{
  Model_t_Lhyquid2DBI tPar;
  
  tPar.RapPRange	= mRapPRange;
  tPar.RapSRange	= mRapSRange;
  tPar.TauI		= mTauI    * kHbarC;
  tPar.TempF		= mThermo->GetTemperature() * 1000.0;
  tPar.MuB		= mThermo->GetMuB() * 1000.0;
  tPar.MuI		= mThermo->GetMuI() * 1000.0;
  tPar.MuS		= mThermo->GetMuS() * 1000.0;
  tPar.MuC		= mThermo->GetMuC() * 1000.0;
  tPar.CollidingEnergy	= mCollidingEnergy;
  tPar.CentralityMin	= mCentralityMin;
  tPar.CentralityMax	= mCentralityMax;
  tPar.ImpactParameter	= mImpactParameter;
  tPar.TempI		= mTempI;
  sprintf(tPar.DeviceName,	"%s",mDeviceName);
  sprintf(tPar.CollidingSystem,	"%s",mCollidingSystem);
  aTree->Branch(_MODEL_T_BRANCH_, &tPar, _MODEL_T_FORMAT_LHYQUID2DBI_)->Fill();
}

void Model_Lhyquid2DBI::ReadParameters()
{
// calculate parameter hash
  ostringstream oss;
  oss << sModel;
  oss << mRapPRange << mRapSRange << mTauI;
  oss << mThermo->GetTemperature() << mThermo->GetMuB() << mThermo->GetMuI() << mThermo->GetMuS() << mThermo->GetMuC();
  oss << mCollidingEnergy << mCentralityMin << mCentralityMax << mImpactParameter << mTempI;
  oss << mDeviceName << mCollidingSystem;
  CalculateHash(TString(oss.str().data()));

// create event subdirectory if needed
  CreateEventSubDir();
}
