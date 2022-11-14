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

#include <sys/stat.h>
#include "Configurator.h"
#include "Model.h"
#include "Crc32.h"
#include "THGlobal.h"

using namespace std;
using namespace TMath;

extern Configurator *sMainConfig;
extern TString	sEventDIR;
extern int	sModel;
extern int	sIntegrateSample;
extern int	sRandomize;

Model::Model()
    : Xt(0.0), Xx(0.0), Xy(0.0), Xz(0.0),
    Pe(0.0), Px(0.0), Py(0.0), Pz(0.0),
    mHyperCube(0.0),
    mRandom(0),
    mSpectralFncs(0), mSpectralFncIntegrals(0),
    mDB(0)
{
  mName="";
  mHash="";
  mDescription="";
}

Model::Model(TRandom2* aRandom)
: Xt(0.0), Xx(0.0), Xy(0.0), Xz(0.0),
Pe(0.0), Px(0.0), Py(0.0), Pz(0.0),
mHyperCube(0.0),
mRandom(aRandom)
{
  mName="";
  mHash="";
  mDescription="";
  mSpectralFncs = new std::map<int, TF1*>;
  mSpectralFncIntegrals = new std::map<int, double>;

}
Model::~Model()
{
  delete mSpectralFncs;
  delete mSpectralFncIntegrals;
}

void Model::AddParameterBranch(TTree* aTree)
{
  Model_t tPar;
  
  tPar.dummy = -1.0;
  aTree->Branch(_MODEL_T_BRANCH_, &tPar, _MODEL_T_FORMAT_);
}

void Model::SetParticlePX(Particle* aParticle)
{
  aParticle->SetParticlePX(Pe,Px,Py,Pz,Xt,Xx,Xy,Xz);
}  

double Model::GetHyperCubeVolume()
{
  return mHyperCube;
}

const char* Model::GetHash()
{
  return mHash.Data();
}

const char* Model::GetName()
{
  return mName.Data();
}

const char* Model::GetDescription()
{
  return mDescription.Data();
}

void Model::CreateEventSubDir()
{
  struct stat tStatus;
  TString tEventDir = sMainConfig->GetParameter("EventDir"); tEventDir.Prepend("./");
  TString tSubDirs  = sEventDIR;
  int     tPos      = 0;

  tSubDirs.ReplaceAll(tEventDir,"");
  while(tPos <= tSubDirs.Length()) {   
    if(tEventDir.EndsWith("/")) {
      if (stat(tEventDir.Data(), &tStatus) == -1) {
	PRINT_DEBUG_3("<Model::CreateEventSubDir>\tDirectory " <<tEventDir<<" does not exist. Tying to create.");
        if (mkdir(tEventDir.Data(), S_IRWXU | S_IXUSR | S_IRWXG | S_IXGRP | S_IROTH | S_IXOTH) == -1) {
          PRINT_MESSAGE("<Model::CreateEventSubDir>\tFailed to create directory " << tEventDir);
          exit(_ERROR_GENERAL_FILE_NOT_FOUND_);
        } else {
          PRINT_DEBUG_2("<Model::CreateEventSubDir>\tDirectory "<<tEventDir<<" created.");
        }
      } else {
	PRINT_DEBUG_3("<Model::CreateEventSubDir>\tDirectory " <<tEventDir<<" does exist.");
      }
    }
    tEventDir += tSubDirs[tPos];
    tPos++;
  }  
}

void Model::CalculateHash(TString aPreHash) {
  Crc32 tHash;
  tHash.Update(aPreHash.Data(), aPreHash.Length());
  tHash.Finish();
  mHash = tHash.GetValueHex();
  PRINT_DEBUG_1("<Model::Hash>\t"<<aPreHash.Data()<<" -> 0x"<<mHash);
}


TH1F * p33() {

    static bool create = true;

    // this is from PML eBW (from his data33 column marked as A)
    static Double_t y4[] = {

	0.000335537, 0.00200371, 0.00464069, 0.00821601,
	0.0127918, 0.0184775, 0.0254215, 0.0338109, 0.0438759,
	0.0558947, 0.0702011, 0.0871924, 0.107337, 0.13118,
	0.159347, 0.192535, 0.231494, 0.276977, 0.329664,
	0.390015, 0.458079, 0.533225, 0.613851, 0.697139,
	0.778991, 0.854333, 0.917805, 0.964779, 0.992344,
	0.999881, 0.988991, 0.962874, 0.925482, 0.880741,
	0.832059, 0.782097, 0.732765, 0.685319, 0.640506,
	0.598705, 0.56004, 0.524474, 0.49187, 0.462039,
	0.434765, 0.40983, 0.387017, 0.366124, 0.346964,
	0.329365, 0.313172, 0.298246, 0.284462, 0.271709,
	0.259888, 0.24891, 0.238697, 0.229178, 0.220292,
	0.21198, 0.204195, 0.19689, 0.190026, 0.183567,
	0.177479, 0.171733, 0.166303, 0.161166, 0.156298,
	0.151681, 0.147296, 0.143128, 0.13916, 0.13538,
	0.131776, 0.128334, 0.125046, 0.121902, 0.118892,
	0.116009, 0.113244, 0.110591, 0.108044, 0.105596,
	0.103243, 0.100977, 0.0987961, 0.0966944, 0.094668,
	0.092713, 0.0908259, 0.0890032, 0.0872417, 0.0855385,
	0.0838908, 0.0822961, 0.0807517, 0.0792555, 0.0778053,
	0.076399, 0.0750347, 0.0737106, 0.0724249, 0.0711761,
	0.0699626, 0.0687831, 0.067636, 0.0665201, 0.0654342,
	0.0643771, 0.0633477, 0.062345, 0.0613679, 0.0604156,
	0.059487, 0.0585814, 0.0576979, 0.0568357, 0.0559942,
	0.0551725, 0.05437, 0.0535861, 0.0528202, 0.0520716,
	0.0513397, 0.0506241, 0.0499242, 0.0492396, 0.0485696,
	0.047914, 0.0472722, 0.0466439, 0.0460285, 0.0454259
    } ;
    // original from PML
    static Double_t y5[] = {
        0.487868830660434  ,  0.974036309173992  ,   1.48796409314226  ,   1.88297294151952  ,   2.28387966214487  ,   2.69710433117348  ,   3.06126905224752  ,   3.48369820687455  ,   3.94681509541577  ,   4.40057706566204  ,
        4.93936605645115  ,   5.48283889953534  ,   6.11314163074160  ,   6.81950751171442  ,   7.55858848427632  ,   8.39291146773919  ,   9.34237968621590  ,   10.3691236214243  ,   11.4890241422109  ,   12.6864224931698  ,
        13.9961542566581  ,   15.3159690791888  ,   16.5829494782608  ,   17.8292131736408  ,   18.9152955911726  ,   19.7159202944624  ,   20.2343075893699  ,   20.3069956416159  ,   20.0125696007672  ,   19.4358933652416  ,
        18.4928474508350  ,   17.3823315861657  ,   16.1716157278873  ,   14.9038674730610  ,   13.6447991751455  ,   12.5056812406078  ,   11.4205163366001  ,   10.4059676207370  ,   9.56824474209078  ,   8.71912443140713  ,
        8.01315763782008  ,   7.39793576023726  ,   6.80293967011530  ,   6.32923543881047  ,   5.89510052808095  ,   5.48821214622965  ,   5.13744465308651  ,   4.85174031833841  ,   4.57035038945562  ,   4.30429423098407  ,
        4.07695313297681  ,   3.89130969244800  ,   3.71215490911060  ,   3.56308036882298  ,   3.40112470196175  ,   3.25304206818625  ,   3.16482535139844  ,   3.04715152970495  ,   2.96881716587616  ,   2.88140869791066  ,
        2.78585022621101  ,   2.72700023847186  ,   2.64608188068563  ,   2.60279320128658  ,   2.54784690812262  ,   2.49439864564337  ,   2.44898753815908  ,   2.37474004226641  ,   2.38889146786364  ,   2.32064636249309  ,
        2.25598151655673  ,   2.26163499621386  ,   2.21317590725296  ,   2.18258711838417  ,   2.14657132226298  ,   2.11184984790125  ,   2.07710318841325  ,   2.04127250297028  ,   2.00975961369809  ,   1.96467048836692  ,
        1.96033089737522  ,   1.92523920165693  ,   1.86071345511037  ,   1.83261799374429  ,   1.80431223657368  ,   1.74073722854448  ,   1.70206527339233  ,   1.70885460224525  ,   1.62966568770821  ,   1.60909660760260  ,
        1.57063349526692  ,   1.53328104700055  ,   1.50929453898890  ,   1.44028961773135  ,   1.41651195253619  ,   1.37802510806219  ,   1.34406780855180  ,   1.30897465984558  ,   1.27409180694395  ,   1.23918522190402  ,
        1.20425635771394  ,   1.16951260420208  ,   1.13384329729916  ,   1.10156913915201  ,   1.05909016156322  ,   1.04733955656242  ,   1.03466485115872  ,  0.994816847474634  ,  0.955153954468695  ,  0.941716527201837  ,
        0.932576069943827  ,  0.882687531334134  ,  0.872806631363069  ,  0.862001630988857  ,  0.814744066283836  ,  0.798214874167238  ,  0.807009326536368  ,  0.792949729784183  ,  0.738143504668653  ,  0.750470752195434  ,
        0.762983110400590  ,  0.707412710433718  ,  0.695986993562356  ,  0.697185321244401  ,  0.703998672833197  ,  0.678025155222950  ,  0.665835263500309  ,  0.640255699384825  ,  0.644986822224574  ,  0.654246037039566  ,
        0.617376174539702  ,  0.639099508356292  ,  0.603315124799318  ,  0.608256543443586  , 
    };

    static TH1F *h5 = new TH1F("hBform","B(m)",134,1.078,1.743); // PML paper
    static TH1F *h4 = new TH1F("hV","V(m)",134,1.073,1.748); // vova

    if (create) {
        for (int b = 1; b <= h5->GetNbinsX(); ++b) {
            h5->SetBinContent(b,y5[b-1]);
            h4->SetBinContent(b,y4[b-1]);
        }
        create = false;
    }

//    return h4;  // eBW
    return h5;   // PML

}


double Model::CalcMass(ParticleType *aPartType, double &statWeight) {
    double m0 = aPartType->GetMass();
    double Gamma = aPartType->GetGamma();

    if (Gamma < 1e-6) {
        statWeight = 1.0;
        return m0;
    }

    double m_min = m0 - 2*Gamma;
    double m_max = m0 + 2*Gamma;
    if (m_min < 0.) {
        m_min = 0;
        m_max = 2*m0;
    }

    double norm = 1.0;  
    TF1 *fnc;
    if (mSpectralFncs->find(aPartType->GetPDGCode()) == mSpectralFncs->end()) {
        fnc = new TF1(Form("spectralFunction_pdg%i",aPartType->GetPDGCode()),"x*[0]*[1]/((x*x-[0]*[0])*(x*x-[0]*[0])+[0]*[0]*[1]*[1])",m_min,m_max);
        fnc->SetParameters(m0,Gamma);
        norm = fnc->Integral(m_min,m_max);
        
        (*mSpectralFncs)[aPartType->GetPDGCode()] = fnc;
        (*mSpectralFncIntegrals)[aPartType->GetPDGCode()] = norm;

       // cout << "mass, Gamma, m_min, m_max " << m0 << " " << Gamma << " " << m_min << " " << m_max << endl;
       // cout << "norm " << norm << endl;
    } else {
        fnc = (*mSpectralFncs)[aPartType->GetPDGCode()];
        norm = (*mSpectralFncIntegrals)[aPartType->GetPDGCode()];
    }
    double m = fnc->GetRandom(m_min,m_max);
    statWeight = fnc->Eval(m) / fnc->Eval(m0) / (m_max-m_min);
    return m;
}


double Model::CalcMassDistr(ParticleType *aPartType, double &statWeight) {

    int pdg = aPartType->GetPDGCode();
    if (Abs(pdg) == 2224 || Abs(pdg) == 2214 || Abs(pdg) == 2114 || Abs(pdg) == 1114) {
        TH1F *distr = p33();
        double m = distr->GetRandom();
        statWeight = distr->Interpolate(m);
        return m;
    } else {
	statWeight = 1;
        return aPartType->GetMass();
       // return CalcMass(aPartType, statWeight);
    }
}

void Model::GetParticleMass(ParticleType *aPartType, bool finiteWidth, double &M, double &spectralFunctionWeight ){
    if (finiteWidth) {
	M = CalcMassDistr(aPartType, spectralFunctionWeight);
    } else {
	M = aPartType->GetMass();
	spectralFunctionWeight = 1.0;
    }
}


