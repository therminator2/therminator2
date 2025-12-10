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

#include <stdlib.h>
#include <sys/stat.h>
#include <math.h>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TProfile.h>
#include <TMath.h>
#include <TDatime.h>
#include "ParticleCoor.h"
#include "StructEvent.h"
#include "Configurator.h"
#include "Parser.h"
#include "THGlobal.h"
#include "Messages.h"

#define PIPID 211
#define KPID 321
#define PPID 2212
#define SPPID 3222
#define LAPID 3122
#define ABSRAP 0.7
#define ETAABS 0.35
#define PTMIN 0.12
#define PTMAX 0.9
#define TMAX 500

#define COULOMBSTEPS 170

#define KAONMASS 0.493677
#define PIONMASS 0.13957

#define FOURPI 1.25663706143591725e+01

using namespace std;

TString sMainINI;
TString sEventDir;
TString sTimeStamp;
int     sParentPID;

static ParticleCoor particle1, particle2, particle1lcms, particle2lcms, particle1prf, particle2prf, buf;
ParticleCoor **evbuf;
ParticleCoor **evbuf2;

Int_t *evtscount;
Int_t *evtscount2;

double mKStarLong, mKStarOut, mKStarSide, mKStarSigned, mRLong, mDTime, mRTrans, mROut, mRSide, mRSidePairCMS, mRLongPairCMS, mDTimePairLCMS, mROutPairCMS, mDTimePairCMS, mRStar, mBetat;
double pairphi;

double mRL, mRS, mRO, mDT;
double mKO, mKS, mKL, mDE;
double mKR, mKP, mKC;

double mROS, mRSS, mRLS, mRSt;

double mKT;
double ktmin, ktmax;
double ptmin1, ptmax1, ptmin2, ptmax2;

double btmin, btmax;
int    uselcms;
int    pairtype;
int    partpid;
int    partpid2;
int    scount;

int    writesourcehists;

double gamov[2000];
long double pionac;
long double oneoveracsq;
long double twopioverac;
long double coulqscpart;     
int  twospin;
long double euler;
long double f0;
long double d0;

TH1D *tROutLCMS;
TH1D *tRSideLCMS;
TH1D *tRLongLCMS;
TH1D *tRTimeLCMS;

TH1D *tROutPRF;
TH1D *tRSidePRF;
TH1D *tRLongPRF;
TH1D *tRTimePRF;

TH3D *tROutSideLongPRF;
TH3D *tROutSideLongLCMS;

TH2D *tROutLCMSQinv;
TH2D *tRSideLCMSQinv;
TH2D *tRLongLCMSQinv;
TH2D *tRTimeLCMSQinv;

TH2D *tROutPRFQinv;
TH2D *tRSidePRFQinv;
TH2D *tRLongPRFQinv;
TH2D *tRTimePRFQinv;

TH1D *tRInvPRF;

TH2D *simage1;
TH2D *simage2;

TH1D *time1;
TH1D *time2;

struct _dcomplex 
{
  long double re;
  long double im;
};

typedef struct _dcomplex dcomplex;

dcomplex d0s;
dcomplex f0s;
dcomplex d0t;
dcomplex f0t;

void		InitializeGamow();
double		Gamow(double arg);
dcomplex	mult(dcomplex arga, dcomplex argb);
dcomplex	conj(dcomplex arg);
dcomplex	mult(dcomplex arga, long double argb);
long double	modl2(dcomplex arg);
long double	modl(dcomplex arg);
dcomplex	invr(dcomplex arg);
void		GetFFsingle(dcomplex *ffp, int sign = 1);
void		GetFFdouble(dcomplex *ffp, dcomplex *ffm);
double		GetQuantumCoulomb();
double		GetCoulomb();
long double	geth(long double eta);
long double	chiim(long double eta);
void		bfunpfun(long double eta, long double rho, long double &bret, long double &pret);
dcomplex	GetG(long double eta, long double rho, long double hfun);
double		GetQuantumCoulombStrong();
void		PairKinematics(ParticleCoor part1, ParticleCoor part2);
void		MessageIntro();
void		MessageHelp();
void		MessageVersion();
void		CopyINIFile();

int main(int argc, char **argv)
{
  Int_t evtomix = 20;
  Double_t tcut = TMAX;
  euler = 0.577215665;
  f0 = 7.77 / 0.197327;
  d0 = 2.77 / 0.197327;
  Int_t weakpid[16] = {3334, -3334, 3312, -3312, 3322, -3322,   3112,  -3112,   3212,  -3212,   3222,  -3222,  3122, -3122,  311, -311};
  Int_t isweak;

  partpid2 = 0;

  int nbin;
  int onlyprim;// = 0;
  int docoulomb = 0;

  TDatime tDate;
  tDate.Set();
  sTimeStamp = tDate.AsSQLString();
// ##############################################################
// # Command line analysis					#
// ##############################################################
  int tEventFiles;
  sParentPID = 0;

    if (argc > 3) 
    {
        nbin	= atoi(argv[1]);
        sEventDir	= argv[2];
        tEventFiles	= atoi(argv[3]);
        if (argc > 4) 
        {
            sMainINI	= argv[4];
            if(sMainINI.IsDigit()) 
            {
                sParentPID = sMainINI.Atoi();
                sMainINI   = "./femto.ini";
            }
        } 
        else
            sMainINI	= "./femto.ini";
        if (argc > 5)
            sParentPID = atoi(argv[5]);
    } 
    else 
    {
        if (argc > 1) 
        {
            sMainINI = argv[1];
            if(sMainINI.Contains("-h") || sMainINI.Contains("--help")) 
            {
                Messages::HelpFmt();
                return 0;
            }
            if(sMainINI.Contains("-v") || sMainINI.Contains("--version")) 
            {
                Messages::Version();
                return 0;
            }    
        }
        Messages::HelpFmt();
        return _ERROR_GENERAL_FILE_NOT_FOUND_;
    }
  
  Messages::Intro();
// ##############################################################
// # Read configuration file (femto.ini)			#
// ##############################################################
  Configurator* tMainConfig	= new Configurator;
  Parser* 	tParser		= new Parser(sMainINI);
  tParser->ReadINI(tMainConfig);
  delete tParser;
  
// ##############################################################
// # ReadParameters						#
// ##############################################################  
  try {
    TString tPairType;
    tPairType = tMainConfig->GetParameter("PairType");
    if (tPairType == "pion-pion")	   pairtype = 0;
    else if (tPairType == "kaon-kaon")	   pairtype = 1;
    else {
      PRINT_MESSAGE("therm2_femto Unknown pair type: " << tPairType);
      PRINT_MESSAGE("Please provide the proper pair name in the main INI file.");
      exit(_ERROR_FEMTO_UNKNOWN_PAIRTYPE_);
    }
    tcut	= tMainConfig->GetParameter("TimeCut").Atof();
    evtomix	= tMainConfig->GetParameter("EventsToMix").Atoi();
    if(tMainConfig->GetParameter("EnableOnlyPrimordial") == "yes")
      onlyprim	= 1;
    else
      onlyprim	= 0;
//    if(tMainConfig->GetParameter("EnableCoulombCorrection") == "yes")
//      docoulomb	= 1;
//    else
//      docoulomb	= 0;
    if(tMainConfig->GetParameter("EnableSourceHistograms") == "yes")
      writesourcehists = 1;
    else
      writesourcehists = 0;
  }
  catch (TString &tError) {
    PRINT_DEBUG_1("therm2_femto - Caught exception " << tError);
    PRINT_MESSAGE("Did not find one of the necessary parameters in the parameters file.");
    exit(_ERROR_CONFIG_PARAMETER_NOT_FOUND_);
  }
  
  CopyINIFile();
  
// ##############################################################
// # Pair and particle cuts					#
// ##############################################################
  btmin = 0.0;
  btmax = 1.0;

  switch(nbin) {
  case 0:
    ktmin = 0.15;
    ktmax = 0.25;
    break;
  case 1:
    ktmin = 0.25;
    ktmax = 0.35;
    break;
  case 2:
    ktmin = 0.35;
    ktmax = 0.45;
    break;
  case 3:
    ktmin = 0.45;
    ktmax = 0.6;
    break;
  case 10:
    ktmin = 0.2;
    ktmax = 0.3;
    break;
  case 11:
    ktmin = 0.3;
    ktmax = 0.36;
    break;
  case 12:
    ktmin = 0.36;
    ktmax = 0.42;
    break;
  case 13:
    ktmin = 0.42;
    ktmax = 0.48;
    break;
  case 14:
    ktmin = 0.48;
    ktmax = 0.54;
    break;
  case 15:
    ktmin = 0.54;
    ktmax = 0.60;
    break;
  case 16:
    ktmin = 0.60;
    ktmax = 0.75;
    break;
  case 17:
    ktmin = 0.75;
    ktmax = 1.0;
    break;
  case 18:
    ktmin = 1.0;
    ktmax = 1.2;
    break;
  case 19:
    ktmin = 0.2;
    ktmax = 2.0;
    break;
  case 20:
    ktmin = 0.15;
    ktmax = 0.35;
    break;
  case 21:
    ktmin = 0.35;
    ktmax = 0.55;
    break;
  case 22:
    ktmin = 0.55;
    ktmax = 0.80;
    break;
  case 30:
    ktmin = 0.0;
    ktmax = 2.0;
    btmin = 0.3;
    btmax = 0.95;
    break;
  case 31:
    ktmin = 0.0;
    ktmax = 2.0;
    btmin = 0.35;
    btmax = 0.5;
    break;
  case 32:
    ktmin = 0.0;
    ktmax = 2.0;
    btmin = 0.5;
    btmax = 0.65;
    break;
  case 33:
    ktmin = 0.0;
    ktmax = 2.0;
    btmin = 0.65;
    btmax = 0.8;
    break;
  case 34:
    ktmin = 0.0;
    ktmax = 2.0;
    btmin = 0.8;
    btmax = 0.95;
    break;
  case 61:
    ktmin = 0.0;
    ktmax = 2.0;
    btmin = 0.5;
    btmax = 0.9;
    break;
  case 71:
    ktmin = 0.5;
    ktmax = 2.5;
    btmin = 0.0;
    btmax = 1.0;
    break;
  case 81:
    ktmin = 0.1;
    ktmax = 2.5;
    btmin = 0.0;
    btmax = 1.0;
    break;
  case 91:
    ktmin = 0.0;
    ktmax = 2.0;
    btmin = 0.65;
    btmax = 0.9;
    break;
  case 101:
    ktmin = 0.0;
    ktmax = 2.0;
    btmin = 0.65;
    btmax = 0.9;
    break;
  }

  switch (pairtype) {
  case 0:
    ptmin1 = PTMIN;
    ptmin2 = PTMIN;
    ptmax1 = PTMAX;
    ptmax2 = PTMAX;
    break;
  case 1:
    ptmin1 = 0.3;
    ptmin2 = 0.3;
    ptmax1 = 0.9;
    ptmax2 = 0.9;
    break;
  case 2:
    ptmin1 = 0.4;
    ptmin2 = 0.4;
    ptmax1 = 1.8;
    ptmax2 = 1.8;
    break;
  case 3:
  case 4:
    ptmin1 = 0.05;
    ptmin2 = 0.05;
    ptmax1 = PTMAX;
    ptmax2 = 1.5;
    break;
  case 5:
  case 6:
    ptmin1 = 0.05;
    ptmin2 = 0.3;
    ptmax1 = PTMAX;
    ptmax2 = 2.0;
    break;
  case 7:
    ptmin1 = 0.05;
    ptmin2 = 0.05;
    ptmax1 = 1.2;
    ptmax2 = 1.5;
    break;
  }

  if (docoulomb) {
    InitializeGamow();
  }
  else {
    twospin = 0;
    if (pairtype == 0) 
      partpid = PIPID;
    else if (pairtype == 1) 
      partpid = KPID;
    else if (pairtype == 2) 
      partpid = PPID;
  }

// ##############################################################
// # Allocate arrays in memory					#
// ##############################################################
  
  evbuf = (ParticleCoor **) malloc(sizeof(ParticleCoor *) * evtomix);
  evtscount = (Int_t *) malloc(sizeof(Int_t) * evtomix);
  if (partpid2 != 0) {
    evbuf2 = (ParticleCoor **) malloc(sizeof(ParticleCoor *) * evtomix);
    evtscount2 = (Int_t *) malloc(sizeof(Int_t) * evtomix);
  }

  for (int imix=0; imix<evtomix; imix++) {
    evbuf[imix] = (ParticleCoor *) malloc(sizeof(ParticleCoor) * 3000);
    evtscount[imix] = 0;
    if (partpid2 != 0) {
      evbuf2[imix] = (ParticleCoor *) malloc(sizeof(ParticleCoor) * 3000);
      evtscount2[imix] = 0;
    }
  }
// ##############################################################
// # Histograms							#
// ##############################################################

  TH3D *cnumas;
  TH3D *cnuma;
  TH3D *cdena;

  TH3D *cnumasph;
  TH3D *cdenasph;

  // 3d correlation function for non-identical particles
  TH3D *cnumanonid;
  TH3D *cdenanonid;
  TH3D *cnumanonidtrue;

  TH1D *num1d;
  TH1D *num1dqsc;
  TH1D *num1dc;
  TH1D *den1d;

  TH1D *num1dp;
  TH1D *den1dp;

  TH1D *num1dn;
  TH1D *den1dn;

  TH1D *num1dptrue;
  TH1D *den1dptrue;

  TH1D *num1dntrue;
  TH1D *den1dntrue;

  TH1D *hbetat;
  TH1D *hkt;

  cnuma = new TH3D("cnuma","cnuma", 121, -0.201666, 0.201666, 121, -0.201666, 0.201666, 121, -0.201666, 0.201666);
  cdena = new TH3D("cdena","cdena", 121, -0.201666, 0.201666, 121, -0.201666, 0.201666, 121, -0.201666, 0.201666);

  cnumasph = new TH3D("cnumasph","cnumasph", 60, 0.0, 0.2, 12, -TMath::Pi(), TMath::Pi(), 12, -1.0, 1.0);
  cdenasph = new TH3D("cdenasph","cdenasph", 60, 0.0, 0.2, 12, -TMath::Pi(), TMath::Pi(), 12, -1.0, 1.0);

  cnumas = new TH3D("cnumas","cnumas", 121, -0.201666, 0.201666, 121, -0.201666, 0.201666, 121, -0.201666, 0.201666);

  cnumanonid = new TH3D("cnumanonid","cnumanonid", 60, -0.2, 0.2, 60, -0.2, 0.2, 60, -0.2, 0.2);
  cdenanonid = new TH3D("cdenanonid","cdenanonid", 60, -0.2, 0.2, 60, -0.2, 0.2, 60, -0.2, 0.2);
  cnumanonidtrue = new TH3D("cnumanonidtrue","cnumanonidtrue", 60, -0.2, 0.2, 60, -0.2, 0.2, 60, -0.2, 0.2);

  num1d = new TH1D("num1d", "num1d;2k* [GeV/c];C(k*)", 200, 0.0, 0.4);
  num1dqsc = new TH1D("num1dqsc", "num1dqsc;2k* [GeV/c];C(k*)", 200, 0.0, 0.4);
  num1dc = new TH1D("num1dc", "num1dc;2k* [GeV/c];C(k*)", 200, 0.0, 0.4);
  den1d = new TH1D("den1d", "den1d;2k* [GeV/c];C(k*)", 200, 0.0, 0.4);

  num1dp = new TH1D("num1dp", "num1dp", 100, 0.0, 0.2);
  den1dp = new TH1D("den1dp", "den1dp", 100, 0.0, 0.2);

  num1dn = new TH1D("num1dn", "num1dn", 100, 0.0, 0.2);
  den1dn = new TH1D("den1dn", "den1dn", 100, 0.0, 0.2);

  num1dptrue = new TH1D("num1dptrue", "num1dptrue", 100, 0.0, 0.2);
  den1dptrue = new TH1D("den1dptrue", "den1dptrue", 100, 0.0, 0.2);

  num1dntrue = new TH1D("num1dntrue", "num1dntrue", 100, 0.0, 0.2);
  den1dntrue = new TH1D("den1dntrue", "den1dntrue", 100, 0.0, 0.2);

  hbetat = new TH1D("hbetat", "hbetat", 100, 0.0, 1.0);
  hkt    = new TH1D("hkt",    "hkt"   , 100, 0.0, 1.2);

  if (writesourcehists) {
    tROutLCMS = new TH1D("ROutLCMS", "ROut in LCMS", 400, -200.0, 200.0);
    tRSideLCMS = new TH1D("RSideLCMS", "RSide in LCMS", 400, -200.0, 200.0);
    tRLongLCMS = new TH1D("RLongLCMS", "RLong in LCMS", 400, -200.0, 200.0);
    tRTimeLCMS = new TH1D("RTimeLCMS", "DTime in LCMS", 400, -200.0, 200.0);
    
    tROutSideLongLCMS = new TH3D("ROutSideLongLCMS", "Space distribution in LCMS", 200, -100.0, 100.0, 200, -100.0, 100.0, 200, -100.0, 100.0);
    
    tROutPRF = new TH1D("ROutPRF", "ROut in PRF", 400, -200.0, 200.0);
    tRSidePRF = new TH1D("RSidePRF", "RSide in PRF", 400, -200.0, 200.0);
    tRLongPRF = new TH1D("RLongPRF", "RLong in PRF", 400, -200.0, 200.0);
    tRTimePRF = new TH1D("RTimePRF", "DTime in PRF", 400, -200.0, 200.0);
    
    tROutSideLongPRF = new TH3D("ROutSideLongPRF", "Space distribution in PRF", 200, -100.0, 100.0, 200, -100.0, 100.0, 200, -100.0, 100.0);
    
    tROutLCMSQinv  = new TH2D("ROutLCMSKStar",  "ROut in LCMS vs qinv", 400, -200.0, 200.0, 60, 0.0, 0.3);
    tRSideLCMSQinv = new TH2D("RSideLCMSKStar", "RSide in LCMS vs qinv", 400, -200.0, 200.0, 60, 0.0, 0.3);
    tRLongLCMSQinv = new TH2D("RLongLCMSKStar", "RLong in LCMS vs qinv", 400, -200.0, 200.0, 60, 0.0, 0.3);
    tRTimeLCMSQinv = new TH2D("RTimeLCMSKStar", "DTime in LCMS vs qinv", 400, -200.0, 200.0, 60, 0.0, 0.3);
    
    tROutPRFQinv  = new TH2D("ROutPRFKStar",  "ROut in PRF vs qinv", 400, -200.0, 200.0, 60, 0.0, 0.3);
    tRSidePRFQinv = new TH2D("RSidePRFKStar", "RSide in PRF vs qinv", 400, -200.0, 200.0, 60, 0.0, 0.3);
    tRLongPRFQinv = new TH2D("RLongPRFKStar", "RLong in PRF vs qinv", 400, -200.0, 200.0, 60, 0.0, 0.3);
    tRTimePRFQinv = new TH2D("RTimePRFKStar", "DTime in PRF vs qinv", 400, -200.0, 200.0, 60, 0.0, 0.3);
    
    tRInvPRF = new TH1D("RInvPRF", "RInvariant", 400, -200.0, 200.0);

    simage1 = new TH2D("simage1", "simage1", 100, -15.0, 15.0, 100, -15.0, 15.0);
    simage2 = new TH2D("simage2", "simage2", 100, -15.0, 15.0, 100, -15.0, 15.0);
    
    time1 = new TH1D("time1", "Emission time particle1;t [fm/c]; dN/dt", 100, 0.0, 50.0);
    time2 = new TH1D("time2", "Emission time particle2;t [fm/c]; dN/dt", 100, 0.0, 50.0);
  }

  num1d->Sumw2();
  num1dqsc->Sumw2();
  num1dc->Sumw2();
  den1d->Sumw2();

  num1dp->Sumw2();
  den1dp->Sumw2();

  num1dn->Sumw2();
  den1dn->Sumw2();

  num1dptrue->Sumw2();
  den1dptrue->Sumw2();

  num1dntrue->Sumw2();
  den1dntrue->Sumw2();

  cnuma->Sumw2();
  cdena->Sumw2();

  cnumasph->Sumw2();
  cdenasph->Sumw2();

  cnumas->Sumw2();

  cnumanonid->Sumw2();
  cdenanonid->Sumw2();
  cnumanonidtrue->Sumw2();

  if (writesourcehists) {
    tROutLCMS->Sumw2();
    tRSideLCMS->Sumw2();
    tRLongLCMS->Sumw2();
    tRTimeLCMS->Sumw2();
    
    tROutSideLongLCMS->Sumw2();
    
    tROutPRF->Sumw2();
    tRSidePRF->Sumw2();
    tRLongPRF->Sumw2();
    tRTimePRF->Sumw2();
    
    tROutSideLongPRF->Sumw2();
    
    tROutLCMSQinv->Sumw2();
    tRSideLCMSQinv->Sumw2();
    tRLongLCMSQinv->Sumw2();
    tRTimeLCMSQinv->Sumw2();
    
    tROutPRFQinv->Sumw2();
    tRSidePRFQinv->Sumw2();
    tRLongPRFQinv->Sumw2();
    tRTimePRFQinv->Sumw2();
    
    tRInvPRF->Sumw2();
  }

// ##############################################################
// # Read event*.root files					#
// ##############################################################

  TChain *chn   = new TChain(_PARTICLES_TREE_);
  TChain *chnEv = new TChain(_EVENTS_TREE_);
  int pcount = 0;
  
  for(int i=0; i<tEventFiles; i++) {
    char Buff[kFileNameMaxChar];
    sprintf(Buff,"%sevent%03i.root",sEventDir.Data(),i);
    PRINT_DEBUG_1("Adding file: " << Buff);
    chn->Add(Buff);
    chnEv->Add(Buff);
  }
  
  Int_t npart		= (int) chn->GetEntries();
  int   tEventIter	= 1;
  int   tNumberOfEvents	= chnEv->GetEntries(); 
  PRINT_DEBUG_1("Particle entries: "<<npart);
  //it must be initialized
  chn->SetBranchAddress(_PARTICLE_BRANCH_,&buf);

  if(tNumberOfEvents<20000) {
    PRINT_MESSAGE("WARNING The number of events read ("<<tNumberOfEvents<<") may be not sufficient, statistics is too low.");
    PRINT_MESSAGE("In order to get more stable results please use at least 20000 - 50000 events.");
  } 
  
// ##############################################################
// # Log file entry 						#
// ##############################################################
  {
    TString tLogName;
    ofstream tLogFile;
    
    try {
      tLogName = tMainConfig->GetParameter("LogFile");
    }
    catch (TString &tError) {
    }   
    if (!tLogName.IsNull()) {
      tDate.Set();
      tLogFile.open(tLogName, ios_base::app);
      if (static_cast<long>(tLogFile.tellp()) == 0) {
        tLogFile << "# THERMINATOR 2 Log File"<<endl;
      }
      tLogFile << '['<<tDate.AsSQLString()<<"]\ttherm2_femto\t"<<sParentPID<<'\t';
      tLogFile << "[input]\t"<<nbin<<'\t'<<sEventDir<<'\t'<<tEventFiles<<'\t'<<sMainINI<<'\t'<<sParentPID<<'\t';
      tLogFile << endl;
      tLogFile << '['<<tDate.AsSQLString()<<"]\ttherm2_femto\t"<<sParentPID<<'\t';
      tLogFile << "[events]\t"<<tEventFiles<<'\t'<<tNumberOfEvents<<'\t'<<evtomix<<'\t';
      tLogFile << "[pair]\t"<<pairtype<<'\t'<<onlyprim<<'\t'<<docoulomb<<'\t';
      tLogFile << endl;
      tLogFile.close();
    }
  }
  
// ##############################################################
// # Particle mixing						#
// ##############################################################
  
  tDate.Set();
  PRINT_MESSAGE("["<<tDate.AsSQLString()<<"]\tMixing particles");   
  
  UInt_t curev = -1;
  Int_t  eviter = -1;
  Double_t pt, rap, peta;
  double coulombweight;
  double quantumweight;
   
  for (int iter=0; iter<npart; iter++) {
    chn->GetEntry(iter);
    if (buf.eventid != curev) {
      if (eviter>-1) {
	if (partpid2 == 0) {
	  // Mix different-event particles;
	  // for identical particle correlations
	  for (int fiter=0; fiter<evtscount[eviter]; fiter++) {
	    for (int mixiter=0; mixiter<evtomix; mixiter++) {
	      if (mixiter != eviter) {
		for (int siter=0; siter<evtscount[mixiter]; siter++) {
		  PairKinematics(evbuf[eviter][fiter], 
				 evbuf[mixiter][siter]);
 		  if ((mKT>ktmin) && (mKT<ktmax)) {
		    hbetat->Fill(mBetat);
		    hkt->Fill(mKT);
		    cdena->Fill(mKO, mKS, mKL, 1.0);
		    den1d->Fill(fabs(mKStarSigned)*2.0, 1.0);

		    cdenasph->Fill(mKR, mKP, mKC);

		    if (twospin == 0) {
		      quantumweight = 1.0+TMath::Cos(-mKO*mRO - mKS*mRS - mKL*mRL + mDE*mDT);
		    }
		    else if (twospin == 1) {
		      if (pcount ==3) {
			quantumweight = 1.0+TMath::Cos(-mKO*mRO - mKS*mRS - mKL*mRL + mDE*mDT);
		      }
		      else {
			quantumweight = 1.0-TMath::Cos(-mKO*mRO - mKS*mRS - mKL*mRL + mDE*mDT);
		      }
		      pcount++;
		      if (pcount == 4) pcount=0;
		    }
		    if (fabs(evbuf[eviter][fiter].t - evbuf[mixiter][siter].t) > 500.0)
		      quantumweight = 1.0;
		    cnuma->Fill(mKO, mKS, mKL, quantumweight);
		    num1d->Fill(fabs(mKStarSigned)*2.0, quantumweight);
		    
		    if (docoulomb) {
		      if (fabs(evbuf[eviter][fiter].t - evbuf[mixiter][siter].t) > 500.0) {
			coulombweight = 1.0;
			coulqscpart = 1.0;
		      }
		      else if (pairtype != 2)
			coulombweight = GetQuantumCoulomb();
		      else 
			coulombweight = GetQuantumCoulombStrong();
		      
		      if (coulombweight > 10.0) {
			cout << "Weigth for " << fabs(mKStarSigned) << " " << fiter << " " << siter << " " << " is " << coulombweight  << endl;
		      }
		      
		      cnumas->Fill(mKO, mKS, mKL, coulombweight);
		      num1dqsc->Fill(fabs(mKStarSigned)*2.0, coulombweight);
		      num1dc->Fill(fabs(mKStarSigned)*2.0, coulqscpart);
		      cnumasph->Fill(mKR, mKP, mKC, coulombweight);
		    }
// 		    if (fabs(mKStarSigned) < 0.1) {
// 		      double rad = hypot(evbuf[eviter][fiter].x, evbuf[eviter][fiter].y);
// 		      double phi = TMath::ATan2(evbuf[eviter][fiter].y, evbuf[eviter][fiter].x);
// 		      simage1->Fill(rad * cos(phi - pairphi), rad * sin(phi - pairphi), 1.0);
// 		      time1->Fill(evbuf[eviter][fiter].t);
// 		    }
 		  }
		}
	      }
	    }
	  }
	}
	//Mix different-event particles
	//for non-identical particle correlations
	else {
	  for (int fiter=0; fiter<evtscount[eviter]; fiter++) {
	    for (int mixiter=0; mixiter<evtomix; mixiter++) {
	      for (int siter=0; siter<evtscount2[mixiter]; siter++) {
		PairKinematics(evbuf[eviter][fiter], 
			       evbuf2[mixiter][siter]);
		if ((mKT>ktmin) && (mKT<ktmax) && (mBetat>btmin) && (mBetat<btmax)) {
		  hbetat->Fill(mBetat);
		  hkt->Fill(mKT);
		  
		  if (mKStarOut > 0.0) 
		    den1dp->Fill(fabs(mKStarSigned), 1.0);
		  else
		    den1dn->Fill(fabs(mKStarSigned), 1.0);
		  
		  cdenanonid->Fill(mKStarOut, mKStarSide, mKStarLong, 1);
		  
		  cdenasph->Fill(mKR, mKP, mKC);

		  if (mKStarOut > 0.0) {
		    if (mixiter != eviter) 
		      den1dptrue->Fill(fabs(mKStarSigned), 1.0);
		  }
		  else {
		    if (mixiter != eviter) 
		      den1dntrue->Fill(fabs(mKStarSigned), 1.0);
		  }
		  if (docoulomb) {
		    if (fabs(evbuf[eviter][fiter].t - evbuf[mixiter][siter].t) > 500.0)
		      coulombweight = 1.0;
		    else
		      coulombweight = GetCoulomb();
		    num1d->Fill(fabs(mKStarSigned), coulombweight); 
		    cnumanonid->Fill(mKStarOut, mKStarSide, mKStarLong, coulombweight);
		    
		    cnumasph->Fill(mKR, mKP, mKC, coulombweight);
		    
		    if (mKStarOut > 0.0)
		      num1dp->Fill(fabs(mKStarSigned), coulombweight);
		    else
		      num1dn->Fill(fabs(mKStarSigned), coulombweight);
		    if (mixiter == eviter) {
		      cnumanonidtrue->Fill(mKStarOut, mKStarSide, mKStarLong, coulombweight);
		      if (mKStarOut > 0.0)
			num1dptrue->Fill(fabs(mKStarSigned), coulombweight);
		      else
			num1dntrue->Fill(fabs(mKStarSigned), coulombweight);
		    }
		  }
		  if (fabs(mKStarSigned) < 0.1) {
		    double rad = hypot(evbuf[eviter][fiter].x, evbuf[eviter][fiter].y);
		    double phi = TMath::ATan2(evbuf[eviter][fiter].y, evbuf[eviter][fiter].x);
		    simage1->Fill(rad * cos(phi - pairphi), rad * sin(phi - pairphi), 1.0);
		    time1->Fill(evbuf[eviter][fiter].t);

		    rad = hypot(evbuf2[mixiter][siter].x, evbuf2[mixiter][siter].y);
		    phi = TMath::ATan2(evbuf2[mixiter][siter].y, evbuf2[mixiter][siter].x);
		    simage2->Fill(rad * cos(phi - pairphi), rad * sin(phi - pairphi), 1.0);
		    time2->Fill(evbuf2[mixiter][siter].t);
		    
		  }
		}
	      }
	    }
	  }
	}
      }
      
      curev = buf.eventid;
      eviter = (eviter + 1) % evtomix;
      
      cout << "\r\tevent " << tEventIter <<"/"<< tNumberOfEvents << ", particle 1 count " << evtscount[eviter];
      evtscount[eviter] = 0;
      if (partpid2 != 0) {
	cout << ", particle 2 count " << evtscount2[eviter];
	evtscount2[eviter] = 0;
      }      
      cout.flush();
      tEventIter++;
    }
    if (((buf.pid == partpid) || ((buf.pid == partpid2) && (partpid2 !=0))) && (buf.t<tcut)) {
      pt = hypot(buf.px, buf.py);
      if ((((pt>ptmin1) && (pt<ptmax1) && (buf.pid == partpid)) ||
	   ((pt>ptmin2) && (pt<ptmax2) && (buf.pid == partpid2) && (partpid2 != 0))) &&
	  ((onlyprim == 0) || ((buf.fatherpid == partpid) || ((buf.fatherpid == partpid2) && (partpid2 != 0))))) {
	rap = 0.5*log((buf.e+buf.pz)/(buf.e-buf.pz));
	peta = -TMath::Log(TMath::Tan(TMath::ATan2(pt, buf.pz)/2.0));
	
	// Check if it is coming from a weak decay
	isweak = 0;
	for (int witer=0; witer<16; witer++) 
	  if (buf.fatherpid == weakpid[witer]) isweak = 1;
	      
	if ((fabs(rap) < ABSRAP) &&
	    ((fabs(peta) < ETAABS) || (nbin < 10) || (nbin>19)) &&
	    ((!isweak) || (onlyprim == 2)))
	  {
	    if (buf.pid == partpid) {
	      evbuf[eviter][evtscount[eviter]].x = buf.x;
	      evbuf[eviter][evtscount[eviter]].y = buf.y;
	      evbuf[eviter][evtscount[eviter]].z = buf.z;
	      evbuf[eviter][evtscount[eviter]].t = buf.t;
	      evbuf[eviter][evtscount[eviter]].px = buf.px;
	      evbuf[eviter][evtscount[eviter]].py = buf.py;
	      evbuf[eviter][evtscount[eviter]].pz = buf.pz;
	      evbuf[eviter][evtscount[eviter]].e = buf.e;
	      evbuf[eviter][evtscount[eviter]].pid = buf.pid;
	      evbuf[eviter][evtscount[eviter]].fatherpid = buf.fatherpid;
	      evbuf[eviter][evtscount[eviter]].rootpid = buf.rootpid;
	      evbuf[eviter][evtscount[eviter]].mass = buf.mass;
	      evbuf[eviter][evtscount[eviter]].eventid = buf.eventid;
	      evtscount[eviter]++;
	    }
	    else if ((buf.pid == partpid2) && (partpid2 != 0)) {
	      evbuf2[eviter][evtscount2[eviter]].x = buf.x;
	      evbuf2[eviter][evtscount2[eviter]].y = buf.y;
	      evbuf2[eviter][evtscount2[eviter]].z = buf.z;
	      evbuf2[eviter][evtscount2[eviter]].t = buf.t;
	      evbuf2[eviter][evtscount2[eviter]].px = buf.px;
	      evbuf2[eviter][evtscount2[eviter]].py = buf.py;
	      evbuf2[eviter][evtscount2[eviter]].pz = buf.pz;
	      evbuf2[eviter][evtscount2[eviter]].e = buf.e;
	      evbuf2[eviter][evtscount2[eviter]].pid = buf.pid;
	      evbuf2[eviter][evtscount2[eviter]].fatherpid = buf.fatherpid;
	      evbuf2[eviter][evtscount2[eviter]].rootpid = buf.rootpid;
	      evbuf2[eviter][evtscount2[eviter]].mass = buf.mass;
	      evbuf2[eviter][evtscount2[eviter]].eventid = buf.eventid;
	      evtscount2[eviter]++;
	    }
	  }
      }
    }
  }
  
// ##############################################################
// # Save Histograms to files					#
// ##############################################################  
  PRINT_MESSAGE("");
  char bufs[kFileNameMaxChar];
  if (pairtype == 0)
    sprintf(bufs, "%sfemtopipi%i%s.root",sEventDir.Data(), nbin, onlyprim ? "p" : "a");
  else if (pairtype == 1)
    sprintf(bufs, "%sfemtokaka%i%s.root",sEventDir.Data(), nbin, onlyprim ? "p" : "a");

  TFile *ofile = new TFile(bufs, "RECREATE");
  ofile->cd();

  if (pairtype == 0) {
    cnuma->Write();
    cdena->Write();
    num1d->Write();
    den1d->Write();

    if (docoulomb) {
      cnumasph->Write();
      cdenasph->Write();    
      cnumas->Write();
      num1dqsc->Write();
      num1dc->Write();
    }
  }

  if (pairtype != 0) {
    cnumanonid->Write();
    cdenanonid->Write();
    cnumanonidtrue->Write();
    
    num1dp->Write();
    den1dp->Write();
    
    num1dn->Write();
    den1dn->Write();
    
    num1dptrue->Write();
    den1dptrue->Write();
    
    num1dntrue->Write();
    den1dntrue->Write();
  }

  if (writesourcehists) {
    tROutLCMS->Write();
    tRSideLCMS->Write();
    tRLongLCMS->Write();
    tRTimeLCMS->Write();
    
    tROutSideLongLCMS->Write();
    
    tROutPRF->Write();
    tRSidePRF->Write();
    tRLongPRF->Write();
    tRTimePRF->Write();
    
    tROutSideLongPRF->Write();
    
    tROutLCMSQinv->Write();
    tRSideLCMSQinv->Write();
    tRLongLCMSQinv->Write();
    tRTimeLCMSQinv->Write();
    
    tROutPRFQinv->Write();
    tRSidePRFQinv->Write();
    tRLongPRFQinv->Write();
    tRTimePRFQinv->Write();
    
    tRInvPRF->Write();
    
    simage1->Write();
    time1->Write();
    if (partpid2 != 0) {
      simage2->Write();
      time2->Write();
    }
  }

  hbetat->Write();
  hkt->Write();
  
  tDate.Set();
  PRINT_MESSAGE("["<<tDate.AsSQLString()<<"]\tFile "<<bufs<<" written.");
   
// ##############################################################
// # Log file entry 						#
// ##############################################################
  {
    TString tLogName;
    TString tCFName;
    ofstream tLogFile;   
    
    try {
      tLogName = tMainConfig->GetParameter("LogFile");
    }
    catch (TString &tError) {
    }   
    if (!tLogName.IsNull()) {
      tDate.Set();
      tCFName = bufs;
      tCFName.ReplaceAll(sEventDir,"");
      
      tLogFile.open(tLogName, ios_base::app);
      if (static_cast<long>(tLogFile.tellp()) == 0) {
        tLogFile << "# THERMINATOR 2 Log File"<<endl;
      }
      tLogFile << '['<<tDate.AsSQLString()<<"]\ttherm2_femto\t"<<sParentPID<<'\t';
      tLogFile << "[output]\t"<<sEventDir<<'\t'<<tCFName<<'\t';
      tLogFile << endl;
      tLogFile.close();
    }
  }

// ##############################################################
// # Temporary file						#
// ##############################################################  
  ofstream tTmpFile;
  char     tTmpFileName[kFileNameMaxChar];
   
  sprintf(tTmpFileName,"./femto_%i.tmp",sParentPID);
  tTmpFile.open(tTmpFileName, ios_base::app);
  if((tTmpFile) && tTmpFile.is_open()) {
    tTmpFile << bufs << endl;
  } else {
    PRINT_MESSAGE("Unable to create temp file "<<tTmpFileName);
    exit(_ERROR_GENERAL_FILE_NOT_FOUND_);
  }
  
  return 0; 
}

void InitializeGamow()
{
  twospin = 0;
  if (pairtype == 0) {
    pionac = 387.5 / 0.197327;
    partpid = PIPID;
    partpid2 = 0;
  }
  else if (pairtype == 1) {
    pionac = 109.55 / 0.197327;
    partpid = KPID;
    partpid2 = 0;
  }
  else if (pairtype == 2) {
    pionac = 57.63975274 / 0.197327;
    partpid = PPID;
    partpid2 = 0;
    twospin = 1;
    d0s.re = 2.77 / 0.197327;
    f0s.re = 7.77 / 0.197327;
    d0t.re = 1.7  / 0.197327;
    f0t.re = -5.4 / 0.197327;
  }
  else if (pairtype == 3) {
    pionac = 248.519 / 0.197327;
    partpid = PIPID;
    partpid2 = KPID;
  }
  else if (pairtype == 4) {
    pionac = -248.519 / 0.197327;
    partpid = -PIPID;
    partpid2 = KPID;
  }
  else if (pairtype == 5) {
    pionac = 222.564 / 0.197327;
    partpid = PIPID;
    partpid2 = PPID;
  }
  else if (pairtype == 6) {
    pionac = -222.564 / 0.197327;
    partpid = -PIPID;
    partpid2 = PPID;
  }
  else if (pairtype == 7) {
    pionac = 83.59432006 / 0.197327;
    partpid = KPID;
    partpid2 = PPID;
  }
  else {
    cout << "Unknown pair type " << pairtype << endl;
    exit(0);
  }

  oneoveracsq = 1.0/(pionac * pionac);
  twopioverac = 2.0*TMath::Pi()/pionac;
  double tpaoverk;
  
  for (int iter=0; iter<2000; iter++) {
    tpaoverk = twopioverac/(iter*0.0002 + 0.0001);
    gamov[iter] = tpaoverk * 1.0 / (exp(tpaoverk) - 1);
  }
}

double Gamow(double arg)
{
  long double eta = twopioverac / arg;
  return (eta) * 1.0/(exp(eta) - 1.0);
}

dcomplex mult(dcomplex arga, dcomplex argb)
{
  dcomplex res;
  
  res.re = arga.re * argb.re - arga.im * argb.im;
  res.im = arga.re * argb.im + argb.re * arga.im;

  return res;
}
dcomplex conj(dcomplex arg)
{
  dcomplex res;
  
  res.re = arg.re;
  res.im = -arg.im;

  return res;
}

dcomplex mult(dcomplex arga, long double argb)
{
  dcomplex res;
  
  res.re = arga.re * argb;
  res.im = arga.im * argb;

  return res;
}

long double modl2(dcomplex arg)
{
  return arg.re*arg.re + arg.im*arg.im;
}

long double modl(dcomplex arg)
{
  return hypot(arg.re, arg.im);
}

dcomplex invr(dcomplex arg)
{
  dcomplex res;
  
  res.re = arg.re/modl2(arg);
  res.im = - arg.im/modl2(arg);

  return res;
}

// Calculates the confluent hypergeometric function F
// from single orientation of cos(theta*)
// For non-symmetrized wave-function (non-identical particles)
void GetFFsingle(dcomplex *ffp, int sign)
{
  double comprep[COULOMBSTEPS];
  double compimp[COULOMBSTEPS];
  double eta, ksip;
  dcomplex alfa, zetp;
  
  int nsteps;

  double kstar = fabs(mKStarSigned);
  double kstrst = mKStarOut*mROS + mKStarSide*mRSS + mKStarLong*mRLS;
  double coskr = sign * kstrst/(fabs(mKStarSigned) * mRSt);

  if (kstar*mRSt*(1+coskr) > 15.0)
    nsteps = 170;
  else if (kstar*mRSt*(1+coskr) > 10.0)
    nsteps = 45;
  else if (kstar*mRSt*(1+coskr) > 5.0)
    nsteps = 35;
  else
    nsteps = 15;

  eta = 1.0/(kstar * pionac);
  alfa.re = 0.0;
  alfa.im = -eta;

  dcomplex fcomp, scompp;
  double tcomp;
  dcomplex sump;
  dcomplex fcmult;

  double rad = mRSt;

  ksip = kstar*rad*(1+coskr);

  zetp.re = 0.0;
  zetp.im = ksip;
      
  fcomp.re = 1.0;
  fcomp.im = 0.0;
  scompp.re = 1.0; 
  scompp.im = 0.0;
  tcomp = 1.0;
      
  for (int istep=0; istep<nsteps; istep++) {
    sump = mult(fcomp, scompp);

    sump = mult(sump, 1.0/(tcomp*tcomp));
	
    if (istep == 0) {
      comprep[istep] = sump.re;
      compimp[istep] = sump.im;
    }
    else {
      comprep[istep] = comprep[istep-1] + sump.re;
      compimp[istep] = compimp[istep-1] + sump.im;
    }
    
    fcmult.re = alfa.re + istep;
    fcmult.im = alfa.im;
	
    fcomp = mult(fcomp, fcmult);
    scompp = mult(scompp, zetp);
    tcomp *= (istep+1);

    if ((sump.re*sump.re+sump.im*sump.im) < 1.0e-14) {
      nsteps = istep;
      break;
    }
  }
  
  ffp->re = comprep[nsteps-1];
  ffp->im = compimp[nsteps-1]; 
}

// Calculates the confluent hypergeometric function
// For two orientations of cos(theta*) 
// For symmetrized wave-function (identical particles)
void GetFFdouble(dcomplex *ffp, dcomplex *ffm)
{
  long double comprep[COULOMBSTEPS];
  long double compimp[COULOMBSTEPS];
  long double comprem[COULOMBSTEPS];
  long double compimm[COULOMBSTEPS];
  long double eta, ksip, ksim;
  dcomplex alfa, zetp, zetm;
  
  int nsteps;

  long double kstar = fabs(mKStarSigned);
  long double kstrst = mKStarOut*mROS + mKStarSide*mRSS + mKStarLong*mRLS;
  long double coskr = kstrst/(fabs(mKStarSigned) * mRSt);

  if ((kstar*mRSt*(1+coskr) < 5.0) &&
      (kstar*mRSt*(1-coskr) < 5.0))
    nsteps = 25;
  else if ((kstar*mRSt*(1+coskr) < 10.0) &&
	   (kstar*mRSt*(1-coskr) < 10.0))
    nsteps = 45;
  else if ((kstar*mRSt*(1+coskr) < 15.0) &&
	   (kstar*mRSt*(1-coskr) < 15.0))
    nsteps = 110;
  else
    nsteps = 150;
  
  eta = 1.0/(kstar * pionac);
  alfa.re = 0.0;
  alfa.im = -eta;

  dcomplex fcomp, scompp, scompm;
  long double tcomp;
  dcomplex sump, summ;
  dcomplex fcmult;

  long double rad = mRSt;

  ksip = kstar*rad*(1+coskr);
  ksim = kstar*rad*(1-coskr);

  zetp.re = 0.0;
  zetp.im = ksip;
      
  zetm.re = 0.0;
  zetm.im = ksim;

  fcomp.re = 1.0;
  fcomp.im = 0.0;
  scompp.re = 1.0; 
  scompp.im = 0.0;
  scompm.re = 1.0; 
  scompm.im = 0.0;
  tcomp = 1.0;
      
  for (int istep=0; istep<nsteps; istep++) {
    sump = mult(fcomp, scompp);
    summ = mult(fcomp, scompm);

    sump = mult(sump, 1.0/(tcomp*tcomp));
    summ = mult(summ, 1.0/(tcomp*tcomp));
	
	
    if (istep == 0) {
      comprep[istep] = sump.re;
      compimp[istep] = sump.im;
      
      comprem[istep] = summ.re;
      compimm[istep] = summ.im;
    }
    else {
      comprep[istep] = comprep[istep-1] + sump.re;
      compimp[istep] = compimp[istep-1] + sump.im;
      
      comprem[istep] = comprem[istep-1] + summ.re;
      compimm[istep] = compimm[istep-1] + summ.im;
    }
    
    fcmult.re = alfa.re + istep;
    fcmult.im = alfa.im;
	
    fcomp = mult(fcomp, fcmult);
    scompp = mult(scompp, zetp);
    scompm = mult(scompm, zetm);
    tcomp *= (istep+1);
  }
  
  ffp->re = comprep[nsteps-1];
  ffp->im = compimp[nsteps-1];

  ffm->re = comprem[nsteps-1];
  ffm->im = compimm[nsteps-1];
}

// Calculate the wave-function modulus sqaured
// for identical bosons (symmetrized)
// with Coulomb interaction included
double GetQuantumCoulomb()
{
  if (mRSt < 0.0000000001)
    return 1.0;

  double kstrst = mKStarOut*mROS + mKStarSide*mRSS + mKStarLong*mRLS;
  int ccase = 0;
  static int pcount = 0;
  int wavesign = 1;

  if (twospin == 1) {
    if (pcount == 3)
      wavesign = 1;
    else
      wavesign = -1;
    pcount++;
    if (pcount == 4) pcount = 0;
  }

  // Classical limit - if distance is larger than Coulomb radius, 
  // the interaction does not matter
  if (fabs(mRSt) > fabs(pionac)) return (1.0 + wavesign*cos(2*kstrst));

  // Classical limit - in the case of large k* we go to 
  // classical coulomb interaction
  long double testp = fabs(mKStarSigned) * mRSt * (1.0 + kstrst/(mRSt*fabs(mKStarSigned)));
  long double testm = fabs(mKStarSigned) * mRSt * (1.0 - kstrst/(mRSt*fabs(mKStarSigned)));

  if ((testp> 15.0) && (testm> 15.0))
    {
      double fasymplus  = (1.0 - 1.0/((mRSt+kstrst)*pionac*mKStarSigned*mKStarSigned));
      double fasymminus = (1.0 - 1.0/((mRSt-kstrst)*pionac*mKStarSigned*mKStarSigned));
      return 0.5 * ((fasymplus + fasymminus)*cos(2*kstrst) + (2.0*sqrt(fasymplus*fasymminus)));
    }
  
  dcomplex ffplus, ffminus;
  // Check for the classical limit in both functions separately
  if (((testp< 15.0) && (testm< 15.0))) // ||
    {
      // Calculate the F function
      GetFFdouble(&ffplus, &ffminus);
      ccase = 1;
    }
  else if (testp< 15.0)
    {
      double asym;
      GetFFsingle(&ffplus);
      asym = (1.0 - 1.0/(mRSt*(1.0 - kstrst/(mRSt*fabs(mKStarSigned))*pionac*mKStarSigned*mKStarSigned)))/Gamow(fabs(mKStarSigned));
      asym = sqrt(asym);
      if (asym < 1.0) 
	ffminus.re = 1.0 + (asym -1.0) *2.0;
      else
	ffminus.re = 1.0 + (asym -1.0) /2.0;
      ffminus.im = sqrt(asym*asym - ffminus.re*ffminus.re);
      ccase = 2;
    }
  else 
    {
      double asym;
      GetFFsingle(&ffminus, -1);
      asym = (1.0 - 1.0/(mRSt*(1.0 + kstrst/(mRSt*fabs(mKStarSigned))*pionac*mKStarSigned*mKStarSigned)))/Gamow(fabs(mKStarSigned));
      asym = sqrt(asym);
      if (asym < 1.0) 
	ffplus.re = 1.0 + (asym -1.0) *2.0;
      else
	ffplus.re = 1.0 + (asym -1.0) /2.0;
      ffplus.im = sqrt(asym*asym - ffplus.re*ffplus.re);
      ccase = 3;
    }

  dcomplex expikr;
  expikr.re = cos(kstrst);
  expikr.im = sin(kstrst);
  dcomplex ffplusc = conj(ffplus);
  dcomplex ffminusc = conj(ffminus);

  dcomplex expikr2 = mult(expikr, expikr);
  dcomplex expikrc2 = conj(expikr2);
  dcomplex sterm = mult(expikr2,  mult(ffplus, ffminusc));
  dcomplex tterm = mult(expikrc2, mult(ffminus, ffplusc));


  if (!isfinite(ffminus.re))
    cout << "FFMinus Re not a number !" << " " << ccase<< endl;
  
  if (!isfinite(ffminus.im))
    cout << "FFMinus Im not a number !" << " " << ccase<< endl;

  if ((ffplus.re > 2.0) || (ffplus.re < -2.0))
    cout << "FFplus Re wild !" << ffplus.re << endl;

  if ((ffplus.im > 2.0) || (ffplus.im < -2.0))
    cout << "FFplus Im wild !" << ffplus.im << endl;

  if ((ffminus.re > 2.0) || (ffminus.re < -2.0))
    cout << "FFminus Re wild !" << ffminus.re << " " << ccase << endl;
  
  if ((ffminus.im > 2.0) || (ffminus.im < -2.0))
    cout << "FFminus Im wild !" << ffminus.im << " " << ccase <<endl;

  coulqscpart = 0.5 * Gamow(fabs(mKStarSigned)) * (modl2(ffplus) + modl2(ffminus));

  return (0.5 * Gamow(fabs(mKStarSigned)) * 
	  (modl2(ffplus) + wavesign*sterm.re + wavesign*tterm.re + modl2(ffminus)));
}

// Calculate the wave-function modulus sqaured
// for non-identical particles
// that is the Coulomb interaction
double GetCoulomb()
{
  double kstrst = mKStarOut*mROS + mKStarSide*mRSS + mKStarLong*mRLS;

  // Classical limit - if distance is larger than Coulomb radius, 
  // the interaction does not matter
  if (fabs(mRSt) > fabs(pionac)) return (1.0);
  if (fabs(mRSt) == 0.0) return (Gamow(fabs(mKStarSigned)));

  // Classical limit - in the case of large k*r* product we go to 
  // classical coulomb interaction
  if (fabs(mKStarSigned) * mRSt * (1.0 + kstrst/(mRSt*fabs(mKStarSigned)))> 15.0)
    return (1.0 - 1.0/(mRSt*pionac*mKStarSigned*mKStarSigned));
  
  // Calculate the F function
  dcomplex ffplus;
  GetFFsingle(&ffplus);

  if (!isfinite(ffplus.re)) {
    cout << "FFPlus Re not a number !" << " " << endl;
    cout << mRSt << " " << mKStarSigned << " " << mROS << " " << mRSS << " " << mRLS << endl;
  }

  if (!isfinite(ffplus.im))
    cout << "FFPlus Im not a number !" << " " << endl;


  return (Gamow(fabs(mKStarSigned)) * modl2(ffplus));
}

// Calculates the h function for the strong interaction
long double geth(long double eta)
{
  long double etasum = log(1.0/eta) - euler;
  long double series = 0.0;
  long double x2inv = (eta*eta);
  long double element;
  long double save;
  for (int iter=1; iter<1000000; iter++) {
    element = ((1.0*iter)*(1.0*iter) + x2inv)*(1.0*iter);
    element = 1.0/element;
    if (iter == 1) save = 1.0e-10 * element;

    series += element;
    if (element < save) break;
  }
  series *= x2inv;
  etasum += series;

  return etasum;
}

long double chiim(long double eta)
{
  return Gamow(1.0/(eta*pionac))/(2.0*eta);
}

void bfunpfun(long double eta, long double rho, long double &bret, long double &pret)
{
  long double b0 = 1;
  long double b1 = eta*rho;
  long double bsum = b0+b1;
  long double bnpu, bn, bnmu;
  long double p0 = 1.0;
  long double p1 = 0.0;
  long double psum = p0+p1;
  long double pnpu, pn, pnmu;

  if (rho > FOURPI) {
    bret = sin(rho)/rho;
    pret = cos(rho);
    return;
  }

  bn = b1; bnmu = b0;
  pn = p1; pnmu = p0;
  for (int iter=1; iter<100000; iter++) {
    bnpu = 2 * eta*rho *bn - rho*rho*bnmu;
    bnpu /= (1.0*iter+1.0)*(1.0*iter+2.0);
    bsum += bnpu;

    pnpu = 2 * eta*rho *pn - rho*rho*pnmu - (2.0*iter+1.0)*2.0*eta*rho*bn;
    pnpu /= (1.0*iter)*(1.0*iter+1.0);
    psum += pnpu;

    bnmu = bn;
    bn   = bnpu;

    pnmu = pn;
    pn   = pnpu;
    if ((fabs(pnmu) + fabs (bnmu) + fabs(bnpu) + fabs(pnpu)) < 1.0e-20) 
      {
	break;
      }
  }

  bret = bsum;
  pret = psum;
}

// Calculate the Gtilde function
dcomplex GetG(long double eta, long double rho, long double hfun)
{
  dcomplex gtemp;
  long double bres, pres;
  long double bmult;

  bfunpfun(eta, rho, bres, pres);
  bmult = 2.0*eta*rho*bres;

  gtemp.re = pres + bmult*(log(fabs(2.0*eta*rho))+2.0*euler-1.0+hfun);
  gtemp.im = bmult * chiim(eta);

  return gtemp;
}

void Getfc(long double kstar, long double eta, long double hfun, dcomplex &fcs, dcomplex &fct)
{
  dcomplex ci;
  dcomplex cia;

  dcomplex fis;
  dcomplex dis;
  dcomplex fcinvs;

  dcomplex fit;
  dcomplex dit;
  dcomplex fcinvt;

  ci.re = hfun;
  ci.im = chiim(eta);
  cia = mult(ci, 2.0/pionac); 

  fis = invr(f0s);
  dis = mult(d0s, 0.5 * kstar * kstar);
  fcinvs.re = fis.re + dis.re - cia.re;
  fcinvs.im = fis.im + dis.im - cia.im;
  fcs = invr(fcinvs);
  
  fit = invr(f0t);
  dit = mult(d0t, 0.5 * kstar * kstar);
  fcinvt.re = fit.re + dit.re - cia.re;
  fcinvt.im = fit.im + dit.im - cia.im;
  fct = invr(fcinvt);
}

// Calculate the wave-function modulus sqaured
// for identical bosons (symmetrized) or fermions (antisymmetrized)
// with Coulomb interaction and Strong interaction included
double GetQuantumCoulombStrong()
{
  if (mRSt < 0.0000000001)
    return 1.0;

  if (mRSt < 1.0/0.197327) 
    return GetQuantumCoulomb();
  
  double tKstRst = mKStarOut*mROS + mKStarSide*mRSS + mKStarLong*mRLS;
  long double kstar = fabs(mKStarSigned);
  long double rho = mRSt * kstar;
  
  int ccase = 0;
  int wavesign = 1;

  // Classical limit - if distance is larger than Coulomb radius, 
  // the interaction does not matter
  //  if (fabs(mRSt) > fabs(pionac)) return (1.0 + wavesign*cos(2*tKstRst));

  // Classical limit - in the case of large k* we go to 
  // classical coulomb interaction
  long double testp = rho * (1.0 + tKstRst/(rho));
  long double testm = rho * (1.0 - tKstRst/(rho));

  dcomplex ffplus, ffminus;
  if ((testp> 15.0) && (testm> 15.0))
    {
      double asym;
      asym = (1.0 - 1.0/(mRSt*(1.0 - tKstRst/rho)*pionac*kstar*kstar))/Gamow(kstar);
      //      cout << "as1 " << asym << endl;
      asym = sqrt(asym);
      if (asym < 1.0) 
	ffminus.re = 1.0 + (asym -1.0) *2.0;
      else
	ffminus.re = 1.0 + (asym -1.0) /2.0;
      ffminus.im = sqrt(asym*asym - ffminus.re*ffminus.re);

      asym = (1.0 - 1.0/(mRSt*(1.0 + tKstRst/rho)*pionac*kstar*kstar))/Gamow(kstar);
      //      cout << "as2 " << asym << endl;
      asym = sqrt(asym);
      if (asym < 1.0) 
	ffplus.re = 1.0 + (asym -1.0) *2.0;
      else
	ffplus.re = 1.0 + (asym -1.0) /2.0;
      ffplus.im = sqrt(asym*asym - ffplus.re*ffplus.re);
      
    }
  
  // Check for the classical limit in both functions separately
  else if (((testp< 15.0) && (testm< 15.0))) // ||
    {
      // Calculate the F function
      GetFFdouble(&ffplus, &ffminus);
      ccase = 1;
    }
  else if (testp< 15.0)
    {
      double asym;
      GetFFsingle(&ffplus,1);
      GetFFsingle(&ffminus, -1);
      if ((fabs(ffminus.re) > 2.0) || fabs(ffminus.im) > 2.0) {
	asym = (1.0 - 1.0/(mRSt*(1.0 - tKstRst/(rho)*pionac*kstar*kstar)))/Gamow(kstar);
	asym = sqrt(asym);
	if (asym < 1.0) 
	  ffminus.re = 1.0 + (asym -1.0) *2.0;
	else
	  ffminus.re = 1.0 + (asym -1.0) /2.0;
	ffminus.im = sqrt(asym*asym - ffminus.re*ffminus.re);
      }
      ccase = 2;
    }
  else 
    {
      double asym;
      GetFFsingle(&ffminus, -1);
      GetFFsingle(&ffplus,   1);
      if ((fabs(ffplus.re) > 2.0) || fabs(ffplus.im) > 2.0) {
	asym = (1.0 - 1.0/(mRSt*(1.0 + tKstRst/(rho)*pionac*kstar*kstar)))/Gamow(kstar);
	asym = sqrt(asym);
	if (asym < 1.0) 
	  ffplus.re = 1.0 + (asym -1.0) *2.0;
	else
	  ffplus.re = 1.0 + (asym -1.0) /2.0;
	ffplus.im = sqrt(asym*asym - ffplus.re*ffplus.re);
      }
      ccase = 3;
    }

  long double eta  = 1.0/(kstar*pionac);
  long double hfun = geth(eta);
  dcomplex gtilde  = GetG(eta, rho, hfun);
  dcomplex gtilor  = mult(gtilde, 1.0/mRSt);

  dcomplex fcouls, fcoult;
  Getfc(kstar, eta, hfun, fcouls, fcoult);
  
  dcomplex fgs       = mult(gtilor, fcouls);
  long double fgmods = modl2(fgs);

  dcomplex fgt       = mult(gtilor, fcoult);
  long double fgmodt = modl2(fgt);

  dcomplex expikr;
  expikr.re = cos(tKstRst);
  expikr.im = -sin(tKstRst);
  dcomplex expikrc = conj(expikr);
  dcomplex ffplusc = conj(ffplus);
  dcomplex ffminusc = conj(ffminus);

  dcomplex expikr2 = mult(expikr, expikr);
  dcomplex expikrc2 = conj(expikr2);
  dcomplex sterm = mult(expikr2,  mult(ffplus, ffminusc));
  dcomplex tterm = mult(expikrc2, mult(ffminus, ffplusc));

  dcomplex epfpc = mult(expikr, ffplus);
  dcomplex emfmc = mult(expikrc, ffminus);

  long double fcgefhs = (fgs.re*emfmc.re + fgs.im*emfmc.im);
  long double fcgefgs = (fgs.re*epfpc.re + fgs.im*epfpc.im);

  long double fcgefht = (fgt.re*emfmc.re + fgt.im*emfmc.im);
  long double fcgefgt = (fgt.re*epfpc.re + fgt.im*epfpc.im);

  long double smult = 1+wavesign;

  if (!isfinite(ffminus.re))
    cout << "FFMinus Re not a number ! " << testp << " " << testm << " " << ccase<< endl;
  
  if (!isfinite(ffminus.im))
    cout << "FFMinus Im not a number !" << testp << " " << testm << " " << ccase<< endl;

  if ((ffplus.re > 2.0) || (ffplus.re < -2.0))
    cout << "FFplus Re wild !" << ffplus.re << " case " << ccase << " " << testp << " " << testm << endl;

  if ((ffplus.im > 2.0) || (ffplus.im < -2.0))
    cout << "FFplus Im wild !" << ffplus.im << " case " << ccase << " " << testp << " " << testm << endl;

  if ((ffminus.re > 2.0) || (ffminus.re < -2.0))
    cout << "FFminus Re wild !" << ffminus.re << " case " << ccase << endl;
  
  if ((ffminus.im > 2.0) || (ffminus.im < -2.0))
    cout << "FFminus Im wild !" << ffminus.im << " case " << ccase <<endl;

  if (twospin == 1) {
    wavesign = 1;
    smult = 2;
    long double singlet = (0.5 * Gamow(kstar) *
			   (2.0 * fgmods * smult + modl2(ffplus) + modl2(ffminus) + 
			    wavesign*sterm.re + wavesign*tterm.re +
			    smult * 2 * (fcgefhs + fcgefgs)));
    wavesign = -1;
    smult = 0;
    long double triplet = (0.5 * Gamow(kstar) *
			   (2.0 * fgmodt * smult + modl2(ffplus) + modl2(ffminus) + 
			    wavesign*sterm.re + wavesign*tterm.re +
			    smult * 2 * (fcgefht + fcgefgt)));
    return (0.25 * singlet + 0.75 * triplet);
  }
  else
    return (0.5 * Gamow(kstar) *
	    (2.0 * fgmods * smult + modl2(ffplus) + modl2(ffminus) + 
	     wavesign*sterm.re + wavesign*tterm.re +
	     smult * 2 * (fcgefhs + fcgefgs)));
}

void PairKinematics(ParticleCoor part1, ParticleCoor part2)
{
  // Calculate pair variables
  double tPx = part1.px+part2.px;
  double tPy = part1.py+part2.py;
  double tPz = part1.pz+part2.pz;
  double tE  = part1.e +part2.e;
  double tPt = tPx*tPx + tPy*tPy;
  double tMt = tE*tE - tPz*tPz;//mCVK;
  double tM  = sqrt(tMt - tPt);
  tMt = sqrt(tMt);
  tPt = sqrt(tPt);
  mKT = tPt/2.0;
  pairphi = TMath::ATan2(tPy, tPx);

  if ((mKT<ktmin) || (mKT>ktmax))
    return;

  mBetat = tPt/tMt;
    
  if ((mBetat<btmin) || (mBetat>btmax)) {
    mKT = -1.0;
    return;
  }

  // Boost to LCMS
  double tBeta = tPz/tE;
  double tGamma = tE/tMt;	    
  mKStarLong = tGamma * (part1.pz - tBeta * part1.e);
  double tE1L = tGamma * (part1.e  - tBeta * part1.pz);
  
  // Transform to LCMS
  particle1lcms.z = tGamma * (part1.z - tBeta * part1.t);
  particle1lcms.t = tGamma * (part1.t - tBeta * part1.z);
  
  particle2lcms.z = tGamma * (part2.z - tBeta * part2.t);
  particle2lcms.t = tGamma * (part2.t - tBeta * part2.z);
  
  particle1prf.pz = particle1lcms.pz = tGamma * (part1.pz - tBeta * part1.e);
  particle1lcms.e  = tGamma * (part1.e  - tBeta * part1.pz);
  
  particle2prf.pz = particle2lcms.pz = tGamma * (part2.pz - tBeta * part2.e);
  particle2lcms.e  = tGamma * (part2.e  - tBeta * part2.pz);
  
  // Rotate in transverse plane
  mKStarOut  = ( part1.px*tPx + part1.py*tPy)/tPt;
  mKStarSide = (-part1.px*tPy + part1.py*tPx)/tPt;
      
  particle1lcms.px = mKStarOut;
  particle1lcms.py = mKStarSide;
  
  particle2lcms.px = (part2.px*tPx + part2.py*tPy)/tPt;
  particle2lcms.py = (part2.py*tPx - part2.px*tPy)/tPt;;
  
  mKO = particle1lcms.px - particle2lcms.px;
  mKS = particle1lcms.py - particle2lcms.py;
  mKL = particle1lcms.pz - particle2lcms.pz;
  mDE = particle1lcms.e  - particle2lcms.e;
  
  // save the rotated coordinates in LCMS variables
  particle1lcms.x = ( part1.x*tPx + part1.y*tPy)/tPt;
  particle1lcms.y = (-part1.x*tPy + part1.y*tPx)/tPt;

  particle2lcms.x = ( part2.x*tPx + part2.y*tPy)/tPt;
  particle2lcms.y = (-part2.x*tPy + part2.y*tPx)/tPt;
  
  // Boost to pair cms
  mKStarOut = tMt/tM * (mKStarOut - tPt/tMt * tE1L);
  
  Double_t tBetat = tPt/tMt;
  Double_t tGammat = 1.0/TMath::Sqrt(1.0-tBetat*tBetat);
  
  particle1prf.x = tGammat*(particle1lcms.x - tBetat*particle1lcms.t);
  particle1prf.t = tGammat*(particle1lcms.t - tBetat*particle1lcms.x);
  
  particle2prf.x = tGammat*(particle2lcms.x - tBetat*particle2lcms.t);
  particle2prf.t = tGammat*(particle2lcms.t - tBetat*particle2lcms.x);
  
  mRO = (particle1lcms.x - particle2lcms.x)/0.197327;
  mRS = (particle1lcms.y - particle2lcms.y)/0.197327;
  mRL = (particle1lcms.z - particle2lcms.z)/0.197327;
  mDT = (particle1lcms.t - particle2lcms.t)/0.197327;

  mKO = particle1lcms.px - particle2lcms.px;

  mKStarSigned = mKStarOut>0.? 1. : -1.;
  mKStarSigned *= sqrt(mKStarSide*mKStarSide + mKStarOut*mKStarOut + mKStarLong*mKStarLong);

  mKR = fabs(mKStarSigned);
  if ( mKR < 1e-10 || fabs(mKR) < 1e-10 ) mKC = 0.0;
  else mKC=mKStarLong/mKR;
  mKP=atan2(mKStarSide,mKStarOut);

  double tDX = part1.x-part2.x;
  double tDY = part1.y-part2.y;
  mRLong = part1.z-part2.z;
  mDTime = part1.t-part2.t;

  mRTrans = tDX>0.? ::sqrt(tDX*tDX+tDY*tDY) : -1.*::sqrt(tDX*tDX+tDY*tDY);
  mROut = (tDX*tPx + tDY*tPy)/tPt;
  mRSide = (-tDX*tPy + tDY*tPx)/tPt;

  mRSidePairCMS = mRSide;
  mRSS = mRSidePairCMS/0.197327;

  mRLongPairCMS = tGamma*(mRLong - tBeta* mDTime);
  mDTimePairLCMS = tGamma*(mDTime - tBeta* mRLong);

  mRLS = mRLongPairCMS/0.197327;
  tBeta = tPt/tMt;
  tGamma = tMt/tM;

  mROutPairCMS = tGamma*(mROut - tBeta* mDTimePairLCMS);
  mROS = mROutPairCMS/0.197327;
  mDTimePairCMS = tGamma*(mDTimePairLCMS - tBeta* mROut);

  mRStar = ::sqrt(mROutPairCMS*mROutPairCMS + mRSidePairCMS*mRSidePairCMS +
		  mRLongPairCMS*mRLongPairCMS);
  mRSt = mRStar/0.197327;

  if (writesourcehists) {
    tROutLCMS ->Fill(mROut);
    tRSideLCMS->Fill(mRSide);
    tRLongLCMS->Fill(mRLong);
    tRTimeLCMS->Fill(mDTimePairLCMS);
    tROutSideLongLCMS->Fill(mROut, mRSide, mRLong);
    
    tROutPRF ->Fill(mROutPairCMS);
    tRSidePRF->Fill(mRSidePairCMS);
    tRLongPRF->Fill(mRLongPairCMS);
    tRTimePRF->Fill(mDTimePairCMS);
    tROutSideLongPRF->Fill(mROutPairCMS, mRSidePairCMS, mRLongPairCMS);
    
    tROutLCMSQinv->Fill(mROut, fabs(mKStarSigned)*2.0);
    tRSideLCMSQinv->Fill(mRSide, fabs(mKStarSigned)*2.0);
    tRLongLCMSQinv->Fill(mRLong, fabs(mKStarSigned)*2.0);
    tRTimeLCMSQinv->Fill(mDTimePairLCMS, fabs(mKStarSigned)*2.0);
    
    tROutPRFQinv->Fill(mROutPairCMS, fabs(mKStarSigned)*2.0);
    tRSidePRFQinv->Fill(mRSidePairCMS, fabs(mKStarSigned)*2.0);
    tRLongPRFQinv->Fill(mRLongPairCMS, fabs(mKStarSigned)*2.0);
    tRTimePRFQinv->Fill(mDTimePairCMS, fabs(mKStarSigned)*2.0);
    
    tRInvPRF->Fill(mRStar, 1.0/(mRStar*mRStar));
  }
}

void CopyINIFile()
{
  TString  tINI;
  ifstream ifs;
  ofstream ofs;
  
  tINI = sEventDir + sMainINI;
  ifs.open(sMainINI, std::ios::binary);
  ofs.open(tINI, std::ios::binary);
  if((ifs) && (ofs) && ifs.is_open() && ofs.is_open()) {
    ofs << ifs.rdbuf();
    ifs.close();
    ofs.close();
  } else
    PRINT_MESSAGE("<therm2_femto::main()>\tUnable to copy "<<sMainINI<<" to "<<tINI); 
}
