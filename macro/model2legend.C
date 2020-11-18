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

#include <cstdio>
#include <iostream>
#include <TLegend.h>
#include <TChain.h>
#include <TString.h>
#include "./../build/include/StructModel.h"

using namespace std;

void model2legend(TString aEventDir, Int_t aEventFiles, TLegend* aLegend)
{
  TChain* tChain;
  int     tIntegrateSample;
  int     tRandomize;
  char    tTimeStamp[21];
  int     tModelID;
  char    tModelHash[9];
  char    tModelName[30];
  char    tModelDescription[2000];
  char    Buff[100];
  
  tChain = new TChain(_PARAMETERS_TREE_);
  for(Int_t i=0; i<aEventFiles; i++) {
    sprintf(Buff,"%sevent%03i.root",aEventDir.Data(),i);
    tChain->AddFile(Buff);
  }
  tChain->SetBranchAddress(_INTEGRATESAMPLE_BRANCH_,	&tIntegrateSample);
  tChain->SetBranchAddress(_RANDOMIZE_BRANCH_,		&tRandomize);
  tChain->SetBranchAddress(_TIMESTAMP_BRANCH_,		&tTimeStamp);
  tChain->SetBranchAddress(_MODELID_BRANCH_,		&tModelID);
  tChain->SetBranchAddress(_MODELNAME_BRANCH_,		tModelName);
  tChain->SetBranchAddress(_MODELHASH_BRANCH_,		tModelHash);
  tChain->SetBranchAddress(_MODELDESCRIPTION_BRANCH_,	tModelDescription);
  tChain->GetEntry(0);

  sprintf(Buff,"%s", tModelName); aLegend->SetHeader(Buff);

  switch (tModelID) {
// HYDRO INSPIRED
  case 0: {
    Model_t_KrakowSFO tPar;
    tChain->SetBranchAddress(_MODEL_T_BRANCH_, &tPar);
    tChain->GetEntry(0);
    sprintf(Buff,"Y^{rng} = %g",		tPar.RapPRange); aLegend->AddEntry("",Buff,"");
    sprintf(Buff,"Y_{S}^{rng} = %g",		tPar.RapSRange); aLegend->AddEntry("",Buff,"");
    sprintf(Buff,"#rho_{max}= %g [fm]",	tPar.RhoMax);	aLegend->AddEntry("",Buff,"");
    sprintf(Buff,"#tau_{C}    = %g [fm]",	tPar.TauC);	aLegend->AddEntry("",Buff,"");
    sprintf(Buff,"T     = %g [MeV]",	tPar.Temp);	aLegend->AddEntry("",Buff,"");
    if(tPar.MuB!=0.0) {
      sprintf(Buff,"#mu_{B}    = %g [Mev]",	tPar.MuB);	aLegend->AddEntry("",Buff,"");
    }
    if(tPar.MuI!=0.0) {
      sprintf(Buff,"#mu_{I_{3}}    = %g [Mev]",	tPar.MuI);	aLegend->AddEntry("",Buff,"");
    }
    if(tPar.MuS!=0.0) {
      sprintf(Buff,"#mu_{S}    = %g [Mev]",	tPar.MuS);	aLegend->AddEntry("",Buff,"");
    }
    if(tPar.MuC!=0.0) {
      sprintf(Buff,"#mu_{C}    = %g [Mev]",	tPar.MuC);	aLegend->AddEntry("",Buff,"");
    }
    break;
  }
  case 1: {
    Model_t_BlastWave tPar;
    tChain->SetBranchAddress(_MODEL_T_BRANCH_, &tPar);
    tChain->GetEntry(0);
    sprintf(Buff,"Y^{rng} = %g",		tPar.RapPRange); aLegend->AddEntry("",Buff,"");
    sprintf(Buff,"Y_{S}^{rng} = %g",		tPar.RapSRange); aLegend->AddEntry("",Buff,"");
    sprintf(Buff,"#rho_{max} = %g [fm]",	tPar.RhoMax);	aLegend->AddEntry("",Buff,"");
    sprintf(Buff,"#tau      = %g [fm]",		tPar.Tau);	aLegend->AddEntry("",Buff,"");
    sprintf(Buff,"v_{T}    = %g [c]",		tPar.Vt);	aLegend->AddEntry("",Buff,"");
    sprintf(Buff,"T     = %g [MeV]",	tPar.Temp);	aLegend->AddEntry("",Buff,"");
    if(tPar.MuB!=0.0) {
      sprintf(Buff,"#mu_{B}    = %g [Mev]",	tPar.MuB);	aLegend->AddEntry("",Buff,"");
    }
    if(tPar.MuI!=0.0) {
      sprintf(Buff,"#mu_{I_{3}}    = %g [Mev]",	tPar.MuI);	aLegend->AddEntry("",Buff,"");
    }
    if(tPar.MuS!=0.0) {
      sprintf(Buff,"#mu_{S}    = %g [Mev]",	tPar.MuS);	aLegend->AddEntry("",Buff,"");
    }
    if(tPar.MuC!=0.0) {
      sprintf(Buff,"#mu_{C}    = %g [Mev]",	tPar.MuC);	aLegend->AddEntry("",Buff,"");
    }
    break;
  }
  case 2:
  case 3:
  case 4:
  case 5:
  case 6: {
    Model_t_BWA tPar;
    tChain->SetBranchAddress(_MODEL_T_BRANCH_, &tPar);
    tChain->GetEntry(0);
    sprintf(Buff,"Y^{rng} = %g",		tPar.RapPRange); aLegend->AddEntry("",Buff,"");
    sprintf(Buff,"Y_{S}^{rng} = %g",		tPar.RapSRange); aLegend->AddEntry("",Buff,"");
    sprintf(Buff,"#rho_{max} = %g [fm]",	tPar.RhoMax);	aLegend->AddEntry("",Buff,"");
    sprintf(Buff,"#tau      = %g [fm]",		tPar.Tau);	aLegend->AddEntry("",Buff,"");
    sprintf(Buff,"v_{T}    = %g [c]",		tPar.Vt);	aLegend->AddEntry("",Buff,"");
    sprintf(Buff,"A     = %g [c]",		tPar.A);	aLegend->AddEntry("",Buff,"");
    sprintf(Buff,"#lambda      = %g [c]",	tPar.Delay);	aLegend->AddEntry("",Buff,"");
    sprintf(Buff,"T    = %g [MeV]",	tPar.Temp);	aLegend->AddEntry("",Buff,"");
    if(tPar.MuB!=0.0) {
      sprintf(Buff,"#mu_{B}    = %g [Mev]",	tPar.MuB);	aLegend->AddEntry("",Buff,"");
    }
    if(tPar.MuI!=0.0) {
      sprintf(Buff,"#mu_{I_{3}}    = %g [Mev]",	tPar.MuI);	aLegend->AddEntry("",Buff,"");
    }
    if(tPar.MuS!=0.0) {
      sprintf(Buff,"#mu_{S}    = %g [Mev]",	tPar.MuS);	aLegend->AddEntry("",Buff,"");
    }
    if(tPar.MuC!=0.0) {
      sprintf(Buff,"#mu_{C}    = %g [Mev]",	tPar.MuC);	aLegend->AddEntry("",Buff,"");
    }
    break;
  }
// HYDRO BASED
  case 10:
  {
    Model_t_Lhyquid3D tPar;
    tChain->SetBranchAddress(_MODEL_T_BRANCH_, &tPar);
    tChain->GetEntry(0);
    sprintf(Buff,"Y^{rng} = %g",		tPar.RapPRange); aLegend->AddEntry("",Buff,"");
    sprintf(Buff,"#Lambda = %g [fm]",		tPar.Lambda);	aLegend->AddEntry("",Buff,"");
    sprintf(Buff,"#tau_{i}      = %g [fm]",	tPar.TauI);	aLegend->AddEntry("",Buff,"");
    sprintf(Buff,"T_{i}     = %g [MeV]",	tPar.TempI);	aLegend->AddEntry("",Buff,"");
    sprintf(Buff,"T     = %g [MeV]",		tPar.TempF);	aLegend->AddEntry("",Buff,"");
    if(tPar.MuB!=0.0) {
      sprintf(Buff,"#mu_{B}    = %g [Mev]",	tPar.MuB);	aLegend->AddEntry("",Buff,"");
    }
    if(tPar.MuI!=0.0) {
      sprintf(Buff,"#mu_{I_{3}}    = %g [Mev]",	tPar.MuI);	aLegend->AddEntry("",Buff,"");
    }
    if(tPar.MuS!=0.0) {
      sprintf(Buff,"#mu_{S}    = %g [Mev]",	tPar.MuS);	aLegend->AddEntry("",Buff,"");
    }
    if(tPar.MuC!=0.0) {
      sprintf(Buff,"#mu_{C}    = %g [Mev]",	tPar.MuC);	aLegend->AddEntry("",Buff,"");
    }
    if(tPar.ImpactParameter >= 0.0) {
      sprintf(Buff,"b   = %g [fm]",		tPar.ImpactParameter);	aLegend->AddEntry("",Buff,"");
    }
    if((tPar.CentralityMin >= 0.0) && (tPar.CentralityMax >= 0.0)) {
      sprintf(Buff,"c      = %g-%g [%%]",		tPar.CentralityMin, tPar.CentralityMax);	aLegend->AddEntry("",Buff,"");
    }
    sprintf(Buff,"#sqrt{s_{NN}} = %g MeV",	tPar.CollidingEnergy);	aLegend->AddEntry("",Buff,"");
    sprintf(Buff,"%s %s",	tPar.DeviceName,tPar.CollidingSystem);	aLegend->AddEntry("",Buff,""); 
    break;
  }
  case 11: {
    Model_t_Lhyquid2DBI tPar;
    tChain->SetBranchAddress(_MODEL_T_BRANCH_, &tPar);
    tChain->GetEntry(0);
    sprintf(Buff,"Y^{rng} = %g",		tPar.RapPRange); aLegend->AddEntry("Y_{rng}",Buff,"");
    sprintf(Buff,"Y_{S}^{rng} = %g",		tPar.RapSRange); aLegend->AddEntry("",Buff,"");
    sprintf(Buff,"#tau_{i}      = %g [fm]",	tPar.TauI);	aLegend->AddEntry("",Buff,"");
    sprintf(Buff,"T_{i}     = %g [MeV]",	tPar.TempI);	aLegend->AddEntry("",Buff,"");
    sprintf(Buff,"T     = %g [MeV]",	tPar.TempF);	aLegend->AddEntry("",Buff,"");
    if(tPar.MuB!=0.0) {
      sprintf(Buff,"#mu_{B}    = %g [Mev]",	tPar.MuB);	aLegend->AddEntry("",Buff,"");
    }
    if(tPar.MuI!=0.0) {
      sprintf(Buff,"#mu_{I_{3}}    = %g [Mev]",	tPar.MuI);	aLegend->AddEntry("",Buff,"");
    }
    if(tPar.MuS!=0.0) {
      sprintf(Buff,"#mu_{S}    = %g [Mev]",	tPar.MuS);	aLegend->AddEntry("",Buff,"");
    }
    if(tPar.MuC!=0.0) {
      sprintf(Buff,"#mu_{C}    = %g [Mev]",	tPar.MuC);	aLegend->AddEntry("",Buff,"");
    }
    if(tPar.ImpactParameter >= 0.0) {
      sprintf(Buff,"b   = %g [fm]",		tPar.ImpactParameter);				aLegend->AddEntry("",Buff,"");
    }
    if((tPar.CentralityMin >= 0.0) && (tPar.CentralityMax >= 0.0)) {
      sprintf(Buff,"c      = %g-%g [%%]",	tPar.CentralityMin, tPar.CentralityMax);	aLegend->AddEntry("",Buff,"");
    }
    sprintf(Buff,"#sqrt{s_{NN}} = %g MeV",	tPar.CollidingEnergy);	aLegend->AddEntry("",Buff,"");
    sprintf(Buff,"%s %s",	tPar.DeviceName,tPar.CollidingSystem);	aLegend->AddEntry("",Buff,""); 
    break;
  }
  default:
    cerr << "Unable to write Model parameter - unknown model type.";
  };
}

/*! @file model2legend.C
 * @brief ROOT macro with model2legend() function.
 */
/*! @fn void model2legend(TString aEventDir, Int_t aEventFiles, TLegend* aLegend)
 * @brief Creates a legend with Model parameters. Function adds entries to an external TLegend object.
 * @param [in] aEventDir directory with event files
 * @param [in] aEventFiles number of event files
 * @param [out] aLegend a pointer to a TLegend object
 */
