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

#include <iostream>
#include <fstream>
#include <sstream>
#include <TH1D.h>
#include <TH3D.h>
#include <TFile.h>
#include <TF1.h>
#include <TF3.h>
#include <TGraph.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TImage.h>
#include <TMath.h>
#include <TDatime.h>
#include <math.h>
#include "Configurator.h"
#include "Parser.h"
#include "THGlobal.h"

using namespace std;

Configurator* sMainConfig;
TString  sMainINI;
TString  sEventDir;
TString  sTimeStamp;
int      sParentPID;

TH3D    *ratq;
TH1D    *hpxc;
TH1D    *hpyc;
TH1D    *hpzc;
TF3     *funqg;
TF3     *funqk;
TGraph  *kcoul;

void		getfitprojc(TH3D *expden, TH3D **projhist);
void		preparepad();
void		preparehist(TH1D *hist, int type);
TH1D*		getproj(TH3D *numq, TH3D *denq, int nproj, int wbin, double norm);
Double_t	fungek(Double_t *x, Double_t *par);
int		GetFitParameter(TString aKeyword, Double_t *parval, Int_t *isfixed, Double_t *parmin, Double_t *parmax);
void		MessageIntro();
void		MessageHelp();
void		MessageVersion();
void		CopyINIFile();

int main(int argc, char **argv)
{
  TDatime  tDate;
  TString  tInRootName;
  Int_t    dofix[5];  
  Double_t pars[5];
  Double_t parmin[5];
  Double_t parmax[5]; 
  Double_t maxrange;
  Double_t maxx=0.0, maxy=0.0, maxz=0.0;
  TString  numname;
  TString  denname;
  
  tDate.Set();
  sTimeStamp = tDate.AsSQLString();
  
// ##############################################################
// # Command line analysis					#
// ##############################################################
  
  sParentPID = 0;
  
  if (argc > 1) {
    tInRootName = argv[1];
    if(tInRootName.Contains("-h") || tInRootName.Contains("--help")) {
      MessageHelp();
      return 0;
    } else if(tInRootName.Contains("-v") || tInRootName.Contains("--version")) {
      MessageVersion();
      return 0;
    }
    if (argc > 2) {
      sMainINI	= argv[2];
      if(sMainINI.IsDigit()) {
	sParentPID = sMainINI.Atoi();
	sMainINI   = "./hbtfit.ini";
      }
    } else
      sMainINI	= "./hbtfit.ini";
    if (argc > 3)
      sParentPID = atoi(argv[3]);
  } else {
    MessageHelp();
    return _ERROR_GENERAL_FILE_NOT_FOUND_;
  }

  MessageIntro();  

  sEventDir = tInRootName;
  {
    int tPos=0;
    for(int i=1; i<sEventDir.Length(); i++) {
      if(sEventDir[i]=='/')
        tPos = i;
    }
    sEventDir.Resize(tPos+1);
  }

// ##############################################################
// # Read configuration file (hbtfit.ini)			#
// ############################################################## 
  sMainConfig = new Configurator;
  Parser* tParser = new Parser(sMainINI);
  tParser->ReadINI(sMainConfig);
  delete tParser;
  
// ##############################################################
// # ReadParameters						#
// ##############################################################  
  try {
    GetFitParameter("Norm"  ,&pars[0], &dofix[0], &parmin[0], &parmax[0]);
    GetFitParameter("Lambda",&pars[1], &dofix[1], &parmin[1], &parmax[1]);
    GetFitParameter("Rout"  ,&pars[2], &dofix[2], &parmin[2], &parmax[2]);
    GetFitParameter("Rside" ,&pars[3], &dofix[3], &parmin[3], &parmax[3]);
    GetFitParameter("Rlong" ,&pars[4], &dofix[4], &parmin[4], &parmax[4]);
    maxrange = sMainConfig->GetParameter("MaxFitRange").Atof();
    numname  = sMainConfig->GetParameter("Numerator");
    denname  = sMainConfig->GetParameter("Denominator");
  }
  catch (TString tError) {
    PRINT_DEBUG_1("therm2_hbtfit - Caught exception " << tError);
    PRINT_MESSAGE("Did not find one of the necessary parameters in the parameters file.");
    exit(_ERROR_CONFIG_PARAMETER_NOT_FOUND_);
  }
  try {
    maxx = sMainConfig->GetParameter("MaxFitRangeX").Atof();
    maxy = sMainConfig->GetParameter("MaxFitRangeY").Atof();
    maxz = sMainConfig->GetParameter("MaxFitRangeZ").Atof();
  }
  catch (TString tError) {
  }
  
  CopyINIFile();
// ##############################################################
// # Log file entry 						#
// ##############################################################
  {
    TString tLogName;
    ofstream tLogFile;   
    
    try {
      tLogName = sMainConfig->GetParameter("LogFile");
    }
    catch (TString tError) {
    }   
    if (!tLogName.IsNull()) {
      tDate.Set();
      
      tLogFile.open(tLogName, ios_base::app);
      if (static_cast<long>(tLogFile.tellp()) == 0) {
        tLogFile << "# THERMINATOR 2 Log File"<<endl;
      }
      tLogFile << '['<<tDate.AsSQLString()<<"]\ttherm2_hbtfit\t"<<sParentPID<<'\t';
      tLogFile << "[input]\t"<<tInRootName<<'\t'<<sMainINI<<'\t'<<sParentPID<<'\t';
      tLogFile << endl;
      tLogFile.close();
    }
  }
  
// ##############################################################
// # Fit							#
// ############################################################## 
  TFile* tInRootFile;
  TH3D*  numq;
  TH3D*  denq;

  PRINT_MESSAGE("["<<sTimeStamp<<"]\tFitting.");
  
  funqg = new TF3("funqg", fungek, -0.15, 0.15, -0.15, 0.15, -0.15, 0.15, 5);
  if (maxx > 0)
    funqg->SetRange(0, 0, 0, maxx, maxy, maxz);
  else
    funqg->SetRange(-maxrange, -maxrange, -maxrange, maxrange, maxrange, maxrange);
  
  funqk = funqg;

  for (int iter=0; iter<5; iter++)
    if (dofix[iter]) {
      if (dofix[iter] == 2) {
	funqk->SetParLimits(iter, parmin[iter], parmax[iter]);
	funqk->SetParameter(iter, pars[iter]);
      }
      else
	funqk->FixParameter(iter, pars[iter]);
    }
    else
      funqk->SetParameter(iter, pars[iter]);

  tInRootFile = new TFile(tInRootName);
  numq = new TH3D(*((TH3D *) tInRootFile->Get(numname)));
  denq = new TH3D(*((TH3D *) tInRootFile->Get(denname)));
  ratq = new TH3D(*numq);
  ratq->Reset("ICE");
  ratq->Divide(numq, denq, 1.0, 1.0);
  ratq->SetName("ratq");
  ratq->SetTitle("ratq");
  ratq->Fit(funqk, "RB");

// ##############################################################
// # Save fit values						#
// ############################################################## 
  TString   tOutTextName;
  ofstream* tOutTextFile;
  
  tOutTextName = sEventDir + "hbtradii.txt";
  tOutTextFile = new ofstream(tOutTextName, ios_base::app);
  (*tOutTextFile) << "["<<sTimeStamp<<"]\t" << tInRootName << endl;
  (*tOutTextFile) << "Norm   " << fabs(funqk->GetParameter(0)) << " +/- " <<fabs(funqk->GetParError(0)) << endl;
  (*tOutTextFile) << "Lambda " << fabs(funqk->GetParameter(1)) << " +/- " <<fabs(funqk->GetParError(1)) << endl;
  (*tOutTextFile) << "Rout   " << fabs(funqk->GetParameter(2)) << " +/- " <<fabs(funqk->GetParError(2)) << endl;
  (*tOutTextFile) << "Rside  " << fabs(funqk->GetParameter(3)) << " +/- " <<fabs(funqk->GetParError(3)) << endl;
  (*tOutTextFile) << "Rlong  " << fabs(funqk->GetParameter(4)) << " +/- " <<fabs(funqk->GetParError(4)) << endl;
  tOutTextFile->close();
  delete tOutTextFile;

  PRINT_DEBUG_1("Fit results:");
  PRINT_DEBUG_1("\tNorm   " << fabs(funqk->GetParameter(0)) << " +/- " <<fabs(funqk->GetParError(0)));
  PRINT_DEBUG_1("\tlambda " << fabs(funqk->GetParameter(1)) << " +/- " <<fabs(funqk->GetParError(1)));
  PRINT_DEBUG_1("\tRout   " << fabs(funqk->GetParameter(2)) << " +/- " <<fabs(funqk->GetParError(2)));
  PRINT_DEBUG_1("\tRside  " << fabs(funqk->GetParameter(3)) << " +/- " <<fabs(funqk->GetParError(3)));
  PRINT_DEBUG_1("\tRlong  " << fabs(funqk->GetParameter(4)) << " +/- " <<fabs(funqk->GetParError(4)));
 
// ##############################################################
// # Make plots							#
// ############################################################## 
  TH3D *fitnq;

  TH1D *hpx1 = getproj(numq, denq, 0, 0, pars[0]);
  TH1D *hpy1 = getproj(numq, denq, 1, 0, pars[0]);
  TH1D *hpz1 = getproj(numq, denq, 2, 0, pars[0]);
  
  TH1D *hpx2 = getproj(numq, denq, 0, 4, pars[0]);
  TH1D *hpy2 = getproj(numq, denq, 1, 4, pars[0]);
  TH1D *hpz2 = getproj(numq, denq, 2, 4, pars[0]);
  
  TH1D *hpx3 = getproj(numq, denq, 0, 10, pars[0]);
  TH1D *hpy3 = getproj(numq, denq, 1, 10, pars[0]);
  TH1D *hpz3 = getproj(numq, denq, 2, 10, pars[0]);
  
  getfitprojc(denq, &fitnq);

  TH1D *fpx1 = getproj(fitnq, denq, 0, 0, pars[0]);
  TH1D *fpy1 = getproj(fitnq, denq, 1, 0, pars[0]);
  TH1D *fpz1 = getproj(fitnq, denq, 2, 0, pars[0]);
 
  TH1D *fpx2 = getproj(fitnq, denq, 0, 4, pars[0]);
  TH1D *fpy2 = getproj(fitnq, denq, 1, 4, pars[0]);
  TH1D *fpz2 = getproj(fitnq, denq, 2, 4, pars[0]);
  
  TH1D *fpx3 = getproj(fitnq, denq, 0, 10, pars[0]);
  TH1D *fpy3 = getproj(fitnq, denq, 1, 10, pars[0]);
  TH1D *fpz3 = getproj(fitnq, denq, 2, 10, pars[0]);
  
  hpx1->SetName("CFproj_out_1");
  hpy1->SetName("CFproj_side_1");
  hpz1->SetName("CFproj_long_1");

  hpx2->SetName("CFproj_out_2");
  hpy2->SetName("CFproj_side_2");
  hpz2->SetName("CFproj_long_2");

  hpx3->SetName("CFproj_out_3");
  hpy3->SetName("CFproj_side_3");
  hpz3->SetName("CFproj_long_3");

  fpx1->SetName("FITproj_out_1");
  fpy1->SetName("FITproj_side_1");
  fpz1->SetName("FITproj_long_1");

  fpx2->SetName("FITproj_out_2");
  fpy2->SetName("FITproj_side_2");
  fpz2->SetName("FITproj_long_2");

  fpx3->SetName("FITproj_out_3");
  fpy3->SetName("FITproj_side_3");
  fpz3->SetName("FITproj_long_3");

  preparehist(hpx1,1);
  preparehist(hpy1,1);
  preparehist(hpz1,1);

  preparehist(hpx2,1);
  preparehist(hpy2,1);
  preparehist(hpz2,1);

  preparehist(hpx3,1);
  preparehist(hpy3,1);
  preparehist(hpz3,1);

  preparehist(fpx1,3);
  preparehist(fpy1,3);
  preparehist(fpz1,3);

  preparehist(fpx2,3);
  preparehist(fpy2,3);
  preparehist(fpz2,3);

  preparehist(fpx3,3);
  preparehist(fpy3,3);
  preparehist(fpz3,3);

  TCanvas *cancf = new TCanvas ("cancf", "cancf", 1200, 1200);
  gPad->SetFillColor(0);
  gPad->SetFillStyle(4000);
  cancf->Divide(3,3,0.0001,0.0001);
  
  hpx1->GetXaxis()->SetTitle("q_{out} [GeV/c]");
  hpy1->GetXaxis()->SetTitle("q_{side} [GeV/c]");
  hpz1->GetXaxis()->SetTitle("q_{long} [GeV/c]");

  hpx1->GetYaxis()->SetTitle("C(q_{out})");
  hpy1->GetYaxis()->SetTitle("C(q_{side})");
  hpz1->GetYaxis()->SetTitle("C(q_{long})");

  hpx2->GetXaxis()->SetTitle("q_{out} [GeV/c]");
  hpy2->GetXaxis()->SetTitle("q_{side} [GeV/c]");
  hpz2->GetXaxis()->SetTitle("q_{long} [GeV/c]");

  hpx2->GetYaxis()->SetTitle("C(q_{out})");
  hpy2->GetYaxis()->SetTitle("C(q_{side})");
  hpz2->GetYaxis()->SetTitle("C(q_{long})");

  hpx3->GetXaxis()->SetTitle("q_{out} [GeV/c]");
  hpy3->GetXaxis()->SetTitle("q_{side} [GeV/c]");
  hpz3->GetXaxis()->SetTitle("q_{long} [GeV/c]");

  hpx3->GetYaxis()->SetTitle("C(q_{out})");
  hpy3->GetYaxis()->SetTitle("C(q_{side})");
  hpz3->GetYaxis()->SetTitle("C(q_{long})");

  cancf->cd(1);
  preparepad();

  hpx1->Draw("HISTPE1");
  fpx1->Draw("SAMEHISTL");

  cancf->cd(2);
  preparepad();

  hpy1->Draw("HISTPE1");
  fpy1->Draw("SAMEHISTL");

  cancf->cd(3);
  preparepad();

  hpz1->Draw("HISTPE1");
  fpz1->Draw("SAMEHISTL");

  cancf->cd(4);
  preparepad();

  hpx2->Draw("HISTPE1");
  fpx2->Draw("SAMEHISTL");

  cancf->cd(5);
  preparepad();

  hpy2->Draw("HISTPE1");
  fpy2->Draw("SAMEHISTL");

  cancf->cd(6);
  preparepad();

  hpz2->Draw("HISTPE1");
  fpz2->Draw("SAMEHISTL");

  cancf->cd(7);
  preparepad();

  hpx3->Draw("HISTPE1");
  fpx3->Draw("SAMEHISTL");

  cancf->cd(8);
  preparepad();

  hpy3->Draw("HISTPE1");
  fpy3->Draw("SAMEHISTL");

  cancf->cd(9);
  preparepad();

  hpz3->Draw("HISTPE1");
  fpz3->Draw("SAMEHISTL");

// ##############################################################
// # Save plots							#
// ##############################################################   
  
  TString   tOutRootName;
  TFile*    tOutRootFile;
  tOutRootName = tInRootName;
  tOutRootName.ReplaceAll("femto","hbtfit");
  tOutRootFile = new TFile(tOutRootName,"RECREATE");
  tOutRootFile->cd();
 
  hpx1->Write();
  hpy1->Write();
  hpz1->Write();

  hpx2->Write();
  hpy2->Write();
  hpz2->Write();

  hpx3->Write();
  hpy3->Write();
  hpz3->Write();

  fpx1->Write();
  fpy1->Write();
  fpz1->Write();

  fpx2->Write();
  fpy2->Write();
  fpz2->Write();

  fpx3->Write();
  fpy3->Write();
  fpz3->Write();

  tOutRootFile->Close();
  tDate.Set(); 
  PRINT_MESSAGE("["<<tDate.AsSQLString()<<"]\tFile "<<tOutTextName<<" and");
  PRINT_MESSAGE("\tfile "<<tOutRootName<<" written.");

// ##############################################################
// # Log file entry 						#
// ##############################################################
  {
    TString  tLogName;
    TString  tORootName;
    TString  tOTextName;
    ofstream tLogFile;   
    
    try {
      tLogName = sMainConfig->GetParameter("LogFile");
    }
    catch (TString tError) {
    }   
    if (!tLogName.IsNull()) {
      tDate.Set();
      tORootName = tOutRootName;
      tOTextName = tOutTextName;
      tORootName.ReplaceAll(sEventDir,"");      
      tOTextName.ReplaceAll(sEventDir,"");
      
      tLogFile.open(tLogName, ios_base::app);
      if (static_cast<long>(tLogFile.tellp()) == 0) {
        tLogFile << "# THERMINATOR 2 Log File"<<endl;
      }
      tLogFile << '['<<tDate.AsSQLString()<<"]\ttherm2_hbtfit\t"<<sParentPID<<'\t';
      tLogFile << "[output]\t"<<sEventDir<<'\t'<<tORootName<<'\t'<<tOTextName<<'\t';
      tLogFile << endl;
      tLogFile.close();
    }    
  }
  
// ##############################################################
// # Temporary file						#
// ##############################################################  
  ofstream tTmpFile;
  char     tTmpFileName[kFileNameMaxChar];
   
  sprintf(tTmpFileName,"./hbtfit_%i.tmp",sParentPID);
  tTmpFile.open(tTmpFileName, ios_base::app);
  if((tTmpFile) && tTmpFile.is_open()) {
    tTmpFile << tOutRootName << endl;
  } else {
    PRINT_MESSAGE("Unable to create temp file "<<tTmpFileName);
    exit(_ERROR_GENERAL_FILE_NOT_FOUND_);
  }
  
  return 0;
}



void getfitprojc(TH3D *expden, TH3D **projhist)
{
  (*projhist) = new TH3D(*expden);

  for (int q1int = 1; q1int <= expden->GetNbinsX(); q1int++)
    for (int q2int = 1; q2int <= expden->GetNbinsY(); q2int++) 
      for (int q3int = 1; q3int <= expden->GetNbinsZ(); q3int++) {
	(*projhist)->SetBinContent(q1int, q2int, q3int,
				   expden->GetBinContent(q1int, q2int, q3int) *
				   funqk->Eval(expden->GetXaxis()->GetBinCenter(q1int),
					       expden->GetYaxis()->GetBinCenter(q2int),
					       expden->GetZaxis()->GetBinCenter(q3int)));
      }
}

void preparepad()
{
  gPad->SetFillStyle(4000);
  gPad->SetFillColor(0);
  gPad->SetRightMargin(0.02);
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.12);
  gPad->SetGridy();
}

void preparehist(TH1D *hist, int type)
{
  if ((type==2) || (type == 3)) {
    if (type == 2)
      hist->SetLineColor(4);
    else
      hist->SetLineColor(8);
    hist->SetLineWidth(2);
    hist->SetLineStyle(2);
  }
  else {
    hist->SetMarkerSize(1.0);
    hist->SetMarkerColor(type+1);
    hist->SetMarkerStyle(type==0?8:24);
    hist->SetTitle("");
    hist->SetMinimum(hist->GetMinimum()*0.9);
    hist->SetMaximum(hist->GetMaximum()*1.1);
    hist->GetXaxis()->SetLabelSize(0.055);
    hist->GetYaxis()->SetLabelSize(0.055);
    hist->GetXaxis()->SetTitleSize(0.055);
    hist->GetYaxis()->SetTitleSize(0.055);
    hist->GetYaxis()->SetTitleOffset(0.8);
  }
}

TH1D *getproj(TH3D *numq, TH3D *denq, int nproj, int wbin, double norm)
{
  Double_t intnexp = norm;
  Double_t intdexp = 1.0;

  Int_t binc = numq->GetXaxis()->FindFixBin(0.0);
  Int_t binm = binc + wbin;

  TH1D *denbuf;
  TH1D *numbuf;

  for (int ia=nproj; ia<=nproj; ia++) {

    numq->GetXaxis()->SetRange(binc, binm);
    numq->GetYaxis()->SetRange(binc, binm);
    numq->GetZaxis()->SetRange(binc, binm);
    denq->GetXaxis()->SetRange(binc, binm);
    denq->GetYaxis()->SetRange(binc, binm);
    denq->GetZaxis()->SetRange(binc, binm);
    switch (ia) {
    case 0:
      denq->GetXaxis()->SetRange(1,ratq->GetNbinsX());
      numq->GetXaxis()->SetRange(1,ratq->GetNbinsX());
      denbuf = new TH1D(*((TH1D *) denq->Project3D("x")));
      numbuf = new TH1D(*((TH1D *) numq->Project3D("x")));
      hpxc = new TH1D(*numbuf);
      hpxc->Sumw2();
      hpxc->Reset("ICE");
      hpxc->Divide(numbuf, denbuf, 1.0, 1.0, "");
      for (int iter=1; iter<hpxc->GetNbinsX(); iter++)
	if (numbuf->GetBinContent(iter)) {
	  Double_t dn = numbuf->GetBinError(iter);
	  Double_t an = numbuf->GetBinContent(iter);
	  Double_t dd = denbuf->GetBinError(iter);
	  Double_t ad = denbuf->GetBinContent(iter);
	  hpxc->SetBinError(iter, TMath::Sqrt((dn*dn*ad*ad + dd*dd*an*an + dd*dd*dn*dn)/(ad*ad*ad*ad)));
	}
	else 
	  hpxc->SetBinContent(iter,0);
      hpxc->Scale(intdexp/intnexp);
      break;
    case 1:
      denq->GetYaxis()->SetRange(1,ratq->GetNbinsY());
      numq->GetYaxis()->SetRange(1,ratq->GetNbinsY());
      denbuf = new TH1D(*((TH1D *) denq->Project3D("y")));
      numbuf = new TH1D(*((TH1D *) numq->Project3D("y")));
      hpyc = new TH1D(*numbuf);
      hpyc->Sumw2();
      hpyc->Reset("ICE");
      hpyc->Divide(numbuf, denbuf, 1.0, 1.0, "");
      for (int iter=1; iter<hpyc->GetNbinsX(); iter++)
	if (numbuf->GetBinContent(iter)) {
	  Double_t dn = numbuf->GetBinError(iter);
	  Double_t an = numbuf->GetBinContent(iter);
	  Double_t dd = denbuf->GetBinError(iter);
	  Double_t ad = denbuf->GetBinContent(iter);
	  hpyc->SetBinError(iter, TMath::Sqrt((dn*dn*ad*ad + dd*dd*an*an + dd*dd*dn*dn)/(ad*ad*ad*ad)));
	}
	else 
	  hpyc->SetBinContent(iter,0);
      hpyc->Scale(intdexp/intnexp);
      break;
    case 2:
      denq->GetZaxis()->SetRange(1,ratq->GetNbinsZ());
      numq->GetZaxis()->SetRange(1,ratq->GetNbinsZ());
      denbuf = new TH1D(*((TH1D *) denq->Project3D("z")));
      numbuf = new TH1D(*((TH1D *) numq->Project3D("z")));
      hpzc = new TH1D(*numbuf);
      hpzc->Sumw2();
      hpzc->Reset("ICE");
      hpzc->Divide(numbuf, denbuf, 1.0, 1.0, "");
      for (int iter=1; iter<hpzc->GetNbinsX(); iter++)
	if (numbuf->GetBinContent(iter)) {
	  Double_t dn = numbuf->GetBinError(iter);
	  Double_t an = numbuf->GetBinContent(iter);
	  Double_t dd = denbuf->GetBinError(iter);
	  Double_t ad = denbuf->GetBinContent(iter);
	  hpzc->SetBinError(iter, TMath::Sqrt((dn*dn*ad*ad + dd*dd*an*an + dd*dd*dn*dn)/(ad*ad*ad*ad)));
	}
	else 
	  hpzc->SetBinContent(iter,0);
      hpzc->Scale(intdexp/intnexp);
      break;
    }
  }
  if (nproj == 0)
    return hpxc;
  else if (nproj == 1)
    return hpyc;
  else
    return hpzc;
}

Double_t fungek(Double_t *x, Double_t *par)
{
  Double_t qosq = x[0]*x[0];
  Double_t qssq = x[1]*x[1];
  Double_t qlsq = x[2]*x[2];

  Double_t paq;

  Double_t lam =  TMath::Abs(par[1]);

  Double_t gpart = exp((-par[2]*par[2]*qosq-par[3]*par[3]*qssq-par[4]*par[4]*qlsq)/0.038938);

  return (par[0] * (1 + lam*gpart));
}

int GetFitParameter(TString aKeyword, Double_t *parval, Int_t *isfixed, Double_t *parmin, Double_t *parmax)
{
  *parval = sMainConfig->GetParameter(aKeyword).Atof(); 
  if (sMainConfig->GetParameter(aKeyword+"FitType") == "fixed") {
    *isfixed = 1;
    PRINT_DEBUG_1("Fitting parameter: "<<aKeyword<<" = "<<(*parval)<<" fixed.");
  } else if (sMainConfig->GetParameter(aKeyword+"FitType") == "limit") {
    *parmin  = sMainConfig->GetParameter(aKeyword+"Min").Atof();
    *parmax  = sMainConfig->GetParameter(aKeyword+"Max").Atof();
    *isfixed = 2;
    PRINT_DEBUG_1("Fitting parameter: "<<aKeyword<<" = "<<(*parval)<<" limited ["<<(*parmin)<<", "<<(*parmax)<<"]");
  } else {
    *isfixed = 0;
    PRINT_DEBUG_1("Fitting parameter: "<<aKeyword<<" = "<<(*parval)<<" free.");
  }
}

void MessageIntro()
{
  PRINT_MESSAGE("  ***********************************************************************"	);
  PRINT_MESSAGE("  *\t\tTHERMINATOR 2 HBTFIT version "<<_THERMINATOR2_VERSION_<<"\t\t\t*"	);
  PRINT_MESSAGE("  *\t\t\t\t\t\t\t\t\t*"							);
  PRINT_MESSAGE("  * authors: M.Chojnacki, A.Kisiel, W.Florkowski, W.Broniowski\t\t*"		);
  PRINT_MESSAGE("  * cite as: arXiv:1102.0273\t\t\t\t\t\t*"					);
  PRINT_MESSAGE("  * webpage: http://therminator2.ifj.edu.pl/\t\t\t\t*"				);
  PRINT_MESSAGE("  ***********************************************************************"	);
}

void MessageHelp()
{
  PRINT_MESSAGE("Usage: ");
  PRINT_MESSAGE("therm2_hbtfit <CORR_FILE> [HBTFIT_INI] [PPID]");
  PRINT_MESSAGE("therm2_hbtfit [OPTION]");
  PRINT_MESSAGE("  <CORR_FILE>\t\tROOT file with the two particle correlation function");
  PRINT_MESSAGE("  [HBTFIT_INI]\t\tmain settings file;\t\tdefault: hbtfit.ini");
  PRINT_MESSAGE("  [PPID]\t\tparent's system process ID;\tdefault: 0");
  PRINT_MESSAGE("  [OPTION]");
  PRINT_MESSAGE("    -h | --help\t\tthis screen");
  PRINT_MESSAGE("    -v | --version\tversion information");
}

void MessageVersion()
{
  PRINT_MESSAGE("version:\tTHERMINATOR 2 HBTFIT version "<<_THERMINATOR2_VERSION_);
  PRINT_MESSAGE("compiled with:\t"<<_CXX_VER_<<", ROOT("<<_ROOT_VER_<<")");
  std::cout <<  "  preprocessor: ";
#ifdef _DEBUG_LEVEL_
  std::cout << "DEBUG="<<_DEBUG_LEVEL_<<" ";
#endif
  std::cout << std::endl;
}

void CopyINIFile()
{
  TString  tINI;
  ifstream ifs;
  ofstream ofs;
  
  tINI = sEventDir + sMainINI;
  ifs.open(sMainINI.Data(), std::ios::binary);
  ofs.open(tINI.Data(), std::ios::binary);
  if((ifs) && (ofs) && ifs.is_open() && ofs.is_open()) {
    ofs << ifs.rdbuf();
    ifs.close();
    ofs.close();
  } else
    PRINT_MESSAGE("<therm2_events::main()>\tUnable to copy "<<sMainINI<<" to "<<tINI);
}
