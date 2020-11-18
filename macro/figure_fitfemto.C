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

#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TFile.h>
#include <TLatex.h>

void PrepCF(TH1D* aH1)
{
  Int_t tMinBin, tMaxBin;
  aH1->SetMarkerColor(46);
  aH1->SetMarkerStyle(24);
  aH1->GetXaxis()->SetTitleSize(0.1);
  aH1->GetXaxis()->CenterTitle(kTRUE);
  aH1->GetXaxis()->SetNdivisions(4,5,0);
  aH1->GetXaxis()->SetLabelSize(0.1);
  aH1->GetYaxis()->SetTitleSize(0.1);
  aH1->GetYaxis()->SetTitleOffset(1.05);
  aH1->GetYaxis()->CenterTitle(kTRUE);
  aH1->GetYaxis()->SetNdivisions(3,5,0);
  aH1->GetYaxis()->SetLabelSize(0.1); 
  tMinBin = aH1->GetMinimumBin();	aH1->SetMinimum((aH1->GetBinContent(tMinBin) - aH1->GetBinError(tMinBin)) * 0.95);
  tMaxBin = aH1->GetMaximumBin();	aH1->SetMaximum((aH1->GetBinContent(tMaxBin) + aH1->GetBinError(tMaxBin)) * 1.05);
}

void PrepFT(TH1D* aH1)
{
  aH1->SetLineColor(9);
  aH1->SetLineStyle(1);
}

void PrePad()
{
  gPad->SetFillColor(0);
  gPad->SetFillStyle(4000);
  gPad->SetTopMargin(0.05);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.25);
  gPad->SetLeftMargin(0.22);
}

void figure_fitfemto(TString aFileHBTFit)
{
  gStyle->SetOptStat("");
  gStyle->SetStatBorderSize(0);
  gStyle->SetTextFont(42);

  TFile*   tFile   = new TFile(aFileHBTFit);
  TCanvas* tCanvas = new TCanvas ("canproj1","Projections", 800, 600);
  PrePad();
  tCanvas->Divide(3,3,0.001, 0.001);

  TH1D* hCF[3][3];
  TH1D* hFT[3][3];

  hCF[0][0] = (TH1D*) tFile->Get("CFproj_out_1");	PrepCF(hCF[0][0]);
  hCF[0][1] = (TH1D*) tFile->Get("CFproj_out_2");	PrepCF(hCF[0][1]);
  hCF[0][2] = (TH1D*) tFile->Get("CFproj_out_3");	PrepCF(hCF[0][2]);
  hCF[1][0] = (TH1D*) tFile->Get("CFproj_side_1");	PrepCF(hCF[1][0]);
  hCF[1][1] = (TH1D*) tFile->Get("CFproj_side_2");	PrepCF(hCF[1][1]);
  hCF[1][2] = (TH1D*) tFile->Get("CFproj_side_3");	PrepCF(hCF[1][2]);
  hCF[2][0] = (TH1D*) tFile->Get("CFproj_long_1");	PrepCF(hCF[2][0]);
  hCF[2][1] = (TH1D*) tFile->Get("CFproj_long_2");	PrepCF(hCF[2][1]);
  hCF[2][2] = (TH1D*) tFile->Get("CFproj_long_3");	PrepCF(hCF[2][2]);

  hFT[0][0] = (TH1D*) tFile->Get("FITproj_out_1");	PrepFT(hFT[0][0]);
  hFT[0][1] = (TH1D*) tFile->Get("FITproj_out_2");	PrepFT(hFT[0][1]);
  hFT[0][2] = (TH1D*) tFile->Get("FITproj_out_3");	PrepFT(hFT[0][2]);
  hFT[1][0] = (TH1D*) tFile->Get("FITproj_side_1");	PrepFT(hFT[1][0]);
  hFT[1][1] = (TH1D*) tFile->Get("FITproj_side_2");	PrepFT(hFT[1][1]);
  hFT[1][2] = (TH1D*) tFile->Get("FITproj_side_3");	PrepFT(hFT[1][2]);
  hFT[2][0] = (TH1D*) tFile->Get("FITproj_long_1");	PrepFT(hFT[2][0]);
  hFT[2][1] = (TH1D*) tFile->Get("FITproj_long_2");	PrepFT(hFT[2][1]);
  hFT[2][2] = (TH1D*) tFile->Get("FITproj_long_3");	PrepFT(hFT[2][2]);

  TLatex tLatex;
  tLatex.SetTextSize(0.095);

// Out projection 1 bin
  tCanvas->cd(1);	PrePad();
  hCF[0][0]->Draw();
  hFT[0][0]->Draw("SAMEHISTL");
  tLatex.DrawLatex(-0.18, 2.0, "Projections     1 bin");
// Out projection 4 bins
  tCanvas->cd(2);	PrePad();
  hCF[0][1]->Draw();
  hFT[0][1]->Draw("SAMEHISTL");
  tLatex.DrawLatex(-0.18, 1.8, "Projections     4 bins");
// Out projection 10 bins
  tCanvas->cd(3);	PrePad();
  hCF[0][2]->Draw();
  hFT[0][2]->Draw("SAMEHISTL");
  tLatex.DrawLatex(-0.18, 1.515, "Projections     10 bins");
// Side projection 1 bin
  tCanvas->cd(4);	PrePad();
  hCF[1][0]->Draw();
  hFT[1][0]->Draw("SAMEHISTL");
// Side projection 4 bins
  tCanvas->cd(5);	PrePad();
  hCF[1][1]->Draw();
  hFT[1][1]->Draw("SAMEHISTL");
// Side projection 10 bins
  tCanvas->cd(6);	PrePad();
  hCF[1][2]->Draw();
  hFT[1][2]->Draw("SAMEHISTL");
// Long projection 1 bin
  tCanvas->cd(7);	PrePad();
  hCF[2][0]->Draw();
  hFT[2][0]->Draw("SAMEHISTL");
// Long projection 4 bins
  tCanvas->cd(8);	PrePad();
  hCF[2][1]->Draw();
  hFT[2][1]->Draw("SAMEHISTL");
// Long projection 10 bins
  tCanvas->cd(9);	PrePad();
  hCF[2][2]->Draw();
  hFT[2][2]->Draw("SAMEHISTL");

  aFileHBTFit.ReplaceAll(".root",".eps");
  aFileHBTFit.ReplaceAll("hbtfit","fig_cfproj");
  tCanvas->SaveAs(aFileHBTFit);
}

