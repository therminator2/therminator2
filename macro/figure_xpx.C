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

#include <TH2D.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TStyle.h>
#include "events2chain.C"
#include "model2legend.C"

#define  _FIGURE_NAME_ "fig_xpx"

void figure_xpx(TString aEventDir = "./events/", Int_t aEventFiles = 1)
{
// ##########################################################################
// # READ ROOT FILES
// ##########################################################################
  static ParticleCoor Particle;
  Int_t   Events;
  TChain* Chain = events2chain(aEventDir, aEventFiles, &Particle, &Events);

// ##########################################################################
// HISTOGRAM
// ##########################################################################
  Int_t   XBins  = 100;
  Float_t XMin   = -20.0;
  Float_t XMax   =  20.0;
  Int_t   YBins  = 100;
  Float_t YMin   = -1.5;
  Float_t YMax   =  1.5;
  TH2D*   H2D;
// Create histograms
  H2D = new TH2D("H0", "x-p_{x} correlation (#pi^{+}, midrapidity)", XBins, XMin, XMax, YBins, YMin, YMax);
  H2D->GetXaxis()->SetTitle("r_{x} [fm]");
  H2D->GetYaxis()->SetTitle("p_{x} [GeV/c]");
// Fill histograms
  for(Int_t i=0; i<Chain->GetEntries(); i++) {
    Chain->GetEntry(i);
    if(Particle.e == Particle.pz)
      continue;
    if( (Particle.pid == 211) && (0.5 * TMath::Log((Particle.e+Particle.pz) / (Particle.e-Particle.pz)) < 1.0))
      H2D->Fill(Particle.x, Particle.px);
  }

// ##########################################################################
// # PLOT HISTOGRAM
// ##########################################################################
  Int_t tMinBin, tMaxBin;
  gStyle->SetOptStat("");
// Canvas
  TCanvas* Canvas	= new TCanvas("canvas", H2D->GetYaxis()->GetTitle(), 800, 600);
  gPad->SetFillColor(0);
  gPad->SetFillStyle(4000);
  gPad->SetMargin(0.14, 0.02, 0.14, 0.08);
// Model legend
  TLegend* LegendModel	= new TLegend(0.70, 0.40, 0.98, 0.92, "model",	"brNDC");
  LegendModel->SetFillColor(0);
  LegendModel->SetFillStyle(4000);
  LegendModel->SetTextSize(0.035);
  model2legend(aEventDir, aEventFiles, LegendModel);
// Histograms  
  H2D->GetXaxis()->SetTitleSize(0.06);
  H2D->GetXaxis()->CenterTitle(kTRUE);
  H2D->GetXaxis()->SetLabelSize(0.05); 
  H2D->GetYaxis()->SetTitleSize(0.06);
  H2D->GetYaxis()->CenterTitle(kTRUE); 
  H2D->GetYaxis()->SetLabelSize(0.05);
  H2D->SetMarkerColor(2);	H2D->SetMarkerStyle(21);
// Plot
  H2D->Draw("COL");
  LegendModel->Draw();
// Save to files
  Canvas->SaveAs(aEventDir + _FIGURE_NAME_ + ".eps");
//  Canvas->SaveAs(aEventDir + _FIGURE_NAME_ + ".png");
//  Canvas->SaveAs(aEventDir + _FIGURE_NAME_ + ".C");
}