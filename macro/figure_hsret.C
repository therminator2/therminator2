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

#include <TH3D.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TStyle.h>
#include "events2chain.C"
#include "model2legend.C"

#define  _FIGURE_NAME_ "fig_hsret"

void figure_hsret(TString aEventDir = "./events/", Int_t aEventFiles = 1)
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
  Int_t   XBins  = 150;
  Float_t XMin   =  0.0;
  Float_t XMax   = 10.0;
  Int_t   YBins  = 150;
  Float_t YMin   = -8.0;
  Float_t YMax   =  8.0;
  Int_t   ZBins  = 150;
  Float_t ZMin   =  0.0;
  Float_t ZMax   = 15.0;
  TH3D*   H3D[2];
// Create histograms
  H3D[0] = new TH3D("H0", "hypersurface in #rho-#eta-#tau", XBins, XMin, XMax, YBins, YMin, YMax, ZBins, ZMin, ZMax);
  H3D[0]->GetXaxis()->SetTitle("#rho [fm]");
  H3D[0]->GetYaxis()->SetTitle("#eta");
  H3D[0]->GetZaxis()->SetTitle("#tau [fm]");
  H3D[1] = (TH3D*) H3D[0]->Clone("H1");  H3D[1]->SetTitle("resonances");
// Fill histograms
  Float_t EtaS, Rho, Phi, Tau;
  for(Int_t i=0; i<Chain->GetEntries(); i++) {
    Chain->GetEntry(i);
    EtaS = 0.5*TMath::Log((Particle.t+Particle.z) / (Particle.t-Particle.z));
    Rho  = TMath::Sqrt(Particle.x*Particle.x + Particle.y*Particle.y);
    Tau  = TMath::Sqrt(Particle.t*Particle.t - Particle.z*Particle.z);
    if(Particle.pid == Particle.fatherpid)
      H3D[0]->Fill(Rho, EtaS, Tau);
    else
      H3D[1]->Fill(Rho, EtaS, Tau);
  }

// ##########################################################################
// # PLOT HISTOGRAM
// ##########################################################################
  TCanvas* Canvas	= new TCanvas("canvas", H3D[0]->GetTitle(), 800, 600);
  TLegend* LegendModel	= new TLegend(0.75, 0.60, 0.99, 0.99, "model",	"brNDC");
  gStyle->SetOptStat("");
  gPad->SetFillColor(0);
  gPad->SetFillStyle(4000);
  H3D[0]->SetMarkerColor(4);
  H3D[1]->SetMarkerColor(3);
// Plot
  H3D[0]->Draw("");
  //H3D[1]->Draw("SAME");
// Save to files
//  Canvas->SaveAs(aEventDir + _FIGURE_NAME_ + ".eps");
  Canvas->SaveAs(aEventDir + _FIGURE_NAME_ + ".png");
//  Canvas->SaveAs(aEventDir + _FIGURE_NAME_ + ".C");
}