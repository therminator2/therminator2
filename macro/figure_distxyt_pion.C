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
#include <TMath.h>
#include <TStyle.h>
#include "events2chain.C"
#include "model2legend.C"
#include "hist2xml.C"

#define  _FIGURE_NAME_ "fig_distxyt_pion"
#define  _N_HISTOGRAMS_ 8

void figure_distxyt_pion(TString aEventDir = "./events/", Int_t aEventFiles = 1)
{
// ##########################################################################
// # READ ROOT FILES
// ##########################################################################
  static ParticleCoor Particle;
  Int_t   Events;
  TChain* Chain = events2chain(aEventDir, aEventFiles, &Particle, &Events);

// ##########################################################################
// HISTOGRAMS
// ##########################################################################
  Int_t   XBins  = 60;
  Float_t XMin   =  0.0;
  Float_t XMax   = 60.0;
  Float_t dX     = (XMax - XMin) / XBins;
  TH1D*   H1D[_N_HISTOGRAMS_];
// Create histograms
  H1D[0] = new TH1D("H0", "#pi^{+} (total)", XBins, XMin, XMax);
  H1D[0]->GetXaxis()->SetTitle("r_{T} [fm]");
  H1D[0]->GetYaxis()->SetTitle("dN / (r_{T} dr_{T} dY)");
  H1D[0]->Sumw2();
  H1D[1] = (TH1D*) H1D[0]->Clone("H1");  H1D[1]->SetTitle("#pi^{+} (from #rho^{0} and #rho^{+})");
  H1D[2] = (TH1D*) H1D[0]->Clone("H2");  H1D[2]->SetTitle("#pi^{+} (from #omega)");
  H1D[3] = (TH1D*) H1D[0]->Clone("H3");  H1D[3]->SetTitle("#pi^{+} (primordial)");
  H1D[4] = new TH1D("H4", "time #pi^{+} (total)", XBins, XMin, XMax);
  H1D[4]->GetXaxis()->SetTitle("t [fm/c]");
  H1D[4]->GetYaxis()->SetTitle("dN / (dt dY)");
  H1D[4]->Sumw2();
  H1D[5] = (TH1D*) H1D[4]->Clone("H5");  H1D[5]->SetTitle("time #pi^{+} (from #rho^{0} and #rho^{+})");
  H1D[6] = (TH1D*) H1D[4]->Clone("H6");  H1D[6]->SetTitle("time #pi^{+} (from #omega)");
  H1D[7] = (TH1D*) H1D[4]->Clone("H7");  H1D[7]->SetTitle("time #pi^{+} (primordial)");

// Fill histograms
  Float_t rT;
  for(Int_t i=0; i<Chain->GetEntries(); i++) {
    Chain->GetEntry(i);
    if(Particle.e == Particle.pz)
      continue;
    if( (TMath::Abs(0.5 * log((Particle.e+Particle.pz) / (Particle.e-Particle.pz))) < 1.0) && (Particle.pid == 211) ) { // |y| < 1.0 and pi+
      rT   = TMath::Sqrt(Particle.x*Particle.x + Particle.y*Particle.y);
      H1D[0]->Fill(rT, 1.0/rT);
      H1D[4]->Fill(Particle.t, 1.0);
      if((Particle.fatherpid == 113) || (Particle.fatherpid == 213)) {
        H1D[1]->Fill(rT, 1.0/rT);
        H1D[5]->Fill(Particle.t, 1.0);
      }
      if(Particle.fatherpid == 223) {
        H1D[2]->Fill(rT, 1.0/rT);
        H1D[6]->Fill(Particle.t, 1.0);
      }
      if(Particle.fatherpid == 211) {
        H1D[3]->Fill(rT, 1.0/rT);
        H1D[7]->Fill(Particle.t, 1.0);
      }
    }
  }
// Rescale histograms
  for(Int_t i=0; i<_N_HISTOGRAMS_; i++)
    H1D[i]->Scale(1.0 / (Events *2*1.0 * dX ));

// ##########################################################################
// # SAVE HISTOGRAMS TO XML FILE
// ##########################################################################
  hist2xml(aEventDir + _FIGURE_NAME_ + ".xml", H1D, _N_HISTOGRAMS_, Events, Chain->GetEntries());

// ##########################################################################
// # PLOT HISTOGRAMS
// ##########################################################################
  TCanvas* Canvas	= new TCanvas("canvas", H1D[0]->GetYaxis()->GetTitle(), 800, 600);
  TLegend* LegendPart	= new TLegend(0.14, 0.28, 0.44, 0.44, "",	"brNDC");
  TLegend* LegendModel	= new TLegend(0.20, 0.60, 0.44, 0.99, "model",	"brNDC");
  gStyle->SetOptStat("");
  gPad->SetFillColor(0);
  gPad->SetFillStyle(4000);
  gPad->SetFillColor(0);
// Markers
  H1D[0]->SetMarkerColor(2);	H1D[0]->SetMarkerStyle(20);
  H1D[1]->SetMarkerColor(4);	H1D[1]->SetMarkerStyle(21);
  H1D[2]->SetMarkerColor(6);	H1D[2]->SetMarkerStyle(23);
  H1D[3]->SetMarkerColor(1);	H1D[3]->SetMarkerStyle(22);
  H1D[4]->SetMarkerColor(2);	H1D[4]->SetMarkerStyle(20);
  H1D[5]->SetMarkerColor(4);	H1D[5]->SetMarkerStyle(21);
  H1D[6]->SetMarkerColor(6);	H1D[6]->SetMarkerStyle(23);
  H1D[7]->SetMarkerColor(1);	H1D[7]->SetMarkerStyle(22);
  LegendPart->SetFillStyle(0);
// Legends
  for(Int_t i=0; i<4; i++)
    LegendPart->AddEntry(H1D[i], H1D[i]->GetTitle());
  model2legend(aEventDir, aEventFiles, LegendModel);
// Plot
  H1D[0]->SetTitle("#pi^{+} r_{T} anatomy");
  H1D[4]->SetTitle("#pi^{+} time anatomy");
  Canvas->Divide(2,2,0.001,0.001);

  Canvas->cd(1);
  H1D[0]->Draw();
  H1D[1]->Draw("SAME");
  H1D[2]->Draw("SAME");
  H1D[3]->Draw("SAME");

  Canvas->cd(2);
  gPad->SetLogy();
  H1D[0]->Draw();
  H1D[1]->Draw("SAME");
  H1D[2]->Draw("SAME");
  H1D[3]->Draw("SAME");

  Canvas->cd(3);
  H1D[4]->Draw();
  H1D[5]->Draw("SAME");
  H1D[6]->Draw("SAME");
  H1D[7]->Draw("SAME");

  Canvas->cd(4);
  gPad->SetLogy();
  H1D[4]->Draw();
  H1D[5]->Draw("SAME");
  H1D[6]->Draw("SAME");
  H1D[7]->Draw("SAME");

  Canvas->cd();
  LegendPart ->Draw();
  LegendModel->Draw();

// Save to files
  Canvas->SaveAs(aEventDir + _FIGURE_NAME_ + ".eps");
//  Canvas->SaveAs(aEventDir + _FIGURE_NAME_ + ".png");
//  Canvas->SaveAs(aEventDir + _FIGURE_NAME_ + ".C");
}