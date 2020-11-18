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
//#include "setStyle.C"
#include "/u/harabasz/HistogramTemplate/drawStyle.C"
#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TStyle.h>
#include "events2chain.C"
#include "model2legend.C"
#include "hist2xml.C"

#include <vector>

using namespace std;

#define  _FIGURE_NAME_ "pionM_midY_T60_muB760_muS0_Yleft01.gif"
#define  _N_HISTOGRAMS_ 13

struct MyParticle {
  Float_t e;
  Float_t px;
  Float_t py;
  Float_t pz;
  Float_t eventid;
  Float_t fathereid;
  Float_t fatherpid;
};

void figure_distminv_pimprot(TString aEventDir = "./", Int_t aEventFiles = 100)
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
//  Int_t   XBins  = 40;
  Int_t   XBins  = 200;
  Float_t XMin   = 1.0;
  Float_t XMax   = 2;
//  Float_t XMax   = 1.0;
  Float_t dX     = (XMax - XMin) / XBins;
  TH1D*   H1D[_N_HISTOGRAMS_];
// Create histograms
  H1D[0] = new TH1D("H0", "#pi^{-} (total)", XBins, XMin, XMax);
//  H1D[0]->GetYaxis()->SetTitle("1/mt * dN/(2 #pi p_{T} dp_{T} dy)");
  H1D[0]->Sumw2();
  H1D[1]  = (TH1D*) H1D[0]->Clone("H1");  H1D[1]->SetTitle("primordial");
  H1D[2]  = (TH1D*) H1D[0]->Clone("H2");  H1D[2]->SetTitle("#Delta(1232), #Delta(1600)");
  H1D[3]  = (TH1D*) H1D[0]->Clone("H3");  H1D[3]->SetTitle("N(1440), N(1520)");
  H1D[4]  = (TH1D*) H1D[0]->Clone("H4");  H1D[4]->SetTitle("#rho^{0} and #rho^{+}");
  H1D[5]  = (TH1D*) H1D[0]->Clone("H5");  H1D[5]->SetTitle("#omega");

  H1D[6]  = (TH1D*) H1D[0]->Clone("H6");  H1D[6]->SetTitle("#Lambda");
  H1D[7]  = (TH1D*) H1D[0]->Clone("H7");  H1D[7]->SetTitle("K^{0}");
  H1D[8]  = (TH1D*) H1D[0]->Clone("H8");  H1D[8]->SetTitle("#Lambda(1520)");
  H1D[9]  = (TH1D*) H1D[0]->Clone("H9");  H1D[9]->SetTitle("#Sigma^{-}");
  H1D[10] = (TH1D*) H1D[0]->Clone("H10"); H1D[10]->SetTitle("#Sigma(1385)");
  H1D[11] = (TH1D*) H1D[0]->Clone("H11"); H1D[11]->SetTitle("#Xi");
  H1D[12] = (TH1D*) H1D[0]->Clone("H12"); H1D[12]->SetTitle("#Xi(1530)");
  Int_t x = kBlue;
  Int_t color[]= {kBlack, kRed, kPink-4, kMagenta-1, kViolet-3, kBlue+2, kBlue+1, kCyan+3, kCyan-2, kTeal+3, kGreen-3,x+4,x+2,x+0,x-4,x-7,x-9,x-8,kGray,x-10,x+3,x+1,x-1,x-3,x-6,x-7,x+4,x+2,x+0,x-4,x-7};


  for(Int_t m=0;m<_N_HISTOGRAMS_;m++){
      if(m!=9){
          setTH1(H1D[m],"m_{inv} [GeV/c^{2}]","dN / dm_{inv}",m+1);
          H1D[m]->SetLineColor(m+1);
          H1D[m]->SetMarkerColor(m+1);
          H1D[m]->SetMarkerStyle(1);
      }
      else setTH1(H1D[m],"m_{T} - m_{0} [GeV/c]","1/(m_{T}^{9} * d^{2}N/dm_{T}dy)",kGreen+2);
  }

  std::vector<MyParticle> pions;
  std::vector<MyParticle> protons;

  // Fill histograms
  Float_t Rap, Pt, Ptot, Mt, M0, Eta;
  Int_t   fatherpid;
  for(Int_t i=0; i<Chain->GetEntries(); i++) {
    Chain->GetEntry(i);
    if(Particle.pid != -211 && Particle.pid != 2212)  continue;

    if(Particle.e == Particle.pz)
      continue;
    fatherpid = Particle.fatherpid;
    if(fatherpid == Particle.pid) {
      continue;
    }
    MyParticle part;
    part.e = Particle.e;
    part.px = Particle.px;
    part.py = Particle.py;
    part.pz = Particle.pz;
    part.eventid = Particle.eventid;
    part.fathereid = Particle.fathereid;
    part.fatherpid = Particle.fatherpid;
    if (Particle.pid == -211) {
      pions.push_back(part);
    } else {
      protons.push_back(part);
    }
  }

  cout << "vector sizes " << protons.size() << " " << pions.size() << endl;

  for (int ii = 0; ii <= protons.size(); ++ii) {
    for (int jj = 0; jj <= pions.size(); ++jj) {
      MyParticle pion = pions[jj];
      MyParticle proton = protons[ii];
      if (pion.eventid != proton.eventid) continue;
      if (pion.fathereid != proton.fathereid) continue;

      fatherpid = pion.fatherpid;
      Float_t deno = 1.;
      Float_t e1 = pion.e;
      Float_t px1 = pion.px;
      Float_t py1 = pion.py;
      Float_t pz1 = pion.pz;
      Float_t e2 = proton.e;
      Float_t px2 = proton.px;
      Float_t py2 = proton.py;
      Float_t pz2 = proton.pz;

      Float_t mass = TMath::Sqrt( TMath::Power( e1 + e2, 2 ) - TMath::Power( px1 + px2, 2 ) - TMath::Power( py1 + py2, 2 ) - TMath::Power( pz1 + pz2, 2 ) );

      H1D[0]->Fill(mass, deno);
      if(fatherpid == -211) H1D[1]->Fill(mass, deno);
      if(fatherpid == 1114 || fatherpid == 2114 || fatherpid ==2214 || fatherpid == 2224||
         fatherpid == 31114 || fatherpid == 32114 ||fatherpid == 32214 || fatherpid == 32224)
      {
  	H1D[2]->Fill(mass, deno);
      }
      if( fatherpid == 12112 || fatherpid == 12212 || fatherpid == 1214 || fatherpid == 2124)
  	H1D[3]->Fill(mass, deno);
      if((fatherpid == 113) || (fatherpid == 213)) H1D[4]->Fill(mass, deno);
      if(fatherpid == 223) H1D[5]->Fill(mass, deno);

      if( fatherpid == 311)  H1D[7]->Fill(mass, deno);
      if( fatherpid == 3122)  H1D[6]->Fill(mass, deno);
      if( fatherpid == 3124)  H1D[8]->Fill(mass, deno);
      if( fatherpid == 3112)  H1D[9]->Fill(mass, deno);
      if( fatherpid == 3114 || fatherpid==3214)  H1D[10]->Fill(mass, deno);
      if( fatherpid == 3321 )  H1D[11]->Fill(mass, deno);
      if( fatherpid == 3324 ) {
          cout << H1D[12] << endl;
          H1D[12]->Fill(mass, deno);
      }

    }
  }

  cout<<"events : "<<Events<<endl;


  for(Int_t i=0; i<_N_HISTOGRAMS_; i++)
      H1D[i]->Scale(1.0 / (Events *0.1* dX ));
  // 0.1 : szerokosc binu w Y
  // dX : szerokosc binu w mt

// ##########################################################################
// # SAVE HISTOGRAMS TO XML FILE
// ##########################################################################
  //hist2xml(aEventDir + _FIGURE_NAME_ + ".xml", H1D, _N_HISTOGRAMS_, Events, Chain->GetEntries());

// ##########################################################################
// # PLOT HISTOGRAMS
// ##########################################################################

  TFile *in_exp = new TFile("/u/sudol/apr12/loopDST/gen8/heidi/pim_hist_4Therminator.root","read");

  TH1F *exp , *exp1;

  exp = (TH1F*)in_exp->Get("Mt_Corr_pim_y6_mult3"); // this is for 10% centrality
//  exp1 = (TH1F*)in_exp->Get("Mt_Corr_pim_y6_mult6");
  Int_t binLow =  exp->FindBin(0);
  Int_t binMax  = exp->FindBin(800);
  Int_t binLow1 = H1D[0]->FindBin(0);
  Int_t binMax1 = H1D[0]->FindBin(800);
//  exp->Add(exp1);

  Float_t int_exp = exp->Integral(binLow,binMax,"width");
  Float_t int_the = H1D[0]->Integral(binLow1,binMax1,"width");

/*
  H1D[0]->Scale(int_exp/int_the);
  H1D[1]->Scale(int_exp/int_the);
  H1D[2]->Scale(int_exp/int_the);
  H1D[3]->Scale(int_exp/int_the);
  H1D[4]->Scale(int_exp/int_the);
  H1D[5]->Scale(int_exp/int_the);
*/
  TH1F *kk;

  kk = (TH1F*)H1D[3]->Clone("kk");
  kk->Add(H1D[4]);


 // exp->Scale(int_the/int_exp);

  cout<<"scaling : "<<int_the<< " / " << int_exp<<endl;
  cout<<"scaling : "<<int_the/int_exp<<endl;


  TFile *ur = new TFile("/u/sudol/apr12/loopDST/gen8/heidi/theory/urqmd/final_state_piminus_spectra_cent.root","read");
  TH1F *h_other = (TH1F*) ur->Get("h_sum_mt_y_8");
  h_other->SetLineColor(8);
  h_other->SetLineWidth(3);

  Int_t tMinBin, tMaxBin;
  TCanvas* Canvas = new TCanvas("Canvas", H1D[0]->GetYaxis()->GetTitle(),10,10,600,600);
  //setOPT_canvas(Canvas);
  Canvas->SetLogy();

  // Histogram legend
  TLegend* LegendPart	= new TLegend(0.42,0.70,0.71,0.99, "",	"brNDC");
  LegendPart->SetFillColor(0);
  LegendPart->SetFillStyle(4000);
  LegendPart->SetTextSize(0.03);
  LegendPart->SetName("LegendPart");
  // for(Int_t i=0; i<_N_HISTOGRAMS_; i++)
  for(Int_t i=0; i<8; i++)
      if (H1D[i]->Integral() > 0) {
          LegendPart->AddEntry(H1D[i], H1D[i]->GetTitle());
      }
  // Model legend
  TLegend* LegendModel	= new TLegend(0.72,0.70,0.99,0.99, "model",	"brNDC");
  LegendModel->SetFillColor(0);
  LegendModel->SetFillStyle(4000);
  LegendModel->SetTextSize(0.03);
  LegendModel->SetName("LegendModel");
  model2legend(aEventDir, aEventFiles, LegendModel);
  // Histograms
//  cout<<"kuku"<<endl;

 // H1D[0]->SetTitle("#pi^{-} anatomy");
  H1D[0]->SetTitle("");
  H1D[0]->Draw("hist");
  H1D[1]->Draw("SAMEhist");
  H1D[2]->Draw("SAMEhist");
  H1D[3]->Draw("SAMEhist");
  H1D[4]->Draw("SAMEhist");
  H1D[5]->Draw("SAMEhist");
  H1D[6]->Draw("SAMEhist");
  H1D[7]->Draw("SAMEhist");
/*
  H1D[8]->Draw("SAMEhist");
  H1D[9]->Draw("SAMEhist");
  H1D[10]->Draw("SAMEhist");
  H1D[11]->Draw("SAMEhist");
*/
//  H1D[12]->Draw("SAMEhist");
//  exp->Draw("same");
  //h_other->Draw("same");
  LegendPart ->Draw();
  LegendModel->Draw();
  //setOPT_pad();
// Save to files
// Canvas->SaveAs(_FIGURE_NAME_);
//  Canvas->SaveAs(TString("plots/") + _FIGURE_NAME_ + ".png");
//  Canvas->SaveAs(TString("plots/") + _FIGURE_NAME_ + ".C");



  TFile *out = new TFile("pimprot.root","recreate");
  out->cd();
  for(Int_t i=0; i<_N_HISTOGRAMS_; i++){
      H1D[i]->Write();
  }
  LegendPart ->Write();
  LegendModel->Write();
  exp->Write();
}
