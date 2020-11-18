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
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <fstream>

#define  _FIGURE_NAME_ "fig_hbtradii"

void preparehist(TH1D *aH1)
{
  Int_t tMinBin, tMaxBin;
  tMinBin = aH1->GetMinimumBin();	aH1->SetMinimum((aH1->GetBinContent(tMinBin) - aH1->GetBinError(tMinBin))*0.95);
  tMaxBin = aH1->GetMaximumBin();	aH1->SetMaximum((aH1->GetBinContent(tMaxBin) + aH1->GetBinError(tMaxBin))*1.05);
}

void preparegraph(TGraphErrors *aGE)
{
  aGE->SetTitle("");
  aGE->GetXaxis()->SetTitleSize(0.1);
  aGE->GetXaxis()->SetLabelSize(0.1);
  aGE->GetXaxis()->CenterTitle(kTRUE);
  aGE->GetXaxis()->SetNdivisions(4,5,0);
  aGE->GetYaxis()->SetTitleSize(0.1);
  aGE->GetYaxis()->SetTitleOffset(0.8);
  aGE->GetYaxis()->SetLabelSize(0.1);
  aGE->GetYaxis()->CenterTitle(kTRUE);
  aGE->GetYaxis()->SetNdivisions(3,5,0);
  aGE->SetMarkerStyle(20);
  aGE->SetMarkerColor(1);
}

void preparepad()
{
  gPad->SetFillStyle(4000);
  gPad->SetFillColor(0);
  gPad->SetRightMargin(0.02);
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.22);
  gPad->SetLeftMargin(0.17);
}

void figure_hbtradii(TString aEventDir)
{
  Double_t xy[4] = { 0.2, 0.3, 0.4, 0.525 };

  Double_t ykp[4];
  Double_t ykpe[4];
  Double_t yko[4];
  Double_t ykoe[4];
  Double_t yks[4];
  Double_t ykse[4];
  Double_t ykl[4];
  Double_t ykle[4];

  ifstream *istr = new ifstream(aEventDir+"hbtradii.txt");
  char name[200];
  char dum[200];
  
  char buf[2000];

  for (int ik=0; ik<4; ik++) {
    istr->getline(name, 2000);
//    cout << "Fit from: "<<name << endl;
    istr->ignore(20000,'\n');
    (*istr) >> name >> ykp[ik] >> dum >> ykpe[ik];
//    cout << name <<" = "<< ykp[ik] <<" +- "<< ykpe[ik] << endl;
    istr->ignore(20000,'\n');
    (*istr) >> name >> yko[ik] >> dum >> ykoe[ik];
//    cout << name <<" = "<< yko[ik] <<" +- "<< ykoe[ik] << endl;
    istr->ignore(20000,'\n');
    (*istr) >> name >> yks[ik] >> dum >> ykse[ik];
//    cout << name <<" = "<< yks[ik] <<" +- "<< ykse[ik] << endl;
    istr->ignore(20000,'\n');
    (*istr) >> name >> ykl[ik] >> dum >> ykle[ik];
//    cout << name <<" = "<< ykl[ik] <<" +- "<< ykle[ik] << endl;
    istr->ignore(20000,'\n');
  }
  
  TGraphErrors *grp = new TGraphErrors(4, xy, ykp, NULL, ykpe);
  TGraphErrors *gro = new TGraphErrors(4, xy, yko, NULL, ykoe);
  TGraphErrors *grs = new TGraphErrors(4, xy, yks, NULL, ykse);
  TGraphErrors *grl = new TGraphErrors(4, xy, ykl, NULL, ykle);

  TCanvas *canproj1 = new TCanvas ("canfitrad","FittedRadii", 800, 400);
  preparepad();
  canproj1->Divide(2,2,0.001, 0.001);
  
  canproj1->cd(1); preparepad();
  preparegraph(grp);
  grp->SetTitle(";k_{T} (GeV/c);\\lambda");
  grp->Draw("ACP");

  canproj1->cd(2); preparepad();
  preparegraph(gro);
  gro->SetTitle(";k_{T} (GeV/c);R_{out} (fm)");
  gro->Draw("ACP");

  canproj1->cd(3); preparepad();
  preparegraph(grl);
  grl->SetTitle(";k_{T} (GeV/c);R_{long} (fm)");
  grl->Draw("ACP");

  canproj1->cd(4); preparepad();
  preparegraph(grs);
  grs->SetTitle(";k_{T} (GeV/c);R_{side} (fm)");
  grs->Draw("ACP");

  canproj1->SaveAs(aEventDir + _FIGURE_NAME_ + ".eps");
}

