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
#include "/u/harabasz/HistogramTemplate/drawStyle.C"
#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TStyle.h>
#include <TObjArray.h>
#include <TGraphErrors.h>
#include "events2chain.C"
#include "model2legend.C"
#include "hist2xml.C"

#include "hphysicsconstants.h"
#include "integrate_mt.h"

#define  _FIGURE_NAME_ "pip_mt2_fig_distmt_pionP_T55_mu808"
#define  _N_HISTOGRAMS_ 12
#define  _N_IDS_ 5
//

int indexOf(int size, int *arr, int elem) {
     for (int i = 0; i < size; ++i) {
         if (arr[i] == elem) return i;
     }
     return -1;
}

int fatherIndex(int pdg, int fatherpid) {
    if (pdg == 211) {
        if(fatherpid == 211)  return 1;
        if(fatherpid == 1114   || fatherpid == 2114  ||fatherpid ==2214   || fatherpid == 2224||
                fatherpid == 31114  || fatherpid == 32114 ||fatherpid == 32214 || fatherpid == 32224) return 2;
        if( fatherpid == 12112 || fatherpid == 12212 || fatherpid == 1214 || fatherpid == 2124)  return 3;
        if(fatherpid == 113    || fatherpid == 213)  return 4;
        if(fatherpid == 223)                         return 5;
        if(fatherpid == 311)                         return 6;
        if(fatherpid == 323)                         return 7;
        if(fatherpid == 3122)                        return 8;
        if(fatherpid == 3214)                        return 9;
        if(fatherpid == 3222)                        return 10;
        if(fatherpid == 3224)                        return 11;
    }
    if (pdg == -211) {
        if(fatherpid == -211) return 1;
        if(fatherpid == 1114 || fatherpid == 2114 || fatherpid ==2214 || fatherpid == 2224||
                fatherpid == 31114 || fatherpid == 32114 ||fatherpid == 32214 || fatherpid == 32224) return 2;
        if( fatherpid == 12112 || fatherpid == 12212 || fatherpid == 1214 || fatherpid == 2124) return 3;
        if((fatherpid == 113) || (fatherpid == 213)) return 4;
        if(fatherpid == 223)                         return 5;

        if( fatherpid == 311)                        return 7;
        if( fatherpid == 3122)                       return 6;
        if( fatherpid == 3124)                       return 8;
        if( fatherpid == 3112)                       return 9;
        if( fatherpid == 3114 || fatherpid==3214)    return 10;
        if( fatherpid == 3321 )                      return 11;
        if( fatherpid == 3324 )                      return 12;
    }
    if (pdg == 2212) {
        if(fatherpid == 2212)
            return 1;
        if((fatherpid == 3212) || (fatherpid == 3222))
            return 2;
        if((fatherpid == 2114) || (fatherpid == 2214) || (fatherpid ==  2224))
            return 3;
        if(fatherpid == 3122)
            return 4;
        if( fatherpid == 12112 || fatherpid == 12212 || fatherpid == 1214 || fatherpid == 2124)
            return 5;
    }
    if (pdg == 321) {
        if(fatherpid == 321)  return 1;
        if(fatherpid == 1114   || fatherpid == 2114  ||fatherpid ==2214   || fatherpid == 2224||
                fatherpid == 31114  || fatherpid == 32114 ||fatherpid == 32214 || fatherpid == 32224) return 2;
        if( fatherpid == 12112 || fatherpid == 12212 || fatherpid == 1214 || fatherpid == 2124)  return 3;
        if(fatherpid == 113    || fatherpid == 213)  return 4;
        if(fatherpid == 223)                         return 5;
    }
    if (pdg == 3122) {
        if(fatherpid == 3122) return 1;
        if(fatherpid == 3212) return 2;
    }
    
    return -1;
}

TH1F *getExp(int pdg) {
    TH1F *exp = NULL;
    if (pdg == 211) {
        TFile *in_exp = new TFile("/u/sudol/apr12/loopDST/gen8/heidi/pip_hist_4Therminator.root","read");
        exp = (TH1F*)in_exp->Get("Mt_Corr_pim_y6_mult3");
    } else if (pdg == -211) {
        TFile *in_exp = new TFile("/u/sudol/apr12/loopDST/gen8/heidi/pim_hist_4Therminator.root","read");
        exp = (TH1F*)in_exp->Get("Mt_Corr_pim_y6_mult3"); // this is for 10% centrality
    } else if (pdg == 2212) {
        TFile *in_exp = new TFile("/lustre/nyx/hades/user/schuldes/apr12/Analyse/Happyparticles/p/Happy_p_jan16_Rk400_MM3_SL_work.root","read");
        exp = (TH1F*)in_exp->Get("Mt_Corr_p_y6_mult3");
        exp->Scale(1./(0.1*25.*1.07235890000000000e+07));
    } else if (pdg == 321) {
        TFile *in_exp = new TFile("/u/harabasz/apr12/theory/therminator2_2.0.3_coul/exp_data/kp_mt2_corr_010centr.root","read");
        exp = (TH1F*)in_exp->Get("Mt_Corr_Kp_y6_sys3");
        exp->Scale(1./25/0.1/5.3e8);
    } else if (pdg == 3122) {
        TFile *in_exp = new TFile("/u/harabasz/apr12/theory/therminator2_2.0.3_coul/exp_data/Lambda_dNdmt.root","read");
        exp = (TH1F*)in_exp->Get("hSignalCur8");
        exp->Scale(1./1e7);
    }
    return exp;
}
vector<TObject *> getExpY(int pdg) {
    vector<TObject *> result;
    if (pdg == -211) {
        TFile *in          = new TFile("/u/sudol/apr12/loopDST/gen8/heidi//test_fitUpFrom100_working.root","read");
        TH1F *data      = (TH1F*)in->Get(Form("hY_2S_mult_3"));
        TH1F *data_ref  = (TH1F*)in->Get(Form("hY_2S_mult_refl_3"));
        result.push_back(data);
        result.push_back(data_ref);
    }
    if (pdg == 211) {
        TFile *in          = new TFile("/u/sudol/pip_y_0_10.root");
        TFile *in_ref      = new TFile("/u/sudol/pip_yRef_0_10.root");
        TH1F *data      = (TH1F*)in->Get(Form("hY_2S_mult_3"));
        TH1F *data_ref  = (TH1F*)in_ref->Get(Form("hY_2S_mult_refl_3"));
        result.push_back(data);
        result.push_back(data_ref);
    }
    if (pdg == 2212) {
        TFile *inD = new TFile("/u/sudol/protons_0_10_Melanie.root","read");
        TH1F *hData = (TH1F*) inD->Get("hY_fSR_0");
        result.push_back(hData);
    }
    if (pdg == 321) {
        TFile *inD = new TFile("/u/harabasz/apr12/theory/therminator2_2.0.3_coul/macro/kp.root","read");
        TH1F *hData = (TH1F*) inD->Get("hY_scale_0");
        result.push_back(hData);
    }
    if (pdg == 3122) {
        TFile *in          = new TFile("/u/sudol/Gosia/Lambda_0_10.root");
        TFile *in_ref      = new TFile("/u/sudol/Gosia/Lambda_0_10_ref.root");
        TH1F *data      = (TH1F*)in->Get("hdNdyCent1");
        TH1F *data_ref  = (TH1F*)in_ref->Get("hdNdyMirrorCent1");
        result.push_back(data);
        result.push_back(data_ref);
    }
    return result;
}

int maxHist(int pdg) {
    if (pdg == 211) return 5;
    if (pdg == -211) return 7;
    if (pdg == 2212) return 5;
    if (pdg == 321) return 1;
    if (pdg == 3122) return 2;
    return 0;
}

TString titles(int pdg, int index) {
    if (pdg == 211) {
        if (index == 0) return "#pi^{+} (total)";
        if (index == 1) return "#pi^{+} primordial";
        if (index == 2) return "#Delta(1232), #Delta(1600)";
        if (index == 3) return "N(1440), N(1520)";
        if (index == 4) return "#rho^{0} and #rho^{+}";
        if (index == 5) return "#omega";
        if (index == 6) return "K^{0}";
        if (index == 7) return "K^{*}(892)";
        if (index == 8) return "#Lambda(1405)";
        if (index == 9) return "#Sigma";
        if (index == 10) return "#Sigma^{+}";
        if (index == 11) return "#Sigma(1385)";
    }
    if (pdg == -211) {
        if (index == 0) return "#pi^{+} (total)";
        if (index == 1) return "#pi^{+} primordial";
        if (index == 2) return "#Delta(1232), #Delta(1600)";
        if (index == 3) return "N(1440), N(1520)";
        if (index == 4) return "#rho^{0} and #rho^{+}";
        if (index == 5) return "#omega";
        if (index == 6) return "#Lambda";
        if (index == 7) return "K^{0}";
        if (index == 8) return "#Lambda(1520)";
        if (index == 9) return "#Sigma^{-}";
        if (index == 10) return "#Sigma(1385)";
        if (index == 11) return "#Xi";
        if (index == 12) return "#Xi(1530)";
    }
    if (pdg == 2212) {
        if (index == 0) return "p (total)";
        if (index == 1) return "p primordial";
        if (index == 2) return "#Sigma";
        if (index == 3) return "#Delta";
        if (index == 4) return "#Lambda^{0}";
        if (index == 5) return "N*";
    }
    if (321) {
        if (index == 0) return "K^{+} (total)";
        if (index == 1) return "K^{+} primordial";
    }
    if (3122) {
        if (index == 0) return "#Lambda (total)";
        if (index == 1) return "#Lambda primordial";
        if (index == 2) return "#Sigma^{0}";
    }
    return "";
}

void figure_distmt(TString aEventDir, Int_t aEventFiles = 100)
{

    TString partNames[] = { "#pi^{-}",  "#pi^{+}",  "p", "K^{+}", "#Lambda" };  
    Int_t ids[_N_IDS_] = { -211, 211, 2212, 321, 3122 };

    TString outDir;
    TObjArray *tokens = aEventDir.Tokenize("/");
    int ntokens = tokens->GetEntries();
    TObjString *ostr = (TObjString*)tokens->At(ntokens-1);
    outDir = ostr->String();
    cout << outDir.Data() << endl;
    gSystem->MakeDirectory(outDir);

    // ##########################################################################
    // # READ ROOT FILES
    // ##########################################################################
    static ParticleCoor Particle;
    Int_t   Events;
    TChain* Chain = events2chain(aEventDir, aEventFiles, &Particle, &Events);

    // ##########################################################################
    // HISTOGRAMS
    // ##########################################################################
    Int_t   XBins  = 80;
    Float_t XMin   = 0.0;
    Float_t XMax   = 800;

    //  Int_t   XBins  = 100;
    //  Float_t XMin   = 0.0;
    //  Float_t XMax   = 1000;
    //  Float_t XMax   = 1.0;
    Float_t dX     = (XMax - XMin) / XBins;


    TH1D*   H1D[_N_IDS_][_N_HISTOGRAMS_];
    TH1D*   Y1D[_N_IDS_][_N_HISTOGRAMS_];
    // Create histograms
    for (int id = 0; id < _N_IDS_; ++id) {
        TString txtid = TString::Itoa(abs(ids[id]),10);
        if (ids[id] < 0) txtid.Prepend("m");
        H1D[id][0] = new TH1D(Form("H0_pid%s",txtid.Data()), titles(ids[id],0), XBins, XMin, XMax);
        H1D[id][0]->Sumw2();
        Y1D[id][0] = new TH1D(Form("Y0_pid%s",txtid.Data()), titles(ids[id],0), 100, -4, 4);
        Y1D[id][0]->Sumw2();
        for(Int_t m=0;m<_N_HISTOGRAMS_;m++){
            if (m > 0) H1D[id][m] = (TH1D*)  H1D[id][0]->Clone(Form("H%i_pid%s",m,txtid.Data()));
            H1D[id][m]->SetTitle(titles(ids[id],m));
            setTH1(H1D[id][m],"m_{T} - m_{0} [MeV/c]","1/m_{t}^{2} d^{2}N/dm_{t}dy#cbar_{y=0} [(MeV/c)^{-3}]",m+1);
            H1D[id][m]->GetXaxis()->SetNdivisions(507);
            H1D[id][m]->SetLineColor(m+1);
            H1D[id][m]->SetMarkerColor(m+1);
            H1D[id][m]->SetMarkerStyle(1);

            if (m > 0) Y1D[id][m] = (TH1D*)  Y1D[id][0]->Clone(Form("Y%i_pid",m,txtid.Data()));
            Y1D[id][m]->SetTitle(titles(ids[id],m));
            setTH1(Y1D[id][m],"y","dN / dy",m+1);
            Y1D[id][m]->GetXaxis()->SetNdivisions(507);
            Y1D[id][m]->SetLineColor(m+1);
            Y1D[id][m]->SetMarkerColor(m+1);
            Y1D[id][m]->SetMarkerStyle(1);

        }
    }

    // Fill histograms
    Float_t Rap, Pt, Mt, M0;
    Int_t   fatherpid;
    for(Int_t i=0; i<Chain->GetEntries(); i++) {
        Chain->GetEntry(i);
        //    if(Particle.pid != -211)
        int index = indexOf(_N_IDS_, ids, Particle.pid);
        if(index == -1)
            continue;
        if(Particle.e == Particle.pz)
            continue;
        Rap = 0.5 * TMath::Log((Particle.e+Particle.pz) / (Particle.e-Particle.pz));
        fatherpid = Particle.fatherpid;
        Pt  = TMath::Sqrt(Particle.px*Particle.px + Particle.py*Particle.py)*1000.;
        Mt  = TMath::Sqrt(Particle.mass*Particle.mass + Particle.px*Particle.px + Particle.py*Particle.py)*1000.;
        M0  = Particle.mass*1000.;
        if (!TMath::Finite(Mt)) continue;

        Float_t deno = 1./(Mt*Mt);
       // Float_t deno = 1.;
        //    Float_t deno = 1./(2*TMath::Pi()*Mt);

        if(TMath::Abs(Rap) <= 0.05 ) {
            H1D[index][0]->Fill(Mt-M0, deno);
            if (fatherIndex(Particle.pid, fatherpid) > 0) {
                H1D[index][fatherIndex(Particle.pid, fatherpid)]->Fill(Mt-M0, deno);
            }
        }
        Y1D[index][0]->Fill(Rap);
        if (fatherIndex(Particle.pid, fatherpid) > 0) {
            Y1D[index][fatherIndex(Particle.pid, fatherpid)]->Fill(Rap);
        }

    }
    cout<<"dX : "<<dX<<endl;
    // Rescale histograms
    TFile *out = new TFile(outDir+"/hists.root","recreate");
    for (int id = 0; id < _N_IDS_; ++id) {
        TString txtid = TString::Itoa(abs(ids[id]),10);
        if (ids[id] < 0) txtid.Prepend("m");

        for(Int_t i=0; i<_N_HISTOGRAMS_; i++)
            H1D[id][i]->Scale(1.0 / (Events *0.1* dX ));
        // 0.1 : szerokosc binu w Y
        // dX : szerokosc binu w mt

        // ##########################################################################
        // # SAVE HISTOGRAMS TO XML FILE
        // ##########################################################################
        hist2xml(outDir + "/pid" + txtid + ".xml", H1D[id], _N_HISTOGRAMS_, Events, Chain->GetEntries());

        // ##########################################################################
        // # PLOT HISTOGRAMS
        // ##########################################################################
        TCanvas* Canvas = new TCanvas(Form("canvas_%s",txtid.Data()), txtid,10,10,600,600);
        Canvas->SetLogy();

        TH1F *exp = getExp(ids[id]);
        cout << id << " " << ids[id] << " " << exp << endl;
        Int_t binLow =  exp->FindBin(70);
        Int_t binMax  = exp->FindBin(800);
        Int_t binLow1 = H1D[id][0]->FindBin(70);
        Int_t binMax1 = H1D[id][0]->FindBin(800);

        Float_t int_exp = integrate_mt(exp,0,800);
        Float_t int_the = integrate_mt(H1D[id][0],0,800);

        cout<<"scaling factor : "<<int_the<<"/"<<int_exp<<"="<<int_the/int_exp<<endl;

        // Histogram legend
        TLegend* LegendPart	= new TLegend(0.42,0.70,0.71,0.99, "",	"brNDC");
        LegendPart->SetFillColor(0);
        LegendPart->SetFillStyle(1001);
        LegendPart->SetTextSize(0.025);
        LegendPart->SetName("LegendPart");
        for(Int_t i=0; i<maxHist(ids[id]); i++)
            LegendPart->AddEntry(H1D[id][i], H1D[id][i]->GetTitle());
        // Model legend
        TLegend* LegendModel	= new TLegend(0.72,0.70,0.99,0.99, "model",	"brNDC");
        LegendModel->SetFillColor(0);
        LegendModel->SetFillStyle(1001);
        LegendModel->SetTextSize(0.025);
        LegendModel->SetName("LegendModel");
        model2legend(aEventDir, aEventFiles, LegendModel);
        // Histograms
        //H1D[0]->SetTitle("#pi^{+} anatomy");
        H1D[id][0]->SetTitle("");
        H1D[id][0]->Draw("hist");
        for (int i = 1; i <= maxHist(ids[id]); ++i) {
            H1D[id][i]->Draw("histSAME");
        }
        exp->SetMarkerColor(1);
        exp->SetMarkerSize(1.2);
        exp->SetLineColor(1);
        exp->Draw("same");
        if (ids[id] == 211 || ids[id] == -211) {
            LegendPart ->Draw();
        }
        LegendModel->Draw();

        TPaveText *partLabel = setOPT_text2(partNames[id].Data(),0.819398,0.5761821,0.9949833,0.676007,kBlack,0.4);

        out->cd();
        for (int i = 0; i <= maxHist(ids[id]); ++i) {
            H1D[id][0]->Write();
            H1D[id][1]->Write();
            H1D[id][2]->Write();
            H1D[id][3]->Write();
            H1D[id][4]->Write();
            H1D[id][5]->Write();
        }
        LegendPart ->Write();
        LegendModel->Write();
        exp->Write();

        Canvas->SaveAs(outDir+"/pid"+txtid+"_mt.gif");
        Canvas->SaveAs(outDir+"/pid"+txtid+"_mt.eps");

        TCanvas *canY = new TCanvas(Form("canY_%s",txtid.Data()),Form("canY_%s",txtid.Data()),10,10,500,500);
        Float_t yWidth = Y1D[id][0]->GetBinWidth(10);
        Float_t scale_norm = 1.0 / (Events *yWidth );

        cout<<"evt : "<<Events<<"  y width  "<<yWidth<<endl;

        Y1D[id][0]->SetTitle("");
        Y1D[id][0]->Draw("hist");
        for (int i = 1; i <= maxHist(ids[id]); ++i) {
            Y1D[id][i]->Draw("histSAME");
        }
        for (int i = 0; i <= maxHist(ids[id]); ++i) {
            Y1D[id][i]->Scale(scale_norm);
        }
        vector<TObject *> data = getExpY(ids[id]);
        for (int i = 0; i < data.size(); ++i) {
            data[i]->Draw("SAME");
        }
        if (data.size() > 0) {
            Y1D[id][0]->GetYaxis()->SetRangeUser(0.01, dynamic_cast<TH1*>(data[0])->GetMaximumStored()*1.2);
            Y1D[id][0]->GetXaxis()->SetRangeUser(-2,2);
        }
        canY->SaveAs(outDir+"/pid"+txtid+"_y.gif");
        canY->SaveAs(outDir+"/pid"+txtid+"_y.eps");

        partLabel->Draw("SAME");

    }
}
