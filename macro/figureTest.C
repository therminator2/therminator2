#include <iostream>
#include <fstream>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH2.h>
#include <TString.h>
#include <TStyle.h>
#include <TColor.h>
#include <TSystem.h>
#include <TLegend.h>
#include <TAxis.h>

using namespace std;

void figureTest()
{
    gStyle->SetOptStat(0);
    
    const int noE = 7;
    const TString folderPath = "/home/jedkol/Downloads/modelSR/MotoMotoTemp";
    const TString epsType[] = {"E6","E5","E4","E3","E2","E1","E0"};
    const TString modelName[] = {"gRinv_5_12","gRout_5_12","gRside_5_12","gRlong_5_12"};
    const TString rName[] = {"R_{inv}","R_{out}","R_{side}","R_{long}"};
    const double epsVal[] = {-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0};

    ifstream ifFile;
    TGraph *gRes[4];
    TH2D *hRes;
    TCanvas *canv1, *canv2;
    TString sName;
    double yVal[4][noE];
    TLegend *tLeg;

    for(int i = 0; i < noE; i++)
    {
        ifFile.open(Form("%s/H165%sD4femto/Q2result.txt",folderPath.Data(),epsType[i].Data()));
        ifFile >> sName >> sName;
        for(int j = 0; j < 4; j++)
            ifFile >> sName >> yVal[j][i];

        ifFile.close();
    }

    hRes = new TH2D("hRes","",7,-0.65,0.05,4,0.5,4.5);
    

    canv1 = new TCanvas("canv1","canv1",800,1000);
    canv1->Divide(1,4);

    canv2 = new TCanvas("canv2","canv2",800,800);

    for(int i = 0; i < 4; i++)
    {
        canv1->cd(i+1);
        gRes[i] = new TGraph(noE,epsVal,&yVal[i][0]);
        //gRes[i]->SetMarkerStyle(20);
        //gRes[i]->SetMarkerColor(kRed);
        gRes[i]->SetLineColor(kRed);
        gRes[i]->SetLineWidth(5);
        gRes[i]->SetTitle(rName[i].Data());
        gRes[i]->GetXaxis()->SetRangeUser(-0.7,0.1);
        
        gRes[i]->Draw("ac");

        for(int j = 0; j < noE; j++)
            hRes->SetBinContent(j+1,i+1,yVal[i][j]);

        hRes->GetYaxis()->SetBinLabel(i+1,rName[i]);
    }

    canv2->cd();
    hRes->Draw("lego2 fb bb");
}