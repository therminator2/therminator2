#include <TPad.h>

#include "HBTFit.h"
#include "THGlobal.h"

HBTFit::HBTFit()
{
}

HBTFit::~HBTFit()
{
}

TH3D* HBTFit::getfitprojc(TH3D *expden, TF3 *funqk)
{
    TH3D *projhist = new TH3D(*expden);

    for (int q1int = 1; q1int <= expden->GetNbinsX(); q1int++)
        for (int q2int = 1; q2int <= expden->GetNbinsY(); q2int++) 
            for (int q3int = 1; q3int <= expden->GetNbinsZ(); q3int++) 
                projhist->SetBinContent(q1int, q2int, q3int,expden->GetBinContent(q1int, q2int, q3int) * funqk->Eval(expden->GetXaxis()->GetBinCenter(q1int),expden->GetYaxis()->GetBinCenter(q2int),expden->GetZaxis()->GetBinCenter(q3int)));

    return projhist;
}

void HBTFit::preparepad(TH1D *hExp, TH1D* hFit)
{
    gPad->SetFillStyle(4000);
    gPad->SetFillColor(0);
    gPad->SetRightMargin(0.02);
    gPad->SetTopMargin(0.03);
    gPad->SetBottomMargin(0.12);
    gPad->SetGridy();

    hExp->Draw("HISTPE1");
    hFit->Draw("SAMEHISTL");
}

void HBTFit::preparehist(TH1D *hist, int projType, int wType, TString type)
{
    if(projType < 0 || projType > 2 || wType < 1 || wType > 3)
    {
        PRINT_MESSAGE("<HBTFit::preparehist>: Unknown projection type or width");
        exit(_ERROR_GENERAL_UNSUPORTED_VALUE_);
    }

    if (!type.CompareTo("FIT")) 
    {
        hist->SetLineColor(8);
        hist->SetLineWidth(2);
        hist->SetLineStyle(2);
    }
    else if(!type.CompareTo("CF"))
    {
        hist->SetMarkerSize(1.0);
        hist->SetMarkerColor(2);
        hist->SetMarkerStyle(kOpenCircle);
        hist->SetTitle(Form(";q_{%s} [GeV/c];C(q_{%s})",sProjNames[projType].Data(),sProjNames[projType].Data()));
        hist->SetMinimum(hist->GetMinimum()*0.9);
        hist->SetMaximum(hist->GetMaximum()*1.1);
        hist->GetXaxis()->SetLabelSize(0.055);
        hist->GetYaxis()->SetLabelSize(0.055);
        hist->GetXaxis()->SetTitleSize(0.055);
        hist->GetYaxis()->SetTitleSize(0.055);
        hist->GetYaxis()->SetTitleOffset(0.8);
    }
    else
    {
        PRINT_MESSAGE("<HBTFit::preparehist>: Unknown histogram type");
        exit(_ERROR_GENERAL_UNSUPORTED_VALUE_);
    }

    hist->SetName(Form("%sproj_%s_%d",type.Data(),sProjNames[projType].Data(),wType));
}

TH1D* HBTFit::getproj(TH3D *numq, TH3D *denq, int nproj, int wbin, double norm)
{
    Double_t intnexp = norm;
    Double_t intdexp = 1.0;
    TString sProj = nullptr;
    TH1D *hproj,*denbuf,*numbuf;

    Int_t binc = numq->GetXaxis()->FindFixBin(0.0);
    Int_t binm = binc + wbin;

    numq->GetXaxis()->SetRange(binc, binm);
    numq->GetYaxis()->SetRange(binc, binm);
    numq->GetZaxis()->SetRange(binc, binm);
    denq->GetXaxis()->SetRange(binc, binm);
    denq->GetYaxis()->SetRange(binc, binm);
    denq->GetZaxis()->SetRange(binc, binm);

    switch (nproj) 
    {
        case 0:
            sProj = "x";
            denq->GetXaxis()->SetRange(1,numq->GetNbinsX());
            numq->GetXaxis()->SetRange(1,numq->GetNbinsX());
            break;
        case 1:
            sProj = "y";
            denq->GetYaxis()->SetRange(1,numq->GetNbinsY());
            numq->GetYaxis()->SetRange(1,numq->GetNbinsY());
            break;
        case 2:
            sProj = "z";
            denq->GetZaxis()->SetRange(1,numq->GetNbinsZ());
            numq->GetZaxis()->SetRange(1,numq->GetNbinsZ());
            break;
        default:
            PRINT_MESSAGE("<HBTFit::getproj>: Unknown projection type");
            exit(_ERROR_GENERAL_UNSUPORTED_VALUE_);
    }
    
    denbuf = new TH1D(*((TH1D *) denq->Project3D(sProj)));
    numbuf = new TH1D(*((TH1D *) numq->Project3D(sProj)));

    hproj = new TH1D(*numbuf);
    hproj->Sumw2();
    hproj->Reset("ICE");
    hproj->Divide(numbuf, denbuf, 1.0, 1.0, "");

    for (int iter=1; iter<hproj->GetNbinsX(); iter++)
        if (numbuf->GetBinContent(iter)) 
        {
            Double_t dn = numbuf->GetBinError(iter);
            Double_t an = numbuf->GetBinContent(iter);
            Double_t dd = denbuf->GetBinError(iter);
            Double_t ad = denbuf->GetBinContent(iter);
            hproj->SetBinError(iter, TMath::Sqrt((dn*dn*ad*ad + dd*dd*an*an + dd*dd*dn*dn)/(ad*ad*ad*ad)));
        }

    hproj->Scale(intdexp/intnexp);

    return hproj;
}

Double_t HBTFit::fungek(Double_t *x, Double_t *par)
{
    Double_t qosq = x[0]*x[0];
    Double_t qssq = x[1]*x[1];
    Double_t qlsq = x[2]*x[2];
    Double_t lam =  TMath::Abs(par[1]);

    Double_t gpart = exp((-par[2]*par[2]*qosq-par[3]*par[3]*qssq-par[4]*par[4]*qlsq)/0.038938);

    return (par[0] * (1 + lam*gpart));
}

Double_t HBTFit::fungek1D(Double_t *x, Double_t *par)
{
    Double_t qinv = x[0]*x[0];
    Double_t lam =  TMath::Abs(par[0]);
    Double_t gpart = exp((-par[1]*par[1]*qinv)/0.038938);

    return (par[2] * (1 + lam*gpart));
}