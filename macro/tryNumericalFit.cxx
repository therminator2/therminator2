#include "Math/Integrator.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/AdaptiveIntegratorMultiDim.h"
#include "Math/WrappedParamFunction.h"
#include "Math/SpecFunc.h"
#include "Math/Functor.h"
#include "Math/ParamFunctor.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"

double Vr(double H, double r)
{
    return tanh(H*r);
}

double gamma(double delt, double H, double theta, double r)
{
    return 1. / sqrt(1 - ((1 + delt * cos(2*theta)) * Vr(H,r) * Vr(H,r)));
}

double Ep(double m, double pt, double Y)
{
    return sqrt(m*m + pt*pt) * cosh(Y);
}

double pz(double m, double pt, double Y)
{
    return sqrt(m*m + pt*pt) * sinh(Y);
}

double integrand(const double *x, const double *par)
{
    return x[0]*x[0]*sin(x[1])*Ep(par[0],par[1],par[2])*exp((- gamma(par[3],par[4],x[1],x[0]) * Ep(par[0],par[1],par[2]))/par[5])*exp((gamma(par[3],par[4],x[1],x[0]) * Vr(par[4],x[0]) * sqrt(1 + par[3]) * pz(par[0],par[1],par[2]) * cos(x[1])) / par[5]) * ROOT::Math::cyl_bessel_i(0,((gamma(par[3],par[4],x[1],x[0]) * Vr(par[4],x[0]))/par[5]) * sqrt(1 - par[3]) * par[1] * sin(x[1]));
}
/*
par:
0 - m
1 - pt
2 - Y
3 - delt
4 - H
5 - T
*/

double ptModSR(double *x, double *par)
{
    double pt = sqrt((par[0] + x[0])*(par[0] + x[0]) - par[0]*par[0]);
    ROOT::Math::IntegratorMultiDimOptions::SetDefaultAbsTolerance(1.e-6);
    double a[] = {0,0};
    double b[] = {par[5]/0.1973269631,M_PI};

    ROOT::Math::WrappedParamFunction<> wf(&integrand,2,6);
    double p[] = {par[0],pt,par[1],par[2],par[3],par[4]};
    wf.SetParameters(p);
    ROOT::Math::AdaptiveIntegratorMultiDim ig;
    ig.SetFunction(wf);
    return sqrt(1+par[6]) * (1-par[6]) * (2*par[7]+1)/(2*M_PI) * pow(par[8],par[9]) * exp((par[10]*par[11] + par[12]*par[13] + par[14]*par[15])/par[4]) *  ig.Integral(a,b) * 1.e3 * (par[0] + x[0]) / (1.e6);
}
/*
par:
0 - m
1 - Y
2 - delt
3 - H
4 - T
5 - R
6 - eps
7 - s
8 - gammaS
9 - NS
10 - muB
11 - B
12 - muS
13 - S
14 - muI3
15 - I3
*/
void tryNumericalFit()
{
    const TString key = "Proton"; // Proton PionPlus PionMinus Deuteron Triton Helium3

    std::map<TString,double> massMap = {{"Proton",0.938272},{"PionPlus",0.1395699},{"PionMinus",0.1395699},{"Deuteron",1.87561},{"Triton",2.80892},{"Helium3",2.80839}};
    std::map<TString,double> sMap = {{"Proton",0.5},{"PionPlus",0.},{"PionMinus",0.},{"Deuteron",0.},{"Triton",0.5},{"Helium3",0.5}};
    std::map<TString,double> BMap = {{"Proton",1.},{"PionPlus",0.},{"PionMinus",0.},{"Deuteron",2.},{"Triton",3.},{"Helium3",3.}};
    std::map<TString,double> I3Map = {{"Proton",0.5},{"PionPlus",1.},{"PionMinus",-1.},{"Deuteron",0.},{"Triton",-0.5},{"Helium3",0.5}};
    std::map<TString,double> scaleMap = {{"Proton",0.9229},{"PionPlus",1.242},{"PionMinus",1.019},{"Deuteron",1.121},{"Triton",0.5504},{"Helium3",0.5001}};
    double m = massMap[key], Y = 0, delt = 0.3, H = 0.016, T = 0.0625, R = 9.215, eps = 0, s = sMap[key], gammaS = 1., NS = 0, muB = 0.792, B = BMap[key], muS = 0.1505, S = 0, muI3 = -0.0252, I3 = I3Map[key]; //
    TF1 *f = new TF1("mt dist of protons",ptModSR,0.,1.,16);
    double p[] = {m,Y,delt,H,T,R,eps,s,gammaS,NS,muB,B,muS,S,muI3,I3};
    TString name[] = {"m0","Y","delta","H","T","R","epsilon","s","gammaS","NS","muB","B","muS","S","muI3","I3"};
    f->SetParameters(p);
    for (int i = 0; i < 16; i++)
    {            
        f->SetParName(i,name[i]);
        if (i == 2)
            f->SetParLimits(i,0.,0.8);
        else if (i == 3)
            f->SetParLimits(i,0.01,0.021);
        else if (i == 5)
            f->SetParLimits(i,7.,12.);
        else
            f->FixParameter(i,p[i]);
    }

    TFile *inpFile = TFile::Open(Form("/home/jedkol/Downloads/HADES/modelSR/NewExpSpectra/%sGeV.root",key.Data()));
    TH1D *hExp = static_cast<TH1D*>(inpFile->Get("hMtYSliceM0Std"));
    hExp->SetMarkerStyle(20);
    hExp->SetMarkerColor(1);
    hExp->Scale(scaleMap[key]);
    hExp->Fit(f,"L");

    TCanvas *c = new TCanvas("c","",800,450);
    c->SetLogy();
    f->Draw();
    hExp->Draw("hist p same");
}