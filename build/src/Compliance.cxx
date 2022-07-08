#include <TMath.h>
#include "Accessibility.h"
#include "Compliance.h"
#include "THGlobal.h"

using namespace std;

Compliance::Compliance():mModName("")
{
}

Compliance::Compliance(TString sModFilePath, TString sModName, TString sExpFilePath, TString sExpName):mModName(sModName)
{
    mEFile = TFile::Open(sExpFilePath);
    if(!mEFile->IsOpen() || mEFile->IsZombie()) 
    {
        PRINT_MESSAGE("<Compliance::Compliance>\tFile "<<sExpFilePath<<" not found.");
        exit(_ERROR_GENERAL_FILE_NOT_FOUND_);
    }

    mMFile = TFile::Open(sModFilePath);
    if(!mMFile->IsOpen() || mMFile->IsZombie()) 
    {
        PRINT_MESSAGE("<Compliance::Compliance>\tFile "<<sModFilePath<<" not found.");
        exit(_ERROR_GENERAL_FILE_NOT_FOUND_);
    }

    mMod = (TGraphAsymmErrors*) mMFile->Get(sModName);
    mExp = (TGraph*) mEFile->Get(sExpName);

    if(mMod == nullptr || mExp == nullptr)
    {
        PRINT_MESSAGE("<Compliance::Compliance>\tGraphs are null");
        exit(_ERROR_GENERAL_UNSUPORTED_VALUE_);
    }

    mEventDir = Accessibility::getEventDir(sModFilePath);
}

Compliance::~Compliance()
{
    if(mEFile->IsOpen())
        mEFile->Close();

    if(mMFile->IsOpen())
        mMFile->Close();
}

double Compliance::printResult(TString testType, bool moreInfo)
{
    double totDif,testRes;
    int ndf;

    if(testType == "Q2test")
    {
        testRes = Q2Test(ndf,totDif);

        if(testRes < 0)
            PRINT_MESSAGE("<Compliance::printResult>\tNDF is zero.")
        else
        {
            PRINT_MESSAGE("================== Q2 Test Result ==================");
            PRINT_MESSAGE(Form("Parameter:\t%s",mModName.Data()));
            PRINT_MESSAGE(Form("Q2:\t%f",testRes));
            if(moreInfo)
                PRINT_MESSAGE(Form("Total difference on x axis per NDF:\t%.2f/%d",totDif,ndf));
            PRINT_MESSAGE("====================================================");
        }
    }
    else if(testType == "Chi2test")
    {
        testRes = Chi2Test(ndf,totDif);

         if(testRes < 0)
            PRINT_MESSAGE("<Compliance::printResult>\tNDF is zero.")
        else
        {
            PRINT_MESSAGE("================= Chi2 Test Result =================");
            PRINT_MESSAGE(Form("Parameter:\t%s",mModName.Data()));
            PRINT_MESSAGE(Form("Chi2:\t%f",testRes));
            if(moreInfo)
                PRINT_MESSAGE(Form("Total difference on x axis per NDF:\t%.2f/%d",totDif,ndf));
            PRINT_MESSAGE("====================================================");
        }
    }
    else
    {
        PRINT_MESSAGE("<Compliance::printResult>\tUnknown test type:" << testType)
        exit(_ERROR_GENERAL_UNSUPORTED_VALUE_);
    }

    return testRes;
}

bool Compliance::getModelName(TString expName, int iter, int minkT, int maxkT, TString &modName)
{
    if(expName == "none")
    {
        modName = "";
        return false;
    }
    if(iter > 1)
        iter++;
    
    modName = TString::Format("g%s_%i_%i",sParNames[iter].Data(),minkT,maxkT);
    return true;
}

double Compliance::Q2Test(int &ndf,double &totDiff)
{
    ndf = mMod->GetN();
    int iter;
    double xMod,yMod,xExp,yExp,Q2 = 0;

    totDiff = 0;

    for(int i = 0; i < ndf; i++)
    {
        mMod->GetPoint(i,xMod,yMod);
        totDiff += getClosest(xMod,iter,xExp,yExp);
        if(yExp > 0)
            Q2 += (yMod*yMod - 2*yMod*yExp + yExp*yExp)/(yExp*yExp);
    }

    if(ndf != 0)
        return TMath::Abs(Q2)/ndf;
    else
        return -1;
}

double Compliance::Chi2Test(int &ndf, double &totDiff)
{
    ndf = mMod->GetN();
    int iter;
    double xMod,yMod,xExp,yExp,yErr,Chi2 = 0;
    totDiff = 0;

    for(int i = 0; i < ndf; i++)
    {
        mMod->GetPoint(i,xMod,yMod);
        yErr = -1;

        totDiff += getClosest(xMod,iter,xExp,yExp);
        if(yMod < yExp)
            yErr = mMod->GetErrorYlow(i);
        else
            yErr = mMod->GetErrorYhigh(i);
        if(yErr > 0)
            Chi2 += (yMod*yMod - 2*yMod*yExp + yExp*yExp)/yErr;
    }

    ndf--;
    if(ndf != 0)
        return TMath::Abs(Chi2)/ndf;
    else
        return -1;
}

double Compliance::getClosest(double xVal,int &i, double &x, double &y)
{
    const int elem = mExp->GetN();
    double xx,yy;
    double dist,distTemp;

    mExp->GetPoint(0,x,y);
    dist = TMath::Abs(xVal-x);
    if(dist < numeric_limits<double>::epsilon()) // baisically means "if(xVal == x)" but for floating point variables
        return dist;

    for(int ii = 1; ii < elem; ii++)
    {
        mExp->GetPoint(ii,xx,yy);
        distTemp = TMath::Abs(xVal-xx);
        if(distTemp < dist)
        {
            dist = distTemp;
            i = ii;
            x = xx;
            y = yy;
        }
    }
    return dist;
}

