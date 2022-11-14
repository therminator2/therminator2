#include <fstream>
#include "Accessibility.h"
#include "THGlobal.h"

using namespace std;

Accessibility::Accessibility()
{
}

Accessibility::~Accessibility()
{
}

void Accessibility::GetParameter(Configurator* aINI,TString aKeyword, Double_t *parval, Int_t *isfixed, Double_t *parmin, Double_t *parmax)
{
    *parval = aINI->GetParameter(aKeyword).Atof(); 
    if (aINI->GetParameter(aKeyword+"FitType") == "fixed") 
    {
        *isfixed = 1;
        PRINT_DEBUG_1("Fitting parameter: "<<aKeyword<<" = "<<(*parval)<<" fixed.");
    } 
    else if (aINI->GetParameter(aKeyword+"FitType") == "limit") 
    {
        *parmin  = aINI->GetParameter(aKeyword+"Min").Atof();
        *parmax  = aINI->GetParameter(aKeyword+"Max").Atof();
        *isfixed = 2;
        PRINT_DEBUG_1("Fitting parameter: "<<aKeyword<<" = "<<(*parval)<<" limited ["<<(*parmin)<<", "<<(*parmax)<<"]");
    } 
    else 
    {
        *isfixed = 0;
        PRINT_DEBUG_1("Fitting parameter: "<<aKeyword<<" = "<<(*parval)<<" free.");
    }
}

void Accessibility::CopyINIFile(TString sMainINI, TString sEventDir)
{
    TString  tINI;
    ifstream ifs;
    ofstream ofs;
    
    tINI = sEventDir + sMainINI;
    ifs.open(sMainINI.Data(), std::ios::binary);
    ofs.open(tINI.Data(), std::ios::binary);
    if((ifs) && (ofs) && ifs.is_open() && ofs.is_open()) 
    {
        ofs << ifs.rdbuf();
        ifs.close();
        ofs.close();
    } 
    else
        PRINT_MESSAGE("<Accessibility::CopyINIFile()>\tUnable to copy "<<sMainINI<<" to "<<tINI);
}

TString Accessibility::getEventDir(TString sFilePath)
{
    int tPos=0;
    for(int i=1; i<sFilePath.Length(); i++) 
        if(sFilePath[i]=='/')
            tPos = i;
    
    sFilePath.Resize(tPos+1);

    return sFilePath;
}
