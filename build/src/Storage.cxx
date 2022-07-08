#include <TDatime.h>
#include "THGlobal.h"
#include "Storage.h"

Storage::Storage()
{
}

Storage::Storage(TString sFileNam,TString ovrrd)
{
    TDatime tDate;
    TString sTime;
    if(ovrrd == "false")
        ofStr = new ofstream(sFileNam,ios_base::app);
    else if(ovrrd == "true")
        ofStr = new ofstream(sFileNam,ios_base::out);
    else
    {
        PRINT_MESSAGE("<Storage::Storage>\tOverride not bool");
        exit(_ERROR_CONFIG_PARAMETER_NOT_FOUND_);
    }

    tDate.Set();
    sTime = tDate.AsSQLString();

    *ofStr << "[" << sTime << "]" << endl;
}

Storage::~Storage()
{
    delete ofStr;
}

void Storage::appendToTxt(TString sValName, double parVal, double parErr)
{
    (*ofStr) << sValName << "\t" << parVal << "\t" << parErr << endl;
}

void Storage::appendToTxt(TString sValName, double parVal)
{
    (*ofStr) << sValName << "\t" << parVal << endl;
}