#include <iostream>
#include <TString.h>
#include "Accessibility.h"
#include "Compliance.h"
#include "Configurator.h"
#include "Messages.h"
#include "Parser.h"
#include "THGlobal.h"
#include "Storage.h"

using namespace std;

int main(int argc, char **argv)
{
    // ##############################################################
    // # Constants definition					
    // ##############################################################

    const TString sParList[] = {"LambdaInv","RadiusInv","LambdaOSL","RadiusOut","RadiusSide","RadiusLong"};
    const int len = (int) sizeof(sParList)/sizeof(sParList[0]);

    // ##############################################################
    // # Variable initialisation					
    // ##############################################################

    int noGraphs,minkT,maxkT;
    TString sFilePath,sMainINI,sEventDir,sExpFilePath,sParam,sModelName,sExpGraphName[len],sOvrrd,sTest,sCompResFile;  
    Configurator *sMainConfig;
    Accessibility *accss;
    Compliance *comp;
    Parser *tParser;
    Storage *store;

    // ##############################################################
    // # Command line analysis					
    // ##############################################################

    sMainINI = "./statTest.ini";

    if(argc > 1)
    {
        sFilePath = argv[1];
        if(sFilePath.Contains("-h") || sFilePath.Contains("--help")) 
        {
            Messages::HelpTst();
            return 0;
        }
        else if(sFilePath.Contains("-v") || sFilePath.Contains("--version")) 
        {
            Messages::Version();
            return 0;
        }
        if(argc > 2)
            sMainINI = argv[2];

        if(argc > 3)
            PRINT_MESSAGE("Warning: Omitting unnecessary arguments");
    }
    else 
    {
        Messages::HelpTst();
        return _ERROR_GENERAL_FILE_NOT_FOUND_;
    }

    Messages::Intro();  

    sEventDir = Accessibility::getEventDir(sFilePath);

    // ##############################################################
    // # Read configuration file (statTest.ini)			
    // ############################################################## 

    sMainConfig = new Configurator;
    tParser = new Parser(sMainINI);
    tParser->ReadINI(sMainConfig);
    delete tParser;
  
    // ##############################################################
    // # ReadParameters						
    // ##############################################################  

    accss = new Accessibility;
    noGraphs = 0;

    try 
    {
        sExpFilePath = sMainConfig->GetParameter("Loc");
        for(int i = 0; i < len; i++)
        {
            sParam = sMainConfig->GetParameter(sParList[i].Data());
            sExpGraphName[i] = sParam;
            if(sParam != "none")
                noGraphs++;
        }
        
        minkT = sMainConfig->GetParameter("MinkT").Atoi();
        maxkT = sMainConfig->GetParameter("MaxkT").Atoi(); 

        sOvrrd = sMainConfig->GetParameter("Override");

        sTest = sMainConfig->GetParameter("Test");
    }
    catch (TString &tError) 
    {
        PRINT_DEBUG_1("Q2test - Caught exception " << tError);
        PRINT_MESSAGE("Did not find one of the necessary parameters in the parameters file.");
        exit(_ERROR_CONFIG_PARAMETER_NOT_FOUND_);
    }
    
    accss->CopyINIFile(sMainINI,sEventDir);

    if(noGraphs == 0)
    {
        PRINT_MESSAGE("All graphs in " + sMainINI + " set to \"none\". Finishing execution.");
        exit(0);
    }

    // ##############################################################
    // # Fit
    // ##############################################################

    sCompResFile = "";

    if(!sTest.CompareTo("Q2test"))
        sCompResFile = "Q2result.txt";
    else if(!sTest.CompareTo("Chi2test"))
        sCompResFile = "Chi2result.txt";
    else
    {
        PRINT_MESSAGE("Unknown compliance test type. Aborting execution.");
        exit(_ERROR_GENERAL_UNSUPORTED_VALUE_);
    }

    store = new Storage(sEventDir + sCompResFile,sOvrrd);
    
    for(int i = 0; i < len; i++)
    {
        if(Compliance::getModelName(sExpGraphName[i],i,minkT,maxkT,sModelName))
        {
            comp = new Compliance(sFilePath,sModelName,sExpFilePath,sExpGraphName[i]);
            store->appendToTxt(sModelName,comp->printResult(sTest,true));
            delete comp;
        }
    }

    delete store;
    return 0;
}
