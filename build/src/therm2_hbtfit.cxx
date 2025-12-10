#include <TH1D.h>
#include <TH3D.h>
#include <TFile.h>
#include <TF1.h>
#include <TF3.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include "Accessibility.h"
#include "Configurator.h"
#include "HBTFit.h"
#include "Messages.h"
#include "Parser.h"
#include "Storage.h"
#include "THGlobal.h"

using namespace std;

int main(int argc, char **argv)
{
    Configurator* sMainConfig;
    HBTFit *hbtFit;
    Storage *store;
    TH3D *ratq;
    TH1D *ratq1,*hExpProj[3][3],*hFitProj[3][3];
    TF3 *funqg,*funqk;
    TF1 *funqg1,*funqk1;
    TDatime tDate;
    TString tInRootName,numname,denname,numname1d,denname1d,partName,sMainINI,sPairType,sEventDir,sTimeStamp,sOvvr;
    bool losl,linv;
    int minkT,maxkT,firstkT,stepkT,dofix[NoParams];
    double lowCut,highCut,lambda,pars[NoParams],parmin[NoParams],parmax[NoParams],maxrange,maxx=0.0, maxy=0.0, maxz=0.0;
    vector<Messages::Summ> outOfBound;
    vector<TString> missFiles;

    Accessibility *accss = new Accessibility();
    
    tDate.Set();
    sTimeStamp = tDate.AsSQLString();

// ##############################################################
// # Command line analysis					
// ##############################################################
  
    minkT = 0;
    maxkT = 0;
    
    if (argc > 1) 
    {
        tInRootName = argv[1];
        if(tInRootName.Contains("-h") || tInRootName.Contains("--help")) 
        {
            Messages::HelpHbt();
            return 0;
        } 
        else if(tInRootName.Contains("-v") || tInRootName.Contains("--version")) 
        {
            Messages::Version();
            return 0;
        }
        if (argc > 2) 
        {
            sPairType = argv[2];
            if(sPairType == "pion-pion" || sPairType == "pionM-pionM")
            {
                if(sPairType == "pion-pion")
                    partName = "pipi";
                else    
                    partName = "piMpiM";

                sMainINI = argv[3];
                if(sMainINI.IsDigit()) 
                {
                    minkT = sMainINI.Atoi();
                    sMainINI = argv[4];
                    if(sMainINI.IsDigit())
                    {
                        maxkT = sMainINI.Atoi();

                        if(minkT > maxkT)
                        {
                            PRINT_DEBUG_1("Command line: wrong kT order");
                            PRINT_MESSAGE("Error: kTmin is greater than kTmax");
                            return _ERROR_FEMTO_WRONG_KT_ORDER_;
                        }
                    }
                    else
                    {
                        PRINT_DEBUG_1("Command line: kTmax = NAN");
                        PRINT_MESSAGE("Error: kTmax is not a number");
                        return _ERROR_FEMTO_WRONG_KT_ORDER_;
                    }
                }
                else
                {
                    PRINT_DEBUG_1("Command line: kTmin = NAN");
                    PRINT_MESSAGE("Error: kTmin is not a number");
                    return _ERROR_FEMTO_WRONG_KT_ORDER_;
                }
                sMainINI = "./hbtfit.ini";
            }
            else
            {
                PRINT_DEBUG_1("Command line: unknown pair type");
                PRINT_MESSAGE("Error: Unknown pair type");
                return _ERROR_FEMTO_UNKNOWN_PAIRTYPE_;
            }
        } 
        else
        {
            PRINT_MESSAGE("Error: unknown argument " + tInRootName);
            return _ERROR_GENERAL_FILE_NOT_FOUND_;
        }
        if (argc > 5)
            sMainINI = argv[5];
    } 
    else 
    {
        Messages::HelpHbt();
        return _ERROR_GENERAL_FILE_NOT_FOUND_;
    }

    Messages::Intro();  

    sEventDir = Accessibility::getEventDir(tInRootName);

// ##############################################################
// # Read configuration file (femtoFit.ini)			
// ############################################################## 

    sMainConfig = new Configurator;
    Parser* tParser = new Parser(sMainINI);
    tParser->ReadINI(sMainConfig);
    delete tParser;
  
// ##############################################################
// # ReadParameters						
// ##############################################################  

    try 
    {
        for(int i = 0; i < NoParams; i++)
            accss->GetParameter(sMainConfig,sParNames[i],&pars[i], &dofix[i], &parmin[i], &parmax[i]);

        maxrange = sMainConfig->GetParameter("MaxFitRange").Atof();
        numname  = sMainConfig->GetParameter("Numerator");
        denname  = sMainConfig->GetParameter("Denominator");
        numname1d  = sMainConfig->GetParameter("Numerator1D");
        denname1d  = sMainConfig->GetParameter("Denominator1D");
        firstkT = sMainConfig->GetParameter("FirstkT").Atoi();
        stepkT = sMainConfig->GetParameter("WidthkT").Atoi();
        sOvvr = sMainConfig->GetParameter("Override");
    }
    catch (TString &tError) 
    {
        PRINT_DEBUG_1("therm2_hbtfit - Caught exception " << tError);
        PRINT_MESSAGE("Did not find one of the necessary parameters in the parameters file.");
        exit(_ERROR_CONFIG_PARAMETER_NOT_FOUND_);
    }
    try 
    {
        maxx = sMainConfig->GetParameter("MaxFitRangeX").Atof();
        maxy = sMainConfig->GetParameter("MaxFitRangeY").Atof();
        maxz = sMainConfig->GetParameter("MaxFitRangeZ").Atof();
    }
    catch (TString &tError) 
    {
    }
    try 
    {
        lowCut = sMainConfig->GetParameter("LambdaCutLow").Atof();
        highCut = sMainConfig->GetParameter("LambdaCutHigh").Atof();
    }
    catch (TString &tError) 
    {
    }
    
    accss->CopyINIFile(sMainINI,sEventDir);

// ##############################################################
// # Log file entry 						
// ##############################################################

    {
        TString tLogName;
        ofstream tLogFile;   
        
        try 
        {
            tLogName = sMainConfig->GetParameter("LogFile");
        }
        catch (TString &tError) 
        {
        }   
        if (!tLogName.IsNull()) 
        {
            tDate.Set();
            
            tLogFile.open(tLogName, ios_base::app);
            if (static_cast<long>(tLogFile.tellp()) == 0) 
                tLogFile << "# THERMINATOR 2 Log File"<<endl;
            
            tLogFile << '['<<tDate.AsSQLString()<<"]\ttherm2_hbtfit\t"<<sPairType<<'\t';
            tLogFile << "[input]\t"<<tInRootName<<'\t'<<sMainINI<<'\t'<<sPairType<<'\t';
            tLogFile << endl;
            tLogFile.close();
        }
    }

// ##############################################################
// # Fit							
// ############################################################## 

    hbtFit = new HBTFit();
    TFile* tInRootFile = new TFile();
    TH3D *numq,*denq;
    TH1D *numq1,*denq1;
    TGraphErrors *gFitRes[NoParams];
    vector<double> parVals[NoParams];
    vector<double> parErrs[NoParams];
    vector<double> kTVal[NoParams];
    int iter2;

    funqg = new TF3("funqg",hbtFit, &HBTFit::fungek, -0.15, 0.15, -0.15, 0.15, -0.15, 0.15, 5);
    if (maxx > 0)
        funqg->SetRange(0, 0, 0, maxx, maxy, maxz);
    else
        funqg->SetRange(-maxrange, -maxrange, -maxrange, maxrange, maxrange, maxrange);
    
    funqk = funqg;
    for (int iter = 2; iter < NoParams; iter++)
    {
        iter2 = iter-2;
        if (dofix[iter]) 
        {
            if (dofix[iter] == 2) 
            {
                funqk->SetParLimits(iter2, parmin[iter], parmax[iter]);
                funqk->SetParameter(iter2, pars[iter]);
            }
            else
                funqk->FixParameter(iter2, pars[iter]);
        }
        else
            funqk->SetParameter(iter2, pars[iter]);

        funqk->SetParName(iter2, sParNames[iter]);
    }

    funqg1 = new TF1("funqg1",hbtFit, &HBTFit::fungek1D, -0.15, 0.15, 3);
    if (maxx > 0)
        funqg1->SetRange(0, 0, 0, maxx, maxy, maxz);
    else
        funqg1->SetRange(-maxrange, -maxrange, -maxrange, maxrange, maxrange, maxrange);
    
    funqk1 = funqg1;
    for (int iter = 0; iter < 3; iter++)
    {
        if (dofix[iter]) 
        {
            if (dofix[iter] == 2) 
            {
                funqk1->SetParLimits(iter, parmin[iter], parmax[iter]);
                funqk1->SetParameter(iter, pars[iter]);
            }
            else
                funqk1->FixParameter(iter, pars[iter]);
        }
        else
            funqk1->SetParameter(iter, pars[iter]);

        funqk1->SetParName(iter, sParNames[iter]);
    }
    
    for(int ii = minkT; ii <= maxkT; ii++)
    {
        lambda = 0;
        PRINT_MESSAGE("["<<sTimeStamp<<"]\tFitting.");

        tInRootName = TString::Format("%sfemto%s%ia.root",sEventDir.Data(),partName.Data(),ii);
        tInRootFile = TFile::Open(tInRootName.Data());
        if(tInRootFile == nullptr)
        {
            missFiles.push_back(tInRootName);
            continue;
        }
            
        numq1 = new TH1D(*((TH1D *) tInRootFile->Get(numname1d)));
        denq1 = new TH1D(*((TH1D *) tInRootFile->Get(denname1d)));
        ratq1 = new TH1D(*numq1);
        ratq1->Reset("ICE");
        ratq1->Divide(numq1, denq1, 1.0, 1.0);
        ratq1->SetName("ratq1");
        ratq1->SetTitle("ratq1");
        ratq1->Fit(funqk1, "RBMPE");

        numq = new TH3D(*((TH3D *) tInRootFile->Get(numname)));
        denq = new TH3D(*((TH3D *) tInRootFile->Get(denname)));
        ratq = new TH3D(*numq);
        ratq->Reset("ICE");
        ratq->Divide(numq, denq, 1.0, 1.0);
        ratq->SetName("ratq");
        ratq->SetTitle("ratq");
        ratq->Fit(funqk, "RBMPE");

// ##############################################################
// # Save fit values						
// ############################################################## 

        TString   tOutTextName;
        double parVal,parErr;
        linv = false;
        losl = false;
        
        tOutTextName = sEventDir + "hbtradii.txt";
        if(ii == minkT)
            store = new Storage(tOutTextName, sOvvr);
        else
            store = new Storage(tOutTextName);

        PRINT_DEBUG_1("Fit results:");
        for (int iter = 0; iter < NoParams; iter++)
        {
            if(iter < 3)
            {
                parVal = fabs(funqk1->GetParameter(iter));
                parErr = fabs(funqk1->GetParError(iter));
                if(iter == 0)
                    lambda = parVal;
            }
            else
            {
                parVal = fabs(funqk->GetParameter(iter-2));
                parErr = fabs(funqk->GetParError(iter-2));
                if(iter == 3)
                    lambda = parVal;
            }

            store->appendToTxt(sParNames[iter],parVal,parErr);
            if(lambda > lowCut && lambda < highCut)
            {
                kTVal[iter].push_back(firstkT+(ii-minkT)*stepkT);
                parVals[iter].push_back(parVal);
                parErrs[iter].push_back(parErr);
            }
            else
            {
                if(iter == 0)
                    linv = true;
                else if(iter == 3)
                    losl = true;
            }
            PRINT_DEBUG_1("\t" + sParNames[iter] + "\t"<< parVal << " +/- " << parErr);
        }
        delete store;    

        if(linv || losl)
            outOfBound.push_back({ii,linv,losl});  
    
// ##############################################################
// # Make plots							
// ############################################################## 

        int padItr = 1;

        TH3D *fitnq = hbtFit->getfitprojc(denq, funqk);

        TCanvas *cancf = new TCanvas ("cancf", "cancf", 1200, 1200);
        gPad->SetFillColor(0);
        gPad->SetFillStyle(4000);
        cancf->Divide(3,3,0.0001,0.0001);

        for(int jj = 0; jj < 3; jj++)
            for(int kk = 0; kk < 3; kk++)
            {
                hExpProj[jj][kk] = hbtFit->getproj(numq, denq, jj, projWidth[kk], pars[0]);
                hFitProj[jj][kk] = hbtFit->getproj(fitnq, denq, jj, projWidth[kk], pars[0]);

                hbtFit->preparehist(hExpProj[jj][kk],jj,kk+1,"CF");
                hbtFit->preparehist(hFitProj[jj][kk],jj,kk+1,"FIT");

                cancf->cd(padItr);
                hbtFit->preparepad(hExpProj[jj][kk],hFitProj[jj][kk]);

                padItr++;
            }

// ##############################################################
// # Save plots							
// ##############################################################   
    
        TString   tOutRootName;
        TFile*    tOutRootFile = new TFile();
        tOutRootName = tInRootName;
        tOutRootName.ReplaceAll("femto","hbtfit");
        tOutRootFile = TFile::Open(tOutRootName,"RECREATE");
        tOutRootFile->cd();

        for(int jj = 0; jj < 3; jj++)
            for(int kk = 0; kk < 3; kk++)
            {
                hExpProj[jj][kk]->Write();
                hFitProj[jj][kk]->Write();
            }

        tOutRootFile->Close();
        tInRootFile->Close();
        tDate.Set(); 
        PRINT_MESSAGE("["<<tDate.AsSQLString()<<"]\tFile "<<tOutTextName<<" and");
        PRINT_MESSAGE("\tfile "<<tOutRootName<<" written.");

        delete cancf;

// ##############################################################
// # Log file entry 						
// ##############################################################

        {
            TString  tLogName;
            TString  tORootName;
            TString  tOTextName;
            ofstream tLogFile;   
            
            try 
            {
                tLogName = sMainConfig->GetParameter("LogFile");
            }
            catch (TString &tError) 
            {
            }   
            if (!tLogName.IsNull()) 
            {
                tDate.Set();
                tORootName = tOutRootName;
                tOTextName = tOutTextName;
                tORootName.ReplaceAll(sEventDir,"");      
                tOTextName.ReplaceAll(sEventDir,"");
                
                tLogFile.open(tLogName, ios_base::app);
                if (static_cast<long>(tLogFile.tellp()) == 0) 
                    tLogFile << "# THERMINATOR 2 Log File"<<endl;

                tLogFile << '['<<tDate.AsSQLString()<<"]\ttherm2_hbtfit\t"<<sPairType<<'\t';
                tLogFile << "[output]\t"<<sEventDir<<'\t'<<tORootName<<'\t'<<tOTextName<<'\t';
                tLogFile << endl;
                tLogFile.close();
            }    
        }
    
// ##############################################################
// # Temporary file						
// ##############################################################  

        ofstream tTmpFile;
        TString tTmpFileName;
        
        tTmpFileName = TString::Format("./hbtfit_%s.tmp",sPairType.Data()); 
        tTmpFile.open(tTmpFileName, ios_base::app);
        if((tTmpFile) && tTmpFile.is_open()) 
            tTmpFile << tOutRootName << endl;
        else
        {
            PRINT_MESSAGE("Unable to create temp file "<<tTmpFileName);
            exit(_ERROR_GENERAL_FILE_NOT_FOUND_);
        }
    }

    TFile *paramFile = new TFile();
    paramFile = TFile::Open(Form("%sparameterFit%s.root",sEventDir.Data(),partName.Data()),"RECREATE");
    for (int iter = 0; iter < NoParams; iter++) 
    {
        gFitRes[iter] = new TGraphErrors(kTVal[iter].size(),&kTVal[iter][0],&parVals[iter][0],0,&parErrs[iter][0]);
        gFitRes[iter]->Write(Form("g%s_%i_%i",sParNames[iter].Data(),minkT,maxkT));
    }
    
    paramFile->Close();

    Messages::Summary(lowCut,highCut,missFiles,outOfBound);

    return 0;
}
