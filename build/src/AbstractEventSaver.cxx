#include "THGlobal.h"
#include "Configurator.h"
#include "StructEvent.h"
#include "AbstractEventSaver.h"
#include "Event.h"
#include <fstream>

using namespace std;

extern void AddLogEntry(const char* aEntry);
extern Configurator *sMainConfig;
extern TString	sEventDIR;
extern int      sParentPID;                                                                                                                                                                
extern int	    sModel;
extern TString	sTimeStamp;

AbstractEventSaver::AbstractEventSaver() : mFileCounter(0), kEventsPerFile(_EVENTS_PER_FILE_), mNumberOfEvents(0), mParameterTree(0) {
    FindPreviousEventFiles();
    ReadParameters();
}
AbstractEventSaver::~AbstractEventSaver() {
  char tBuff[2*kFileNameMaxChar];
  int  tFileCounter;
  
  if(mNumberOfEvents) {
    tFileCounter = mFileCounter - (mNumberOfEvents / kEventsPerFile) - (((mNumberOfEvents % kEventsPerFile) > 0) ? (1) : (0));
    if(mFileCounter-tFileCounter-1 == 0)
      sprintf(tBuff,"[output]\t%s\t\"event%03i.root\"\t[events]\t%i\t%i\t%i",sEventDIR.Data(),tFileCounter,(mFileCounter-tFileCounter),mNumberOfEvents,kEventsPerFile);  
    else
      sprintf(tBuff,"[output]\t%s\t\"event%03i.root-event%03i.root\"\t[events]\t%i\t%i\t%i",sEventDIR.Data(),tFileCounter,mFileCounter-1,(mFileCounter-tFileCounter),mNumberOfEvents,kEventsPerFile);
    AddLogEntry(tBuff);
  }
}

void AbstractEventSaver::FindPreviousEventFiles()
{
  // if previous files in this directory are not in a sequence then they may be overwritten. 
  fstream tFile;
  char	tFileName[kFileNameMaxChar];
  
  mFileCounter = 0;
  do {
    tFile.clear(std::ios::failbit);
    sprintf(tFileName,"%sevent%03i.root",sEventDIR.Data(),mFileCounter);
    tFile.open(tFileName);
    tFile.close();
    if(!tFile.fail())
      mFileCounter++ ;
  } while(!tFile.fail());
  
  if(mFileCounter) {
    sprintf(tFileName,"event%03i.root",mFileCounter);
    PRINT_DEBUG_1("<EventGenerator::FindPreviousEventFiles>\tFound "<<mFileCounter<<" previous event file(s) in "<<sEventDIR.Data()<<".");
    PRINT_MESSAGE("\tEvent files names from this run will start with \""<<tFileName<<"\"");
  }
}
  
void AbstractEventSaver::SetEventsTemp() {
  ofstream tFile;
  char tFileName[kFileNameMaxChar];
  
  sprintf(tFileName,"./events_%i.tmp",sParentPID);
  tFile.open(tFileName);
  if((tFile) && tFile.is_open()) {
    tFile << sEventDIR    << endl;
    tFile << mFileCounter << endl;
    tFile.close();
  } else {
    PRINT_MESSAGE("<EventGenerator::SetEventsINI>\tUnable to create file "<<tFileName);
    exit(_ERROR_GENERAL_FILE_NOT_FOUND_);
  }
}

void AbstractEventSaver::ReadParameters()
{
  try {
    mNumberOfEvents	= (sMainConfig->GetParameter("NumberOfEvents")).Atoi();
    sEventDIR	= sMainConfig->GetParameter("EventDir"); sEventDIR.Prepend("./");
  }
  catch (TString &tError) {
    PRINT_MESSAGE("<EventGenerator::ReadParameters>\tCaught exception " << tError);
    PRINT_MESSAGE("\tDid not find one of the necessary parameters in the parameters file.");
    exit(_ERROR_CONFIG_PARAMETER_NOT_FOUND_);
  }
}

void AbstractEventSaver::SaveParameters(Model *tModel)
{
  // Events read from files (therm2_coulomb) come without a model
  if (tModel == nullptr) {
    PRINT_MESSAGE("<AbstractEventSaver::SaveParameters>\tNo model, the parameters tree is not written");
    return;
  }
  Configurator *aMainConfig = tModel->GetMainConfig();
        mParameterTree = new TTree(_PARAMETERS_TREE_,"parameters and model description tree");
char tTimeStamp[21];
        sprintf(tTimeStamp,"%s",sTimeStamp.Data());
        
        if (aMainConfig->HasParameter("MaxIntegrateSamples")) {
            mMaxIntegrationSamples = aMainConfig->GetParameter("MaxIntegrateSamples").Atoi();
        }
        if (aMainConfig->HasParameter("IntegTolerance")) {
            mIntegTolerance = aMainConfig->GetParameter("IntegTolerance").Atof();
        }
        if (aMainConfig->HasParameter("IntegToleranceInterval")) {
            mIntegToleranceInterval = aMainConfig->GetParameter("IntegToleranceInterval").Atoi();
        }
        if (aMainConfig->HasParameter("IntegToleranceNSuccessive")) {
            mIntegToleranceNSuccessive = aMainConfig->GetParameter("IntegToleranceNSuccessive").Atoi();
        }
        if (aMainConfig->HasParameter("Randomize")) {
            mRandomize = aMainConfig->GetParameter("Randomize");
        }

        mParameterTree->Branch(_MAX_INTEGRATESAMPLE_BRANCH_, &mMaxIntegrationSamples, "i");
        mParameterTree->Branch(_INTEG_TOLERANCE_BRANCH_, &mIntegTolerance, "f");
        mParameterTree->Branch(_INTEG_TOLERANCE_INTERVAL_BRANCH_, &mIntegToleranceInterval, "i");
        mParameterTree->Branch(_INTEG_TOLERANCE_N_SUCCESSIVE_BRANCH_, &mIntegToleranceNSuccessive, "i");

        mParameterTree->Branch(_RANDOMIZE_BRANCH_,		(UInt_t*) &mRandomize,						 "i"			);
        mParameterTree->Branch(_TIMESTAMP_BRANCH_,		(Char_t*) tTimeStamp,						 _TIMESTAMP_FORMAT_	);
        mParameterTree->Branch(_MODELID_BRANCH_,		(UInt_t*) &sModel,						 "i"			);
        mParameterTree->Branch(_MODELNAME_BRANCH_,		(Char_t*) tModel->GetName(),	 _MODELNAME_FORMAT_	);
        mParameterTree->Branch(_MODELHASH_BRANCH_,		(Char_t*) tModel->GetHash(),	 _MODELHASH_FORMAT_	);
        mParameterTree->Branch(_MODELDESCRIPTION_BRANCH_,	(Char_t*) tModel->GetDescription(), _MODELDESCRIPTION_FORMAT_);    
        tModel->AddParameterBranch(mParameterTree);
        mParameterTree->Fill();    

}

