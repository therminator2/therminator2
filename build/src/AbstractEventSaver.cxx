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

AbstractEventSaver::AbstractEventSaver() : mFileCounter(0), kEventsPerFile(_EVENTS_PER_FILE_), mNumberOfEvents(0) {
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
  catch (TString tError) {
    PRINT_MESSAGE("<EventGenerator::ReadParameters>\tCaught exception " << tError);
    PRINT_MESSAGE("\tDid not find one of the necessary parameters in the parameters file.");
    exit(_ERROR_CONFIG_PARAMETER_NOT_FOUND_);
  }
}

