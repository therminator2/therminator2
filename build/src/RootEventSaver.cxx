#include "THGlobal.h"
#include "RootEventSaver.h"
#include "StructEvent.h"
#include "TDatime.h"

#include <list>

using namespace std;

extern TString	sEventDIR;
extern TString	sTimeStamp;
extern int	sIntegrateSample;
extern int	sModel;
extern int	sRandomize;

RootEventSaver::RootEventSaver() : mFile(0), mEventTree(0), mParameterTree(0), mParticleTree(0) {
}

void RootEventSaver::Save(Event *tEvent, Model *tModel, int tEventCounter)
{
  TDatime tDate; 
  
  char	  tEventFName[kFileNameMaxChar];
  char	  tTempFName[kFileNameMaxChar];
  static ParticleCoor	tPartCoor;
  static StructEvent	tStructEvent;
  list<Particle>::iterator	tParticle_i;
  
// open new file every _EVENTS_PER_FILE_ events
  if( !((tEventCounter - 1) % kEventsPerFile)) { 
    sprintf(tTempFName,"%sevent%03i.tmp", sEventDIR.Data(), mFileCounter);
    mFile = new TFile(tTempFName,"RECREATE");
    mFile->cd();
    mParticleTree  = new TTree(_PARTICLES_TREE_, "particle tree");    
    mEventTree     = new TTree(_EVENTS_TREE_,    "event tree");
    mParameterTree = new TTree(_PARAMETERS_TREE_,"parameters and model description tree");
    // (void*) cast to avoid some compilation errors on older ROOT versions
    char tTimeStamp[21];
    sprintf(tTimeStamp,"%s",sTimeStamp.Data());
    mParticleTree ->Branch(_PARTICLE_BRANCH_,		&tPartCoor,							 _PARTICLE_FORMAT_	);
    mEventTree    ->Branch(_EVENTS_BRANCH_,		&tStructEvent,							 _EVENTS_FORMAT_	);   
    mParameterTree->Branch(_INTEGRATESAMPLE_BRANCH_,	(UInt_t*) &sIntegrateSample,					 "i"			);
    mParameterTree->Branch(_RANDOMIZE_BRANCH_,		(UInt_t*) &sRandomize,						 "i"			);
    mParameterTree->Branch(_TIMESTAMP_BRANCH_,		(Char_t*) tTimeStamp,						 _TIMESTAMP_FORMAT_	);
    mParameterTree->Branch(_MODELID_BRANCH_,		(UInt_t*) &sModel,						 "i"			);
    if (tModel != NULL) {
      mParameterTree->Branch(_MODELNAME_BRANCH_,		(Char_t*) tModel->GetName(),	 _MODELNAME_FORMAT_	);
      mParameterTree->Branch(_MODELHASH_BRANCH_,		(Char_t*) tModel->GetHash(),	 _MODELHASH_FORMAT_	);
      mParameterTree->Branch(_MODELDESCRIPTION_BRANCH_,	(Char_t*) tModel->GetDescription(), _MODELDESCRIPTION_FORMAT_);    
      tModel->AddParameterBranch(mParameterTree);
    }
    mParameterTree->Fill();    
    PRINT_DEBUG_2("<EventGenerator::SaveAsRoot>\tCreated file "<<tTempFName);
  }
// add all Particle entries to file and Event information
  if(mFile) {
    for(tParticle_i = tEvent->GetParticleList()->begin(); tParticle_i != tEvent->GetParticleList()->end(); tParticle_i++) {
      tPartCoor = (*tParticle_i);
      tPartCoor.t *= kHbarC;
      tPartCoor.x *= kHbarC;
      tPartCoor.y *= kHbarC;
      tPartCoor.z *= kHbarC;
      tPartCoor.eventid = tEvent->GetEventID();
      mParticleTree->Fill();
    }
    tStructEvent.eventID     = tEvent->GetEventID();
    tStructEvent.entries     = tEvent->GetParticleList()->size();
    tStructEvent.entriesprev = 0; // not used here
    mEventTree->Fill();
    PRINT_DEBUG_2("<EventGenerator::SaveAsRoot>\tEvent "<<tEventCounter<<" saved.");
  }
// close & rename file   
  if( !(tEventCounter % kEventsPerFile) || (tEventCounter == mNumberOfEvents) ) {
    mFile->Write();
    mFile->Close();
    mFile = 0;
    // rename temporary file to ROOT
    sprintf(tEventFName,"%sevent%03i.root",sEventDIR.Data(),mFileCounter);
    sprintf(tTempFName, "%sevent%03i.tmp", sEventDIR.Data(),mFileCounter);
    if(rename(tTempFName,tEventFName) != 0) {
      PRINT_MESSAGE("<EventGenerator::SaveAsRoot>\tUnable to rename temp file to " << tTempFName);
      exit(_ERROR_GENERAL_FILE_NOT_FOUND_);
    } else {
      tDate.Set();
      PRINT_MESSAGE("\n["<<tDate.AsSQLString() << "]\tFile "<<tEventFName<<" written");
    }
    mFileCounter++;
  }
}


