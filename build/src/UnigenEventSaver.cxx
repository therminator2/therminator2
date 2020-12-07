#include "THGlobal.h"
#include "UnigenEventSaver.h"
#include "StructEvent.h"
#include "TDatime.h"

#include <list>

using namespace std;

extern TString	sEventDIR;
extern TString	sTimeStamp;
extern int	sIntegrateSample;
extern int	sModel;
extern int	sRandomize;

UnigenEventSaver::UnigenEventSaver() : mFile(0), mEventTree(0), mParameterTree(0) {
}

void UnigenEventSaver::Save(Event *tEvent, int tEventCounter)
{
    TDatime tDate; 

    char	  tEventFName[kFileNameMaxChar];
    char	  tTempFName[kFileNameMaxChar];
    static ParticleCoor	tPartCoor;
    static StructEvent	tStructEvent;
    list<Particle>::iterator	tParticle_i;

    UEvent *ev = new UEvent;
    ev->Clear();
    ev->SetEventNr(tEventCounter);
    ev->SetB(0);    // no impact parameter in Therminator
    ev->SetPhi(0);  // no event plane (at the moment) in Therminator
    ev->SetNes(0);
    mDecayIndex = 0;

    // open new file every _EVENTS_PER_FILE_ events
    if( !((tEventCounter - 1) % kEventsPerFile)) { 
        sprintf(tTempFName,"%sunigen%03i.tmp", sEventDIR.Data(), mFileCounter);
        mFile = new TFile(tTempFName,"RECREATE");
        mFile->cd();
        mEventTree     = new TTree(_EVENTS_TREE_,    "event tree");
        mParameterTree = new TTree(_PARAMETERS_TREE_,"parameters and model description tree");

        // (void*) cast to avoid some compilation errors on older ROOT versions
        char tTimeStamp[21];
        sprintf(tTimeStamp,"%s",sTimeStamp.Data());
        mEventTree    ->Branch(_EVENTS_BRANCH_,		"UEvent",                                                        &ev                    );   
        mParameterTree->Branch(_INTEGRATESAMPLE_BRANCH_,	(UInt_t*) &sIntegrateSample,					 "i"			);
        mParameterTree->Branch(_RANDOMIZE_BRANCH_,		(UInt_t*) &sRandomize,						 "i"			);
        mParameterTree->Branch(_TIMESTAMP_BRANCH_,		(Char_t*) tTimeStamp,						 _TIMESTAMP_FORMAT_	);
        mParameterTree->Branch(_MODELID_BRANCH_,		(UInt_t*) &sModel,						 "i"			);
        mParameterTree->Branch(_MODELNAME_BRANCH_,		(Char_t*) tEvent->GetIntegrator()->GetModel()->GetName(),	 _MODELNAME_FORMAT_	);
        mParameterTree->Branch(_MODELHASH_BRANCH_,		(Char_t*) tEvent->GetIntegrator()->GetModel()->GetHash(),	 _MODELHASH_FORMAT_	);
        mParameterTree->Branch(_MODELDESCRIPTION_BRANCH_,	(Char_t*) tEvent->GetIntegrator()->GetModel()->GetDescription(), _MODELDESCRIPTION_FORMAT_);    
        tEvent        ->GetIntegrator()->GetModel()->AddParameterBranch(mParameterTree);
        mParameterTree->Fill();    
        PRINT_DEBUG_2("<UnigenEventSaver::Save>\tCreated file "<<tTempFName);
    }
    // add all Particle entries to file and Event information
    if(mFile) {
        int i = 0;

        std::shared_ptr<DecayNode> decay_chains(new DecayNode(NULL));
        for(tParticle_i = tEvent->GetParticleList()->begin(); tParticle_i != tEvent->GetParticleList()->end(); tParticle_i++) {

            decay_chains = decay_chains->AddParticle(&(*tParticle_i));

            tPartCoor = (*tParticle_i);
            tPartCoor.t *= kHbarC;
            tPartCoor.x *= kHbarC;
            tPartCoor.y *= kHbarC;
            tPartCoor.z *= kHbarC;
            tPartCoor.eventid = tEvent->GetEventID();
        }
        std::vector<std::shared_ptr<DecayNode> >::const_iterator it = decay_chains->children.begin();
        for ( ; it != decay_chains->children.end(); ++it) {
            AddDecayToEvent(ev, *it, -1); // Primary particles have decay index = -1
        }

        tStructEvent.eventID     = tEvent->GetEventID();
        tStructEvent.entries     = tEvent->GetParticleList()->size();
        tStructEvent.entriesprev = 0; // not used here
        mEventTree->Fill();
        PRINT_DEBUG_2("<EventGenerator::SaveAsUnigen>\tEvent "<<tEventCounter<<" saved.");
    }
    // close & rename file   
    if( !(tEventCounter % kEventsPerFile) || (tEventCounter == mNumberOfEvents) ) {
        mFile->Write();
        mFile->Close();
        mFile = 0;
        // rename temporary file to ROOT
        sprintf(tEventFName,"%sunigen%03i.root",sEventDIR.Data(),mFileCounter);
        sprintf(tTempFName, "%sunigen%03i.tmp", sEventDIR.Data(),mFileCounter);
        if(rename(tTempFName,tEventFName) != 0) {
            PRINT_MESSAGE("<EventGenerator::SaveAsUnigen>\tUnable to rename temp file to " << tTempFName);
            exit(_ERROR_GENERAL_FILE_NOT_FOUND_);
        } else {
            tDate.Set();
            PRINT_MESSAGE("\n["<<tDate.AsSQLString() << "]\tFile "<<tEventFName<<" written");
        }
        mFileCounter++;
    }
}

void UnigenEventSaver::AddDecayToEvent(UEvent *ev, std::shared_ptr<DecayNode> node, int parent_decay) {
    int status = 1;
    int decay = -1;
    int mate = -1;
    int child[2];
    if (node->children.size() > 0) {
         child[0] = node->children.front()->particle->eid;
         child[1] = node->children.back()->particle->eid;
         decay = mDecayIndex;
         ++mDecayIndex;
    } else {
        child[0] = child[1] = -1;
    }
    ev->AddParticle(node->particle->eid, node->particle->pid, status, node->particle->fathereid,
            parent_decay,mate-1, decay,child,
            node->particle->px, node->particle->py, node->particle->pz, node->particle->e, node->particle->x, node->particle->y, node->particle->z, node->particle->t, 1.);
    std::vector<std::shared_ptr<DecayNode> >::const_iterator it = node->children.begin();
    for ( ; it != node->children.end(); ++it) {
        AddDecayToEvent(ev, *it, decay);
    }
}

