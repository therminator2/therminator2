#include "THGlobal.h"
#include "UnigenEventSaver.h"
#include "StructEvent.h"
#include "TDatime.h"

#include <list>

using namespace std;

extern TString	sEventDIR;

UnigenEventSaver::UnigenEventSaver() : mFile(0), mEventTree(0), mUEvent(new UEvent) {
}

UnigenEventSaver::~UnigenEventSaver() {
    delete mUEvent;
}

void UnigenEventSaver::Save(Event *tEvent, Model *tModel, int tEventCounter)
{
    TDatime tDate; 

    char	  tEventFName[kFileNameMaxChar];
    char	  tTempFName[kFileNameMaxChar];
    static ParticleCoor	tPartCoor;
    list<Particle>::iterator	tParticle_i;
    mUnsortedParticles.clear();

    UEvent *ev = mUEvent;
    ev->Clear();
    ev->SetEventNr(tEventCounter);
    ev->SetB(0);    // no impact parameter in Therminator
    ev->SetPhi(0);  // no event plane (at the moment) in Therminator
    ev->SetNes(0);
    mDecayIndex = 0;

    // open new file every _EVENTS_PER_FILE_ events

    Configurator *aMainConfig = tModel->GetMainConfig();

    if( !((tEventCounter - 1) % kEventsPerFile)) { 
        sprintf(tTempFName,"%sunigen%03i.tmp", sEventDIR.Data(), mFileCounter);
        mFile = new TFile(tTempFName,"RECREATE");
        mFile->cd();
        mEventTree     = new TTree(_EVENTS_TREE_,    "event tree");

        SaveParameters(tModel);

        // (void*) cast to avoid some compilation errors on older ROOT versions
        mEventTree    ->Branch(_EVENTS_BRANCH_, "UEvent", &mUEvent);   
        PRINT_DEBUG_2("<UnigenEventSaver::Save>\tCreated file "<<tTempFName);
    }
    // add all Particle entries to file and Event information
    if(mFile) {
        std::shared_ptr<DecayNode> decay_chains(new DecayNode(NULL));
        for(tParticle_i = tEvent->GetParticleList()->begin(); tParticle_i != tEvent->GetParticleList()->end(); tParticle_i++) {

            decay_chains = decay_chains->AddParticle(&(*tParticle_i));

            tPartCoor = (*tParticle_i);
            tPartCoor.SetTXYZ(
                    tPartCoor.T() * kHbarC,
                    tPartCoor.X() * kHbarC,
                    tPartCoor.Y() * kHbarC,
                    tPartCoor.Z() * kHbarC
                    );
            tPartCoor.eventid = tEvent->GetEventID();
        }
        std::vector<std::shared_ptr<DecayNode> >::const_iterator it = decay_chains->children.begin();
        for ( ; it != decay_chains->children.end(); ++it) {
            AddDecayToEvent(ev, *it, -1); // Primary particles have decay index = -1
        }

        std::sort(mUnsortedParticles.begin(), mUnsortedParticles.end(), [](const UParticle &a, const UParticle &b) {
                return a.GetIndex() < b.GetIndex();
                });
        for (auto aParticle : mUnsortedParticles) {
            ev->AddParticle(aParticle);
        }

        mEventTree->Fill();


        PRINT_DEBUG_2("<EventGenerator::SaveAsUnigen>\tEvent "<<tEventCounter<<" saved.");
    }
    // close & rename file   
    if( !(tEventCounter % kEventsPerFile) || (tEventCounter == mNumberOfEvents) ) {
        mFile->Write();
        mFile->Close();
        delete mFile;
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
    recur++;
    int status = 1;
    int decay = -1;
    int mate = -1;
    int child[2];
/*
    for (int r = 0; r < recur; ++r) {
        cout << "-";
    }
    cout << " " << node->particle->eid << " " << node->particle->pid << endl;;
*/
    if (node->children.size() > 0) {
         child[0] = node->children.front()->particle->eid;
         child[1] = node->children.back()->particle->eid;
         decay = mDecayIndex;
         ++mDecayIndex;
    } else {
        child[0] = child[1] = -1;
    }
    mUnsortedParticles.push_back(UParticle(node->particle->eid, node->particle->pid, status, node->particle->fathereid,
            parent_decay,mate-1, decay,child,
            node->particle->Px(), node->particle->Py(), node->particle->Pz(), node->particle->E(), node->particle->X(), node->particle->Y(), node->particle->Z(), node->particle->T(), 1.));
//    ev->AddParticle(node->particle->eid, node->particle->pid, status, node->particle->fathereid,
//            parent_decay,mate-1, decay,child,
//            node->particle->px, node->particle->py, node->particle->pz, node->particle->e, node->particle->x, node->particle->y, node->particle->z, node->particle->T(), 1.);
    std::vector<std::shared_ptr<DecayNode> >::const_iterator it = node->children.begin();
    for ( ; it != node->children.end(); ++it) {
        AddDecayToEvent(ev, *it, decay);
    }
    recur--;
}

