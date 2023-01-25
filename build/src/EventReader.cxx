/********************************************************************************
 *                                                                              *
 *             THERMINATOR 2: THERMal heavy-IoN generATOR 2                     *
 *                                                                              *
 * Version:                                                                     *
 *      Release, 2.0.3, 1 February 2011                                         *
 *                                                                              *
 * Authors:                                                                     *
 *      Mikolaj Chojnacki   (Mikolaj.Chojnacki@ifj.edu.pl)                      *
 *      Adam Kisiel         (kisiel@if.pw.edu.pl)                               *
 *      Wojciech Broniowski (Wojciech.Broniowski@ifj.edu.pl)                    *
 *      Wojciech Florkowski (Wojciech.Florkowski@ifj.edu.pl)                    *
 *                                                                              *
 * Project homepage:                                                            *
 *      http://therminator2.ifj.edu.pl/                                         *
 *                                                                              *
 * For the detailed description of the program and further references           *
 * to the description of the model please refer to                              *
 * http://arxiv.org/abs/1102.0273                                               *
 *                                                                              *
 * This code can be freely used and redistributed. However if you decide to     *
 * make modifications to the code, please, inform the authors.                  *
 * Any publication of results obtained using this code must include the         *
 * reference to arXiv:1102.0273 and the published version of it, when           *
 * available.                                                                   *
 *                                                                              *
 ********************************************************************************/

#include <fstream>
#include <TDatime.h>
#include "EventReader.h"
#include "AbstractEventSaver.h"
#include "Event.h"
#include "Configurator.h"
#include "Parser.h"
#include "THGlobal.h"

extern void AddLogEntry(const char* aEntry);

extern Configurator *sMainConfig;
extern TString	sMainINI;
extern TString	sModelINI;
extern int	sModel;
extern int	sRandomize;
extern int	sIntegrateSample;
extern int	sParentPID;

using namespace std;

EventReader::EventReader()
: mDB(0), mEventSaver(0), mPlugins(0), mEvent(0),
  mEventCounter(0), mNumberOfEvents(0)
{
}

EventReader::EventReader(ParticleDB* aDB, AbstractEventSaver *aES, ListPlugin *aLA)
: mDB(aDB), mEventSaver(aES), mPlugins(aLA),
  mEventCounter(0), mNumberOfEvents(0)
{
  ReadParameters();
  mEvent = new Event();
  mEventSaver->FindPreviousEventFiles();
}

EventReader::~EventReader()
{
  delete mEvent;
}

void EventReader::ReadEvents()
{
  TDatime tDate;
  
  tDate.Set();
  PRINT_MESSAGE("["<<tDate.AsSQLString()<<"]\tReading events");

  ParticleCoor tParticleCoor;
  StructEvent tStructEvent;

  TString tEventDir = sMainConfig->GetParameter("EventDir");
  cout << "EventDir " << tEventDir.Data() << " EventSubDir " << mEventSubDir << endl;

  Int_t tEvents;
  std::pair<TChain*, TChain*> tChains = GetChains(tEventDir + mEventSubDir, 1, &tParticleCoor, &tStructEvent, &tEvents);
  TChain *tChainParts = tChains.first;
  TChain *tChainEvent = tChains.second;

  int tAllParticleIter = 0;
  for (int tEventIter=0; tEventIter<TMath::Min(mNumberOfEvents, tEvents); tEventIter++) {

    tChainEvent->GetEntry(tEventIter);

    mEventCounter = tEventIter+1;
    mEvent->Reset(tEventIter);

      for (int tParticleIter = 0; tParticleIter < tStructEvent.entries; ++tParticleIter, ++tAllParticleIter) {
        tChainParts->GetEntry(tAllParticleIter);
        if (tParticleCoor.eventid != tStructEvent.eventID) {
          PRINT_MESSAGE("<EventReader::ReadEvents>\tError");
          PRINT_MESSAGE("\tEvent ID for particle doesn't match.");
          exit(-1);        
        }
        ParticleType *tParticleType = mDB->FindByPID(tParticleCoor.pid);
        Particle tParticle(&tParticleCoor, tParticleType);
        mEvent->GetParticleList()->push_back(tParticle);
      }

    cout << "\r\tevent " << tEventIter+1 <<"/"<< mNumberOfEvents;
    cout.flush();
    mPlugins->Apply(mEvent);
    // TODO: FInd better solution to pass here a correct Model* parameter
    mEventSaver->Save(mEvent,NULL,mEventCounter);
  }
}

std::pair<TChain*, TChain*> EventReader::GetChains(TString aEventDir, Int_t aEventFiles, ParticleCoor* aParticle, StructEvent *aStructEvent, Int_t* aEvents)
{
// ##########################################################################
// WHICH ROOT FILES TO READ
// ##########################################################################
  TChain* tChainParts = new TChain(_PARTICLES_TREE_);
  TChain* tChainEvent = new TChain(_EVENTS_TREE_);

  tChainParts->SetBranchAddress(_PARTICLE_BRANCH_, aParticle);
  tChainEvent->SetBranchAddress(_EVENTS_BRANCH_,   aStructEvent);
  
  for(Int_t i=0; i<aEventFiles; i++) {
    char Buff[100];
    sprintf(Buff,"%sevent%03i.root",aEventDir.Data(),i);
    cerr << "Adding file: " << Buff << endl;
    tChainParts->Add(Buff);
    tChainEvent->Add(Buff);
  }
  (*aEvents) = tChainEvent->GetEntries();

  cerr << "Total number of events: " << (*aEvents) << endl;
  return std::make_pair(tChainParts, tChainEvent);
}

void EventReader::ReadParameters()
{
  Configurator* tModelParam;
  Parser*       tParser;

  tModelParam = new Configurator;
  tParser     = new Parser(sModelINI.Data());                                                                                                                                                                                   
  tParser->ReadINI(tModelParam);
  delete tParser;

  try {
    mNumberOfEvents	= (sMainConfig->GetParameter("NumberOfEvents")).Atoi();
    sIntegrateSample	= (sMainConfig->GetParameter("IntegrateSamples")).Atoi();

    mEventSubDir = tModelParam->GetParameter("EventSubDir");
  }
  catch (TString tError) {
    PRINT_MESSAGE("<EventReader::ReadParameters>\tCaught exception " << tError);
    PRINT_MESSAGE("\tDid not find one of the necessary parameters in the parameters file.");
    exit(_ERROR_CONFIG_PARAMETER_NOT_FOUND_);
  }
  delete tModelParam;

}
