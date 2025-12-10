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
#include "CoulombAfterburner.h"
#include "AbstractEventReader.h"
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

AbstractEventReader::AbstractEventReader()
: mDB(0), mEventSaver(0), mAfterburners(0), mEvent(0),
  mEventCounter(0), mNumberOfEvents(0), mEventIter(0)
{
}

AbstractEventReader::AbstractEventReader(ParticleDB* aDB, AbstractEventSaver *aES, ListAfterburner *aLA)
: mDB(aDB), mEventSaver(aES), mAfterburners(aLA),
  mEventCounter(0), mNumberOfEvents(0), mEventIter(0)
{
  ReadParameters();
  mEvent = new Event();
  mEventSaver->FindPreviousEventFiles();
}

AbstractEventReader::~AbstractEventReader()
{
  delete mEvent;
}

void AbstractEventReader::ReadEvents()
{
  TDatime tDate;
  
  tDate.Set();
  PRINT_MESSAGE("["<<tDate.AsSQLString()<<"]\tReading events");

  TString tEventDir = sMainConfig->GetParameter("EventDir");
  cout << "EventDir " << tEventDir.Data() << " EventSubDir " << mEventSubDir << endl;

  Int_t tEvents;
  SetupInput(tEventDir + mEventSubDir, 1);

  int tAllParticleIter = 0;
  while (NextEvent()) {

    mEventCounter = mEventIter+1;
    mEvent->Reset(mEventIter);

    for (unsigned long tParticleIter = 0; tParticleIter < mStructEvent.entries; ++tParticleIter, ++tAllParticleIter) {
       // std::cerr << "tParticleIter " << tParticleIter << " of " << mStructEvent.entries << ", ";
        GetParticle(tAllParticleIter);
        if (mParticleCoor.eventid != mStructEvent.eventID) {
            PRINT_MESSAGE("<AbstractEventReader::ReadEvents>\tError");
            PRINT_MESSAGE("\tEvent ID for particle doesn't match.");
            exit(-1);        
        }
        ParticleType *tParticleType = mDB->FindByPID(mParticleCoor.pid);
        Particle tParticle(&mParticleCoor, tParticleType);
        mEvent->GetParticleList()->push_back(tParticle);
    }

    Float_t dummy;
    ParticleType *tSpectatorType = CoulombAfterburner::SpectatorType(mEvent->GetParticleList(), dummy);
    for (auto &tParticle : *mEvent->GetParticleList()) {
        if (tParticle.GetParticleType() == nullptr) {
            tParticle.SetParticleType(tSpectatorType);
//            Info("AbstractEventReader::ReadEvents", "No particle type found in DB, using spectator type with mass %f", tSpectatorType->GetMass());
        }
    } 

    cout << "\r\tevent " << mEventCounter <<"/"<< mNumberOfEvents;
    cout.flush();
    mAfterburners->Apply(mEvent);
    // TODO: FInd better solution to pass here a correct Model* parameter
    mEventSaver->Save(mEvent,NULL,mEventCounter);
  }
}

void AbstractEventReader::ReadParameters()
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
  catch (TString &tError) {
    PRINT_MESSAGE("<AbstractEventReader::ReadParameters>\tCaught exception " << tError);
    PRINT_MESSAGE("\tDid not find one of the necessary parameters in the parameters file.");
    exit(_ERROR_CONFIG_PARAMETER_NOT_FOUND_);
  }
  delete tModelParam;

}
