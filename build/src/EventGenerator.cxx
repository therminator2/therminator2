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
#include "EventGenerator.h"
#include "AbstractEventSaver.h"
#include "Event.h"
#include "Configurator.h"
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

EventGenerator::EventGenerator()
: mDB(0), mEventSaver(0), mInteg(0), mEvent(0),
  mEventCounter(0), mNumberOfEvents(0)
{
}

EventGenerator::EventGenerator(ParticleDB* aDB, AbstractEventSaver *aES)
: mDB(aDB), mEventSaver(aES),
  mEventCounter(0), mNumberOfEvents(0)
{
  ReadParameters();
  mInteg = new Integrator(sIntegrateSample);
  mInteg->SetMultiplicities(mDB);
  mEvent = new Event(mDB, mInteg);
  if (sRandomize) {
    mInteg->Randomize();
    mEvent->Randomize();
  }
  mEventSaver->FindPreviousEventFiles();
}

EventGenerator::~EventGenerator()
{
  delete mInteg;
  delete mEvent;
}

void EventGenerator::GenerateEvents()
{
  TDatime tDate;
  
  tDate.Set();
  PRINT_MESSAGE("["<<tDate.AsSQLString()<<"]\tGenerating events");
  for (int tIter=0; tIter<mNumberOfEvents; tIter++) {
    mEventCounter = tIter+1;
    mEvent->Reset(tIter);
    if (sRandomize) {
      mEvent->GeneratePrimordials();
      mEvent->DecayParticles();
    } else {
      mEvent->GeneratePrimordials (43212 - tIter * 2);
      mEvent->DecayParticles	  (43212 - tIter * 2);
    }
    cout << "\r\tevent " << tIter+1 <<"/"<< mNumberOfEvents;
    cout.flush();
    mEventSaver->Save(mEvent,mEventCounter);
  }
}

void EventGenerator::ReadParameters()
{
  try {
    mNumberOfEvents	= (sMainConfig->GetParameter("NumberOfEvents")).Atoi();
    sIntegrateSample	= (sMainConfig->GetParameter("IntegrateSamples")).Atoi();
  }
  catch (TString tError) {
    PRINT_MESSAGE("<EventGenerator::ReadParameters>\tCaught exception " << tError);
    PRINT_MESSAGE("\tDid not find one of the necessary parameters in the parameters file.");
    exit(_ERROR_CONFIG_PARAMETER_NOT_FOUND_);
  }
}
