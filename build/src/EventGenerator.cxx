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
#include "ParticleDecayer.h"
#include "Event.h"
#include "Configurator.h"
#include "THGlobal.h"

extern void AddLogEntry(const char* aEntry);

extern Configurator *sMainConfig;
extern TString	sMainINI;
extern TString	sModelINI;

extern int	sSeed;
extern int	sModel;
extern int	sRandomize;
extern int	sIntegrateSample;
extern int	sParentPID;

using namespace std;

EventGenerator::EventGenerator()
: mDB(0), mEventSaver(0), mAfterburners(0), mInteg(0), mEvent(0),
  mEventCounter(0), mNumberOfEvents(0), mDistribution(0)
{
  mMultiplicities.clear();
}

EventGenerator::EventGenerator(ParticleDB* aDB, AbstractEventSaver *aES, ListAfterburner *aLA)
: mDB(aDB), mEventSaver(aES), mAfterburners(aLA),
  mEventCounter(0), mNumberOfEvents(0), mDistribution(0)
{
  mMultiplicities.clear();
  mMultiplicities.resize(mDB->GetParticleTypeCount());
  ReadParameters();
  mInteg = new Integrator(sIntegrateSample);
  mInteg->SetMultiplicities(mDB);
  mEvent = new Event();
  if (sRandomize) {
    mInteg->Randomize();
    Randomize();
  }
  mRandom = new TRandom2();
#ifdef _ROOT_4_
  mRandom->SetSeed2(31851+sSeed, 14327);
#else
  mRandom->SetSeed(31851+sSeed);
#endif
  mEventSaver->FindPreviousEventFiles();
}

EventGenerator::~EventGenerator()
{
  delete mInteg;
  delete mEvent;
  delete mRandom;
  mMultiplicities.clear();
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
      GeneratePrimordials();
      DecayParticles();
    } else {
      GeneratePrimordials (43212 - tIter * 2);
      DecayParticles	  (43212 - tIter * 2);
    }
    cout << "\r\tevent " << tIter+1 <<"/"<< mNumberOfEvents;
    cout.flush();
    mAfterburners->Apply(mEvent);
    mEventSaver->Save(mEvent,mInteg->GetModel(),mEventCounter);
  }
}

void EventGenerator::GeneratePrimordials(int aSeed)
{ 
#ifdef _ROOT_4_
  if (aSeed) mRandom->SetSeed2(aSeed, (aSeed*2) % (7*11*23*31));
#else
  if (aSeed) mRandom->SetSeed(aSeed);
#endif

  GenerateMultiplicities();
  for (int tIter=0; tIter<mDB->GetParticleTypeCount(); tIter++) {
    if(! strstr(mDB->GetParticleType(tIter)->GetName(),"gam000zer")) { // disable primordial photons production
      mInteg->GenerateParticles(mDB->GetParticleType(tIter), mMultiplicities[tIter], mEvent->GetParticleList());
    }
  }
}

void EventGenerator::GenerateMultiplicities()
{
  if(mDistribution == 0) { // Poisson
    for (int tIter=0; tIter<mDB->GetParticleTypeCount(); tIter++) {
      mMultiplicities[tIter] = mRandom->Poisson(mDB->GetParticleType(tIter)->GetMultiplicity());
    }
  } else if(mDistribution == 1) { // Negative Binomial
    for (int tIter=0; tIter<mDB->GetParticleTypeCount(); tIter++)
      mMultiplicities[tIter] = 0; // HOW?
  }
}

void EventGenerator::DecayParticles(int aSeed)
{
  list<Particle>::iterator tIter;
  ParticleType*    tFatherType;
  ParticleDecayer* tDecayer;
  
  tDecayer = new ParticleDecayer(mDB, mEvent->GetParticleList(),mInteg->GetModel());

  if (sRandomize)
    tDecayer->Randomize();
  else
    tDecayer->SeedSet(aSeed);
 

  tIter = mEvent->GetParticleList()->begin();
// as new particles are added from decays the end() of the list moves until all particles had decayed
  while (tIter != mEvent->GetParticleList()->end()) {
    tFatherType = tIter->GetParticleType();
    // if not stable or stable but has a decay table with at least one decay channel
    if((tFatherType->GetGamma() >= 0.0) && (tFatherType->GetTable()) && (tFatherType->GetTable()->GetChannelCount() + 1 > 0))
      tDecayer->DecayParticle( &(*tIter) );
    tIter++;
  }
  delete tDecayer;
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

  TString tDistribution; 
  try {
    tDistribution	= sMainConfig->GetParameter("MultiplicityDistribution");
    if (tDistribution.Contains("NegativeBinomial"))
      mDistribution = 1;
  }
  catch (TString tError) {
    PRINT_DEBUG_1("<Event::ReadParameters>\tUsing default multiplicity distribution: Poissonian");
  }
}

void EventGenerator::Randomize()
{
  TDatime tDate;

#ifdef _ROOT_4_
  mRandom->SetSeed2(tDate.Get() / 2 * 3, tDate.Get() / 11 * 9);
#else
  mRandom->SetSeed(tDate.Get() / 2 * 3);
#endif
}
