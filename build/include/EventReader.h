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

#ifndef _TH2_EVENT_READER_H_
  #define _TH2_EVENT_READER_H_
  
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TString.h>
#include "ParticleDB.h"
#include "Integrator.h"
#include "AbstractEventSaver.h"
#include "ListAfterburner.h"
#include "Event.h"
#include "StructEvent.h"

class EventReader {
  public:
    EventReader();
    EventReader(ParticleDB* aDB, AbstractEventSaver *aES, ListAfterburner *aLA);
    ~EventReader();

    void   ReadEvents();

  private:
    void ReadParameters();  
    void SaveEvent();
    std::pair<TChain*, TChain*> GetChains(TString aEventDir, Int_t aEventFiles, ParticleCoor* aParticle, StructEvent* aStructEvent, Int_t* aEvents);

    ParticleDB*	        mDB;
    AbstractEventSaver* mEventSaver;
    ListAfterburner*    mAfterburners;
    Event*      	mEvent;
    int		        mEventCounter;
    int		        mNumberOfEvents;
    TString       mEventSubDir;
};

#endif
