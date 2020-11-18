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

#include <iostream>
#include <cstdio>
#include <TChain.h>
#include "./../build/include/ParticleCoor.h"
#include "./../build/include/StructEvent.h"

using namespace std;

TChain* events2chain(TString aEventDir, Int_t aEventFiles, ParticleCoor* aParticle, Int_t* aEvents)
{
// ##########################################################################
// WHICH ROOT FILES TO READ
// ##########################################################################
  static StructEvent tStructEvents;
  TChain* tChainParts = new TChain(_PARTICLES_TREE_);
  TChain* tChainEvent = new TChain(_EVENTS_TREE_);

  tChainParts->SetBranchAddress(_PARTICLE_BRANCH_, aParticle);
  tChainEvent->SetBranchAddress(_EVENTS_BRANCH_,   &tStructEvents);
  
  for(Int_t i=0; i<aEventFiles; i++) {
    char Buff[100];
    sprintf(Buff,"%sevent%03i.root",aEventDir.Data(),i);
    cerr << "Adding file: " << Buff << endl;
    tChainParts->Add(Buff);
    tChainEvent->Add(Buff);
  }
  (*aEvents) = tChainEvent->GetEntries();

  cerr << "Total number of events: " << (*aEvents) << endl;
  return tChainParts;
}

/*! @file events2chain.C
 * @brief ROOT macro with events2chain() function.
 */
/*! @fn TChain* events2chain(TString aEventDir, Int_t aEventFiles, ParticleCoor* aParticle, Int_t* aEvents)
 * @brief ROOT macro. Function adds a number of <b>eventXXX.root</b> files, links the given Particle_t object to the TChain.
 * @param [in] aEventDir directory with event files
 * @param [in] aEventFiles number of event files to add to TChain
 * @param [out] aParticle pointer to ParticleCoor object  
 * @param [out] aEvents total number of events added
 * @retval TChain* pointer to TChain
 */
