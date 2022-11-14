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

#include <sstream>
#include <TDatime.h>
#include "Crc32.h"
#include "Configurator.h"
#include "ParticleDecayer.h"
#include "Event.h"
#include "THGlobal.h"

extern Configurator* sMainConfig;
extern TString	sTimeStamp;

using namespace std;

Event::Event()
{
  Reset();
}

Event::~Event()
{
  mParticles.clear();
}

void Event::Reset(int aEventIter)
{
  ostringstream oss;
  Crc32 tEventID;
  
  mParticles.clear();
  Particle::ZeroEID();
  
  oss << sTimeStamp.Data() << "Event: " << aEventIter;
  tEventID.Update(oss.str().data(), oss.str().length());
  tEventID.Finish(); 
  mEventID = tEventID.GetValue();
}

list<Particle>* Event::GetParticleList()
{
  return &mParticles;
}

unsigned int Event::GetEventID() const
{
  return mEventID;
}

