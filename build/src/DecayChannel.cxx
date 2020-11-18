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

#include "DecayChannel.h"

DecayChannel::DecayChannel()
: mParticleType1(0), mParticleType2(0), mParticleType3(-1), mBranchRatio(0.0)
{
}

DecayChannel::DecayChannel(const DecayChannel& aChannel) 
{
  mBranchRatio   = aChannel.GetBranchingRatio();
  mParticleType1 = aChannel.GetParticle1();
  mParticleType2 = aChannel.GetParticle2();
  mParticleType3 = aChannel.GetParticle3();
}

DecayChannel::DecayChannel(double aBranchRatio, int aPartType1, int aPartType2, int aPartType3)
: mParticleType1(aPartType1), mParticleType2(aPartType2), mParticleType3(aPartType3), mBranchRatio(aBranchRatio)
{
}

DecayChannel::~DecayChannel()
{
}

int DecayChannel::GetParticle1() const
{
  return mParticleType1;
}

int DecayChannel::GetParticle2() const
{
  return mParticleType2;
}

int DecayChannel::GetParticle3() const
{
  return mParticleType3;
}

double DecayChannel::GetBranchingRatio() const
{
  return mBranchRatio;
}

int DecayChannel::Is3Particle() const
{
  return (mParticleType3 != -1);
}

void DecayChannel::SetParticle1(int aPartType1)
{
  mParticleType1 = aPartType1;
}

void DecayChannel::SetParticle2(int aPartType2)
{
  mParticleType2 = aPartType2;
}

void DecayChannel::SetParticle3(int aPartType3)
{
  mParticleType3 = aPartType3;
}

void DecayChannel::SetBranchingRatio(double aRatio)
{
  mBranchRatio = aRatio;
}
