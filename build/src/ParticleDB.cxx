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

#include "ParticleDB.h"

ParticleDB::ParticleDB()
{
  mParticleTable.clear();
  mParticleNames.clear();
}

ParticleDB::~ParticleDB()
{
  mParticleTable.clear();
  mParticleNames.clear();
}

int ParticleDB::AddParticleType(ParticleType* aPartType)
{
  mParticleTable.push_back(*aPartType);
  mParticleNames[aPartType->GetName()] = mParticleTable.size() - 1;
  return  mParticleTable.size() - 1;
}

ParticleType* ParticleDB::GetParticleType(int aIndex)
{
  return &(mParticleTable[aIndex]);
}

ParticleType* ParticleDB::GetParticleType(TString aName)
{
  return &(mParticleTable[mParticleNames[aName]]);
}

int ParticleDB::GetParticleTypeIndex(TString aName)
{
  return mParticleNames[aName];
}

int ParticleDB::GetParticleTypeCount()
{
  return mParticleTable.size();
}

int ParticleDB::ExistsParticleType(TString aName)
{
  return mParticleNames.count(aName);
}
