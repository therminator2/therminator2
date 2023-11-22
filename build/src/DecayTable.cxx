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

#include "DecayTable.h"
#include <iostream>
using namespace std;

DecayTable::DecayTable()
{
  mDecayChannels.clear();
  mBranchingRatios.clear();
}

DecayTable::DecayTable(const DecayTable& aTable)
{
  mDecayChannels.clear();
  mBranchingRatios.clear();
  for (int tIter=0; tIter<aTable.GetChannelCount(); tIter++)
    AddDecayChannel( *(aTable.GetDecayChannel(tIter)) );
}

DecayTable::~DecayTable()
{
  mDecayChannels.clear();
  mBranchingRatios.clear();
}

void DecayTable::AddDecayChannel(DecayChannel aChannel)
{
  mDecayChannels.push_back(aChannel);
  RecalculateBranchingRatios();
}

void DecayTable::RecalculateBranchingRatios()
{
  float tSumRatio = 0.0;
  float tCurRatio = 0.0;

  for (int tIter=0; tIter<mDecayChannels.size(); tIter++)
    tSumRatio += mDecayChannels[tIter].GetBranchingRatio();

  for (int tIter=0; tIter<mDecayChannels.size(); tIter++) {
    tCurRatio += mDecayChannels[tIter].GetBranchingRatio() / tSumRatio;
    if (mBranchingRatios.size() <= tIter)
      mBranchingRatios.push_back(tCurRatio);
    else
      mBranchingRatios[tIter] = tCurRatio;
  }
}

void DecayTable::RecalculateBranchingRatiosForMass(Float_t aM) {
  for (int tIter=0; tIter<mDecayChannels.size(); tIter++) 
    mBranchingRatios[tIter] = mDecayChannels[tIter].GetBranchingRatioFunc()->Interpolate(aM);


  float tSumRatio = 0.0;
  float tCurRatio = 0.0;

  for (int tIter=0; tIter<mDecayChannels.size(); tIter++)
    tSumRatio += mBranchingRatios[tIter];

  for (int tIter=0; tIter<mDecayChannels.size(); tIter++) {
    tCurRatio += mBranchingRatios[tIter] / tSumRatio;
    if (mBranchingRatios.size() <= tIter)
      mBranchingRatios.push_back(tCurRatio);
    else
      mBranchingRatios[tIter] = tCurRatio;
  }
    

}

int DecayTable::GetChannelCount() const
{
  return mDecayChannels.size() - 1;
}


const DecayChannel* DecayTable::GetDecayChannel(int aIndex) const
{
  return &(mDecayChannels[aIndex]);
}

float DecayTable::GetDecayStep(int aIndex)
{
  return mBranchingRatios[aIndex];
}

int DecayTable::ChooseDecayChannel(double aProb,Float_t aM)
{
  std::vector<float> branchingRatiosCopy;
  if(mDecayChannels[0].IsBRFunc()) {
    branchingRatiosCopy = mBranchingRatios;
    RecalculateBranchingRatiosForMass(aM);
  }

  int tChanIndex = 0;
  while ((mBranchingRatios[tChanIndex] < aProb) && (tChanIndex < mDecayChannels.size()))
    tChanIndex++;

  if(mDecayChannels[0].IsBRFunc())
    mBranchingRatios = branchingRatiosCopy;
  

  return tChanIndex;
}

int DecayTable::ChooseDecayChannelOrNot(double aProb,Float_t aM)
{

  std::vector<float> branchingRatiosCopy;
  if(mDecayChannels[0].IsBRFunc()) {
    branchingRatiosCopy = mBranchingRatios;
    RecalculateBranchingRatiosForMass(aM);
  }

  float tSumRatio = 0.0;

  for (int tIter=0; tIter<mDecayChannels.size(); tIter++) {
    if ((aProb > tSumRatio) && (aProb <= tSumRatio+mDecayChannels[tIter].GetBranchingRatio())){
      if(mDecayChannels[0].IsBRFunc()) {
        mBranchingRatios = branchingRatiosCopy;
      return tIter;
      }
    }
    tSumRatio += mDecayChannels[tIter].GetBranchingRatio();
  }

  if(mDecayChannels[0].IsBRFunc())
    mBranchingRatios = branchingRatiosCopy;

  return -1;
}
