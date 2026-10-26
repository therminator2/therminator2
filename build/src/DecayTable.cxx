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

  for (unsigned long tIter=0; tIter<mDecayChannels.size(); tIter++)
    tSumRatio += mDecayChannels[tIter].GetBranchingRatio();

  for (unsigned long tIter=0; tIter<mDecayChannels.size(); tIter++) {
    tCurRatio += mDecayChannels[tIter].GetBranchingRatio() / tSumRatio;
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

int DecayTable::ChooseDecayChannel(double aProb)
{
  unsigned long tChanIndex = 0;
  while ((mBranchingRatios[tChanIndex] < aProb) && (tChanIndex < mDecayChannels.size()))
    tChanIndex++;

  return tChanIndex;
}

int DecayTable::ChooseDecayChannelOrNot(double aProb)
{
  float tSumRatio = 0.0;

  for (unsigned long tIter=0; tIter<mDecayChannels.size(); tIter++) {
    if ((aProb > tSumRatio) && (aProb <= tSumRatio+mDecayChannels[tIter].GetBranchingRatio()))
      return tIter;
    tSumRatio += mDecayChannels[tIter].GetBranchingRatio();
  }
  return -1;
}

// Like ChooseDecayChannel(), but only among the channels with aOpen true (e.g. kinematically allowed), with their
// branching ratios renormalised. Returns -1 if no channel is open.
int DecayTable::ChooseDecayChannel(double aProb, const std::vector<bool> &aOpen)
{
  double tSumOpen = 0.0;
  for (unsigned long tIter=0; tIter<mDecayChannels.size(); tIter++)
    if (aOpen[tIter])
      tSumOpen += mDecayChannels[tIter].GetBranchingRatio();
  if (tSumOpen <= 0.0)
    return -1;

  double tCurRatio = 0.0;
  int tLastOpen = -1;
  for (unsigned long tIter=0; tIter<mDecayChannels.size(); tIter++) {
    if (!aOpen[tIter])
      continue;
    tCurRatio += mDecayChannels[tIter].GetBranchingRatio() / tSumOpen;
    tLastOpen = tIter;
    if (aProb <= tCurRatio)
      return tIter;
  }
  return tLastOpen;
}

// Like ChooseDecayChannelOrNot(), but only among the channels with aOpen true: the total probability to decay stays
// the sum of all branching ratios and is shared by the open channels in proportion to their branching ratios.
// Returns -1 for no decay, or if no channel is open.
int DecayTable::ChooseDecayChannelOrNot(double aProb, const std::vector<bool> &aOpen)
{
  double tSumAll = 0.0;
  double tSumOpen = 0.0;
  for (unsigned long tIter=0; tIter<mDecayChannels.size(); tIter++) {
    tSumAll += mDecayChannels[tIter].GetBranchingRatio();
    if (aOpen[tIter])
      tSumOpen += mDecayChannels[tIter].GetBranchingRatio();
  }
  if (tSumOpen <= 0.0)
    return -1;

  double tScale = tSumAll / tSumOpen;
  double tSumRatio = 0.0;
  for (unsigned long tIter=0; tIter<mDecayChannels.size(); tIter++) {
    if (!aOpen[tIter])
      continue;
    double tRatio = mDecayChannels[tIter].GetBranchingRatio() * tScale;
    if ((aProb > tSumRatio) && (aProb <= tSumRatio + tRatio))
      return tIter;
    tSumRatio += tRatio;
  }
  return -1;
}
