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

#ifndef _TH2_DECAY_TABLE_H_
  #define _TH2_DECAY_TABLE_H_

#include <vector>
#include "DecayChannel.h"

class DecayTable {
  public:
    DecayTable();
    DecayTable(const DecayTable& aTable);
    ~DecayTable();

    void		AddDecayChannel(DecayChannel aChannel);
    int			GetChannelCount() const;
    const DecayChannel*	GetDecayChannel(int aIndex) const;
    float		GetDecayStep(int aIndex);
    int			ChooseDecayChannel(double aProb);
    int			ChooseDecayChannelOrNot(double aProb);

  private:  
    void RecalculateBranchingRatios();

    std::vector<DecayChannel>	mDecayChannels;
    std::vector<float>		mBranchingRatios;
};

#endif

/*! @file DecayTable.h
 * @brief Definition of DecayTable class. Stores available DecayChannel list.
 */
/*! @class DecayTable
 * @brief Stores DecayChannel objects in a vector-type container class.
 * 
 * Every particle type has its own DecayTable.
 *
 * @fn DecayTable::DecayTable()
 * @brief Default constructor.
 *
 * @fn DecayTable::DecayTable(const DecayTable& aTable)
 * @brief Copying constructor.
 * @param [in] aTable reference to DecayChannel
 *
 * @fn DecayTable::~DecayTable()
 * @brief Destructor
 *
 * @fn void DecayTable::AddDecayChannel(DecayChannel aChannel)
 * @brief Adds a new DecayChannel to the DecayTable and recalculates branching ratios.
 * @param [in] aChannel decay channel
 *
 * @fn int DecayTable::GetChannelCount() const
 * @brief Returns number of decay channel present.
 *
 * @fn const DecayChannel* DecayTable::GetDecayChannel(int aIndex) const
 * @brief Returns a DecayChannel pointer of the given index.
 * @param [in] aIndex index in DecayTable
 *
 * @fn float DecayTable::GetDecayStep(int aIndex)
 * @brief Returns a branching ration of the indexed channel.
 * @param [in] aIndex index in DecayTable
 *
 * @fn int DecayTable::ChooseDecayChannel(double aProb)
 * @brief Returns the last decay channel index, that has the branching less than the one given.
 * @param [in] aProb decay probability
 *
 * @fn int DecayTable::ChooseDecayChannelOrNot(double aProb)
 * @brief Chooses (or not) a decay channel from all available based on probability of decay.
 * @param [in] aProb decay probability
 * @retval >=0 index of the decay channel in DecayTable
 * @retval -1 decay channel not chosen
 */
