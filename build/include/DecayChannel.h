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

#ifndef _TH2_DECAY_CHANNEL_H_
  #define _TH2_DECAY_CHANNEL_H_
#include "TH1D.h"

class DecayChannel {
  public:
    DecayChannel();
    DecayChannel(const DecayChannel& aChannel);
    DecayChannel(double aBranchRatio, int aPartType1, int aPartType2, int aPartType3 = -1);
    ~DecayChannel();

    int    GetParticle1() const;
    int    GetParticle2() const;
    int    GetParticle3() const;
    double GetBranchingRatio() const;
    TH1D* GetBranchingRatioFunc() const;
    int    Is3Particle() const;
    bool    IsBRFunc() const;

    void   SetParticle1(int aPartType1);
    void   SetParticle2(int aPartType2);
    void   SetParticle3(int aPartType3);
    void   SetBranchingRatio(double aRatio);
    void   SetBranchingRatioFunc(TH1D* aHist);
    

  private:
    int    mParticleType1;
    int    mParticleType2;
    int    mParticleType3;
    double mBranchRatio;
    TH1D *mBRFunc;
    bool  mHasBRFunc;
};

#endif

/*! @file DecayChannel.h
 * @brief Definition of DecayChannel class. Stores the information on the particle decay channel.
 */
/*! @class DecayChannel
 * @brief Contains information on the types of child particles that are produced in a two or three 
 * particle decay and the decay channel branching ratio.
 *
 * @fn DecayChannel::DecayChannel()
 * @brief Default constructor.
 *
 * @fn DecayChannel::DecayChannel(const DecayChannel& aChannel)
 * @brief Copying constructor.
 * @param [in] aChannel reference to DecayChannel
 *
 * @fn DecayChannel::DecayChannel(double aBranchRatio, int aPartType1, int aPartType2, int aPartType3 = -1)
 * @brief Create a new DecayChannel with branching ratio and child particle types indexes in ParticleDB.
 * @param [in] aBranchRatio decay channel branching ratio
 * @param [in] aPartType1 first child particle type index (two or three particle decay)
 * @param [in] aPartType2 second child particle type index (two or three particle decay)
 * @param [in] aPartType3 third child particle type index (three particle decay only)
 *
 * @fn DecayChannel::~DecayChannel()
 * @brief Destructor.
 *
 * @fn int DecayChannel::GetParticle1() const
 * @brief Returns the index in ParticleDB of the first child particle.
 *
 * @fn int DecayChannel::GetParticle2() const
 * @brief Returns the index in ParticleDB of the second child particle.
 *
 * @fn int DecayChannel::GetParticle3() const
 * @brief Returns the index in ParticleDB of the third child particle.
 *
 * @fn double DecayChannel::GetBranchingRatio() const
 * @brief Returns the branching ratio of this decay channel.
 *
 * @fn int DecayChannel::Is3Particle() const
 * @brief If it is a three particle decay channel.
 * @retval 0 two particle decay channel
 * @retval 1 three particle decay channel
 *
 * @fn void DecayChannel::SetParticle1(int aPartType1)
 * @brief Sets the first child particle index.
 * @param [in] aPartType1 particle index in ParticleDB
 *
 * @fn void DecayChannel::SetParticle2(int aPartType2)
 * @brief Sets the second child particle index.
 * @param [in] aPartType2 particle index in ParticleDB
 *
 * @fn void DecayChannel::SetParticle3(int aPartType3)
 * @brief Sets the third child particle index.
 * @param [in] aPartType3 particle index in ParticleDB
 *
 * @fn void DecayChannel::SetBranchingRatio(double aRatio)
 * @brief Sets the decay channel branching ratio
 * @param [in] aRatio branching ratio
 */
