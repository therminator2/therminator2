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

#ifndef _TH2_PARTICLE_TYPE_H_
  #define _TH2_PARTICLE_TYPE_H_

#include <TString.h>
#include "DecayChannel.h"
#include "DecayTable.h"

class ParticleType {
  public:
    ParticleType();
    ParticleType(const ParticleType& aParticleType);
    ~ParticleType();

    int    GetNumber()		const;
    const char* GetName()	const;
    float  GetMass()		const;
    float  GetGamma()		const;
    float  GetSpin()		const;
    int    GetBarionN()		const;
    float  GetI()		const;
    float  GetI3()		const;
    int    GetStrangeN()	const;
    int    GetCharmN()		const;
    int    GetNumberQ()		const;
    int    GetNumberAQ()	const;
    int    GetNumberS()		const;
    int    GetNumberAS()	const;
    int    GetNumberC()		const;
    int    GetNumberAC()	const;
    int    GetCharge()		const;
    int    GetPDGCode()		const;
    double GetMaxIntegrand()	const;
    double GetMultiplicity()	const;
    int    GetDecayChannelCount2() const;
    int    GetDecayChannelCount3() const;
    
    void  SetNumber(int aNumber);
    void  SetName(const char *aName);
    void  SetMass(float aMass);
    void  SetGamma(float aGamma);
    void  SetSpin(float aSpin);
    void  SetBarionN(int aBarionN);
    void  SetI(float aI);
    void  SetI3(float aI3);
    void  SetStrangeN(int aStrangeN);
    void  SetCharmN(int aCharmN);
    void  SetNumberQ(int aNq);
    void  SetNumberAQ(int aNaq);
    void  SetNumberS(int aNs);
    void  SetNumberAS(int aNas);
    void  SetNumberC(int aNc);
    void  SetNumberAC(int aNac);
    void  SetPDGCode(int aCode);
    void  SetMaxIntegrand(double aMaxInt);
    void  SetMultiplicity(double aMulti);
    void  SetDecayChannelCount2(int aDCCount2);
    void  SetDecayChannelCount3(int aDCCount3);

    DecayTable* GetTable() const;
    void        AddDecayChannel(DecayChannel aChannel);

  private:
    int		mNumber;	// particle type number
    TString	mName;		// particle name
    float	mMass;		// mass
    float	mGamma;		// gamma
    float	mSpin;		// spin
    int		mBarionN;	// baryon number
    float	mI;		// isospin
    float	mI3;		// isospin 3rd component
    int		mStrangeN;	// strangeness
    int		mCharmN;	// charm
    int		mNq;		// number of Q  quarks
    int		mNaq;		// number of AQ quarks
    int		mNs;		// number of S  quarks
    int		mNas;		// number of AS quarks
    int		mNc;		// number of C  quarks
    int		mNac;   	// number of AC quarks
    int		mPDGCode; 
    double	mMaxIntegrand;	// maximal value of the integrand
    double	mMultiplicity;	// average multiplicity
    int		mDecayChannelCount2;      //number of channels in this case
    int		mDecayChannelCount3;      //number of channels in this case
    DecayTable*	mTable;
};

#endif

/*! @file ParticleType.h
 * @brief Definition of ParticleType class. %Particle type properties.
 */
/*! @class ParticleType
 * @brief Contains information on particle properties like: mass, spin, quark content, quantum numbers etc.
 *
 * @fn ParticleType::ParticleType()
 * @brief Default constructor.
 *
 * @fn ParticleType::ParticleType(const ParticleType& aParticleType)
 * @brief Copying constructor.
 * @param [in] aParticleType reference to ParticleType
 *
 * @fn ParticleType::~ParticleType()
 * @brief Destructor.
 *
 * @fn int ParticleType::GetNumber() const
 * @brief Returns particle index number in ParticleDB
 *
 * @fn const char* ParticleType::GetName() const
 * @brief Returns particle name
 *
 * @fn float ParticleType::GetMass() const
 * @brief Returns particle mass
 *
 * @fn float ParticleType::GetGamma() const
 * @brief Returns particle width
 *
 * @fn float ParticleType::GetSpin() const
 * @brief Returns particle spin
 *
 * @fn int ParticleType::GetBarionN() const
 * @brief Returns particle baryon number
 *
 * @fn float ParticleType::GetI() const
 * @brief Returns particle isospin
 *
 * @fn float ParticleType::GetI3() const
 * @brief Returns particle 3rd component of the isospin
 *
 * @fn int ParticleType::GetStrangeN() const
 * @brief Returns particle strangeness
 *
 * @fn int ParticleType::GetCharmN() const
 * @brief Returns particle charm number
 *
 * @fn int ParticleType::GetNumberQ() const
 * @brief Returns number of light quarks
 *
 * @fn int ParticleType::GetNumberAQ() const
 * @brief Returns number of anti-light quarks
 *
 * @fn int ParticleType::GetNumberS() const
 * @brief Returns number of strange quarks
 *
 * @fn int ParticleType::GetNumberAS() const
 * @brief Returns number of anti-strange quarks
 *
 * @fn int ParticleType::GetNumberC() const
 * @brief Returns number of charm quarks
 *
 * @fn int ParticleType::GetNumberAC() const
 * @brief Returns number of anti-charm quarks
 *
 * @fn int ParticleType::GetCharge() const
 * @brief Returns particle charge
 *
 * @fn int ParticleType::GetPDGCode() const
 * @brief Returns particle code in ParticleDataGroup
 *
 * @fn double ParticleType::GetMaxIntegrand() const
 * @brief  Returns particle maximal value of the integrand
 *
 * @fn double ParticleType::GetMultiplicity() const
 * @brief Returns particle multiplicity
 *
 * @fn int ParticleType::GetDecayChannelCount2() const
 * @brief Returns number of two-particle decay channels
 *
 * @fn int ParticleType::GetDecayChannelCount3() const
 * @brief Returns number of three-particle decay channels
 *
 * @fn DecayTable* ParticleType::GetTable() const
 * @brief Return a pointer to the DecayTable of that particle
 *
 * @fn void ParticleType::SetNumber(int aNumber)
 * @brief Set particle index number in ParticleDB
 * @param [in] aNumber index number in ParticleDB
 *
 * @fn void ParticleType::SetName(const char *aName)
 * @brief Set particle name
 * @param [in] aName particle name
 *
 * @fn void ParticleType::SetMass(float aMass)
 * @brief Set particle mass
 * @param [in] aMass particle mass
 *
 * @fn void ParticleType::SetGamma(float aGamma)
 * @brief Set particle width
 * @param [in] aGamma particle width
 *
 * @fn void ParticleType::SetSpin(float aSpin)
 * @brief Set particle spin
 * @param [in] aSpin particle spin
 *
 * @fn void ParticleType::SetBarionN(int aBarionN)
 * @brief Set particle baryon number
 * @param [in] aBarionN particle baryon number
 *
 * @fn void ParticleType::SetI(float aI)
 * @brief Set particle isospin
 * @param [in] aI particle isospin
 *
 * @fn void ParticleType::SetI3(float aI3)
 * @brief Set particle 3rd component of the isospin 
 * @param [in] aI3 3rd component of the isospin
 *
 * @fn void ParticleType::SetStrangeN(int aStrangeN)
 * @brief Set particle strangeness
 * @param [in] aStrangeN particle strangeness
 *
 * @fn void ParticleType::SetCharmN(int aCharmN)
 * @brief Set particle charm number
 * @param [in] aCharmN particle charm number
 *
 * @fn void ParticleType::SetNumberQ(int aNq)
 * @brief Set number of light quarks
 * @param [in] aNq number of light quarks
 *
 * @fn void ParticleType::SetNumberAQ(int aNaq)
 * @brief Set number of anti-light quarks
 * @param [in] aNaq number of anti-light quarks
 *
 * @fn void ParticleType::SetNumberS(int aNs)
 * @brief Set number of strange quarks
 * @param [in] aNs number of strange quarks
 *
 * @fn void ParticleType::SetNumberAS(int aNas)
 * @brief Set number of anti-strange quarks
 * @param [in] aNas number of anti-strange quarks
 *
 * @fn void ParticleType::SetNumberC(int aNc)
 * @brief Set number of charm quarks
 * @param [in] aNc number of charm quarks
 *
 * @fn void ParticleType::SetNumberAC(int aNac)
 * @brief Set number of anti-strange quarks
 * @param [in] aNac number of anti-strange quarks
 *
 * @fn void ParticleType::SetPDGCode(int aCode)
 * @brief Set particle code in ParticleDataGroup
 * @param [in] aCode particle code in ParticleDataGroup
 *
 * @fn void ParticleType::SetMaxIntegrand(double aMaxInt)
 * @brief Set particle maximal value of the integrand
 * @param [in] aMaxInt particle maximal value of the integrand
 *
 * @fn void ParticleType::SetMultiplicity(double aMulti)
 * @brief Set particle multiplicity
 * @param [in] aMulti particle multiplicity
 *
 * @fn void ParticleType::SetDecayChannelCount2(int aDCCount2)
 * @brief Set number of two-particle decay channels
 * @param [in] aDCCount2 number of two-particle decay channels
 *
 * @fn void ParticleType::SetDecayChannelCount3(int aDCCount3)
 * @brief Set number of three-particle decay channels
 * @param [in] aDCCount3 number of three-particle decay channels
 *
 * @fn void ParticleType::AddDecayChannel(DecayChannel aChannel)
 * @brief Adds a decay channel to DecayTable of the particle
 * @param [in] aChannel decay channel
 */
