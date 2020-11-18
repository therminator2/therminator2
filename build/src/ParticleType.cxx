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

#include "ParticleType.h"

ParticleType::ParticleType()
: mNumber(0), mName(""),
  mMass(-1.0), mGamma(-1.0), mSpin(-1.0),
  mBarionN(-1), mI(-1.0), mI3(-1.0), mStrangeN(-1), mCharmN(-1),
  mNq(-1), mNaq(-1), mNs(-1), mNas(-1), mNc(-1), mNac(-1),
  mPDGCode(0),
  mMaxIntegrand(0.0), mMultiplicity(0.0),
  mDecayChannelCount2(0), mDecayChannelCount3(0)
{
  mTable	= new DecayTable();
}

ParticleType::ParticleType(const ParticleType& aParticleType)
{
  mNumber	= aParticleType.GetNumber();
  mName		= aParticleType.GetName();
  mMass		= aParticleType.GetMass();
  mGamma	= aParticleType.GetGamma();
  mSpin		= aParticleType.GetSpin();
  mBarionN	= aParticleType.GetBarionN();
  mI		= aParticleType.GetI();
  mI3		= aParticleType.GetI3();
  mStrangeN	= aParticleType.GetStrangeN();
  mCharmN	= aParticleType.GetCharmN();
  mNq		= aParticleType.GetNumberQ();
  mNaq		= aParticleType.GetNumberAQ();
  mNs		= aParticleType.GetNumberS();
  mNas		= aParticleType.GetNumberAS();
  mNc		= aParticleType.GetNumberC();
  mNac		= aParticleType.GetNumberAC();
  mPDGCode	= aParticleType.GetPDGCode();
  mMaxIntegrand	= aParticleType.GetMaxIntegrand();
  mMultiplicity = aParticleType.GetMultiplicity();
  mDecayChannelCount2 = aParticleType.GetDecayChannelCount2();
  mDecayChannelCount3 = aParticleType.GetDecayChannelCount3();
  mTable	= new DecayTable(*(aParticleType.GetTable()));
}

ParticleType::~ParticleType()
{
  if (mTable)
    delete mTable;
}

int	ParticleType::GetNumber() const		{ return mNumber; }
const char*	ParticleType::GetName() const	{ return mName.Data(); }
float	ParticleType::GetMass() const		{ return mMass; }
float	ParticleType::GetGamma() const		{ return mGamma; }
float	ParticleType::GetSpin() const		{ return mSpin; }
int	ParticleType::GetBarionN() const	{ return mBarionN; }
float	ParticleType::GetI() const		{ return mI; }
float	ParticleType::GetI3() const		{ return mI3; }
int	ParticleType::GetNumberQ() const	{ return mNq; }
int	ParticleType::GetNumberAQ() const	{ return mNaq; }
int	ParticleType::GetNumberS() const	{ return mNs; }
int	ParticleType::GetNumberAS() const	{ return mNas; }
int	ParticleType::GetNumberC() const	{ return mNc; }
int	ParticleType::GetNumberAC() const	{ return mNac; }
int	ParticleType::GetPDGCode() const	{ return mPDGCode; }
int	ParticleType::GetStrangeN() const	{ return mStrangeN; }
int	ParticleType::GetCharmN() const		{ return mCharmN; }
double	ParticleType::GetMaxIntegrand() const	{ return mMaxIntegrand; }
double	ParticleType::GetMultiplicity() const	{ return mMultiplicity; }
int	ParticleType::GetDecayChannelCount2() const { return mDecayChannelCount2; }
int	ParticleType::GetDecayChannelCount3() const { return mDecayChannelCount3; }
int	ParticleType::GetCharge() const		{ return static_cast<int> (mI3+(mBarionN+mStrangeN)/2.); }

void	ParticleType::SetNumber(int aNumber)		{ mNumber = aNumber; }
void	ParticleType::SetName(const char *aName)	{ mName = aName; }
void	ParticleType::SetMass(float aMass)		{ mMass = aMass; }
void	ParticleType::SetGamma(float aGamma)		{ mGamma = aGamma; }
void	ParticleType::SetSpin(float aSpin)		{ mSpin = aSpin; }
void	ParticleType::SetBarionN(int aBarionN)		{ mBarionN = aBarionN; }
void	ParticleType::SetI(float aI)			{ mI = aI; }
void	ParticleType::SetI3(float aI3)			{ mI3 = aI3; }
void	ParticleType::SetStrangeN(int aStrangeN)	{ mStrangeN = aStrangeN; }
void	ParticleType::SetCharmN(int aCharmN)		{ mCharmN = aCharmN; }
void	ParticleType::SetNumberQ( int aNq)		{ mNq = aNq; }
void	ParticleType::SetNumberAQ( int aNaq)		{ mNaq = aNaq; }
void	ParticleType::SetNumberS( int aNs)		{ mNs = aNs; }
void	ParticleType::SetNumberAS( int aNas)		{ mNas = aNas; }
void	ParticleType::SetNumberC( int aNc)		{ mNc = aNc; }
void	ParticleType::SetNumberAC( int aNac)		{ mNac = aNac; }
void	ParticleType::SetPDGCode(int aCode)		{ mPDGCode = aCode; }
void	ParticleType::SetMaxIntegrand(double aMaxInt)	{ mMaxIntegrand = aMaxInt; }
void	ParticleType::SetMultiplicity(double aMulti)	{ mMultiplicity = aMulti; }
void	ParticleType::SetDecayChannelCount2(int aDCCount2) { mDecayChannelCount2 = aDCCount2; }
void	ParticleType::SetDecayChannelCount3(int aDCCount3) { mDecayChannelCount3 = aDCCount3; }

DecayTable* ParticleType::GetTable() const
{
  if (mTable)
    return mTable;
  else
    return 0;
}

void ParticleType::AddDecayChannel(DecayChannel aChannel)
{
  if (!mTable)
    mTable = new DecayTable();
  mTable->AddDecayChannel(aChannel);
}
