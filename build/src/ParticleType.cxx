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
#include <iostream>
#include <set>

static std::set<TString> sUnknownTyps;

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

  InitUrQMDCodes();
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

void ParticleType::InitUrQMDCodes()
{
  mPrefixToUrQMDCode["pr0"] = 1;
  mPrefixToUrQMDCode["ne0"] = 1;
  mPrefixToUrQMDCode["Ns1440"] = 2;
  mPrefixToUrQMDCode["Ns1520"] = 3;
  mPrefixToUrQMDCode["Ns1535"] = 4;
  mPrefixToUrQMDCode["Ns1650"] = 5;
  mPrefixToUrQMDCode["Ns1675"] = 6;
  mPrefixToUrQMDCode["Ns1680"] = 7;
  mPrefixToUrQMDCode["Ns1700"] = 8;
  mPrefixToUrQMDCode["Ns1710"] = 9;
  mPrefixToUrQMDCode["Ns1720"] = 10;
  mPrefixToUrQMDCode["Ns1900"] = 11;
  mPrefixToUrQMDCode["Ns1990"] = 12;
  mPrefixToUrQMDCode["Ns2080"] = 13;
  mPrefixToUrQMDCode["Ns2190"] = 14;
  mPrefixToUrQMDCode["Ns2200"] = 15;
  mPrefixToUrQMDCode["Ns2250"] = 16;
  mPrefixToUrQMDCode["Dl1232"] = 17;
  mPrefixToUrQMDCode["Dl1600"] = 18;
  mPrefixToUrQMDCode["Dl1620"] = 19;
  mPrefixToUrQMDCode["Dl1700"] = 20;
  mPrefixToUrQMDCode["Dl1900"] = 21;
  mPrefixToUrQMDCode["Dl1905"] = 22;
  mPrefixToUrQMDCode["Dl1910"] = 23;
  mPrefixToUrQMDCode["Dl1920"] = 24;
  mPrefixToUrQMDCode["Dl1930"] = 25;
  mPrefixToUrQMDCode["Dl1950"] = 26;
  mPrefixToUrQMDCode["Lm1115"] = 27;
  mPrefixToUrQMDCode["Lm1405"] = 28;
  mPrefixToUrQMDCode["Lm1520"] = 29;
  mPrefixToUrQMDCode["Lm1600"] = 30;
  mPrefixToUrQMDCode["Lm1670"] = 31;
  mPrefixToUrQMDCode["Lm1690"] = 32;
  mPrefixToUrQMDCode["Lm1800"] = 33;
  mPrefixToUrQMDCode["Lm1810"] = 34;
  mPrefixToUrQMDCode["Lm1820"] = 35;
  mPrefixToUrQMDCode["Lm1830"] = 36;
  mPrefixToUrQMDCode["Lm1890"] = 37;
  mPrefixToUrQMDCode["Lm2100"] = 38;
  mPrefixToUrQMDCode["Lm2110"] = 39;
  mPrefixToUrQMDCode["Sg1192"] = 40;
  mPrefixToUrQMDCode["Sg1385"] = 41;
  mPrefixToUrQMDCode["Sg1660"] = 42;
  mPrefixToUrQMDCode["Sg1670"] = 43;
  mPrefixToUrQMDCode["Sg1775"] = 44;
  mPrefixToUrQMDCode["Sg1790"] = 45;
  mPrefixToUrQMDCode["Sg1915"] = 46;
  mPrefixToUrQMDCode["Sg1940"] = 47;
  mPrefixToUrQMDCode["Sg2030"] = 48;
  mPrefixToUrQMDCode["Xi1321"] = 49;
  mPrefixToUrQMDCode["Xi1530"] = 50;
  mPrefixToUrQMDCode["Xi1690"] = 51;
  mPrefixToUrQMDCode["Xi1820"] = 52;
  mPrefixToUrQMDCode["Xi1950"] = 53;
  mPrefixToUrQMDCode["Xi2025"] = 54;
  mPrefixToUrQMDCode["Om1672"] = 55;
  
  mPrefixToUrQMDCode["gam"] = 100;
  mPrefixToUrQMDCode["pi"] = 101;
  mPrefixToUrQMDCode["eta5"] = 102;  
  mPrefixToUrQMDCode["om0"] = 103;
  mPrefixToUrQMDCode["rho"] = 104;  
  mPrefixToUrQMDCode["f00"] = 105;
  mPrefixToUrQMDCode["Ka04"] = 106;
  mPrefixToUrQMDCode["eta0prime"] = 107;
  mPrefixToUrQMDCode["Ka08"] = 108;
  mPrefixToUrQMDCode["ph1020"] = 109;  
  mPrefixToUrQMDCode["Ka1430"] = 110;
  mPrefixToUrQMDCode["a0"] = 111;
  mPrefixToUrQMDCode["f0137"] = 112;
  mPrefixToUrQMDCode["Ka12"] = 113;
  mPrefixToUrQMDCode["a1"] = 114;
  mPrefixToUrQMDCode["f112"] = 115;
  mPrefixToUrQMDCode["f114"] = 116;
  mPrefixToUrQMDCode["a2"] = 118;
  mPrefixToUrQMDCode["f212"] = 119;
  mPrefixToUrQMDCode["f215"] = 120;
  mPrefixToUrQMDCode["Ka1400"] = 121;
  mPrefixToUrQMDCode["b1"] = 122;
  mPrefixToUrQMDCode["h1"] = 123;
  mPrefixToUrQMDCode["Ka1410"] = 125;
  mPrefixToUrQMDCode["rh14"] = 126;
  mPrefixToUrQMDCode["om14"] = 127;
  mPrefixToUrQMDCode["ph16"] = 128;
  mPrefixToUrQMDCode["Ka16"] = 129;
  mPrefixToUrQMDCode["rh17"] = 130;
  mPrefixToUrQMDCode["om16"] = 131;
  mPrefixToUrQMDCode["ph18"] = 132;
  mPrefixToUrQMDCode["Jp"] = 135;

}

int ParticleType::GetUrQMDCode() const
{
  int result = 0;
  bool found = false;
  for (auto it = mPrefixToUrQMDCode.begin(); it != mPrefixToUrQMDCode.end(); ++it)
  {
    if (mName.BeginsWith(it->first)) {
      if (found) // already before
      {
        std::cout << "Found ambiguous UrQMD ityp for name " << mName.Data() << ": " << result << " and " << it->second << std::endl;
      }
      else {
        result = it->second;
        found = true;
      }
    }
  }
  if (!found) 
  {
    if (sUnknownTyps.find(mName) == sUnknownTyps.end()) {
      sUnknownTyps.insert(mName);
      std::cout << "Unknown UrQMD ityp for particle with name " << mName.Data() << " (it is printed only once)" << std::endl;
    }
  }
  return result;
}
