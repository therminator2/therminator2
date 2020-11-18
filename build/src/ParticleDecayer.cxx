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

#include "THGlobal.h"
#include "ParticleType.h"
#include "ParticleDecayer.h"
#include "DecayChannel.h"
#include "DecayTable.h"
#include <TMath.h>
#include <TDatime.h>
<<<<<<< HEAD
#include <iostream>

using namespace std;

extern int sSeed;
=======
>>>>>>> 6eef3d1fc1f3ac17aedbca11e4a34ae9c4017af9

ParticleDecayer::ParticleDecayer()
: mParticles(0), mDB(0), mTypeF(0), mFather(0), mChannel(0), mRandom(0)
{
}

ParticleDecayer::ParticleDecayer(ParticleDB *aDB, std::list<Particle>* aParticles)
: mDB(aDB), mParticles(aParticles)
{
  mRandom = new TRandom2();
#ifdef _ROOT_4_
<<<<<<< HEAD
  mRandom->SetSeed2(60356 + sSeed, 18233);
#else
  mRandom->SetSeed(60356 + sSeed);
=======
  mRandom->SetSeed2(60356, 18233);
#else
  mRandom->SetSeed(60356);
>>>>>>> 6eef3d1fc1f3ac17aedbca11e4a34ae9c4017af9
#endif
}

ParticleDecayer::~ParticleDecayer()
{
  delete mRandom;
}

int ParticleDecayer::DecayParticle(Particle* aFather)
{
  DecayTable*	tTable;
  int		tChannelIndex;
  double	tProb;
  
  mFather = aFather;
  mTypeF  = mFather->GetParticleType();
  tTable  = mTypeF->GetTable();
  
  tProb = mRandom->Rndm();
#ifdef _PARTICLE_DECAYER_RESCALE_CHANNELS_
  tChannelIndex = tTable->ChooseDecayChannel(tProb);
#else
  tChannelIndex = tTable->ChooseDecayChannelOrNot(tProb);
  if (tChannelIndex == -1) {
    PRINT_DEBUG_3("<ParticleDecayer::DecayParticle>\tNot decaying " << mTypeF->GetName() << " for prob " << tProb);
    for (int tIter=0; tIter<=tTable->GetChannelCount(); tIter++) {
      if (tTable->GetDecayChannel(tIter)->GetParticle3() > -1) {
	PRINT_DEBUG_3("\tD1 "<<mDB->GetParticleType(tTable->GetDecayChannel(tIter)->GetParticle1())->GetName() << " D2 "<<mDB->GetParticleType(tTable->GetDecayChannel(tIter)->GetParticle2())->GetName() << " BR "<<tTable->GetDecayChannel(tIter)->GetBranchingRatio());
      } else {
	PRINT_DEBUG_3("\tD1 "<<mDB->GetParticleType(tTable->GetDecayChannel(tIter)->GetParticle1())->GetName() << " D2 "<<mDB->GetParticleType(tTable->GetDecayChannel(tIter)->GetParticle2())->GetName() << " D3 "<<mDB->GetParticleType(tTable->GetDecayChannel(tIter)->GetParticle3())->GetName() << " BR "<<tTable->GetDecayChannel(tIter)->GetBranchingRatio());
      }     
    }
    return 0;
  }
#endif
  mChannel = tTable->GetDecayChannel(tChannelIndex);

  if (mChannel->Is3Particle()) {
#ifndef _PARTICLE_DECAYER_DISABLE_THREE_BODY_DECAYS_
    ThreeBodyDecay();
#endif
  } else {
#ifndef _PARTICLE_DECAYER_DISABLE_TWO_BODY_DECAYS_
    TwoBodyDecay();
#endif
  }
}

int ParticleDecayer::TwoBodyDecay()
{
  ParticleType*	tType1;
  ParticleType*	tType2;
  Particle*	tDaughter1;
  Particle*	tDaughter2;
  double	Xt, Xx, Xy, Xz;
  double	Pe, Px, Py, Pz;

<<<<<<< HEAD
  int pdg = mFather->pid;


  mFather->GetPosition(&Xt, &Xx, &Xy, &Xz);
  mFather->GetMomentum(&Pe, &Px, &Py, &Pz);

//  Float_t masa = sqrt(pow(mFather->e,2) - pow(mFather->px,2) - pow(mFather->py,2) - pow(mFather->pz,2));
//  if (abs(pdg) == 2224 || abs(pdg) == 2214 || abs(pdg) == 2114 || abs(pdg) == 1114)  cout<<"masa delty : "<<masa<<"  ";



=======
  mFather->GetPosition(&Xt, &Xx, &Xy, &Xz);
  mFather->GetMomentum(&Pe, &Px, &Py, &Pz);
  
>>>>>>> 6eef3d1fc1f3ac17aedbca11e4a34ae9c4017af9
  tType1 = mDB->GetParticleType(mChannel->GetParticle1());
  tType2 = mDB->GetParticleType(mChannel->GetParticle2());
  tDaughter1 = new Particle(tType1);
  tDaughter2 = new Particle(tType2);

  double tE  = Pe;
<<<<<<< HEAD
  double tM  = TMath::Sqrt(Pe*Pe - Px*Px - Py*Py - Pz*Pz);

  //double tM  = mTypeF->GetMass();
=======
  double tM  = mTypeF->GetMass();
>>>>>>> 6eef3d1fc1f3ac17aedbca11e4a34ae9c4017af9
  double tM1 = tType1->GetMass();
  double tM2 = tType2->GetMass();
  
  double tTime;
  if (mTypeF->GetGamma() == 0.0)
    tTime = 1.0e10;
  else {
    double tTau0 = tE / (mTypeF->GetMass() * mTypeF->GetGamma());
    // When it decays
    tTime = -tTau0 * TMath::Log(mRandom->Rndm());
  }
  // Decay coordinates
  double rxr = Xx + (Px / tE) * tTime;
  double ryr = Xy + (Py / tE) * tTime;
  double rzr = Xz + (Pz / tE) * tTime;
  double rtr = Xt + tTime;

  // Decay energy
  double tMC1 = (tM * tM - (tM1 + tM2) * (tM1 + tM2) );
  double tMC2 = (tM * tM - (tM1 - tM2) * (tM1 - tM2) );
  double tMom = TMath::Sqrt(tMC1 * tMC2) / (2 * tM);
  double tPhi = mRandom->Rndm() * 2 * TMath::Pi();
  double tCosTh = 2.0 * mRandom->Rndm() - 1.0;

  double tPtr  = tMom * TMath::Sqrt(1 - tCosTh * tCosTh);
  double tPxr1 = tPtr * TMath::Cos(tPhi);
  double tPyr1 = tPtr * TMath::Sin(tPhi);
  double tPzr1 = tMom * tCosTh;

  double tVx = Px / Pe;
  double tVy = Py / Pe;
  double tVz = Pz / Pe;

  double tES1 = TMath::Sqrt(tM1*tM1 + tPxr1*tPxr1 + tPyr1*tPyr1 + tPzr1*tPzr1);
  double tES2 = TMath::Sqrt(tM2*tM2 + tPxr1*tPxr1 + tPyr1*tPyr1 + tPzr1*tPzr1);

  double tV2 = tVx*tVx + tVy*tVy + tVz*tVz;
  double tGamma = TMath::Power(1 - tV2, -0.5);
  double tVP = tVx*tPxr1 + tVy*tPyr1 + tVz*tPzr1;
  double tgvp = (tGamma - 1.0) * (1.0/tV2) * tVP;

  double tPx1 =  tPxr1 + (tgvp + tGamma * tES1) * tVx;
  double tPy1 =  tPyr1 + (tgvp + tGamma * tES1) * tVy;
  double tPz1 =  tPzr1 + (tgvp + tGamma * tES1) * tVz;
  double tPe1 = TMath::Sqrt(tM1*tM1 + tPx1*tPx1 + tPy1*tPy1 + tPz1*tPz1);
    
  double tPx2 = -tPxr1 + (-tgvp + tGamma * tES2) * tVx;
  double tPy2 = -tPyr1 + (-tgvp + tGamma * tES2) * tVy;
  double tPz2 = -tPzr1 + (-tgvp + tGamma * tES2) * tVz;
  double tPe2 = TMath::Sqrt(tM2*tM2 + tPx2*tPx2 + tPy2*tPy2 + tPz2*tPz2);

  mFather->SetDecayed();
  tDaughter1->SetParticlePX(tPe1, tPx1, tPy1, tPz1, rtr, rxr, ryr, rzr, mFather);
  tDaughter2->SetParticlePX(tPe2, tPx2, tPy2, tPz2, rtr, rxr, ryr, rzr, mFather);
  mParticles->push_back(*tDaughter1);
  mParticles->push_back(*tDaughter2);

  delete tDaughter1;
  delete tDaughter2;
}

int ParticleDecayer::ThreeBodyDecay()
{
  ParticleType*	tType1;
  ParticleType*	tType2;
  ParticleType*	tType3;
  Particle*	tDaughter1;
  Particle*	tDaughter2;
  Particle*	tDaughter3;
  double	Xt, Xx, Xy, Xz;
  double	Pe, Px, Py, Pz;

  mFather->GetPosition(&Xt, &Xx, &Xy, &Xz);
  mFather->GetMomentum(&Pe, &Px, &Py, &Pz);
  
  tType1 = mDB->GetParticleType(mChannel->GetParticle1());
  tType2 = mDB->GetParticleType(mChannel->GetParticle2());
  tType3 = mDB->GetParticleType(mChannel->GetParticle3());
  tDaughter1 = new Particle(tType1);
  tDaughter2 = new Particle(tType2);
  tDaughter3 = new Particle(tType3);
  
  double tE  = Pe;
<<<<<<< HEAD
  double tM  = TMath::Sqrt(Pe*Pe - Px*Px - Py*Py - Pz*Pz);

  //double tM  = mTypeF->GetMass();
=======
  double tM  = mTypeF->GetMass();
>>>>>>> 6eef3d1fc1f3ac17aedbca11e4a34ae9c4017af9
  double tM1 = tType1->GetMass();
  double tM2 = tType2->GetMass();
  double tM3 = tType3->GetMass();

  double tES1, tES2, tP1, tP2, tCos12, tZ;
  do {
    // Generate E1 and E2 with the Monte-Carlo method
    do {
      tES1 = mRandom->Rndm() * (tM - tM2 - tM3 - tM1) + tM1;
      tES2 = mRandom->Rndm() * (tM - tM1 - tM3 - tM2) + tM2;
    } while (tES1+tES2 > tM); // The sum of both energies must be smaller than the resonance mass
    tP1  = TMath::Sqrt(tES1*tES1 - tM1*tM1);
    tP2  = TMath::Sqrt(tES2*tES2 - tM2*tM2);
    tZ = tM - tES1 - tES2;
    tZ *= tZ;
    tCos12 = (tZ - tP1*tP1 - tP2*tP2 - tM3*tM3)/(2*tP1*tP2);
  } while ((tCos12 < -1.0) || (tCos12 > 1.0)); // Cos Theta must exist (be within -1.0 to 1.0 )

  double tTime;
  if (mTypeF->GetGamma() == 0.0)
    tTime = 1.0e10;
  else {
    double tTau0 = tE / (mTypeF->GetMass() * mTypeF->GetGamma());
    // When it decays
    tTime = -tTau0 * TMath::Log(mRandom->Rndm());
  }

// Decay coordinates
  double rxr = Xx + (Px/tE)*tTime;
  double ryr = Xy + (Py/tE)*tTime;
  double rzr = Xz + (Pz/tE)*tTime;
  double rtr = Xt + tTime;

  double tPxr2 = tP2 * TMath::Sqrt(1-tCos12*tCos12);
  double tPzr2 = tP2*tCos12;
  double tPxr3 = - tPxr2;
  double tPzr3 = - (tP1 + tPzr2);
  double tP3 = TMath::Hypot(tPxr3, tPzr3);
  double tES3 = TMath::Hypot(tM3, tP3);

  // Generating Euler angles
  double tPhi = mRandom->Rndm() * 2 * TMath::Pi();
  double tKsi = mRandom->Rndm() * 2 * TMath::Pi();
  double tCosTh = mRandom->Rndm() * 2.0 - 1.0;

  double sp = TMath::Sin(tPhi);
  double cp = TMath::Cos(tPhi);
  double sk = TMath::Sin(tKsi);
  double ck = TMath::Cos(tKsi);
  double st = TMath::Sqrt(1.0-tCosTh*tCosTh);
  double ct = tCosTh;

  // Rotating the whole system
  double tPxp1 = - st*ck * tP1;
  double tPyp1 = st*sk * tP1;
  double tPzp1 = ct * tP1;

  double tPxp2 = (cp*ct*ck - sp*sk)  * tPxr2 + (-st*ck) * tPzr2;
  double tPyp2 = (-cp*ct*sk - sp*ck) * tPxr2 + (st*sk)  * tPzr2;
  double tPzp2 = cp*st               * tPxr2 + ct       * tPzr2;

  double tPxp3 = (cp*ct*ck - sp*sk)  * tPxr3 + (-st*ck) * tPzr3;
  double tPyp3 = (-cp*ct*sk - sp*ck) * tPxr3 + (st*sk)  * tPzr3;
  double tPzp3 = cp*st               * tPxr3 + ct       * tPzr3;

  double tVx = Px / Pe;
  double tVy = Py / Pe;
  double tVz = Pz / Pe;

  tES1 = TMath::Sqrt(tM1*tM1 + tPxp1*tPxp1 + tPyp1*tPyp1 + tPzp1*tPzp1);
  tES2 = TMath::Sqrt(tM2*tM2 + tPxp2*tPxp2 + tPyp2*tPyp2 + tPzp2*tPzp2);
  tES3 = TMath::Sqrt(tM3*tM3 + tPxp3*tPxp3 + tPyp3*tPyp3 + tPzp3*tPzp3);

  double tV2 = tVx*tVx + tVy*tVy + tVz*tVz;
  double tGamma = TMath::Power(1-tV2,-0.5);

  // Boosting by the parent velocity
  double tVP = tVx*tPxp1 + tVy*tPyp1 + tVz*tPzp1;
  double tgvp = (tGamma - 1.0) * (1.0/tV2) * tVP;

  double tPx1 = tPxp1 + (tgvp + tGamma * tES1) * tVx;
  double tPy1 = tPyp1 + (tgvp + tGamma * tES1) * tVy;
  double tPz1 = tPzp1 + (tgvp + tGamma * tES1) * tVz;
  double tPe1 = TMath::Sqrt(tM1*tM1 + tPx1*tPx1 + tPy1*tPy1 + tPz1*tPz1);

  tVP = tVx*tPxp2 + tVy*tPyp2 + tVz*tPzp2;
  tgvp = (tGamma - 1.0) * (1.0/tV2) * tVP;

  double tPx2 = tPxp2 + (tgvp + tGamma * tES2) * tVx;
  double tPy2 = tPyp2 + (tgvp + tGamma * tES2) * tVy;
  double tPz2 = tPzp2 + (tgvp + tGamma * tES2) * tVz;
  double tPe2 = TMath::Sqrt(tM2*tM2 + tPx2*tPx2 + tPy2*tPy2 + tPz2*tPz2);

  tVP = tVx*tPxp3 + tVy*tPyp3 + tVz*tPzp3;
  tgvp = (tGamma - 1.0) * (1.0/tV2) * tVP;

  double tPx3 = tPxp3 + (tgvp + tGamma * tES3) * tVx;
  double tPy3 = tPyp3 + (tgvp + tGamma * tES3) * tVy;
  double tPz3 = tPzp3 + (tgvp + tGamma * tES3) * tVz;
  double tPe3 = TMath::Sqrt(tM3*tM3 + tPx3*tPx3 + tPy3*tPy3 + tPz3*tPz3);

  mFather->SetDecayed();
  tDaughter1->SetParticlePX(tPe1, tPx1, tPy1, tPz1, rtr, rxr, ryr, rzr, mFather);
  tDaughter2->SetParticlePX(tPe2, tPx2, tPy2, tPz2, rtr, rxr, ryr, rzr, mFather);
  tDaughter3->SetParticlePX(tPe3, tPx3, tPy3, tPz3, rtr, rxr, ryr, rzr, mFather);
  mParticles->push_back(*tDaughter1);
  mParticles->push_back(*tDaughter2);
  mParticles->push_back(*tDaughter3);

  delete tDaughter1;
  delete tDaughter2;
  delete tDaughter3;
}

void ParticleDecayer::Randomize()
{
  TDatime tDate; 
  
#ifdef _ROOT_4_
  mRandom->SetSeed2(tDate.Get(), (tDate.Get() % 11) * 7 + (tDate.Get() / 7));
#else
  mRandom->SetSeed(tDate.Get() / 2 * 3);
#endif

}

void ParticleDecayer::SeedSet(int aSeed)
{
  
#ifdef _ROOT_4_
  mRandom->SetSeed2(aSeed, aSeed * 11 % 9);
#else
  mRandom->SetSeed(aSeed);
#endif
}

inline double ParticleDecayer::BreitWigner(double Mass, double Gamma) const
{
  double x,y;

  y = mRandom->Rndm();
  x = Mass + Gamma/2 * TMath::Tan( TMath::Pi() * (y - 0.5) );

  return x;
}
