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

#include <TMath.h>
#include "THGlobal.h"
#include "ParticleCoor.h"

#include <iostream>

ParticleCoor& ParticleCoor::operator= (ParticleCoor& aPartCoor)
{
  mass	= aPartCoor.mass;
  t	= aPartCoor.t;
  x	= aPartCoor.x;
  y	= aPartCoor.y;
  z	= aPartCoor.z;
  e	= aPartCoor.e;
  px	= aPartCoor.px;
  py	= aPartCoor.py;
  pz	= aPartCoor.pz;
  pid	= aPartCoor.pid;
  fatherpid	= aPartCoor.fatherpid;
  rootpid	= aPartCoor.rootpid;
  eid		= aPartCoor.eid;
  fathereid	= aPartCoor.fathereid;
  decayed	= aPartCoor.decayed;
  eventid	= aPartCoor.eventid;
  return *this;
}

ParticleCoor& ParticleCoor::operator= (ParticleCoor* aPartCoor)
{
  mass	= aPartCoor->mass;
  t	= aPartCoor->t;
  x	= aPartCoor->x;
  y	= aPartCoor->y;
  z	= aPartCoor->z;
  e	= aPartCoor->e;
  px	= aPartCoor->px;
  py	= aPartCoor->py;
  pz	= aPartCoor->pz;
  pid	= aPartCoor->pid;
  fatherpid	= aPartCoor->fatherpid;
  rootpid	= aPartCoor->rootpid;
  eid		= aPartCoor->eid;
  fathereid	= aPartCoor->fathereid;
  decayed	= aPartCoor->decayed;
  eventid	= aPartCoor->eventid;
  return *this;
}

void ParticleCoor::GetMomentum(double* aPe, double* aPx, double* aPy, double* aPz) const
{
  *aPe = e; *aPx = px; *aPy = py; *aPz = pz;
}

void ParticleCoor::GetPosition(double* aXt, double* aXx, double* aXy, double* aXz) const
{
  *aXt = t; *aXx = x; *aXy = y; *aXz = z;
}

Float_t ParticleCoor::GetTau() const
{
  return TMath::Sqrt(t*t - z*z)*kHbarC;
}

Float_t ParticleCoor::GetR() const
{
  return TMath::Sqrt(x*x + y*y + z*z)*kHbarC;
}


Float_t ParticleCoor::GetRho() const
{
  return TMath::Sqrt(x*x + y*y)*kHbarC;
}

Float_t ParticleCoor::GetPhiS() const
{
  return TMath::ATan2(y,x);
}

Float_t ParticleCoor::GetRapidityS() const
{
  return 0.5 * TMath::Log( (t + z) / (t - z) );
}

Float_t ParticleCoor::GetP() const
{
  return TMath::Sqrt(px*px + py*py + pz*pz);
}

Float_t ParticleCoor::GetPt() const
{
  return TMath::Sqrt(px*px + py*py);
}

Float_t ParticleCoor::GetMt() const
{
  return TMath::Sqrt(mass*mass + px*px + py*py);
}

Float_t ParticleCoor::GetPhiP() const
{
  return TMath::ATan2(py,px);
}

Float_t ParticleCoor::GetRapidityP() const
{
  return 0.5 * TMath::Log( (e + pz) / (e - pz) );
}   

Float_t ParticleCoor::GetEtaP() const
{
  Float_t tP = GetP();
  return 0.5 * TMath::Log( (tP + pz) / (tP - pz) );
}

Int_t ParticleCoor::GetDecayed() const
{
  return decayed;
}

void ParticleCoor::SetDecayed()
{
  decayed = 1;
}

void ParticleCoor::TransformToLCMS(Float_t aBetaZ)
{
  Float_t tmp;
  Float_t tGammaZ = 1.0 / TMath::Sqrt(1.0 - aBetaZ*aBetaZ);
  
  tmp = tGammaZ * (e  - aBetaZ * pz);
  pz  = tGammaZ * (pz - aBetaZ * e );
  e   = tmp;
  
  tmp = tGammaZ * (t  - aBetaZ * z );
  z   = tGammaZ * (z  - aBetaZ * t );
  t   = tmp;
}

void ParticleCoor::TransformRotateZ(Float_t aPhi)
{
  Float_t tmp;
  Float_t tCosPhi = TMath::Cos(aPhi);
  Float_t tSinPhi = TMath::Sin(aPhi);
  
  tmp = ( px * tCosPhi + py * tSinPhi);
  py  = ( py * tCosPhi - px * tSinPhi);
  px  = tmp;
  
  tmp = ( x  * tCosPhi + y  * tSinPhi);
  y   = ( y  * tCosPhi - x  * tSinPhi);
  x   = tmp;
}

void ParticleCoor::TransformToPRF(Float_t aBetaT)
{
  Float_t tmp;
  Float_t tGammaT = 1.0 / TMath::Sqrt(1.0 - aBetaT*aBetaT);
  
  tmp = tGammaT * (e  - aBetaT * px);
  px  = tGammaT * (px - aBetaT * e);
  e   = tmp;
  
  tmp = tGammaT * (t  - aBetaT * x);
  x   = tGammaT * (x  - aBetaT * t);
  t   = tmp;
}
