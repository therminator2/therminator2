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

#include <sstream>
#include "THGlobal.h"
#include "Particle.h"

using namespace std;

long int Particle::EIDi = 0;

Particle::Particle()
: ParticleCoor(), mPartType(0) 
{
}

Particle::Particle(ParticleType *aType)
: ParticleCoor(), mPartType(aType)
{
  mass      = mPartType->GetMass();
  pid       = mPartType->GetPDGCode();
  fatherpid = pid;  
  rootpid   = pid;    
  eid       = EIDi;
  fathereid = -1;
}


Particle::Particle(ParticleCoor *aCoor, ParticleType *aType)
: ParticleCoor(*aCoor), mPartType(aType)
{ }

Particle::~Particle()
{
}

ParticleCoor& Particle::operator= (Particle& aParticle)
{
  mass	= aParticle.mass;
  t	= aParticle.t;
  x	= aParticle.x;
  y	= aParticle.y;
  z	= aParticle.z;
  e	= aParticle.e;
  px	= aParticle.px;
  py	= aParticle.py;
  pz	= aParticle.pz;
  pid	= aParticle.pid;
  fatherpid	= aParticle.fatherpid;
  rootpid	= aParticle.rootpid;
  eid		= aParticle.eid;
  fathereid	= aParticle.fathereid;
  decayed	= aParticle.decayed;
  eventid	= aParticle.eventid;
  return *this;
}

ParticleType* Particle::GetParticleType() const
{
  return mPartType;
}

void Particle::SetParticlePX(double aPe, double aPx, double aPy, double aPz,
			     double aXt, double aXx, double aXy, double aXz,
			     Particle* aFather)
{
  e = aPe; px = aPx; py = aPy; pz = aPz;
  t = aXt; x = aXx;  y = aXy;  z = aXz;
  if (aFather) {
    fatherpid = aFather->pid;
    rootpid   = aFather->rootpid;
    fathereid = aFather->eid;
  }  
  eid = EIDi++;
}

const std::string Particle::MakeTEXTEntry()
{
  ostringstream oss;

  oss.flags(std::ios::right);
  oss.precision(6);
  oss.width(6);
  oss << eid <<'\t'<< fathereid <<'\t';
  oss << pid <<'\t'<< fatherpid <<'\t'<< rootpid <<'\t';
  oss << decayed <<'\t';
  oss.flags(std::ios::left | std::ios::scientific);
  oss.width(15);
  oss << mass <<'\t';
  oss << e        <<'\t'<< px       <<'\t'<< py       <<'\t'<< pz       <<'\t';
  oss << t*kHbarC <<'\t'<< x*kHbarC <<'\t'<< y*kHbarC <<'\t'<< z*kHbarC <<'\n';
  return oss.str();
}

void Particle::ZeroEID()
{
  EIDi = 0; 
}
