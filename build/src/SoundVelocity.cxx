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

#include "SoundVelocity.h"

SoundVelocity::SoundVelocity()
: mCsConst(0.0), mCsVar(0), mCsType(-1)
{
}

SoundVelocity::~SoundVelocity()
{
  if(mCsVar)
    delete mCsVar;
}

double SoundVelocity::GetCs() const
{
  return mCsConst;
}

double SoundVelocity::GetCs(double aX, double aY, double aZ) const
{
  return mCsVar->Interpolate(aX, aY, aZ);
}

double SoundVelocity::GetCs2() const
{
  return mCsConst * mCsConst;
}

double SoundVelocity::GetCs2(double aX, double aY, double aZ) const
{
  double tCs = mCsVar->Interpolate(aX, aY, aZ);
  return tCs * tCs;
}

int SoundVelocity::GetCsType() const
{
  return mCsType;
}

void SoundVelocity::SetCs(double aCs)
{
  mCsType  = 0;
  mCsConst = aCs;
}

void SoundVelocity::SetCs(Vector3D* aCs)
{
  mCsType = 1;
  mCsVar  = aCs;
}
