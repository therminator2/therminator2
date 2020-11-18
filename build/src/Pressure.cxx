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

#include "Pressure.h"

Pressure::Pressure()
: mPressConst(0.0), mPressVar(0), mPressType(-1)
{
}

Pressure::~Pressure()
{
  if(mPressVar)
    delete mPressVar;
}

double Pressure::GetPressure() const
{
  return mPressConst;
}

double Pressure::GetPressure(double aX, double aY, double aZ) const
{
  return mPressVar->Interpolate(aX, aY, aZ);
}

int Pressure::GetPressureType() const
{
  return mPressType;
}

void Pressure::SetPressure(double aPress)
{
  mPressType = 0;
  mPressConst = aPress;
}

void Pressure::SetPressure(Vector3D* aPress)
{
  mPressType = 1;
  mPressVar = aPress;
}
