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

#include "Viscosity.h"

Viscosity::Viscosity()
: mViscosityType(0), mBulk(0), mTau(0.0)
{
  for(int i=0; i<10; i++)
    mShear[i] = 0;
}

Viscosity::~Viscosity()
{
  delete mBulk;
  for(int i=0; i<10; i++)
    delete mShear[i];  
}

int Viscosity::GetViscosityType() const
{
  return mViscosityType;
}

double Viscosity::GetBulkTau() const
{
  return mTau;
}

double Viscosity::GetBulkPi(double aX, double aY, double aZ) const
{
  return mBulk->Interpolate(aX,aY,aZ);
}

double Viscosity::GetShearPPdotPi(double aX, double aY, double aZ, double aPe, double aPx, double aPy, double aPz) const
{
  return (
      aPe*aPe * mShear[i00]->Interpolate(aX,aY,aZ)
    + aPx*aPx * mShear[iXX]->Interpolate(aX,aY,aZ)
    + aPy*aPy * mShear[iYY]->Interpolate(aX,aY,aZ)
    + aPz*aPz * mShear[iZZ]->Interpolate(aX,aY,aZ)
    + 2.0 * ( 
        aPx*aPy * mShear[iXY]->Interpolate(aX,aY,aZ)
      + aPx*aPz * mShear[iXZ]->Interpolate(aX,aY,aZ)
      + aPy*aPz * mShear[iYZ]->Interpolate(aX,aY,aZ)
    ) - 2.0 * aPe * (
        aPx * mShear[i0X]->Interpolate(aX,aY,aZ)
      + aPy * mShear[i0Y]->Interpolate(aX,aY,aZ)
      + aPz * mShear[i0Z]->Interpolate(aX,aY,aZ)
    )
  );
}
   
void Viscosity::SetBulkTau(double aTau)
{
  mTau = aTau;
}

void Viscosity::SetBulkPi(Vector3D* aBulk)
{
  mViscosityType = mViscosityType | 0x01;
  mBulk = aBulk;
}

void Viscosity::SetShearPi(eViscPi aIndex, Vector3D* aShear)
{
  mViscosityType = mViscosityType | 0x02;
  mShear[aIndex] = aShear;
}
