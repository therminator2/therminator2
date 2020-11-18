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
#include "Thermodynamics.h"

using namespace TMath;

Thermodynamics::Thermodynamics()
: Chemistry(), Energy(), Entropy(), Pressure(), SoundVelocity(), Temperature(), Viscosity()
{
}

Thermodynamics::~Thermodynamics()
{
}

double Thermodynamics::GetInvFugacityFactor(ParticleType* aPartType) const
{
  switch (GetChemistryType()) {
    case 0:  return Exp( - GetChemicalPotential(aPartType) / GetTemperature() );
    case 2:  return 1.0 / GetFugacity(aPartType);
    default: return 1.0;
  }
}

double Thermodynamics::GetInvFugacityFactor(ParticleType* aPartType, double aX, double aY, double aZ) const
{
  switch (GetChemistryType()) {
    case 1:  return Exp( - GetChemicalPotential(aPartType, aX, aY, aZ) / GetTemperature(aX, aY, aZ) );
    case 3:  return 1.0 / GetFugacity(aPartType, aX, aY, aZ);
    default: return 1.0;
  }
}

double Thermodynamics::GetViscosityCorrection(double aX, double aY, double aZ, double aPe, double aPx, double aPy, double aPz, double aPdotU) const
{
  double tMass2;
  double tTemp;
  double tCs2;
  double tBulkPart  = 0.0;
  double tShearPart = 0.0;
  
  tMass2 = aPe*aPe - aPx*aPx - aPy*aPy - aPz*aPz;
  
  switch (GetTemperatureType()) {
    case 0:
      tTemp = GetTemperature();
      break;
    case 1:
      tTemp = GetTemperature(aX, aY, aZ);
      break;
  }
  switch (GetCsType()) {
    case 0:
      tCs2 = GetCs2();
      break;
    case 1:
      tCs2 = GetCs2(aX, aY, aZ);
      break;
  }
    
  switch (GetViscosityType()) {
    case 0:
      tShearPart  = 1.0 / (2.0*tTemp*tTemp) * GetShearPPdotPi(aX,aY,aZ,aPe,aPx,aPy,aPz);
      break;
    case 1:
      tBulkPart = ( tCs2*aPdotU + (tMass2 - aPdotU*aPdotU) / (3.0*aPdotU) ) * GetBulkTau() * GetBulkPi(aX,aY,aZ) / tTemp;
      break;
    case 2:
      tShearPart  = 1.0 / (2.0*tTemp*tTemp) * GetShearPPdotPi(aX,aY,aZ,aPe,aPx,aPy,aPz);
      tBulkPart = ( tCs2*aPdotU + (tMass2 - aPdotU*aPdotU) / (3.0*aPdotU) ) * GetBulkTau() * GetBulkPi(aX,aY,aZ) / tTemp;
      break;      
  }
  return (tBulkPart + tShearPart);
}
