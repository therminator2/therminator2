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

#ifndef _TH2_THERMODYNAMICS_H_
  #define _TH2_THERMODYNAMICS_H_

#include "Temperature.h"
#include "Chemistry.h"
#include "Energy.h"
#include "Entropy.h"
#include "Pressure.h"
#include "SoundVelocity.h"
#include "Viscosity.h"

class Thermodynamics
: public Chemistry,
  public Energy,
  public Entropy,
  public Pressure,
  public SoundVelocity,
  public Temperature,
  public Viscosity
{
  public:
    Thermodynamics();
    ~Thermodynamics();
    
    double GetInvFugacityFactor(ParticleType* aPartType) const;
    double GetInvFugacityFactor(ParticleType* aPartType, double aX, double aY, double aZ) const;
    double GetViscosityCorrection(double aX, double aY, double aZ, double aPe, double aPx, double aPy, double aPz, double aPdotU) const;
};

#endif

/*! @file Thermodynamics.h
 * @brief Definition of Thermodynamics class. Contains Temperature, Chemistry classes.
 */
/*! @class Thermodynamics
 * @brief Class inherits from Temperature and Chemistry classes.
 *
 * @fn Thermodynamics::Thermodynamics()
 * @brief Default constructor.
 *
 * @fn Thermodynamics::~Thermodynamics()
 * @brief Destructor.
 *
 * @fn double Thermodynamics::GetInvFugacityFactor(ParticleType* aPartType) const
 * @brief Returns the inverse fugacity factor @f$ \Upsilon^{-1} @f$ for a given ParticleType that depends on the GetChemistryType().
 * 
 * In case of local chemical equilibrium @f$ \Upsilon^{-1} = \exp \left[ - \frac{\mu_B B + \mu_{I_3} I_3 + \mu_S S + \mu_C C}{T}\right] @f$.
 * Otherwise the returned value is just simply @f$ \Upsilon^{-1} = \frac{1}{\Upsilon} @f$, see GetFugacity().
 * @param [in] aPartType type of particle
 *
 * @fn double Thermodynamics::GetInvFugacityFactor(ParticleType* aPartType, double aX, double aY, double aZ) const
 * @brief Returns the inverse fugacity factor @f$ \Upsilon^{-1} @f$ for a given ParticleType that depends on the GetChemistryType().
 * 
 * In case of local chemical equilibrium @f$ \Upsilon^{-1} = \exp \left[ - \frac{\mu_B B + \mu_{I_3} I_3 + \mu_S S + \mu_C C}{T}\right] @f$.
 * Otherwise the returned value is just simply @f$ \Upsilon^{-1} = \frac{1}{\Upsilon} @f$, see GetFugacity().
 * @param [in] aPartType type of particle
 * @param [in] aX x-coordinate
 * @param [in] aY y-coordinate
 * @param [in] aZ z-coordinate
 * 
 * @fn double Thermodynamics::GetViscosityCorrection(double aX, double aY, double aZ, double aPe, double aPx, double aPy, double aPz, double aPdotU) const
 * @brief Returns the viscosity correction @f$ \Delta_{vis} @f$ to the distribution function. Type of viscosity correction terms depend on the GetViscosityType().
 * 
 * Value of viscosity correction:
 * @f{equation}{
 *   \Delta_{vis} = \Delta_{bulk} + \Delta_{shear} \nonumber
 * @f}
 * where the bulk viscosity correction term equals
 * @f{equation}{
 *   \Delta_{bulk} = \left( c_s^2 \,\, p_\mu u^\mu + \frac{m^2 - (p_\mu u^\mu)^2}{3 \, p_\mu u^\mu}\right) \frac{\tau_{bulk}}{T} \frac{\Pi}{\zeta} \nonumber
 * @f}
 * and the shear viscosity correction term is equal to
 * @f{equation}{
 *   \Delta_{shear} = \frac{1}{2 T^2}\, p_\mu p_\nu \frac{\pi^{\mu\nu}}{s T} \nonumber
 * @f}
 * The distribution function @f$ f_0 = \left[ \exp \left( \frac{p_\mu u^\mu - \mu_{ch}}{T}\right) \pm 1 \right]^{-1} @f$ with the viscosity correction is equal to
 * @f{equation}{
 *   f = f_0 \left[ 1 + \left(1 \pm f_0\right) \, \Delta_{vis} \right] \nonumber
 * @f}
 * where the (+) sign is for the Fermi-Dirac distribution and (-) sign is for the Bose-Einstein distribution.
 */
