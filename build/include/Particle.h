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

#ifndef _TH2_PARTICLE_H_
  #define _TH2_PARTICLE_H_

#include "ParticleType.h"
#include "ParticleCoor.h"


class Particle
: public ParticleCoor
{
  public:
    Particle();
    Particle(ParticleType* aType);
    ~Particle();

    ParticleCoor& operator= (Particle& aParticle);
    ParticleType* GetParticleType() const;
    void	  SetParticlePX(double  aPe, double  aPx, double  aPy, double  aPz,
				double  aXt, double  aXx, double  aXy, double  aXz,
				Particle* aFather=0);
    const char*	  MakeTEXTEntry();
    static void   ZeroEID();
    
  private:    
    static long int EIDi;
    ParticleType*   mPartType;
};

#endif

/*! @file Particle.h
 * @brief Definition of Particle class. %Particle space-time and energy-momentum coordinates.
 */
/*! @class Particle
 * @brief Class containing particle of ParticleType space-time and energy-momentum coordinates. 
 * 
 * Class has information about the PDG Monte-Carlo number of the ParticleType, father particle and the original primordial particle
 * To back trace the whole decay cascade Particle stores the unique particle number of the particle in event and its fathers.
 * Particle also holds information if the particle has decayed.
 *
 * @fn Particle::Particle()
 * @brief Default constructor.
 *
 * @fn Particle::Particle(ParticleType* aType)
 * @brief Initiate particle of ParticleType.
 * @param [in] aType type of particle
 *
 * @fn Particle::~Particle()
 * @brief Destructor.
 *
 * @fn ParticleCoor& Particle::operator= (Particle& aParticle)
 * @brief Assignment operator.
 * @param [in] aParticle referance to Particle
 * @retval ParticleCoor& referance to ParticleCoor
 * 
 * @fn ParticleType* Particle::GetParticleType() const
 * @brief Returns a pointer to the ParticleType.
 *
 * @fn void Particle::SetParticlePX(double aPe, double aPx, double aPy, double aPz, double aXt, double aXx, double aXy, double aXz, Particle* aFather=0)
 * @brief Set particle energy, momentum, time, position and father IDs if applicable
 * @param [in] aPe particle energy
 * @param [in] aPx particle px-coordinate
 * @param [in] aPy particle py-coordinate
 * @param [in] aPz particle pz-coordinate
 * @param [in] aXt particle time
 * @param [in] aXx particle x-coordinate
 * @param [in] aXy particle y-coordinate
 * @param [in] aXz particle z-coordinate
 * @param [in] aFather pointer to the father Particle
 * 
 * @fn const char* Particle::MakeTEXTEntry()
 * @brief Returns a text string with all Particle properties.
 * 
 * @fn static void Particle::ZeroEID()
 * @brief Resets the internal particle counter to zero. Used in ParticleCoor::eid.
 */
