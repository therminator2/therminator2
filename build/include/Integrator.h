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

#ifndef _TH2_INTEGRATOR_H_
  #define _TH2_INTEGRATOR_H_

#include <list>
#include <TRandom2.h>
#include "ParticleDB.h"
#include "ParticleType.h"
#include "Particle.h"
#include "Model.h"

class Integrator {
  public:
    Integrator();
    Integrator(int aNSamples);
    ~Integrator();
    
    Model* GetModel();
    void   GenerateParticles(ParticleType* aPartType, int aPartCount, std::list<Particle>* aParticles);
    void   SetMultiplicities(ParticleDB* aDB);
    void   Randomize();
    
  private:
    double Integrate(ParticleType* aPartType);
    
    Model*    mFOModel;
    TRandom2* mRandom;
    int       mNSamples;
};

#endif

/*! @file Integrator.h
 * @brief Definition of Integrator class. Integrates the Cooper-Frye formula and generates primordial particles.
 */
/*! @class Integrator
 * @brief Integrates the Cooper-Frye formula using a Monte-Carlo method 
 * and randomly generates primordial particles on a given freeze-out hypersurface.
 *
 * @fn Integrator::Integrator()
 * @brief Default constructor.
 *
 * @fn Integrator::Integrator(int aNSamples)
 * @brief Sets the number of Monte-Carlo samples to determine the maximal value of the integrand and average multiplicity
 * of the primordial particles.
 * @param [in] aNSamples number of Monte-Carlo samples
 *
 * @fn Integrator::~Integrator()
 * @brief Destructor
 *
 * @fn Model* Integrator::GetModel()
 * @brief Returns a pointer to the currently used Model object.
 *
 * @fn void Integrator::GenerateParticles(ParticleType* aPartType, int aPartCount, std::list<Particle>* aParticles)
 * @brief Generates a given number of primordial particles of a given type and stores them in the list-type container.
 * @param [in] aPartType type of particle
 * @param [in] aPartCount number of particles to be generated
 * @param [out] aParticles std::list container object that stores all generated particles
 *
 * @fn void Integrator::SetMultiplicities(ParticleDB* aDB)
 * @brief Set the particle database with appropriate average multiplicity and maximum value of the integrand.
 * 
 * Based on the model parameter a CRC32 hash is created. If a file with such hash exists then the multiplicities and 
 * max integrand data is taken from the file. Otherwise these numbers are calculate with Monte-Carlo with the number of 
 * samples defined.
 * @param [in] aDB particle data-base. See ParticleDB class.
 *
 * @fn void Integrator::Randomize()
 * @brief Sets a new seed for the random number generator based on the system time.
 */
