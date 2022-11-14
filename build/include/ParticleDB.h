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

#ifndef _TH2_PARTICLE_DB_H_
  #define _TH2_PARTICLE_DB_H_

#include <map>
#include <vector>
#include <TString.h>
#include "ParticleType.h"

class ParticleDB {
  public:
    ParticleDB();
    ~ParticleDB();

    int           AddParticleType(ParticleType* aPartType);
    ParticleType* GetParticleType(int aIndex);
    ParticleType* GetParticleType(TString aName);
    ParticleType* FindByPID(int aPID);
    int           GetParticleTypeIndex(TString aName);
    int           GetParticleTypeCount();
    int           ExistsParticleType(TString aName);

  private:
    std::vector<ParticleType> mParticleTable;
    std::map<TString, int>    mParticleNames;
};

#endif

/*! @file ParticleDB.h
 * @brief Definition of ParticleDB class. Database with ParticleType.
 */
/*! @class ParticleDB
 * @brief Contains a database of particle types. The database can be accessed by giving an index number of desired particle or its name.
 *
 * @fn ParticleDB::ParticleDB()
 * @brief Default constructor.
 *
 * @fn ParticleDB::~ParticleDB()
 * @brief Destructor.
 *
 * @fn int ParticleDB::AddParticleType(ParticleType* aPartType)
 * @brief Adds a particle type to the databse.
 * @param [in] aPartType pointer to particle type
 *
 * @fn ParticleType* ParticleDB::GetParticleType(int aIndex)
 * @brief Returns a pointer to ParticleType of a given index in the database
 * @param [in] aIndex particle type index in the database
 *
 * @fn ParticleType* ParticleDB::GetParticleType(TString aName)
 * @brief Returns a pointer to ParticleType of a given name in the database
 * @param [in] aName particle type name in the database
 *
 * @fn int ParticleDB::GetParticleTypeIndex(TString aName)
 * @brief Returns the index of the particle type of the given name
 * @param [in] aName particle type name
 *
 * @fn int ParticleDB::GetParticleTypeCount()
 * @brief Returns the number of particle types in the database
 *
 * @fn int ParticleDB::ExistsParticleType(TString aName)
 * @brief Checks if a particle type is in the database
 * @param [in] aName particle type name
 * @retval 0 particle type not found
 * @retval >0 number of database entries for that particle type
 */
