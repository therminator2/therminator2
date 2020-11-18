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

#ifndef _TH2_PARSER_H_
  #define _TH2_PARSER_H_

#include <fstream>
#include <TString.h>
#include "Configurator.h"
#include "ParticleDB.h"
#include "Hypersurface_Library.h"

class Parser {
  public:
    Parser();
    Parser(const char* aFileName);
    ~Parser();

    void   ReadINI(Configurator* aINI);
    void   ReadSHAREParticles(ParticleDB* aDB);
    void   ReadSHAREDecays(ParticleDB* aDB);
    void   ReadXML(Hypersurface_Library* aXML);
    
  private:
    double SHAREDeltaJ(double aJot1, double aJot2, double aJot);
    double SHAREClebschGordan(double aJot, double aEm, double aJot1, double aEm1, double aJot2, double aEm2);

    XMLtag* XMLTagNew(TString& aBuff);
    
    TString	  mFileName;
    std::ifstream mFile;
};

#endif

/*! @file Parser.h
 * @brief Definition of Parser class. %Parser filling ParticleDB with particle types.
 */
/*! @class Parser
 * @brief A parser for the particle database files.
 *
 * @fn Parser::Parser()
 * @brief Default constructor.
 *
 * @fn Parser::Parser(const char* aFile)
 * @brief %Parser will process this text file.
 * @param [in] aFile text file
 *
 * @fn Parser::~Parser()
 * @brief Destructor.
 *
 * @fn void Parser::ReadINI(Configurator* aINI)
 * @brief Reads the INI file and fills the Configurator.
 * @param [in,out] aINI configuration
 * 
 * @fn void Parser::ReadSHAREParticles(ParticleDB* aDB)
 * @brief Reads <c><b>SHARE</b></c> particles database and fills the ParticleDB. 
 * @param [in,out] aDB particles database
 * 
 * @fn void Parser::ReadSHAREDecays(ParticleDB* aDB)
 * @brief Reads <c><b>SHARE</b></c> decay database and fills the ParticleDB. 
 * @param [in,out] aDB particles database
 * 
 * @fn void Parser::ReadXML(Hypersurface_Library* aXML)
 * @brief Reads the XML file containing the definition of the hypersurface and fills the Hypersurface_Library.
 * @param [in,out] aXML hypersurface library
 */
