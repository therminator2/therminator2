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

#ifndef _TH2_CONFIGURATOR_H_
  #define _TH2_CONFIGURATOR_H_

#include <vector>
#include <TString.h>

struct Parameter
{
  TString keyword;
  TString value;
};

class Configurator {

  public:
    Configurator();
    ~Configurator();

    TString GetParameter(const char* aKeyword) noexcept(false);
    bool    HasParameter(const char* aKeyword) noexcept(true);
    void    AddParameter(Parameter* aPar);
    int     PrintParameters();

  private:    
    std::vector<Parameter> mParameters;	
};

#endif

/*! @file Configurator.h
 * @brief Definition of Configurator class.
 */
/*!@struct Parameter
 * @brief Parameter structure
 * 
 * @var Parameter::keyword
 * @brief parameter keyword
 * 
 * @var Parameter::value
 * @brief parameter value
 */
/*! @class Configurator
 * @brief Class contains parameters for <c><b>THERMINATOR 2</b></c> that are stored in the <em>ini</em> type files.
 * 
 * @ref file_type_ini
 * and extracts settings information. 
 * Returns the parameters value that matches the keyword.
 *
 * @fn Configurator::Configurator()
 * @brief Default constructor.
 *
 * @fn Configurator::~Configurator()
 * @brief Destructor.
 *
 * @fn TString Configurator::GetParameter(const char* aKeyword) throw(TString)
 * @brief Returns the value of the given keyword.
 * 
 * If parameter was not found in the file the parameter keyword will be thrown.
 * @param [in] aKeyword keyword
 * @exception TString requested keyword
 *
 * @fn void Configurator::AddParameter(Parameter* aPar)
 * @brief Adds a new parameter to the list
 * @param [in] aPar pointer to Parameter structure
 * 
 * @fn int Configurator::PrintParameters()
 * @brief Prints all of the stored parameters keywords and its values.
 * Requires compilation variable <c>DEBUG_LEVEL_</c> set to 3.
 */
