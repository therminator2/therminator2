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

#ifndef _TH2_STRUCT_XML_H_
  #define _TH2_STRUCT_XML_H_

#include <list>
#include "TString.h"

struct XMLattribute {
  TString name;
  TString value;
};

struct XMLtag {
  TString  name;
  std::list<XMLattribute> attributes;
  long int begin;
  long int end;
  XMLtag*  father;
  XMLtag*  child;
  XMLtag*  prev;
  XMLtag*  next;
};

#endif

/*! @file StructXML.h
 * @brief Definition of XMLtag structure. XML type entry.
 */
/*! @struct XMLtag
 * @brief Structure of an XML tag entry.
 *
 * @var TString XMLtag::name
 * @brief XML tag name
 *
 * @var std::list<XMLattribute> XMLtag::attributes
 * @brief list of XML tag attributes
 *
 * @var XMLtag* XMLtag::father
 * @brief pointer to father XML tag 
 *
 * @var XMLtag* XMLtag::child
 * @brief pointer to child XML tag
 *
 * @var XMLtag* XMLtag::next
 * @brief pointer to next sibling XML tag
 *
 * @var XMLtag* XMLtag::prev
 * @brief pointer to previous sibling XML tag
 *
 * @var long int XMLtag::begin
 * @brief begin of XML tag content - in file
 *
 * @var long int XMLtag::end
 * @brief end of XML tag content - in file
 */
/*! @struct XMLattribute
 * @brief Structure of an XML tag attribute entry.
 *
 * @var TString XMLattribute::name
 * @brief XML tag attribute name
 *
 * @var TString XMLattribute::value
 * @brief XML tag attribute value
 */
