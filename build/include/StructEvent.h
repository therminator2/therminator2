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

#ifndef _TH2_STRUCT_EVENTS_H_
  #define _TH2_STRUCT_EVENTS_H_

#include <Rtypes.h>

#define _EVENTS_PER_FILE_ 50000
#define _EVENTS_TREE_   "events"
#define _EVENTS_BRANCH_ "event"
#define _EVENTS_FORMAT_ "eventid/i:entries:entriesprev"

struct StructEvent
{
  UInt_t eventID;
  UInt_t entries;
  UInt_t entriesprev;
};

#endif

/*! @file StructEvent.h
 * @brief Definition of StructEvent structure. Entry to the <b>events</b> TTree.
 */
/*! @def _EVENTS_PER_FILE_
 * @brief Macro with the number of events stored per file.
 * 
 * @def _EVENTS_TREE_
 * @brief Macro with the name of the <a href="http://root.cern.ch/root/html/TTree.html">TTree</a> with event information.
 * 
 * @def _EVENTS_BRANCH_
 * @brief Macro with the name of a <a href="http://root.cern.ch/root/html/TBranch.html">TBranch</a> in the TTree with StructEvent.
 * 
 * @def _EVENTS_FORMAT_
 * @brief Macro with the format of StructEvent.
 */
/*! @struct StructEvent
 * @brief Structure of the entry to the <b>events</b> TTree.
 * 
 * @var UInt_t StructEvent::eventID
 * @brief ID of the event, CRC-32
 * 
 * @var UInt_t StructEvent::entries
 * @brief number of particle entries in this event.
 * 
 * @var UInt_t StructEvent::entriesprev
 * @brief total number of previous particle entries.
 */
