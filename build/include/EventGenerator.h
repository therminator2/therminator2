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

#ifndef _TH2_EVENT_GENERATOR_H_
  #define _TH2_EVENT_GENERATOR_H_
  
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include "ParticleDB.h"
#include "Integrator.h"
#include "AbstractEventSaver.h"
#include "Event.h"

class EventGenerator {
  public:
    EventGenerator();
    EventGenerator(ParticleDB* aDB, AbstractEventSaver *aES);
    ~EventGenerator();

    void   GenerateEvents();

  private:
    void ReadParameters();  
    void SaveEvent();

    AbstractEventSaver* mEventSaver;
    ParticleDB*	        mDB;
    Integrator*	        mInteg;
    Event*      	mEvent;
    int		        mEventCounter;
    int		        mNumberOfEvents;
};

#endif

/*! @file EventGenerator.h
 * @brief Definition of EventGenerator class. Generates events and saves them to files.
 */
/*!
 * @def _EVENTS_PER_FILE_
 * @brief Macro with the number of events to be stored per event*.root file.
 */
/*! @class EventGenerator
 * @brief Generates number of Event and saves them to ROOT type or text files.
 * 
 * The <c><b>THERMINATOR 2</b></c> main settings file <b>therminator.ini</b> options used here:
 * <table>
 *   <tr><th>Keyword</th>		<th>Description</th></tr>
 *   <tr><td>NumberOfEvents</td>	<td>number of events to be generated</td></tr>
 *   <tr><td>IntegrateSamples</td>	<td>number of Monte-Carlo samples used by Integrator</td></tr>
 *   <tr><td>EventDir</td>		<td>general director to store events, subdirectory is added by Model</td></tr>
 *   <tr><td>EventFileType</td>		<td>type of files to save events 
 * 					    <table>
 * 					      <tr><td><b>"root"</b></td><td>ROOT file : <b>eventXXX.root [default]</b></td></tr>
 * 					      <tr><td>"text"</td>	<td>TEXT file : <b>event.txt</b></td></tr>
 * 					      <tr><td>"root&text"</td>	<td>both ROOT and TEXT files</td></tr>
 * 					    </table>
 * 					</td></tr>
 *   <tr><td>LogFile</td>		<td>name and location of the LOG file</td></tr>
 * </table>
 * 
 * EventGenerator creates an Event object. For a given number of events the Event is called to generate primordial particles 
 * and to decay resonances. The Particle list is than saved to file, either a ROOT type file or a TEXT file or both.
 * 
 * In case the event directory already has some <b>eventXXX.root</b> EventGenerator save new files with the name counter 
 * set to continue the number sequence.<br />
 * 
 * EventGenerator also create a log entry to a <c>THERMINATOR 2</c> log file [default: <em>therminator.log</em>]. The log entry 
 * consists of:
 * <ul>
 *   <li> current time and date,
 *   <li> Model number (sModel)
 *   <li> directory the events are stored
 *   <li> initial event file number in this program run
 *   <li> last event file number in this program run  
 * </ul>
 *
 * @fn EventGenerator::EventGenerator()
 * @brief Default constructor.
 *
 * @fn EventGenerator::EventGenerator(ParticleDB* aDB)
 * @brief Generates events with a particles form a given database.
 * @param [in] aDB pointer to ParticleDB
 *
 * @fn EventGenerator::~EventGenerator()
 * @brief Destructor.
 *
 * @fn void EventGenerator::GenerateEvents()
 * @brief Creates a given number of events and saves them to a set type of file.
 *
 * @fn void EventGenerator::SetEventsTemp()
 * @brief Creates a temporary file in the main <c><b>THERMINATOR 2</b></c> directory called <b>event_<aPID>.tmp</b>.
 * 
 * Default name is <b>event_0.tmp</b>.The file contains two lines:
 * <ul>
 *   <li> current directory where the event files are stored,
 *   <li> number of event files located at that directory, including previously generate files.
 * </ul>
 */
