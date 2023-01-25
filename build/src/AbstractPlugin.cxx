#include "THGlobal.h"
#include "Configurator.h"
#include "StructEvent.h"
#include "AbstractPlugin.h"
#include "Event.h"
#include <fstream>

using namespace std;

extern void AddLogEntry(const char* aEntry);
extern Configurator *sMainConfig;
extern TString	sEventDIR;
extern int      sParentPID;                                                                                                                                                                

AbstractPlugin::AbstractPlugin() : m_bSkip(false) {
    ReadParameters();
}
AbstractPlugin::~AbstractPlugin() {}

void AbstractPlugin::ReadParameters() {}

