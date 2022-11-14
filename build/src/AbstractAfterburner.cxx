#include "THGlobal.h"
#include "Configurator.h"
#include "StructEvent.h"
#include "AbstractAfterburner.h"
#include "Event.h"
#include <fstream>

using namespace std;

extern void AddLogEntry(const char* aEntry);
extern Configurator *sMainConfig;
extern TString	sEventDIR;
extern int      sParentPID;                                                                                                                                                                

AbstractAfterburner::AbstractAfterburner() {
    ReadParameters();
}
AbstractAfterburner::~AbstractAfterburner() {}

void AbstractAfterburner::ReadParameters() {}

