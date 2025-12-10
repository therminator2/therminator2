#ifndef URQMD_EVENT_SAVER
#define URQMD_EVENT_SAVER

#include "AbstractEventSaver.h"

class UrQMDEventSaver : public AbstractEventSaver {
  public:
    virtual void Save(Event *, Model *, int);
};

#endif
