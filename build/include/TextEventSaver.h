#ifndef TEXT_EVENT_SAVER
#define TEXT_EVENT_SAVER

#include "AbstractEventSaver.h"

class TextEventSaver : public AbstractEventSaver {
  public:
    virtual void Save(Event *, Model *, int);
};

#endif
