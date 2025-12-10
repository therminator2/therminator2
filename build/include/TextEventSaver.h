#ifndef TEXT_EVENT_SAVER
#define TEXT_EVENT_SAVER

#include "AbstractEventSaver.h"

class TextEventSaver : public AbstractEventSaver {
  public:
    TextEventSaver(std::string aFileName = "event.txt") : mFileName(aFileName) { }
    virtual void Save(Event *, Model *, int);
  private:
    std::string mFileName;
};

#endif
