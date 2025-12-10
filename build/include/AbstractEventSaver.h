#ifndef ABSTRACT_EVENT_SAVER
#define ABSTRACT_EVENT_SAVER

#include "Event.h"

class AbstractEventSaver {
  public:
    AbstractEventSaver();
    virtual ~AbstractEventSaver();
    virtual void Save(Event *, Model *, int) = 0;
    void FindPreviousEventFiles();
    void SetEventsTemp();

  protected:
    void ReadParameters();  
    int	        mFileCounter;
    const int   kEventsPerFile;
    int	        mNumberOfEvents;

};

#endif
