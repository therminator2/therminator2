#ifndef ABSTRACT_EVENT_SAVER
#define ABSTRACT_EVENT_SAVER

#include "Event.h"

class AbstractEventSaver {
  public:
    AbstractEventSaver();
    virtual ~AbstractEventSaver();
    virtual void Save(Event *, int) = 0;
    void FindPreviousEventFiles();
    void SetEventsTemp();

  protected:
    void ReadParameters();  
    int	        mFileCounter;
    int	        mNumberOfEvents;
    const int   kEventsPerFile;

};

#endif
