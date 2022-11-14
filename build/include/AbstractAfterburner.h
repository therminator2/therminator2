#ifndef ABSTRACT_AFTERBURNER
#define ABSTRACT_AFTERBURNER

#include "Event.h"

class AbstractAfterburner {
  public:
    AbstractAfterburner();
    virtual ~AbstractAfterburner();
    virtual void Apply(Event *) = 0;

  protected:
    void ReadParameters();  

};

#endif
