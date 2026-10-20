#ifndef ABSTRACT_AFTERBURNER
#define ABSTRACT_AFTERBURNER

#include "Event.h"

class AbstractAfterburner {
  public:
    AbstractAfterburner();
    virtual ~AbstractAfterburner();
    virtual void Init() = 0;
    virtual void Apply(Event *) = 0;
    virtual void Finalize() = 0;

  protected:
    void ReadParameters();  
};

#endif
