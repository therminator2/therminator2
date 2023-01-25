#ifndef ABSTRACT_PLUGIN
#define ABSTRACT_PLUGIN

#include "Event.h"

class AbstractPlugin {
  public:
    AbstractPlugin();
    virtual ~AbstractPlugin();
    virtual void Apply(Event *) = 0;

  protected:
    void ReadParameters();
    bool m_bSkip;

};

#endif
