#ifndef LIST_PLUGIN
#define LIST_PLUGIN

#include "AbstractPlugin.h"
#include <vector>

class ListPlugin : public AbstractPlugin {
  public:
    virtual ~ListPlugin();
    virtual void Apply(Event *);
    void Add(AbstractPlugin *Plugin);
  private:
    std::vector<AbstractPlugin *> Plugins;
};

#endif
