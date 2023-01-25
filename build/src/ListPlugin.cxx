#include "ListPlugin.h"

using namespace std;

ListPlugin::~ListPlugin() {
  vector<AbstractPlugin*>::iterator it = Plugins.begin();
  for ( ; it != Plugins.end(); ++it) { delete *it; }
}
void ListPlugin::Apply(Event *tEvent)
{
  vector<AbstractPlugin*>::iterator it = Plugins.begin();
  for ( ; it != Plugins.end(); ++it) { (*it)->Apply(tEvent); }
}
void ListPlugin::Add(AbstractPlugin *Plugin) { Plugins.push_back(Plugin); }
