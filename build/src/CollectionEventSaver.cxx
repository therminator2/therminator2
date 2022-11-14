#include "CollectionEventSaver.h"

using namespace std;

CollectionEventSaver::~CollectionEventSaver() {
  vector<AbstractEventSaver*>::iterator it = savers.begin();
  for ( ; it != savers.end(); ++it) { delete *it; }
}
void CollectionEventSaver::Save(Event *tEvent, Model *tModel, int tEventCounter)
{
  vector<AbstractEventSaver*>::iterator it = savers.begin();
  for ( ; it != savers.end(); ++it) { (*it)->Save(tEvent, tModel, tEventCounter); }
}
void CollectionEventSaver::Add(AbstractEventSaver *saver) { savers.push_back(saver); }
