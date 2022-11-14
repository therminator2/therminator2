#include "ListAfterburner.h"

using namespace std;

ListAfterburner::~ListAfterburner() {
  vector<AbstractAfterburner*>::iterator it = afterburners.begin();
  for ( ; it != afterburners.end(); ++it) { delete *it; }
}
void ListAfterburner::Apply(Event *tEvent)
{
  vector<AbstractAfterburner*>::iterator it = afterburners.begin();
  for ( ; it != afterburners.end(); ++it) { (*it)->Apply(tEvent); }
}
void ListAfterburner::Add(AbstractAfterburner *afterburner) { afterburners.push_back(afterburner); }
