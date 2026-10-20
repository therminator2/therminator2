#include "ListAfterburner.h"

using namespace std;

ListAfterburner::~ListAfterburner() {
  vector<AbstractAfterburner*>::iterator it = afterburners.begin();
  for ( ; it != afterburners.end(); ++it) { delete *it; }
}
void ListAfterburner::Init()
{
    for (auto afterburner : afterburners) {
        afterburner->Init();
    }
}
void ListAfterburner::Apply(Event *tEvent)
{
    for (auto afterburner : afterburners) {
        afterburner->Apply(tEvent);
    }
    vector<AbstractAfterburner*>::iterator it = afterburners.begin();
    for ( ; it != afterburners.end(); ++it) { (*it)->Apply(tEvent); }
}
void ListAfterburner::Finalize()
{
    for (auto afterburner : afterburners) {
        afterburner->Finalize();
    }
}
void ListAfterburner::Add(AbstractAfterburner *afterburner) { afterburners.push_back(afterburner); }
