#ifndef COLLECTION_EVENT_SAVER
#define COLLECTION_EVENT_SAVER

#include "AbstractEventSaver.h"
#include <vector>

class CollectionEventSaver : public AbstractEventSaver {
  public:
    virtual ~CollectionEventSaver();
    virtual void Save(Event *, int);
    void Add(AbstractEventSaver *saver);
  private:
    std::vector<AbstractEventSaver *> savers;
};

#endif
