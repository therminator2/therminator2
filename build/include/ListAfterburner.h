#ifndef LIST_AFTERBURNER
#define LIST_AFTERBURNER

#include "AbstractAfterburner.h"
#include <vector>

class ListAfterburner : public AbstractAfterburner {
  public:
    virtual ~ListAfterburner();
    virtual void Apply(Event *);
    void Add(AbstractAfterburner *afterburner);
  private:
    std::vector<AbstractAfterburner *> afterburners;
};

#endif
