#ifndef LIST_AFTERBURNER
#define LIST_AFTERBURNER

#include "AbstractAfterburner.h"
#include <vector>

class ListAfterburner : public AbstractAfterburner {
  public:
    virtual ~ListAfterburner();
    virtual void Init();
    virtual void Apply(Event *);
    virtual void Finalize();
    void Add(AbstractAfterburner *afterburner);
  private:
    std::vector<AbstractAfterburner *> afterburners;
};

#endif
