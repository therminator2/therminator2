#ifndef ROOT_EVENT_SAVER
#define ROOT_EVENT_SAVER

#include "AbstractEventSaver.h"
#include "TFile.h"
#include "TTree.h"

class RootEventSaver : public AbstractEventSaver {
  public:
    RootEventSaver();
    virtual void Save(Event *, Model *, int);
  private:
    TFile*	        mFile;
    TTree*	        mParameterTree;
    TTree*        	mEventTree;
    TTree*	        mParticleTree;    
};

#endif
