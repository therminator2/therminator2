#ifndef ROOT_EVENT_READER
#define ROOT_EVENT_READER

#include "AbstractEventReader.h"

class RootEventReader : public AbstractEventReader {

public:
    RootEventReader(ParticleDB* aDB, AbstractEventSaver *aES, ListAfterburner *aLA) : AbstractEventReader(aDB, aES, aLA) {}

    private:
        virtual void SetupInput(TString aEventDir, Int_t aEventFiles);
        virtual bool NextEvent();
        virtual void GetParticle(int aIter);

        TChain* mChainParts;
        TChain* mChainEvent;



};

#endif //ROOT_EVENT_READER
