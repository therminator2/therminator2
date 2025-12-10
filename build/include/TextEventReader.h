#ifndef TEXT_EVENT_READER
#define TEXT_EVENT_READER

#include "AbstractEventReader.h"
#include <fstream>

class TextEventReader : public AbstractEventReader {

public:
    TextEventReader(ParticleDB* aDB, AbstractEventSaver *aES, ListAfterburner *aLA) : AbstractEventReader(aDB, aES, aLA) {}

    private:
        virtual void SetupInput(TString aEventDir, Int_t aEventFiles);
        virtual bool NextEvent();
        virtual void GetParticle(int aIter);

        std::ifstream mInput;



};

#endif //TEXT_EVENT_READER
