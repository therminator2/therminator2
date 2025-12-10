#include "THGlobal.h"
#include "TextEventReader.h"

void TextEventReader::SetupInput(TString aEventDir, Int_t) {

    char Buff[100];
    sprintf(Buff,"%sevent.txt",aEventDir.Data());
    Info("TextEventReader::SetupInput", "Reading file: %s", Buff);
    mInput.open(Buff);
    std::string tLine;

    // Skipping two header lines
    std::getline(mInput, tLine);
    std::getline(mInput, tLine);

    Info("TextEventReader::SetupInput", "Total number of events unknown upfront for a text file.");

}

bool TextEventReader::NextEvent() {
    std::string tLine;
    // if there are more events, there should be another tLine after reading all the particles
    if (std::getline(mInput, tLine)) {
        // Needed in both lines that will be parsed
        std::string tOpenTag, tCloseTag;

        std::stringstream tLineStream(tLine);
        int tEventID;
        tLineStream >> tOpenTag >> tEventID >> tCloseTag;

        // Read another tLine to access the particle number
        std::getline(mInput, tLine);
        tLineStream = std::stringstream(tLine);
        int tNumberOfParticles;
        tLineStream >> tOpenTag >> tNumberOfParticles >> tCloseTag;

        int tEntriesPrev = mStructEvent.entries;
        mStructEvent = StructEvent();
        mStructEvent.entries = tNumberOfParticles;
        mStructEvent.entriesprev = tEntriesPrev;
        mStructEvent.eventID = tEventID;
        mEventIter++;
        return true;
    } else {
        return false;
    }
}

void TextEventReader::GetParticle(int aIter) {
    static long long int i = 0;
    std::string tLine;
    std::getline(mInput, tLine);
//    std::cerr << "tLine:" << tLine << "< " << i++ << std::endl;
/*
<EVENT_ENTRY>  EID     fatherEID       PID     fatherPID       rootPID decayed mass    E       p_x     p_y     p_z     t       x       y       z       </EVENT_ENTRY>
     0  -1      80000   80000   80000   0       1.875613e+00    2.109259e+00    -1.497327e-01   -4.404987e-01   -8.454174e-01   1.000000e+01    2.664469e+00    -3.389931e+00   -2.901068e+00
*/
    std::stringstream tLineStream(tLine);
    tLineStream
        >> mParticleCoor.eid
        >> mParticleCoor.fathereid
        >> mParticleCoor.pid
        >> mParticleCoor.fatherpid
        >> mParticleCoor.rootpid
        >> mParticleCoor.decayed
        >> mParticleCoor.mass
        >> mParticleCoor.e
        >> mParticleCoor.px
        >> mParticleCoor.py
        >> mParticleCoor.pz
        >> mParticleCoor.t
        >> mParticleCoor.x
        >> mParticleCoor.y
        >> mParticleCoor.z;
        mParticleCoor.t /= kHbarC;
        mParticleCoor.x /= kHbarC;
        mParticleCoor.y /= kHbarC;
        mParticleCoor.z /= kHbarC;
    mParticleCoor.eventid = mStructEvent.eventID;    

}
