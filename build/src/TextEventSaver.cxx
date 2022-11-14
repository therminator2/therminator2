#include "THGlobal.h"
#include "TextEventSaver.h"
#include "TDatime.h"
#include <fstream>

using namespace std;

extern TString	sEventDIR;

void TextEventSaver::Save(Event *tEvent, Model *, int)
{
  TDatime  tDate;
  ofstream tFile;
  list<Particle>::iterator tParticle_i;
  
  tFile.open((sEventDIR + "event.txt").Data(), ios_base::app); 
  if((tFile) && (tFile.is_open())) {
    if (static_cast<long>(tFile.tellp()) == 0) {
      tDate.Set();
      PRINT_MESSAGE("\n["<<tDate.AsSQLString() << "]\tFile "<<(sEventDIR + "event.txt").Data()<<" created.");
      tFile << "# THERMINATOR 2 text output" << endl;
      tFile << "#<EVENT_ENTRY>\tEID\tfatherEID\tPID\tfatherPID\trootPID\tdecayed\tmass\tE\tp_x\tp_y\tp_z\tt\tx\ty\tz\t</EVENT_ENTRY>" << endl;
    }
    tFile << "#<EVENT_ID>\t0x"<< hex << uppercase << tEvent->GetEventID() << nouppercase << dec <<"\t</EVENT_ID>"<<endl;
    tFile << "#<NO_OF_PARTICLES>\t"<<tEvent->GetParticleList()->size() <<"\t</NO_OF_PARTICLES>"<<endl;
    for(tParticle_i = tEvent->GetParticleList()->begin(); tParticle_i != tEvent->GetParticleList()->end(); tParticle_i++)
      tFile << tParticle_i->MakeTEXTEntry();
    PRINT_DEBUG_2("<EventGenerator::SaveAsText>\tEventID "<<tEvent->GetEventID()<<" saved.");
    tFile.close();
  } else {
    PRINT_MESSAGE("<EventGenerator::SaveAsText>\tUnable to create file "<<(sEventDIR + "event.txt").Data());
    exit(_ERROR_GENERAL_FILE_NOT_FOUND_);
  }
}


