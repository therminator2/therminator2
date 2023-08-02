#include "THGlobal.h"
#include "HGeantEventSaver.h"
#include "TDatime.h"
#include <fstream>

using namespace std;

extern TString	sEventDIR;

HGeantEventSaver::HGeantEventSaver() : AbstractEventSaver()
{
    // there is probably a better way of doing this, but it will stay like this for now - JJ
    PIDmap[22] = 1; // gamma
    PIDmap[-11] = 1; // positron
    PIDmap[11] = 1; // electron
    PIDmap[12] = 1; // neutrino
    PIDmap[-13] = 1; // muon +
    PIDmap[13] = 1; // muon -
    PIDmap[111] = 1; // pion 0
    PIDmap[211] = 1; // pion +
    PIDmap[-211] = 1; // pion -
    PIDmap[130] = 1; // kaon 0 long
    PIDmap[321] = 1; // kaon +
    PIDmap[-321] = 12; // kaon +
    PIDmap[2112] = 13; // neutron
    PIDmap[2212] = 14; // proton
    PIDmap[-2212] = 15; // antiproton
    PIDmap[310] = 16; // kaon 0 short
    PIDmap[221] = 17; // eta
    PIDmap[3122] = 18; // lambda
    PIDmap[3222] = 19; // sigma +
    PIDmap[3212] = 20; // sigma 0
    PIDmap[3112] = 21; // sigma -
    PIDmap[3322] = 22; // xi 0
    PIDmap[3312] = 23; // xi -
    PIDmap[3334] = 24; // omega -
    PIDmap[-2112] = 25; // antineutron
    PIDmap[-3122] = 26; // antilambda
    PIDmap[-3112] = 27; // antisigma -
    PIDmap[-3212] = 28; // antisigma 0
    PIDmap[-3222] = 29; // antisigma +
    PIDmap[-3322] = 30; // antixi 0
    PIDmap[-3312] = 31; // antixi +
    PIDmap[-3334] = 32; // antiomega +
    PIDmap[80000] = 45; // deuteron
}

void HGeantEventSaver::Save(Event *tEvent, Model *, int eventCounter)
{
    TDatime  tDate;
    TString tName = "event.evt";
    ofstream tFile;
    list<Particle>::iterator tParticle_i;
    std::vector<TString> tTracks;
    int tGeantPID;
    float tE, tPx, tPy, tPz;
    
    tFile.open((sEventDIR + tName).Data(), ios_base::app); 
    if((tFile) && (tFile.is_open())) 
    {
        if (static_cast<long>(tFile.tellp()) == 0) 
        {
            tDate.Set();
            PRINT_MESSAGE("\n["<<tDate.AsSQLString() << "]\tFile "<<(sEventDIR + tName).Data()<<" created.");
        }

        for(tParticle_i = tEvent->GetParticleList()->begin(); tParticle_i != tEvent->GetParticleList()->end(); tParticle_i++)
        {
            tGeantPID = GetGeantID((*tParticle_i).pid);
            tE = (*tParticle_i).e;
            tPx = (*tParticle_i).px;
            tPy = (*tParticle_i).py;
            tPz = (*tParticle_i).pz;
            
            if (tGeantPID == 0 || std::isnan(tE) || std::isnan(tPx) || std::isnan(tPy) || std::isnan(tPz))
                continue;

            tTracks.push_back(Form("%f\t%f\t%f\t%f\t%d\t1\t1\t1", tE, tPx, tPy, tPz, tGeantPID));
        }

        tFile << eventCounter << "\t" << tTracks.size() << "\t2.4\t0.0\t2" << endl; // fixed value of beam energy, b parameter and iflag respectively - JJ
        for (auto trackIter : tTracks)
            tFile << trackIter << endl;
        
        PRINT_DEBUG_2("<EventGenerator::SaveAsGeant>\tEventID "<<tEvent->GetEventID()<<" saved.");
        tFile.close();
    } 
    else
    {
        PRINT_MESSAGE("<EventGenerator::SaveAsGeant>\tUnable to create file "<<(sEventDIR + tName).Data());
        exit(_ERROR_GENERAL_FILE_NOT_FOUND_);
    }

    tTracks.clear();
    tTracks.resize(0);
}

int HGeantEventSaver::GetGeantID(int therminatorPID)
{
    if (! PIDmap.count(therminatorPID))
        return 0;
    else
        return PIDmap[therminatorPID];
}


