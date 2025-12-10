#include "THGlobal.h"
#include "UrQMDEventSaver.h"
#include "TDatime.h"
#include <fstream>
#include "Configurator.h"
#include "Model_SR.h"

using namespace std;

extern TString	sEventDIR;
extern Configurator* sMainConfig;

void UrQMDEventSaver::Save(Event *tEvent, Model *tModel, int tEventCounter)
{
  TDatime  tDate;
  ofstream tFile;
  
  tFile.open((sEventDIR + "events.f10").Data(), ios_base::app); 
  if((tFile) && (tFile.is_open())) {
    tDate.Set();
    char buff[256];
    double tT0;
    Model_SR *tModel_SR = dynamic_cast<Model_SR*>(tModel);
    if (tModel_SR != nullptr) {
      tT0 = tModel_SR->GetT0();
    }
    else {
      tT0 = tEvent->GetParticleList()->begin()->t;
    }
    tFile << "UQMD   version:       30400   1000  30400  output_file  14" << endl;
    tFile << "projectile:  (mass, char)  197  79   target:  (mass, char)  197  79" << endl;
    tFile << "transformation betas (NN,lab,pro)     0.0000000  0.6292914 -0.6292914" << endl;
    tFile << "impact_parameter_real/min/max(fm):    0.00  0.00  0.00  total_cross_section(mbarn):    3141.59" << endl;
    tFile << "equation_of_state:    0  E_lab(GeV/u): 0.1230E+01  sqrt(s)(GeV): 0.2414E+01  p_lab(GeV/u): 0.1955E+01" << endl;
    sprintf(buff, "event#%10i random seed:%12i (fixed)  total_time(fm/c):%9i Delta(t)_O(fm/c):%13.3f",
      tEventCounter, sMainConfig->GetParameter("Seed").Atoi(), (int)tT0, tT0);
    tFile << buff << endl;
      
    tFile << "op  0    0    0    0    1    0    0    0    0    0    0    0    0    0    0" << endl;
    tFile << "op  0    0    0    0    0    0    1    0    1 *  0    0    0    0    2    1" << endl;
    tFile << "op  0    0    0    1    1    0    0    0    0    0    1 *  0    0    1    0" << endl;
    tFile << "op  0    5    0    0    0    0    0    0    0    0    0    1    0    1    0" << endl;
    tFile << "pa 0.1000E+01   0.5200E+00   0.2000E+01   0.3000E+00   0.0000E+00   0.3700E+00   0.0000E+00   0.9300E-01   0.3500E+00   0.2500E+00   0.0000E+00   0.5000E+00" << endl;
    tFile << "pa 0.2700E+00   0.4900E+00   0.2700E+00   0.1000E+01   0.1600E+01   0.8500E+00   0.1550E+01   0.0000E+00   0.0000E+00   0.0000E+00   0.0000E+00   0.0000E+00" << endl;
    tFile << "pa 0.9000E+00   0.5000E+02   0.1000E+01   0.1000E+01   0.1000E+01   0.1500E+01   0.1600E+01   0.0000E+00   0.2500E+01   0.1000E+00   0.3000E+01   0.2750E+00" << endl;
    tFile << "pa 0.4200E+00   0.1080E+01   0.8000E+00   0.5000E+00   0.0000E+00   0.5500E+00   0.5000E+01   0.8000E+00   0.5000E+00   0.8000E+06   0.1000E+01   0.2000E+01" << endl;
    tFile << "pa 0.5000E+00   0.1000E+01   0.1000E+01   0.1000E+01   0.1000E+01   0.1000E+01   0.1000E+01   0.1000E+01   0.1000E+01   0.1000E+01   0.7000E+00   0.3000E+01" << endl;
    tFile << "pa 0.2000E+00   0.2000E+03   0.1000E+01   0.5000E+01   0.1000E+01   0.1000E+11   0.1000E+01   0.1000E+01   0.0000E+00   0.1000E+11   0.2000E+01   0.5500E+00" << endl;

    tFile << "pvec: r0              rx              ry              rz              p0              px              py              pz              m          ityp 2i3 chg lcl#  ncl or " << endl;

    std::vector<TString> tBaryonEntries;
    std::vector<TString> aMesonEntries;

    int tNGoodParticles = 0;
    int tNDecaysTotal = 0;
    int tNResonances = 0;
    int tNResonancesFromResonances = 0;

    for(Particle tParticle : *(tEvent->GetParticleList())) {
      int aUrQMDCode = tParticle.GetParticleType()->GetUrQMDCode();
      bool tIsBaryonResonance = TMath::Abs(aUrQMDCode) > 1 && TMath::Abs(aUrQMDCode) < 100;
      if (tIsBaryonResonance) {
        ++tNResonances;
        if (tParticle.eid != tParticle.fathereid) {
          ++tNResonancesFromResonances;
        }
      }
      if (tParticle.decayed) {
        ++tNDecaysTotal;
      }
      TString tEntry = tParticle.MakeUrQMDEntry();
      if (!tParticle.decayed && aUrQMDCode != 0 && !tEntry.Contains("NAN")) {
        if (aUrQMDCode < 100) {
          tBaryonEntries.push_back(tEntry);
        }
        else
        {
          aMesonEntries.push_back(tEntry);
        }
        ++tNGoodParticles;
      }
    }
    sprintf(buff, "%12i%12i", tNGoodParticles, (int)tT0);
    tFile << buff << endl;
    sprintf(buff, "%8i%8i%8i%8i%8i%8i%8i%8i",
      0, 0, 0, 0, tNDecaysTotal, tNResonances, 0, tNResonancesFromResonances);
    tFile << buff << endl;  

    for (TString aEntry : tBaryonEntries) {
      tFile << aEntry.Data() << endl;
    }
    for (TString aEntry : aMesonEntries) {
      tFile << aEntry.Data() << endl;
    }

    tFile.close();
  } else {
    PRINT_MESSAGE("<EventGenerator::SaveAsText>\tUnable to create file "<<(sEventDIR + "event.txt").Data());
    exit(_ERROR_GENERAL_FILE_NOT_FOUND_);
  }
}


