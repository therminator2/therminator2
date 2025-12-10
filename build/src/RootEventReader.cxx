#include "RootEventReader.h"

void RootEventReader::SetupInput(TString aEventDir, Int_t aEventFiles) {
{
// ##########################################################################
// WHICH ROOT FILES TO READ
// ##########################################################################
  mChainParts = new TChain(_PARTICLES_TREE_);
  mChainEvent = new TChain(_EVENTS_TREE_);

  mChainParts->SetBranchAddress(_PARTICLE_BRANCH_, &mParticleCoor);
  mChainEvent->SetBranchAddress(_EVENTS_BRANCH_,   &mStructEvent);
  
  for(Int_t i=0; i<aEventFiles; i++) {
    char Buff[100];
    sprintf(Buff,"%sevent%03i.root",aEventDir.Data(),i);
    Info("RootEventReader::SetupInput", "Adding file: %s", Buff);
    mChainParts->Add(Buff);
    mChainEvent->Add(Buff);
  }
  Long64_t tEntries = mChainEvent->GetEntries();
  if (tEntries < mNumberOfEvents) {
      Info("RootEventReader::SetupInput", "Replacing the configured number of events (%i) with actusally available (%lli)", mNumberOfEvents, tEntries);
      mNumberOfEvents = tEntries;
  }
  Info("RootEventReader::SetupInput", "Total number of events: %i", mNumberOfEvents);
}


}
bool RootEventReader::NextEvent() {
    if (mEventIter < mNumberOfEvents) {
        mChainEvent->GetEntry(mEventIter);
        mEventIter++;
        return true;
    } else {
        return false;
    }
}
void RootEventReader::GetParticle(int aIter) {
  mChainParts->GetEntry(aIter);
}


