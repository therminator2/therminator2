#include <iostream>
#include "Math/LorentzVector.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TStyle.h"

#include "../build/include/ParticleCoor.h"
#include "events2chain.C"
#include "/home/jedkol/Downloads/indicators/single_include/indicators/indicators.hpp"

using Vector3D = ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<float>>;

struct ChildInfo
{
    Vector3D momVec;
    int parenteid;
};

struct ParentInfo
{
    Vector3D momVec;
    int eid;
};

float GetAngle(Vector3D p1, Vector3D p2)
{
    return acos(p1.Dot(p2) / ((p1.Mag2()) * sqrt(p2.Mag2())));
}

bool CheckIfInList(const ParticleCoor &part, const std::vector<int> &partIdList, std::vector<ParentInfo> &partList)
{
    for (auto &id : partIdList)
        if (part.pid == id)
        {
            ParentInfo parent{{part.px,part.py,part.pz},part.eid};
            partList.emplace_back(parent);
            return true;
        }

    return false;
}

bool CheckIfInList(const ParticleCoor &part, const std::vector<int> &partIdList, std::vector<ChildInfo> &partList)
{
    for (auto &id : partIdList)
        if (part.fatherpid == id)
        {
            ChildInfo child{{part.px,part.py,part.pz},part.fathereid};
            partList.emplace_back(child);
            return true;
        }

    return false;
}

void FillAndClear(TH1D *hinp, std::vector<ChildInfo> &protonVec, std::vector<ParentInfo> &parentVec)
{
    if (protonVec.size() == 0)
        return;

    // assumption: parentVec >= protonVec
    for (auto &proton : protonVec)
        for (auto &parent : parentVec)
        {
            if (proton.parenteid == parent.eid)
            {
                hinp->Fill(GetAngle(proton.momVec,parent.momVec) * 180./ROOT::Math::Pi());
                break;
            }
        }

    protonVec.clear();
    protonVec.resize(0);
    parentVec.clear();
    parentVec.resize(0);
}

void figure_openangle(TString aEventDir = "./events/", Int_t aEventFiles = 1)
{
    gStyle->SetOptStat(0);

    const std::vector<int> fHe4IdList = {1000020140,1000020240,1000020340,1000020440,1000020540,1000020640,1000020740,1000020840,1000020940,1000021040,1000021140,1000021240,1000021340,1000021440,1000021540};
    const std::vector<int> fLi4IdList = {1000030040,1000030140,1000030240,1000030340};
    const std::vector<int> fLi5IdList = {1000030050,1000030150,1000030250};
    const int fProtId = 2212;

    std::vector<ChildInfo> fHe4ProtonVec, fLi4ProtonVec, fLi5ProtonVec;
    std::vector<ParentInfo> fHe4Vec, fLi4Vec, fLi5Vec;

    TH1D *hHe4Angle = new TH1D("hHe4Angle","Opening angle between {}^{4}He and p; #alpha (#circ)",180,0,180);
    TH1D *hLi4Angle = new TH1D("hLi4Angle","Opening angle between {}^{4}Li and p; #alpha (#circ)",180,0,180);
    TH1D *hLi5Angle = new TH1D("hLi5Angle","Opening angle between {}^{5}Li and p; #alpha (#circ)",180,0,180);

    static ParticleCoor fParticle;
    Int_t   fEvents;
    TChain* fChain = events2chain(aEventDir, aEventFiles, &fParticle, &fEvents);
    long int entries = fChain->GetEntries();
    int countdown = fEvents;

    indicators::show_console_cursor(false);
    indicators::BlockProgressBar bar{indicators::option::BarWidth{80}, 
                                    indicators::option::FontStyles{std::vector<indicators::FontStyle>{indicators::FontStyle::bold}}, 
                                    indicators::option::ShowElapsedTime{true}, 
                                    indicators::option::ShowRemainingTime{true}, 
                                    indicators::option::MaxProgress{entries/fEvents}};

    std::cout << "Looping over " << entries << " entries" <<  std::endl;
    for (long int i = 0; i < entries; i++)
    {
        fChain->GetEntry(i);
        if (--countdown == 0) // change progress every fEvents
        {
            bar.tick();
            countdown = fEvents;
        }

        if (fParticle.eid == 0) // if new event
        {
            FillAndClear(hHe4Angle,fHe4ProtonVec,fHe4Vec);
            FillAndClear(hLi4Angle,fLi4ProtonVec,fLi4Vec);
            FillAndClear(hLi5Angle,fLi5ProtonVec,fLi5Vec);
        }
        
        // first check if parent exists
        if (CheckIfInList(fParticle,fHe4IdList,fHe4Vec))
            continue;
        if (CheckIfInList(fParticle,fLi4IdList,fLi4Vec))
            continue;
        if (CheckIfInList(fParticle,fLi5IdList,fLi5Vec))
            continue;

        if (fParticle.pid != fProtId || fParticle.fatherpid == fProtId) // now we only need non-primordial protons
            continue;

        // check if proton comes from He and Li
        if (CheckIfInList(fParticle,fHe4IdList,fHe4ProtonVec))
            continue;
        if (CheckIfInList(fParticle,fLi4IdList,fLi4ProtonVec))
            continue;
        CheckIfInList(fParticle,fLi5IdList,fLi5ProtonVec);
    }
    bar.mark_as_completed();
    // Show cursor
    indicators::show_console_cursor(true);

    TCanvas *canv = new TCanvas("canv","",1800,600);
    canv->Divide(3,1);

    canv->cd(1);
    hHe4Angle->Draw();
    canv->cd(2);
    hLi4Angle->Draw();
    canv->cd(3);
    hLi5Angle->Draw();
}