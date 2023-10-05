#include <vector>

#include "Math/Math.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TStyle.h"

#include "../build/include/ParticleCoor.h"
#include "events2chain.C"
#include "/home/jedkol/Downloads/indicators/single_include/indicators/indicators.hpp"

enum class ProtonType{Undefined,Primordial,LightNuclei};

bool IsInList(ParticleCoor &part, const std::vector<int> &partList)
{
    for (auto &father : partList)
        if (part.fatherpid == father)
            return true;

    return false;
}

int GetRadiusIter(float radius,const std::vector<std::tuple<int,int> > radList)
{
    for (int iter = 0; iter < radList.size(); ++iter)
        if (radius > std::get<0>(radList.at(iter)) && radius <= std::get<1>(radList.at(iter)))
            return iter;

    return -1;
}

void figure_ptfreeze(TString aEventDir = "/home/jedkol/lustre/hades/user/kjedrzej/CaseE/H160E0D3femto/", Int_t aEventFiles = 10)
{
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kPastel);

    const std::vector<int> fHe4IdList = {1000020140,1000020240,1000020340,1000020440,1000020540,1000020640,1000020740,1000020840,1000020940,1000021040,1000021140,1000021240,1000021340,1000021440,1000021540};
    const std::vector<int> fLi4IdList = {1000030040,1000030140,1000030240,1000030340};
    const std::vector<int> fLi5IdList = {1000030050,1000030150,1000030250};
    const std::vector<std::tuple<int,int> > fRlist = {{0,5},{5,10},{10,15},{15,20},{20,25},{25,30},{30,35},{35,40},{40,45},{45,50}};
    const size_t fRbins = fRlist.size();
    const int fProtId = 2212, nBins = 60;
    const float fRapWidth = 0.1, binMin = 0, binMax = 3.;
    const float fPtWidth = (binMax-binMin)/nBins;

    static ParticleCoor fParticle;
    Int_t   fEvents;
    TChain* fChain = events2chain(aEventDir, aEventFiles, &fParticle, &fEvents);
    long int entries = fChain->GetEntries();
    int fRadIter, countdown = fEvents;
    float fPt, fRap, fR;
    ProtonType fDecayMode;
    std::vector<TH1D*> hPrimVec, hDecayVec;
    TLegend *fLegend;

    for (size_t iter = 0; iter < fRlist.size(); ++iter)
    {
        hPrimVec.emplace_back(new TH1D(Form("hPrim%d",iter),Form("p_{T} of primordial protons for R #in (%d,%d) fm;p_{T} (GeV/c); dN/(2 #pi p_{T} dp_{T} dy)",std::get<0>(fRlist.at(iter)),std::get<1>(fRlist.at(iter))),nBins,binMin,binMax));
        hDecayVec.emplace_back(new TH1D(Form("hDecay%d",iter),Form("p_{T} of dacay protons for R #in (%d,%d) fm;p_{T} (GeV/c); dN/(2 #pi p_{T} dp_{T} dy)",std::get<0>(fRlist.at(iter)),std::get<1>(fRlist.at(iter))),nBins,binMin,binMax));
    }

    fLegend = new TLegend(0.6,0.5,0.89,0.89);
    fLegend->SetBorderSize(0);

    indicators::show_console_cursor(false);
    indicators::BlockProgressBar bar{indicators::option::BarWidth{80}, 
                                    indicators::option::FontStyles{std::vector<indicators::FontStyle>{indicators::FontStyle::bold}}, 
                                    indicators::option::ShowElapsedTime{true}, 
                                    indicators::option::ShowRemainingTime{true}, 
                                    indicators::option::MaxProgress{entries/fEvents}};

    std::cout << "Looping over " << entries << " entries" <<  std::endl;
    for (long int i = 0; i < entries; i++)
    {
        fDecayMode = ProtonType::Undefined;

        fChain->GetEntry(i);
        if (--countdown == 0) // change progress every fEvents
        {
            bar.tick();
            countdown = fEvents;
        }

        if (fParticle.pid != fProtId) // if not proton
            continue;
        if(fParticle.e == fParticle.pz) // if inf rapidity
            continue;

        fRap = 0.5 * log((fParticle.e + fParticle.pz) / (fParticle.e - fParticle.pz));
        if( abs(fRap) >= (fRapWidth/2.) ) // if not midrap
            continue;

        if (fParticle.fatherpid == fProtId)
            fDecayMode = ProtonType::Primordial;
        else if (IsInList(fParticle,fHe4IdList) || IsInList(fParticle,fLi4IdList) || IsInList(fParticle,fLi5IdList))
            fDecayMode = ProtonType::LightNuclei;
        else // if not desired decay
            continue;

        fR = sqrt(fParticle.x*fParticle.x + fParticle.y*fParticle.y + fParticle.z*fParticle.x);
        fRadIter = GetRadiusIter(fR,fRlist);
        if (fRadIter == -1) // if radius not in bounds
            continue;

        fPt  = sqrt(fParticle.px * fParticle.px + fParticle.py * fParticle.py);

        switch (fDecayMode)
        {
            case ProtonType::Primordial :
                hPrimVec.at(fRadIter)->Fill(fPt,1./fPt);
                break;
            case ProtonType::LightNuclei :
                hDecayVec.at(fRadIter)->Fill(fPt,1./fPt);
                break;
            
            default:
                break;
        }
    }
    bar.mark_as_completed();
    // Show cursor
    indicators::show_console_cursor(true);

    TCanvas *canv = new TCanvas("canv","",1200,600);
    canv->Divide(2,1);

    for (size_t iter = 0; iter < fRlist.size(); ++iter)
    {
        canv->cd(1)->SetLogy();
        hPrimVec.at(iter)->Scale(1./(fEvents * 2*ROOT::Math::Pi() * fRapWidth * fPtWidth));
        hPrimVec.at(iter)->SetLineWidth(2);
        fLegend->AddEntry(hPrimVec.at(iter),Form("R #in (%d,%d)",std::get<0>(fRlist.at(iter)),std::get<1>(fRlist.at(iter))),"l");
        if (iter == 0)
        {
            hPrimVec.at(iter)->SetTitle("");
            hPrimVec.at(iter)->Draw("hist c plc");
        }
        else
            hPrimVec.at(iter)->Draw("hist c plc same");

        canv->cd(2)->SetLogy();
        hDecayVec.at(iter)->Scale(1./(fEvents * 2*ROOT::Math::Pi() * fRapWidth * fPtWidth));
        hDecayVec.at(iter)->SetLineWidth(2);
        if (iter == 0)
        {
            hDecayVec.at(iter)->SetTitle("");
            hDecayVec.at(iter)->Draw("hist c plc");
        }
        else
            hDecayVec.at(iter)->Draw("hist c plc same");
    }
    fLegend->Draw("same");
}