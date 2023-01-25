#include "CanonicalSuppressionPlugin.h"
#include "Configurator.h"
#include "Parser.h"
#include "TVector3.h"
#include "TMath.h"
#include <limits>
#include <iostream>
#include "THGlobal.h"

using namespace std;
extern TString sModelINI;

CanonicalSuppressionPlugin::CanonicalSuppressionPlugin()
{
    ReadParameters();
}

CanonicalSuppressionPlugin::~CanonicalSuppressionPlugin() {}

void CanonicalSuppressionPlugin::Apply(Event *tEvent)
{
    if (m_bSkip) return;
    Float_t Rc = m_Rc;
    std::list<Particle> *tParticles = tEvent->GetParticleList();
    Int_t N = tParticles->size();
    int i = 0;
    for (const Particle &partI : *tParticles)
    {
        int ns = partI.GetParticleType()->GetNumberS();
        int j = 0;
        for (int ks = 0; ks < ns; ++ks)
        {
            for (const Particle &partJ : *tParticles)
            {
                int nas = partJ.GetParticleType()->GetNumberAS();
                for (int kas = 0; kas < nas; ++kas)
                {
                    Float_t rel_beta_x = partI.px / partI.e - partJ.px / partJ.e;
                    Float_t rel_beta_y = partI.py / partI.e - partJ.py / partJ.e;
                    Float_t rel_beta_z = partI.pz / partI.e - partJ.pz / partJ.e;
                    Float_t rel_beta_2 = TMath::Power(rel_beta_x, 2) + TMath::Power(rel_beta_y, 2) + TMath::Power(rel_beta_z, 2);
                    Float_t rel_gamma = 1 / TMath::Sqrt(1 - rel_beta_2);

                    Float_t dist = ( kHbarC / rel_gamma) * TMath::Sqrt(
                                                TMath::Power(partI.x - partJ.x, 2) +
                                                TMath::Power(partI.y - partJ.y, 2) +
                                                TMath::Power(partI.z - partJ.z, 2));
                    std::cout << "dist = " << dist << std::endl;
                    if (dist < Rc)
                    {
                        grid.at(i, j) = dist;
                    }
                    if (i == 0)
                    { // correct?
                        Int_t eidJ = partJ.eid;
                        if (eidToIs.count(eidJ) == 0)
                        {
                            eidToIs[eidJ] = std::vector<UInt_t>();
                        }
                        eidToJs[eidJ].push_back(j);
                    }
                    ++j;
                }
            }
            Int_t eidI = partI.eid;
            if (eidToIs.count(eidI) == 0)
            {
                eidToIs[eidI] = std::vector<UInt_t>();
            }
            eidToIs[eidI].push_back(i);
            ++i;
        }
    }

    while (!grid.canReshuffleToDiagonal())
    {
        auto rowIs = grid.rowIs();
        for (UInt_t rowI : rowIs)
        {
            std::map<UInt_t, Float_t> row = grid.row(rowI);
            if (row.size() > 1)
            {
                grid.clearRowExceptOne(rowI);
            }
        }

        auto colJs = grid.colJs();
        for (UInt_t colJ : colJs)
        {
            std::map<UInt_t, Float_t> col = grid.col(colJ);
            if (col.size() > 1)
            {
                grid.clearColExceptOne(colJ);
            }
        }

        for (auto &aPart : *tParticles) {
            if (brokenMultistrange(aPart)) {
                Int_t aEid = aPart.eid;
                for (UInt_t i : eidToIs[aEid]) {
                    grid.clearRow(i);
                }
                for (UInt_t j : eidToJs[aEid]) {
                    grid.clearCol(j);
                }
            }
        }
    }
    tParticles->remove_if([this](Particle &part) -> bool
                          {
                              Int_t eid = part.eid;
                              Int_t ns = part.GetParticleType()->GetNumberS();
                              Int_t nas = part.GetParticleType()->GetNumberAS();
                              Int_t netS = ns - nas;
                              return netS != 0 && !(
                                                      std::all_of(eidToIs[eid].begin(), eidToIs[eid].end(), [this](UInt_t i)
                                                                  { return this->grid.hasI(i); }) &&
                                                      std::all_of(eidToJs[eid].begin(), eidToJs[eid].end(), [this](UInt_t j)
                                                                  { return this->grid.hasJ(j); }));
                          });
}
/*
// Maybe still will be used if the lambda-in-lambda approach will not work, but then eidToIs and eidToJs
// still have to be added to the predicate class:
bool ParticleBreaksConservationLaw::operator()(Particle &part) {
    Int_t eid = part.eid;
    return !(
        std::all_of(eidToIs[eid].begin(), eidToIs[eid].end(), [this](UInt_t i) { return this->grid.hasI(i); }) &&
        std::all_of(eidToJs[eid].begin(), eidToJs[eid].end(), [this](UInt_t j) { return this->grid.hasJ(j); })
    );
}
*/

void CanonicalSuppressionPlugin::ReadParameters()
{
    Configurator *tModelParam;
    Parser *tParser;

    tModelParam = new Configurator;
    tParser = new Parser(sModelINI.Data());
    tParser->ReadINI(tModelParam);
    delete tParser;

    try
    {
        m_Rc = tModelParam->GetParameter("Rc").Atof();
    }
    catch (TString &str)
    {
        cout << "Parameter " << str.Data() << " is not known" << endl;
        cout << "CanonicalSuppression Plugin will be skipped" << endl;
        m_bSkip = true;
    }
    delete tModelParam;
}

const Utils::Option<Particle> &&CanonicalSuppressionPlugin::findParticleByEid(Event *tEvent, Int_t tEid)
{
    list<Particle> *aParticles = tEvent->GetParticleList();
    auto result = std::find_if(aParticles->begin(), aParticles->end(), [tEid](Particle particle) -> bool
                               { return particle.eid == tEid; });
    if (result == aParticles->end())
    {
        return std::move(Utils::None<Particle>());
    }
    else
    {
        return std::move(Utils::Some(*result));
    }
}

bool CanonicalSuppressionPlugin::brokenMultistrange(Event *tEvent, Int_t tEid)
{
    const Utils::Option<Particle> &&aMaybeParticle = findParticleByEid(tEvent, tEid);
    if (aMaybeParticle.has_value())
    {
        Particle aPart = (dynamic_cast<const Utils::Some<Particle> &&>(aMaybeParticle)).val;
        return brokenMultistrange(aPart);
    }
    else
    {
        return false;
    }
}

bool CanonicalSuppressionPlugin::brokenMultistrange(Particle &tPart)
{
    Int_t ns = tPart.GetParticleType()->GetNumberS();
    Int_t nas = tPart.GetParticleType()->GetNumberAS();
    Int_t aEid = tPart.eid;
    return !(ns == eidToIs[aEid].size() && nas == eidToJs[aEid].size());
}