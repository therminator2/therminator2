#ifndef UNIGEN_EVENT_SAVER
#define UNIGEN_EVENT_SAVER

#include "utils/Option.h"
#include "AbstractEventSaver.h"
#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <memory>

#include "UEvent.h"

class UnigenEventSaver : public AbstractEventSaver {
    public:
        UnigenEventSaver();
        virtual void Save(Event *, Model *, int);
    private:
        struct DecayNode {
           // DecayNode() : particle(nullptr) {}
            DecayNode(const Particle *p) : particle(p) {}
            DecayNode(const Particle *p, std::vector<std::shared_ptr<DecayNode> > c) : particle(p), children(c) {}
            const Particle *particle;
            const std::vector<std::shared_ptr<DecayNode> > children;

            const Utils::Option<DecayNode> *FindByEid(const Int_t &eid) const {
                if (particle->eid == eid) {
                    return new Utils::Some<DecayNode>(*this);
                }
                std::vector<std::shared_ptr<DecayNode> >::const_iterator it = children.begin();
                for ( ; it != children.end(); ++it) {
                    const Utils::Option<DecayNode> *descendant = (*it)->FindByEid(eid);
                        if (*descendant != Utils::None<DecayNode>()) {
                            return descendant;
                        }

                }
                return &Utils::None<DecayNode>();
            }
            bool operator==(const DecayNode &other) const {
                return particle->eid == other.particle->eid;
            }
            bool operator!=(const DecayNode &other) const {
                return particle->eid != other.particle->eid;
            }
            void Print(unsigned nesting_level) const {
                if (nesting_level > 0) {
                    for (unsigned i = 1; i < nesting_level; ++i) {
                        std::cout << "---";
                    }
                    std::cout << " " << particle->eid << " " << particle->pid << " " << particle->fathereid << " " << particle->fatherpid << std::endl;
                }
                for (unsigned i = 0; i < children.size(); ++i) {
                    children[i]->Print(nesting_level + 1);
                }
            }
            std::shared_ptr<DecayNode> AddParticle(Particle *p) const {
                std::vector<std::shared_ptr<DecayNode> > new_children;
                std::vector<std::shared_ptr<DecayNode> >::const_iterator it = children.begin();
                bool added = false;
                for ( ; it != children.end(); ++it) {
                    if ((*it)->FindByEid(p->fathereid) == &Utils::None<DecayNode>()) {
                        if ((*it)->particle->fathereid == p->eid) {
                            std::vector<std::shared_ptr<DecayNode> > sub_children;
                            sub_children.push_back(*it);
                            new_children.push_back(std::shared_ptr<DecayNode>(new DecayNode(p, sub_children)));
                            added = true;

                        } else {
                            new_children.push_back(*it);
                        }
                    } else {
                        new_children.push_back((*it)->AddParticle(p));
                        added = true;
                    }
                }
                if (!added) {
                    new_children.push_back(std::shared_ptr<DecayNode>(new DecayNode(p)));
                }
                return std::shared_ptr<DecayNode>(new DecayNode(particle, new_children));
            }
        };

        void AddDecayToEvent(UEvent *ev, std::shared_ptr<DecayNode> node, int parent_decay);
        TFile*	        mFile;
        TTree*	        mParameterTree;
        TTree*      	mEventTree;
        Int_t           mDecayIndex;
};

#endif
