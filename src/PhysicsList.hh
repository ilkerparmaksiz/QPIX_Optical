//
// Created by ilker on 6/14/25.
//

#ifndef QPIXG4_PHYSICSLIST_HH
#define QPIXG4_PHYSICSLIST_HH

#include "G4VModularPhysicsList.hh"
#include "FTFP_BERT_HP.hh"
class PhysicsList : public FTFP_BERT_HP {
public:
    PhysicsList();
    virtual ~PhysicsList();
protected:
    virtual void ConstructProcess() override;

};


#endif //QPIXG4_PHYSICSLIST_HH
