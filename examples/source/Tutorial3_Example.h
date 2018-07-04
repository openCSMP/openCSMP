#ifndef TUTORIAL3_EXAMPLE_H
#define TUTORIAL3_EXAMPLE_H

#include "Example.h"

#include "Model.h"
#include "TwoPhaseModel.h"

namespace csmp {

class  Tutorial3_Example : public Example
{

public:

    virtual void Run();
    virtual void Specifications();

private:

    // a simple function that computes the total mobility from the chosen relative permeability model
    void computeTotalMobility( Model<2U>& mdl, TwoPhaseModel<2U>& relperm );

};

} // csmp

#endif // TUTORIAL3_EXAMPLE_H
