#ifndef TUTORIAL3_EXAMPLE_H
#define TUTORIAL3_EXAMPLE_H

#include "Example.h"

#include "Model.h"
#include "TwoPhaseModel.h"

namespace csmp {

/// incompressible two-phase flow of immiscible fluids
class  Tutorial3_Example : public Example
{

public:

    virtual void Run();
    virtual void Specifications();

private:

    /// function computes the total mobility lambda_t  from the relative permeability curves defined in the chosen relative permeability model
    void computeTotalMobility( Model<2U>& mdl, TwoPhaseModel<2U>& relperm );

};

} // csmp

#endif // TUTORIAL3_EXAMPLE_H
