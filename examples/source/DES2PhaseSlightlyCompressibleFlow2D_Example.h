#ifndef DES_2PHASE_SLIGHTLY_COMPRESSBILE_FLOW_2D_EXAMPLE_H
#define DES_2PHASE_SLIGHTLY_COMPRESSBILE_FLOW_2D_EXAMPLE_H

#include "Example.h"

#include "Model.h"
#include "TwoPhaseModel.h"
#include "FlowFunctionsModule.h"

namespace csmp {

class  DES2PhaseSlightlyCompressibleFlow2D_Example : public Example
{

public:
    virtual void Run();
    virtual void Specifications();

private:

    // a simple function that computes the total mobility from the chosen relative permeability model
    void computeTotalMobility( Model<2U>& mdl, FlowFunctionsModule1<2U>& flowfunctions );

};

} // csmp

#endif // DES_2PHASE_SLIGHTLY_COMPRESSBILE_FLOW_2D_EXAMPLE_H
