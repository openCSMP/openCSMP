#ifndef DES_TWOPHASEFLOW2D_EXAMPLE_H
#define DES_TWOPHASEFLOW2D_EXAMPLE_H

#include "Example.h"

#include "Model.h"
#include "TwoPhaseModel.h"
#include "FlowFunctions1.h"

namespace csmp {

class  DESTwoPhaseFlow2D_Example : public Example
{

public:
    virtual void Run();
    virtual void Specifications();

private:

    // a simple function that computes the total mobility from the chosen relative permeability model
    void computeTotalMobility( Model<2U>& mdl, FlowFunctions2<2U>& flowfunctions );

};

} // csmp

#endif // DES_TWOPHASEFLOW2D_EXAMPLE_H
