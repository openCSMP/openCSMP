#ifndef DES_TWOPHASEFLOW3D_EXAMPLE_H
#define DES_TWOPHASEFLOW3D_EXAMPLE_H

#include "Example.h"

#include "Model.h"
#include "TwoPhaseModel.h"
#include "FlowFunctionsModule.h"

namespace csmp {

class  DESTwoPhaseFlow3D_Example : public Example
{

public:
    virtual void Run();
    virtual void Specifications();

private:

    // a simple function that computes the total mobility from the chosen relative permeability model
    void computeTotalMobility( Model<3U>& mdl, FlowFunctionsModule1<3U>& flowfunctions );

};

} // csmp

#endif // DES_TWOPHASEFLOW3D_EXAMPLE_H
