//
//  EffectiveStressDilatation2D_Example.h
//  CSMP_API_library2014
//
//  Created by Stephan Matthai on 1/24/14.
//  Copyright (c) 2014 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_API_library2014__EffectiveStressDilatation2D_Example_h
#define CSMP_API_library2014__EffectiveStressDilatation2D_Example_h

#include <iostream>
#include "Example.h"

namespace csmp {

template<size_t> class Model;
template<size_t> class Region;
template<size_t,template<size_t> class> class PDE_Integrator;

/// gas flow into a horizontal well drilled into overpressured "shale" layer
class  EffectiveStressDilatation2D_Example : public Example {
  public:
    EffectiveStressDilatation2D_Example();
    virtual void Run();
    virtual void Specifications();
  
  private:
    void ComputeTransientFluidPressure( Model<2U>&, double64 time_increment );
    void ComputeTransientFluidPressure( Model<2U>&, const char* split_boundary_name, double64 time_increment );

    PDE_Integrator<2U,Region>*  fluid_pressure_;
    bool  verbose_;
};

} // csmp

#endif /* defined(CSMP_API_library2014__EffectiveStressDilatation2D_Example_h) */
