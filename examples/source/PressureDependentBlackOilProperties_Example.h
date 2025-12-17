//
//  PressureDependentBlackOilProperties_Example.h
//  shaho_sim
//
//  Created by Lukas Mosser on 7/7/14.
//  Copyright (c) 2014 LukasMosser. All rights reserved.
//

#ifndef CSMP_PRESSURE_DEPENDENT_BLACKOIL_PROPERTIES_EXAMPLE_H
#define CSMP_PRESSURE_DEPENDENT_BLACKOIL_PROPERTIES_EXAMPLE_H

#include "Example.h"

namespace csmp {
  
/** @brief PressureDependentBlackOilProperties_Example
 
  File set provided by Lukas Mosser to
  get a black-oil fluid property calculation done.
  
- fitting will have been done for a single temperature only
  (a compositional model is needed to deal with non-isothermal conditions)

- Correlations are contained in:
   *.ptb   file which is the output from a PVT analysis of the black-oil samples done
         with the software PetroleumExperts.  It contains the correlation
         coefficients that will be used by the 
         
        'FluidPropertyLookupTable' 
        
        to compute the density, viscosity and composition of the oil phase
        as well as the gas phase, including all the necessary solubility 
        information.
 */
class PressureDependentBlackOilProperties_Example : public Example
  {
  public:
    virtual void Run();
    virtual void Specifications();
  };
  
} // csmp


#endif // CSMP_PRESSURE_DEPENDENT_BLACKOIL_PROPERTIES_EXAMPLE_H

