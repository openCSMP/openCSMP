//
//  fault_permeability_modeling.h
//  CSMP_ReservoirSimulator
//
//  Created by Stephan Matthai on 12/19/12.
//  Copyright (c) 2012 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_FAULT_PERMEABILITY_MODELLING_H
#define CSMP_FAULT_PERMEABILITY_MODELLING_H

#include <iostream>

namespace csmp {

/**

@todo  identify the rocks that the fault juxtaposes against one-another:
       How?
        - have regions and boundary of the fault at the same time. One for property modeling the other for flow simulation
        - 1. use integer code to transfer region names to variables, storing mapping in a map
        - 2. find elements adjacent to fault surfaces via their nodes: put pointers into map
        - 3. extract these elements that share face with surface
        - 4. group fault elements on the basis which materials are juxtaposed against one-another along them: granite-granite,  granite-sediments etc.
        - 5. form corresponding boundaries
        - 6. make permeability model dependent on which rocks are juxtaposed against eachother along fault

@todo  !!! deal with the case where fault planes have holes because the surfaces are shared with other faults
          (consider only the longest closed loop of tipline segments, drop Dirichlet conditions from interior loops)
          
@todo  ! move stress tensor inside the bounding box of model so that one can see it

@todo  !! SKUA_Interface   for the output of properties to pointcloud

@todo  !!! vertical relative permeability modeling for fractured porous media (make experiments with Shaho's simulator)

*/
void faultPermeabilityModelling( const char* model_name );

}

#endif /* CSMP_FAULT_PERMEABILITY_MODELLING_H */
