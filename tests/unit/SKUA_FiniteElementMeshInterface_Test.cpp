//
//  SKUA_FiniteElementMeshInterface_Test.cpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 2/12/20.
//  Copyright © 2020 Stephan Matthai. All rights reserved.
//

#include "SKUA_FiniteElementMeshInterface_Test.h"
#include "SKUA_Model.h"
#include "VTK_Interface.h"

using namespace std;

namespace csmp {

void SKUA_FiniteElementMeshInterface_Test::run()
 {
    // 0. creating box-shaped test model with nodes on the side surfaces
    SKUA_Model  model( "V6_SOLID2",
                       "CSMP-variables.txt",
                       false,   ///< true = binary, false = ascii */
                       true,    ///< true = reduce regions according to regions file, false = does not redure regions
                       false,   ///< true = creates boundaries around model, false = does not create boundaries 
                       false,   ///< true = creates splitboundaries around model, false = does not create splitboundaries 
                       true );  ///< isoparametric is the only option if model does not only consist of simplex elements

    printModelDimensions( model, true );
    model.RegionsOut();

  //  InputDataManager<3U>  model_configuration;
  //  model_configuration.ConfigureFromFile( model, model_name.c_str(),
  //                                         false, true, true, true, false );

    // checking model boundary flagging
    boxFlagsToVariable( model, "node variable", "element variable" );
    
    VTK_Interface<3U>  vtk_output;
    vtk_output.OutputDataToVTK( model, "porosity", "porosity", 1 );
    vtk_output.OutputDataToVTK( model, "node-flag", "node variable", 1 );
    vtk_output.OutputDataToVTK( model, "element-flag", "element variable", 1 );
    
    // 1. basic testing
    /// consecutive numbering of nodes and elements
    /// mesh consistency (Jacobians, element volume range, non-manifold vertices, triangle boxes
    /// pfverts, disambiguated neighbors of lower-dimensional elements
    /// consistent facing directions of surface elements
    /// consecutive numbering of series of line elements
    /// material IDs
    /// test assignment of box boundary flags (sides, edges, corners)
    /// property assignment
    /// use of regions file (is this still needed?)
    // 2. advanced testing
    /// generation of boundaries  
    /// generation of split boundaries                      

 } // end run


} // end csmp
