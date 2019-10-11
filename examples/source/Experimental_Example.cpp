//
//  Experimental_Example.cpp
//  CSMP_API_library2014
//
//  Created by Stephan Matthai on 1/27/14.
//  Copyright (c) 2014 Stephan Matthai. All rights reserved.
//

#include "CSMP_definitions.h"
#include "Experimental_Example.h"
#include "compareFloats.h"
#include "ANSYS_Model3D.h"
#include "ANSYS_Model2D.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_NT_op_N_dV.h"
#include "PDE_Integrator.h"
#include "PDE_IntegratorExperimental.h"
#include "Face.h"
#include "NumIntegral_NT_op_N_dS.h"
#ifdef CSMP_WITH_SAMG_SOLVER
#include "SAMG_Solver.h"
#include "SAMG_Settings.h"
#endif
#ifdef CSMP_WITH_MESCHACH
#include "Gauss_Solver.h"
#endif

#include <iostream>
#include "InputDataManager.h"
#include "VTK_Interface.h"
#include "VTU_Interface.h"

// include any header files that you need here...
#include "SteadyStateDiffusor.h"
#include "CSMP_highLevelUtilities.h"
#include "ComputationalSettings.h"
#include "IAPWS_H2OPropertiesVisitor.h"
#include "NumIntegral_NT_op_dNi_dV.h"
#include "PT_op.h"
#include "NumIntegral_BT_D_B_dV.h"
#include "NumIntegral_PT_op_dS.h"
#include "NumIntegral_PT_op_dV.h"
#include "NumIntegral_BT_D_op_dV.h"
#include "NumIntegral_BT_op_dV.h"
#include "StressesAndStrains.h"
#include "StressesAndStrains2.h"
#include "NodeCenteredFiniteVolumeTransport.h"
#include "SinglePhaseVelocityVisitor.h"
#include "ExtractTensorVariableComponent.h"

using namespace std;

namespace csmp {

const size_t DIM(3U);

void Experimental_Example::Specifications()
  {
     SetTitle( "Experimental_Example" );
     SetDifficulty( 1 );
     SetCategory( "Software Functionality" );
     AddAuthor( "You!" );
     AddDescription( "source in: Experimental_Example.cpp" );
     AddDescription( "Empty example for the user to experiment with" );
     AddRequirement( "none" );
     AddRequirement( "no predefined model or variables file" );
  }

// OTHER RELATIONSHIPS



/**
     Put CSMP code that you would like to test here and run it as part of the 
     example suite.
     
     Systems Modelling and Design model that calculates the effective stress in a dam, 
     using the gravitational loading and the plane stress assumption.
     
     input models:
     - slope_model1
*/
void Experimental_Example::Run()
{
  // 0. converting ANSYS input file set into a CSMP model and save it to binary file
  // -------------------------------------------------------------------------------
  const bool irregular_mesh(false);         /* true = non-box shaped model, false = box shaped model */
  const bool binary_file(true);             /* true = binary, false = ascii */
  const bool use_regions_file(true);        /* true = reduce regions according to regions file, false = does not redure regions */
  const bool create_boundaries(true);       /* true = creates boundaries around model, false = does not create boundaries */
  const bool create_splitboundaries(false); // FAIL /* true = creates splitboundaries around model, false = does not create splitboundaries */
  // Ansys model
  const string model_name("BoxHalfs2D");
  ANSYS_Model2D ansys_model("BoxHalfs2D", "BoxHalfs2D", "THMC_shear_zone-variables.txt",
                             irregular_mesh, binary_file, use_regions_file, create_boundaries, create_splitboundaries );
                             
  // testing region insertion here
  Region<2> region1_before(ansys_model.Region("MATRIX_RIGHT"));
//   ansys_model.InsertSplitBoundary( "MATRIX_RIGHT", "MATRIX_LEFT" );

  // saving model into CSMP native file format
  ansys_model.OutputToBinaryFile( model_name.c_str() );
   

  // 1. starting the simulation with the creation of a SplitBoundary
  // -------------------------------------------------------------------------------
  // read model model from file and get started with SplitBoundary code
  Model<2> model( model_name );
  model.RegionsOut();
  model.BoundariesOut();
  // checking the regions of the model (visualising their perimeter)
  VTU_Interface<2>  vtk_out( model );
  model.InputPropertyValue( "test variable", makeScalar(ANY,0.) );
  model.Region("MATRIX_LEFT").InputPropertyValue( "test variable", makeScalar(ANY,1.), PERIMETER );
  model.Region("MATRIX_RIGHT").InputPropertyValue( "test variable", makeScalar(ANY,1.), PERIMETER );
  vtk_out.OutputDataToVTU( "region-flag", "test variable", "MATRIX_LEFT", 1 );
  vtk_out.OutputDataToVTU( "region-flag", "test variable", "MATRIX_RIGHT", 2 );

  // create SplitBoundary between the model regions
  Region<2> region1_after(model.Region("MATRIX_RIGHT"));
  model.InsertSplitBoundary( "MATRIX_RIGHT", "MATRIX_LEFT" );
  // putting a lower dimensional region inside of all split boundaries
  set<string>  newly_created_regions;
  model.RegionsFromSplitBoundaries( newly_created_regions );
  assert( !newly_created_regions.empty() );
  model.RegionsOut();
  // getting a reference to the split boundary
  SplitBoundary<2>  interface( model.SplitBoundary("SPLITBOUNDARY_MATRIX_RIGHT_MATRIX_LEFT") );
  // getting a reference to the newly created region
  Region<2>  detached_surface( model.Region( (*newly_created_regions.begin()) ) );

  cout <<"\nExperimental_Example: That's it!\n";
  
} // end Run





/* TODO: get 3D approach to work in 2D, see below
 // create Boundary
 InsertBoundary( const char* region1, const char* region2,
 const bool remove_dim_minus1_region(true);
 model.CreateInternalBoundaryFrom( "STANDARD", remove_dim_minus1_region );
 set<string> strings_in_bundary_name({"STANDARD"});
 string boundary_name = model.FindBoundaryName( strings_in_bundary_name );
 model.BoundariesOut();
 // create SplitBoundary from boundary
 Boundary<2>& boundary_domain(model.Boundary(boundary_name));
 model.CreateSplitBoundaryFrom( boundary_domain );
 model.SplitBoundariesOut();

*/

} // csmp
