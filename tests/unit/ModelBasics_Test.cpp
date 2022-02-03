/*
 *  ModelBasics_Test.h
 *  csmp_core
 *
 *  Created by SKM 2/2/2022.
 *
 */

#include "ModelBasics_Test.h"

#include "vsetMakers.h"
#include "Model.h"
#include "Region.h"
#include "Boundary.h"
#include "VSet.h"
#include "VTK_Interface.h"

using namespace std;

namespace csmp {

ModelBasics_Test::ModelBasics_Test( bool verbose )
    : verbose_(verbose)
    {
    }
    

ModelBasics_Test::~ModelBasics_Test()
    {
    }
    

/**
    Reading of all basic variable types and placements.
*/
void ModelBasics_Test::run()
  {
     VSet<3U>      vset;
     ModelTopology topology;
     test_Create_FracBox( vset, topology );

     // creating a model
     // Model( ModelTopology&, VSet<dim>&, const char* var_file,
     //        bool create_boundary_objects = false, bool box_shaped = true );

     Model<3U>  model( topology, vset, "CSMP-1phase-variables.txt", true, false );
     model.Name("FracBox");
     
     // assigning some dummy values to verify functionality
     Region<3U>& model_domain = model.Region("Model");
     Region<3U>& matrix_domain = model.Region("MATRIX");
     model_domain.InputPropertyValue( "fluid pressure", makeScalar(ANY,1e5) );
     model_domain.InputPropertyValue( "permeability", makeScalar(ANY,1.0e-15) );
     matrix_domain.InputPropertyValue( "permeability", makeScalar(ANY,1.0e-12) );
     //model.Out(); // crashes when trying to print normals to FV Stencil for ISO_LIN_HEX
     model.RegionsOut(); // if boundaries are created, programme crashes here when trying to calculate volume of an element
     
     // saving model to disk
     model.OutputToBinaryFile( "ModelBasics_Test" );
     
     // bringing model back
     set<string>  subset_variables; // all variables
     Model<3U>    restored_model( string{"ModelBasics_Test"}, subset_variables );
     double       pmin, pmax;
     restored_model.MinMaxOf( "permeability", pmin, pmax );
     
     const double tolerance( 1.0e-17 );
     _equal( pmin, 1.0e-15, tolerance );
     _equal( pmax, 1.0e-12, tolerance );
     
     _test( std::numeric_limits<size_t>::max() == UINT_MAX );
     
     // some visual QC, using VTK
     restored_model.RegionsOut();
     restored_model.BoundariesOut();
     //restored_model.Out(); // fails when trying to print normals to FV Stencil for ISO_LIN_HEX
  
   } // end run

  
} // end csmp



