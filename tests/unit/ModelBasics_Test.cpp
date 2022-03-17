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
#include "PropertyConstraints.h"
#include "VSet.h"
#include "VTK_Interface.h"
#include "ANSYS_Model3D.h"
#include "MeshManagementUtilities.h"
#include "compareFloats.h"
#include "ModelComparator.h"

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
     TestModelConstructionFromVSet();
  
     bool create_boundaries_from_faces{false};
     _test( TestWriteModelToDiskAndReadBack( create_boundaries_from_faces ) );
     create_boundaries_from_faces = true;
     _test( TestWriteModelToDiskAndReadBack( create_boundaries_from_faces ) );
 
   } // end run



bool ModelBasics_Test::TestModelConstructionFromVSet()
 {
    const  bool bSkewed{false};
    VSet<3U>    vset, vset1;
    test_Create_Pyramid_Hexa_VSet( vset, bSkewed );
    
    // does the VSet write/reads correctly?
    double time0{3600.123}, time1;
    vset.OutputTo( "ModelBasics_Test", time0 );
    vset1.InputFrom( "ModelBasics_Test", time1 );
    // testing
    _test( vset1 == vset );
    _test( approximatelyEqual(time0,time1) );
    
    // test: basic constructor
    Model<3U>   model( vset, "CSMP-1phase-variables.txt" );
    // save to native binary
    model.OutputToBinaryFile( "ModelBasics_Test" );
    // bring back from binary
    set<string>  subset_variables; // all variables
    Model<3U>    restored_model( string{"ModelBasics_Test"}, subset_variables );
     
    // compare = test
    return true;
    
 } // end TestModelConstructionFromVSet




bool ModelBasics_Test::TestWriteModelToDiskAndReadBack( bool create_boundaries_from_faces )
 {
     VSet<3U>      vset;
     ModelTopology topology;
     test_Create_FracBox( vset, topology );

     // creating model with boundaries
     const bool only_with_regions(true);
     Model<3U>  model( topology, vset, "CSMP-1phase-variables.txt", only_with_regions );
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
     
     _test( approximatelyEqual( pmin, 1.0e-15 ) );
     _test( approximatelyEqual( pmax, 1.0e-12 ) );
     
     // testing fundamental assumption made working with default initialisations of 'size_t'
     size_t default_uint = std::numeric_limits<uint32_t>::max();
     _test( default_uint != UINT_MAX ); // should be false because UINT_MAX is not for size_t
     _test( hasDefaultValueForUnassignedInteger( default_uint ) );
     
     // some visual QC, using VTK
     restored_model.RegionsOut();
     restored_model.BoundariesOut();
     //restored_model.Out(); // fails when trying to print normals to FV Stencil for ISO_LIN_HEX
     
     return true;
 }
 





  
} // end csmp



