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

ModelBasics_Test::ModelBasics_Test()
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
     
     // builds model from regions, converting lower-dimensional elements into faces
     _test( TestWriteModelToDiskAndReadBack() );
     
     // 2D model with split boundaries
     _test( TestWriteModelToDiskAndReadBackWithInterfaces() );
 
   } // end run



bool ModelBasics_Test::TestModelConstructionFromVSet()
 {
    const  bool bSkewed{false};
    VSet<3U>    vset, vset1;
    create_Pyramid_Hexa_VSet( vset, bSkewed );
    
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




bool ModelBasics_Test::TestWriteModelToDiskAndReadBack()
 {
     VSet<3U>      vset;
     ModelTopology topology;
     create_FracBox( topology, vset );

     // creating model
     const bool do_not_use_regions_file{true};
     Model<3U>  model( topology, vset, "CSMP-1phase-variables.txt", do_not_use_regions_file );
     model.Name("FracBox");
     
     // assigning some dummy values to verify functionality
     Region<3U>& model_domain = model.Region("Model");
     Region<3U>& matrix_domain = model.Region("MATRIX");
     model_domain.InputPropertyValue( "fluid pressure", makeScalar(ANY,1e5) );
     model_domain.InputPropertyValue( "permeability", makeScalar(ANY,1.0e-15) );
     matrix_domain.InputPropertyValue( "permeability", makeScalar(ANY,1.0e-12) );

     if ( verbose_ ) model.Out(); // crashes when trying to print normals to FV Stencil for ISO_LIN_HEX
       
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
     uint32_t default_uint = std::numeric_limits<uint32_t>::max();
     _test( std::numeric_limits<size_t>::max() != UINT_MAX ); // false because UINT_MAX is not for size_t
     _test( isUninitialisedInteger( default_uint ) );
     _test( default_uint == UINT_MAX );
     
     // some visual QC, using VTK
     if ( verbose_ ) restored_model.Out();
     
     return true;
 }




// 2D model SPLIT_22 with crossing split boundaries
bool ModelBasics_Test::TestWriteModelToDiskAndReadBackWithInterfaces()
 {
     VSet<2U> vset;
     ModelTopology topology = create_BoundarySplitBoundaryPatch( vset );

     // creating model
     Model<2U>  model( topology, vset, "CSMP-1phase-variables.txt", false );
     model.Name("FracBox");
     
     // assigning some dummy values to verify functionality
     Region<2U>& model_domain = model.Region("Model");
     Region<2U>& matrix_domain = model.Region("lower");
     model_domain.InputPropertyValue( "fluid pressure", makeScalar(ANY,1e5) );
     model_domain.InputPropertyValue( "permeability", makeScalar(ANY,1.0e-15) );
     matrix_domain.InputPropertyValue( "permeability", makeScalar(ANY,1.0e-12) );

     if ( verbose_ ) model.Out(); // crashes when trying to print normals to FV Stencil for ISO_LIN_HEX
       
     // saving model to disk
     model.OutputToBinaryFile( "ModelBasics_Test" );
     
     // bringing model back
     set<string>  subset_variables; // all variables
     Model<2U>    restored_model( string{"ModelBasics_Test"}, subset_variables );
     double       pmin, pmax;
     restored_model.MinMaxOf( "permeability", pmin, pmax );
     
     _test( approximatelyEqual( pmin, 1.0e-15 ) );
     _test( approximatelyEqual( pmax, 1.0e-12 ) );
     
     // testing fundamental assumption made working with default initialisations of 'size_t'
     uint32_t default_uint = std::numeric_limits<uint32_t>::max();
     _test( std::numeric_limits<size_t>::max() != UINT_MAX ); // false because UINT_MAX is not for size_t
     _test( isUninitialisedInteger( default_uint ) );
     _test( default_uint == UINT_MAX );
     
     // some visual QC, using VTK
     if ( verbose_ ) restored_model.Out();
     
     return true;
     
 } // TestWriteModelToDiskAndReadBackWithInterfaces


  
} // end csmp



