/*
 *  ModelBasics_Test.h
 *  csmp_core
 *
 *  Created by SKM 2/2/2022.
 *
 */

#include <type_traits>
#include "ModelBasics_Test.h"

#include "vsetMakers.h"
#include "Model.h"
#include "Region.h"
#include "Boundary.h"
#include "PropertyConstraints.h"
#include "VSet.h"
#include "VTK_Interface.h"
#include "compareFloats.h"

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
    bool create_boundaries_from_faces{false};
    _test( TestWriteModelToDiskAndReadBack( create_boundaries_from_faces ) );
    create_boundaries_from_faces = true;
     _test( TestWriteModelToDiskAndReadBack( create_boundaries_from_faces ) );

    // testing the repair of non-unique regions built using property constraints
    TestRebuiltRegionsFromPropertyConstraints();
    
    cout <<"\nModelBasics_Test::run: raw size of local variable storage of Model: ";
    cout << sizeof(LocalVariableStorage<3,Model>) << endl;
    // around 100-bytes for empty storage!
 
   } // end run




bool ModelBasics_Test::TestWriteModelToDiskAndReadBack( bool create_boundaries_from_faces )
 {
     VSet<3U>      vset;
     ModelTopology topology;
     test_Create_FracBox( vset, topology );

     // creating model with boundaries
     Model<3U>  model( topology, vset, "CSMP-1phase-variables.txt", create_boundaries_from_faces, false );
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
     size_t default_uint = std::numeric_limits<size_t>::max();
     _test( default_uint != UINT_MAX ); // should be false because UINT_MAX is not for size_t
     _test( hasDefaultValueForUnassignedInteger( default_uint ) );
     
     // some visual QC, using VTK
     restored_model.RegionsOut();
     restored_model.BoundariesOut();
     //restored_model.Out(); // fails when trying to print normals to FV Stencil for ISO_LIN_HEX
     
     return true;
 }
 
 
 
 
bool ModelBasics_Test::TestRebuiltRegionsFromPropertyConstraints()
 {
     const bool    using_isoparametric_elements{true};
     ModelTopology topology( "FracBox", using_isoparametric_elements );
     VSet<3U>      vset;
     test_Create_FracBox( vset, topology );

     const bool create_boundaries_from_faces{true}, box_shaped{false};
     Model<3U>  model( topology, vset, "CSMP-1phase-variables.txt", create_boundaries_from_faces, box_shaped );
     
     // assigning some dummy values to verify functionality
     Region<3U>& model_domain = model.Region("Model");
     Region<3U>& matrix_domain = model.Region("MATRIX");
     const csmp::Index phi_key = model.Database().StorageKey("porosity");
     model_domain.InputPropertyValue( "permeability", makeScalar(ANY,1.0e-12) );
     model_domain.InputPropertyValue( "porosity", makeScalar(ANY,1.0) ); // fracture
     matrix_domain.InputPropertyValue( "porosity", makeScalar(ANY,0.4) ); // matrix
     // elevating the porosity of some extra elements
     for ( size_t eidx{50}; eidx<150; ++eidx )
       matrix_domain.E(eidx)->Store( phi_key, makeScalar(ANY,0.8) );
       
     // building new region from PropertyConstraints
     PropertyConstraints porosity_constraints( "porosity", 0.75, 0.85 );
     porosity_constraints.AddConstraint( "permeability", 1.0e-12, 1.0e-12 );
     
     model.FormRegionFrom( "medium porosity", porosity_constraints );
     const Region<3U>& medium_porosity_domain = model.Region("medium porosity");
     _test( medium_porosity_domain.Elements() == 100 );

     return true;
 }


  
} // end csmp



