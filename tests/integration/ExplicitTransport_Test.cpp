//
//  ExplicitTransport_Test.cpp
//  CSMP_API_examples
//
//  Created by Stephan Matthai on 2/09/2015.
//  Copyright (c) 2015 Stephan Matthai. All rights reserved.
//

#include "Model.h"
#include "vset_makers.h"
#include "ExplicitTransport_Test.h"
#include "VTK_Interface.h"
#include "Index.h"

// is being tested
#include "finiteVolumeUniversalFunctions.h"

using namespace std;

namespace csmp {

ExplicitTransport_Test::ExplicitTransport_Test( const char* test_model, const char* test_variables )
 : test_model_(test_model),
   test_variable_file_(test_variables)
 {
 }

ExplicitTransport_Test::~ExplicitTransport_Test() { }

/**
    Tests performed:
    
    1. Verify that we have indeed the right porevolumes facet areas and facet normals
    
    use BoxHalves3D - tetra model with 2 regions for tetra testing
    For poly-elemement mesh testing use model b25
*/
void ExplicitTransport_Test::run()
 {
    const size_t    test_dim(3U);
    VSet<test_dim>  vset;
   
    // create test model
    const bool bSkewed(false);
    test_Create_Hexahedra_VSet( vset, bSkewed );
    //                     singlePhase_advection-variables.txt
    const bool isoparametric(true);
    Model<test_dim> model( vset, test_variable_file_.c_str(), isoparametric );
    printModelDimensions( model );

    // property assigment
    model.InputPropertyValue( "permeability", makeScalar(ANY,1.0e-12) );
    model.InputPropertyValue( "porosity", makeScalar(ANY,0.25) );
    VectorVariable<test_dim> velo(ANY,1.0e-5);
    model.InputPropertyValue( "velocity", velo );
   
VTK_Interface<test_dim> vtk_output;
vtk_output.OutputDataToVTK( model, "test_output1", "permeability", 0 );
   
    // tests initialisation of stencils and property assigment
    model.InstantiateFiniteVolumes();
    Region<test_dim>  gref_model    = model.Region("Model");
    //Region<test_dim>  gref_interior = model.Region("interior");
    // call of testee
    initializeFiniteVolumeProperties( model, gref_model );
   
    // test 1: is the volume correct
    csmp::Index  vol_key(model.Database().StorageKey("finite volume"));
    double64     total_volume(0.);
    for ( vector<Node<test_dim>*>::iterator nit=gref_model.NodesBegin(); nit!=gref_model.NodesEnd(); ++nit ) {
         total_volume += (*nit)->Read( vol_key );
      }

    cout <<"\nrun: model volume vs. finite volume integrated: "<< gref_model.Volume() <<" vs "<< total_volume << endl;
    _equal( total_volume, gref_model.Volume(), numeric_limits<double64>::epsilon() );
   
    // void test_Create_Prism_VSet(VSet<3U>& vset, bool bSkewed=false );
    

    // void test_Create_Prism_Hexa_VSet(VSet<3U>& vset, bool bSkewed=false );
 }


} // end csmp
