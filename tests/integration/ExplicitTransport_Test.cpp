//
//  ExplicitTransport_Test.cpp
//  CSMP_API_examples
//
//  Created by Stephan Matthai on 2/09/2015.
//  Copyright (c) 2015 Stephan Matthai. All rights reserved.
//

#include "Model.h"
#include "vsetMakers.h"
#include "ExplicitTransport_Test.h"
#include "VTK_Interface.h"
#include "Index.h"

// is being tested
#include "finiteVolumeAuxiliaryFunctions.h"

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
 
    2. Flux balance test
 
    3. TVD test
 
    4. Boundary interaction tests
*/
void ExplicitTransport_Test::run()
 {
    // BUILD VSET MAKER PATCH POLY-ELEMENT MODEL
    // =========================================
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

    VTK_Interface<test_dim> vtk_output;
    vtk_output.OutputDataToVTK( model, "test_output1", "permeability", 0 );
   
    // tests initialisation of stencils and property assigment
    model.InstantiateFiniteVolumes();
    Region<test_dim>  model_domain = model.Region("Model");
    // call of testee
    // here the facet fluxes as precomputed are used
    VectorVariable<test_dim> velo(ANY,ANY,ANY, 1., 1., 1. );
    velo.EuclideanNormalize(); // to 1.
    //cout <<"\nvelocity magnitude: "<< velo.Length() <<"\n";
    model.InputPropertyValue( "total velocity", velo );
    const bool initialize_flux( true );
    // uses total velocity
    initializeFiniteVolumeProperties( model, model_domain, initialize_flux );
    printRangeOfVariable( model, "total velocity" );
    printRangeOfVariable( model, "facet flux" );
    printRangeOfVariable( model, "facet normal" );
 
 
    // test 1: is the finite volume equal to the element volume
    // --------------------------------------------------------
    const csmp::Index  vol_key(model.Database().StorageKey("finite volume"));
    double64           total_volume(0.);
    for ( vector<Node<test_dim>*>::iterator nit=model_domain.NodesBegin(); nit!=model_domain.NodesEnd(); ++nit ) {
         total_volume += (*nit)->Read( vol_key );
      }
    cout <<"\nrun: model volume vs. finite volume integrated: "<< model_domain.Volume() <<" vs "<< total_volume << endl;
    // TODO: test fails for skewed hexahedra, fix:
    _equal( total_volume, model_domain.Volume(), numeric_limits<double64>::epsilon() );
 
 
    // test 2: is the flux conserved in the interior of the model (prescribed velocity case)
    // -------------------------------------------------------------------------------------
    const csmp::Index  velo_key(model.Database().StorageKey("velocity"));
    const csmp::Index  flux_key(model.Database().StorageKey("facet flux"));
    model_domain.UpdateMemberIndexes();
    vector<double64>  flux_balance( model_domain.Nodes(), 0. );
    size_t inside_node, outside_node;
    for ( vector<Element<test_dim>*>::iterator it=model_domain.ElementsBegin(); it!=model_domain.PerimeterElementsBegin(); ++it ) {
         (*it)->Read( velo_key, velo );
         for ( size_t i=0U; i<(*it)->Facets(); ++i ) {
              // finding the orientation of the facets
              const size_t facet_ip(0U);
              (*it)->FV()->FacetEdgeNodes( i, inside_node, outside_node );
              flux_balance[ (*it)->N(inside_node)->Idx() ]  += (*it)->Read( i, facet_ip, flux_key );
              flux_balance[ (*it)->N(outside_node)->Idx() ] -= (*it)->Read( i, facet_ip, flux_key );
           }
      }
    auto min_value = *min_element(flux_balance.begin(),flux_balance.end());
    auto max_value = *max_element(flux_balance.begin(),flux_balance.end());
    _equal( fabs(min_value), 0., numeric_limits<double64>::epsilon() );
    _equal( fabs(max_value), 0., numeric_limits<double64>::epsilon() );
    cout <<"\nrun: velocity vs. finite volume flux balances (min/max) for total velocity of 1.: "<< velo <<" vs "<< min_value <<"-"<< max_value << endl;


    // test 3: is the flux conserved along no-flow boundaries
    // ------------------------------------------------------

    // void test_Create_Prism_VSet(VSet<3U>& vset, bool bSkewed=false );
    

    // void test_Create_Prism_Hexa_VSet(VSet<3U>& vset, bool bSkewed=false );
 }


} // end csmp
