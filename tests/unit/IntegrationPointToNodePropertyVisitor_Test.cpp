#include "IntegrationPointToNodePropertyVisitor_Test.h"
#include "ANSYS_Model3D.h"
#include "Model.h"
#include "Region.h"
#include "PropertyHandle.h"
#include "IntegrationPointToNodePropertyVisitor.h"
#include "compareFloats.h"
#include "geometricCalculations.h"

using namespace std;

namespace csmp {

constexpr uint32_t DIM{3U};

void IntegrationPointToNodePropertyVisitor_Test::run()
{
  const bool restrict_regions_with_regions_file{true}; // min variable set, the ones needed will be created
  ANSYS_Model3D  model( "prism_test", "prism_test", "CSMP-variables.txt", restrict_regions_with_regions_file );

  // 1. Test FE interpolation from nodes to FE integration points
  // ------------------------------------------------------------
  Test_NodeToIntegrationPointInterpolation( model );
  
  // 2. Test integration point to node interpolation
  // -----------------------------------------------
  Test_IntegrationPointToNodePropertyVisitor( model );
  
  // 3. Testing the interpolation of a field to the element barycentres
  // ------------------------------------------------------------------
  // (this is important because single-value interpolation may seem fine even when some interpolation function derivatives are wrong)
  Test_CoordinateInterpolationToBaryCenter( model );

} // run



/**
    Test interpolation from integration points to nodes
*/
void IntegrationPointToNodePropertyVisitor_Test::Test_NodeToIntegrationPointInterpolation( Model<3U>& model )
 {
    // Assuring nodes numbered from 0..n-1
    Region<DIM>& model_domain = model.Region( "Model" );
    model_domain.RenumberNodes();

    // Creating test properties, initializing values
    PropertyHandle<DIM> cpointProperty( model, "cpoint property", SCALAR, ELEMENT_INTEGRATION_POINT );
    PropertyHandle<DIM> nodeProperty( model, "node property", SCALAR, NODE );
//    cpointProperty = 0.;
    nodeProperty = 1.;

    // interpolating the node property values to the element integration points
    // ------------------------------------------------------------------------
    const csmp::Index nd_key = model.Database().StorageKey("node property");
    const csmp::Index ip_key = model.Database().StorageKey("cpoint property");
    ScalarVariable    sc;
    // interpolation
    for ( auto& it : model_domain.CellVector() )
      {
         for ( uint32_t i{0U}; i<it->IntegrationPoints(); ++i ) {
              it->PropertyValueAtIntegrationPoint<ScalarVariable>( nd_key, i, sc );
              it->Store( i, ip_key, sc );
           }
      }
    // testing
    for ( auto& it : model_domain.CellVector() )
      for ( uint32_t i{0U}; i<it->IntegrationPoints(); ++i ) {
           it->Read( i, ip_key, sc );
          _test( approximatelyEqual( 1., sc() ) );
        }

    // interpolating the node property values to the element barycentre
    // ----------------------------------------------------------------
    const csmp::Index el_key = model.Database().StorageKey("element variable"); // already exists
    // interpolation
    for ( auto& it : model_domain.CellVector() )
      {
         it->PropertyValueAtBaryCenter<ScalarVariable>( nd_key, sc );
         it->Store( el_key, sc );
      }
    // testing
    for ( auto& it : model_domain.CellVector() )
      _test( approximatelyEqual( 1., it->Read(el_key) ) );


 } // end Test_NodeToIntegrationPointInterpolation






/**
    Test interpolation from integration points to nodes
*/
void IntegrationPointToNodePropertyVisitor_Test::Test_IntegrationPointToNodePropertyVisitor( Model<3U>& model )
 {
    // Assuring nodes numbered from 0..n-1
    Region<DIM>& model_domain = model.Region( "Model" );
    model_domain.RenumberNodes();

    // Creating test properties, initializing values
    PropertyHandle<DIM> cpointProperty( model, "cpoint property", SCALAR, ELEMENT_INTEGRATION_POINT );
    PropertyHandle<DIM> nodeProperty( model, "node property", SCALAR, NODE );
    cpointProperty = 1.;
//    nodeProperty   = 0.;

    // Applying visitor to be tested
    IntegrationPointToNodePropertyVisitor<ScalarVariable,3U> test( model.Database(), "cpoint property",
                                                                   "node property", model.Region("Model").Nodes());
    test.ApplyWeightingToExtrapolatedValues();
    model_domain.ExtrapolateIntegrationPointToNodeProperty("cpoint property","node property");
    model_domain.Accept(test);

    const csmp::Index nkey = model.Database().StorageKey("node property");
    for (auto nit=model_domain.NodesBegin(); nit!=model_domain.NodesEnd(); nit++ )
      {
        //std::cout<<"Value: "<<(*it)->Read(model_.Database().StorageKey("node property"))<<std::endl;
        double nodevalue =(*nit)->Read( nkey );
        _test( approximatelyEqual( nodevalue, 1. ) );
      }

 } // end Test_IntegrationPointToNodePropertyVisitor



/**
      The node coordinate values are interpolated to the element barycentres and then compared with the x,y,z locations of the barycentre points.
*/
void IntegrationPointToNodePropertyVisitor_Test::Test_CoordinateInterpolationToBaryCenter( Model<3U>& model )
 {
    // Assuring nodes numbered from 0..n-1
    Region<DIM>& model_domain = model.Region( "Model" );
    model_domain.RenumberNodes();

    // interpolating the node property values to the element integration points
    // ------------------------------------------------------------------------
    const csmp::Index nd_key = model.Database().StorageKey("nodal vector");
    const csmp::Index el_key = model.Database().StorageKey("element vector"); // already exists
    VectorVariable<DIM> vc;
    // assigning the node coordinates to the 'nodal variable'
    for ( auto& nit : model_domain.NodeVector() ) {
         Point<DIM> nd_coord = nit->Coordinate();
         for ( uint32_t i{0u}; i<DIM; ++i ) vc(i) = nd_coord[i];
         nit->Store( nd_key, vc );
      }
    // interpolation to barycentre
    for ( auto& it : model_domain.CellVector() ) {
         it->PropertyValueAtBaryCenter( nd_key, vc );
         it->Store( el_key, vc );
      }
    // testing
    Point<DIM> xyz_min, xyz_max;
    model.MinMaxCoordinates(xyz_min, xyz_max );
    const double epsilon{ distance( xyz_min, xyz_max ) * numeric_limits<double>::epsilon() };
    
    for ( auto& it : model_domain.CellVector() ) {
          Point<DIM> barycentre = it->BaryCenter();
          it->Read( el_key, vc );
          bool coords_approx_equal{ true };
          // pyramid is known not to be that accurate
          if ( !isPyramid(it->FE_Type() ) ) {
              for ( uint32_t i{0u}; i<DIM; ++i )
                coords_approx_equal = approximatelyEqual( barycentre[i], vc[i], epsilon );
            }
          else {
              // restricts the checks to pyramids with a permissible skewness
              if ( isValidElement( it ) );
                for ( uint32_t i{0u}; i<DIM; ++i )
                  // coords_approx_equal = approximatelyEqual( barycentre[i], vc[i], epsilon * 1.0e13 ); // passes always
                  coords_approx_equal = approximatelyEqual( barycentre[i], vc[i], epsilon * 5.0e12 );
            }
          _test( coords_approx_equal );
          if ( !coords_approx_equal ) {
              cerr <<"\n"<< barycentre <<" vs "<< vc[0] <<" "<< vc[1] <<" "<< vc[2];
              cerr <<" "<< parseAbbreviated_FE_Type(it->FE_Type()) <<": "<< it->Idx() << endl;
           }
        }

 } // end Test_CoordinateInterpolationToBaryCenter



} // end csmp
