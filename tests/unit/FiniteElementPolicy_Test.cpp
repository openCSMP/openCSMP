#include "FiniteElementPolicy_Test.h"
#include "vsetMakers.h"
#include "Model.h"
#include "Element.h"
#include "compareFloats.h"
#include "Exception.h"

using namespace std;

namespace csmp{

FiniteElementPolicy_Test::FiniteElementPolicy_Test()
{
  // debugging
  assert( model_ptr_ != NULL );
}

FiniteElementPolicy_Test::~FiniteElementPolicy_Test()
  {
     delete model_ptr_;
  }


void FiniteElementPolicy_Test::InitialiseModel()
 {
    VSet<3U> vset;
    test_Create_Prism_Hexa_VSet( vset, false );
    model_ptr_ = new Model<3U>( vset, "Variables_Test.txt" );
    e_ptr_     = model_ptr_->Region("Model").E(0);
 }
    


/// using 'test_Create_Prism_Hexa_VSet' dataset, tests fem operations against provided data in text file
void FiniteElementPolicy_Test::run()
{
   // getting an idea about element that is tested
   e_ptr_->Out();
   
  _test( Test_RstToXYZ() );

    /// integrates a scalar property over the element and returns this value
  _test( Test_PropertyIntegral() );

    /// interpolates node property values to the point of interest (in global coordinates); @attention costly for isoparametric elements
  _test( Test_PropertyValueAt() );
  
    /// returns the value of any property interpolated to the element's center of gravity
  _test( Test_PropertyValueAtBaryCenter() );
  
    /// returns the value of any property interpolated to the integration point of interest
  _test( Test_PropertyValueAtIntegrationPoint() );

    /// returns property values at the integration points
  _test( Test_IntegrationPointPropertyVector() );

    /// uses linear extrapolation of property values from the integration points to the nodes
  _test( Test_ExtrapolateIntegrationPointVariableToNodes() );

    /// retrieve barycenter of element face
  _test( Test_FaceBaryCenter() );

} // run




bool FiniteElementPolicy_Test::Test_RstToXYZ()
 {
    Point<3U> rst( 0.5,0.5,0.5 );
    Point<3U> xyz = e_ptr_->RstToXYZ( rst );
    _test( xyz == Point<3U>(1.5,2.1,0.7) );
    
    return xyz == Point<3U>(1.5,2.1,0.7);
 }



    /// integrates a scalar property over the element and returns this value
bool FiniteElementPolicy_Test::Test_PropertyIntegral()
 {
    const csmp::Index key = model_ptr_->Database().StorageKey("nodal variable 1");
    
    // for uniform value
    const double value{ 9.99e+7 };
    model_ptr_->InputPropertyValue( "nodal variable 1", makeScalar(ANY, value ) );
    double integral_over_p = e_ptr_->PropertyIntegral( key );
    _test( approximatelyEqual( integral_over_p, value * e_ptr_->Volume() ) );
  
    return approximatelyEqual( integral_over_p, value * e_ptr_->Volume() );
 }



    /// interpolates node property values to the point of interest (in global coordinates); @attention costly for isoparametric elements
bool FiniteElementPolicy_Test::Test_PropertyValueAt()
 {
    bool passed_all_tests{ true };
 
    // for uniform scalar value
    const double value{ 9.99e+7 };
    model_ptr_->InputPropertyValue( "nodal variable 1", makeScalar(ANY, value ) );
    vector<double> xyz{ 1., 1., 1. };
    ScalarVariable sc;

    const csmp::Index key = model_ptr_->Database().StorageKey("nodal variable 1");
    e_ptr_->PropertyValueAt( key, xyz, sc );
    _test( approximatelyEqual( sc(), value ) );
    passed_all_tests = approximatelyEqual( sc(), value );

    // for uniform vector value
    const VectorVariable<3U> vec(ANY,ANY,ANY, 3., 2., 1. );
    model_ptr_->InputPropertyValue( "nodal vector 1", vec );

    const csmp::Index vkey = model_ptr_->Database().StorageKey("nodal vector 1");
    VectorVariable<3U> test_vec;
    e_ptr_->PropertyValueAt( vkey, xyz, test_vec );
    _test( vec == test_vec );
    passed_all_tests = (vec == test_vec);

    // for uniform tensor value
    TensorVariable<3U> ts;
    ts      = 0.;
    ts(0,0) = 1.;
    ts(1,1) = 2.;
    ts(2,2) = 3.;
    model_ptr_->InputPropertyValue( "nodal tensor 1", ts );

    const csmp::Index tkey = model_ptr_->Database().StorageKey("nodal tensor 1");
    TensorVariable<3U> test_tensor;
    e_ptr_->PropertyValueAt( vkey, xyz, test_tensor );
    _test( ts == test_tensor );
    passed_all_tests = ( ts == test_tensor );
    
    // TODO: test for spatially variable values (x,y,z gradients or other)

    return passed_all_tests;
 }
  
  
  
    /// returns the value of any property interpolated to the element's center of gravity
bool FiniteElementPolicy_Test::Test_PropertyValueAtBaryCenter()
 {
    bool passed_all_tests{ true };

    const csmp::Index skey = model_ptr_->Database().StorageKey("nodal variable 1");
    const csmp::Index vkey = model_ptr_->Database().StorageKey("nodal vector 1");
    const csmp::Index tkey = model_ptr_->Database().StorageKey("nodal tensor 1");

    ScalarVariable     sc(ANY,2.), sc_test;
    VectorVariable<3U> vc(ANY,ANY,ANY,1.,2.,3.), vc_test;
    TensorVariable<3U> ts(ANY,ANY,ANY,1.,0.,0.,0.,2.,0.,0.,0.,3.), ts_test;
    model_ptr_->InputPropertyValue( "nodal variable 1", sc );
    model_ptr_->InputPropertyValue( "nodal vector 1", vc );
    model_ptr_->InputPropertyValue( "nodal tensor 1", ts );
    
    // scalars
    e_ptr_->PropertyValueAtBaryCenter( skey, sc_test );
    _test( approximatelyEqual( sc(), sc_test() ) );
    passed_all_tests = approximatelyEqual( sc(), sc_test() );
    _test( approximatelyEqual( sc(), e_ptr_->PropertyValueAtBaryCenter( skey ) ) );
    passed_all_tests = approximatelyEqual( sc(), e_ptr_->PropertyValueAtBaryCenter( skey ) );

    // vectors
    e_ptr_->PropertyValueAtBaryCenter( vkey, vc_test );
    _test( vc == vc_test );
    passed_all_tests = ( vc == vc_test );

    // tensors
    TensorVariable<3U> tensor, test_tensor;
    e_ptr_->PropertyValueAtBaryCenter( tkey, test_tensor );
    //e_ptr_->Read( tkey, tensor );
    _test( tensor == test_tensor );
    passed_all_tests = ( tensor == test_tensor );
    
    // TODO: make same test for gradient fields

    return passed_all_tests;
 }
  
  
  
    /// returns the value of any property interpolated to the integration point of interest
bool FiniteElementPolicy_Test::Test_PropertyValueAtIntegrationPoint()
 {
    bool passed_all_tests{ true };

    const csmp::Index skey = model_ptr_->Database().StorageKey("nodal variable 1");
    const csmp::Index vkey = model_ptr_->Database().StorageKey("nodal vector 1");
    const csmp::Index tkey = model_ptr_->Database().StorageKey("nodal tensor 1");

    ScalarVariable     sc(ANY,2.), sc_test;
    VectorVariable<3U> vc(ANY,ANY,ANY,1.,2.,3.), vc_test;
    TensorVariable<3U> ts(ANY,ANY,ANY,1.,0.,0.,0.,2.,0.,0.,0.,3.), ts_test;
    model_ptr_->InputPropertyValue( "nodal variable 1", sc );
    model_ptr_->InputPropertyValue( "nodal vector 1", vc );
    model_ptr_->InputPropertyValue( "nodal tensor 1", ts );
    
    // scalars
    e_ptr_->PropertyValueAtIntegrationPoint( skey, 2U, sc_test );
    _test( approximatelyEqual( sc(), sc_test() ) );
    passed_all_tests = approximatelyEqual( sc(), sc_test() );

    // vectors
    e_ptr_->PropertyValueAtIntegrationPoint( vkey, 1U, vc_test );
    _test( vc == vc_test );
    passed_all_tests = ( vc == vc_test );

    // tensors
    TensorVariable<3U> tensor, test_tensor;
    e_ptr_->PropertyValueAtIntegrationPoint( tkey, 3U, test_tensor );
    _test( ts == test_tensor );
    passed_all_tests = ( ts == test_tensor );
    
    // TODO: make same test for gradient fields

    return passed_all_tests;
 }



    /// returns property values at the integration points
bool FiniteElementPolicy_Test::Test_IntegrationPointPropertyVector()
 {
    bool passed_all_tests{ true };

    const csmp::Index skey = model_ptr_->Database().StorageKey("eip scalar 1");
    ScalarVariable sc(ANY,2.);
    model_ptr_->InputPropertyValue( "eip scalar 1", sc );

    // scalar variable
    vector<ScalarVariable>  sc_test_vec;
    e_ptr_->IntegrationPointPropertyVector( skey, sc_test_vec );
    for ( const auto& it : sc_test_vec ) {
         _test( it == sc );
         passed_all_tests = ( it == sc );
      }

    return passed_all_tests;
 }



    /// uses linear extrapolation of property values from the integration points to the nodes
bool FiniteElementPolicy_Test::Test_ExtrapolateIntegrationPointVariableToNodes()
 {
    bool passed_all_tests{ true };

    const csmp::Index skey = model_ptr_->Database().StorageKey("eip scalar 1");
    ScalarVariable sc(ANY,2.);
    model_ptr_->InputPropertyValue( "eip scalar 1", sc );

    const double   var_val{ 4.5 };
    vector<double> ipoint_var_vec( e_ptr_->IntegrationPoints(), var_val );
    vector<double> node_var_vec( e_ptr_->Nodes() );
    uint32_t nvars{ 1U };
    e_ptr_->ExtrapolateIntegrationPointVariableToNodes( nvars, ipoint_var_vec, node_var_vec );
    
    for ( const auto& nit : node_var_vec ) {
         _test( nit == var_val );
         passed_all_tests = ( nit == var_val );
      }

    return passed_all_tests;
 }



    /// retrieve barycenter of element face
bool FiniteElementPolicy_Test::Test_FaceBaryCenter()
 {
    throw csmp::Exception( ERROR, "FiniteElementPolicy_Test::Test_FaceBaryCenter", "test not fully implemented yet");
 
    for ( uint32_t i{0U}; i<e_ptr_->Faces(); ++i ) {
          Point<3U> pt = e_ptr_->FaceBaryCenter(i);
      }

    bool passed_all_tests{ true };

    return passed_all_tests;
 }



} // csmp
