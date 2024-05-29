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
    InitialiseModel();
    // debugging
    assert( model_ptr_ != NULL );
    CreateLinearScalarNodePropertyVariation();
    // testing
    const csmp::Index key = model_ptr_->Database().StorageKey("nodal variable 2");
    vector<ScalarVariable> nvar;
    e_ptr_->NodePropertyVector( key, nvar );
    if ( verbose_ ) {
         // nodal variable values in rows
         cout <<"\nValues of 'nodal variable' varying linearly among the nodes:\n";
         for ( auto i{0U}; i<nvar.size(); ++i )
           cout <<"\n\t"<<"node "<< i <<": "<< nvar[i]();
         e_ptr_->OutputPropertyToVTK( key, "FiniteElementPolicy_Test_prop_variation", "nodal variable 2" );
      }
  }



FiniteElementPolicy_Test::~FiniteElementPolicy_Test()
  {
     delete model_ptr_;
  }



void FiniteElementPolicy_Test::InitialiseModel()
 {
    VSet<3U> vset;
    create_Prism_Hexa_VSet( vset, false );
    model_ptr_ = new Model<3U>( vset, "Variables_Test.txt" );
    e_ptr_     = model_ptr_->Region("Model").E(0);
 }
    


/**
    Creates linear variation of 'node variable 2' in element with a value of 1 at origin of local coordinate system
*/
void FiniteElementPolicy_Test::CreateLinearScalarNodePropertyVariation()
 {
    const csmp::Index key = model_ptr_->Database().StorageKey("nodal variable 2");
    
    // since not all elements have local coordinates and because there is no access to the local coordinates of the element
    // the global coordinates are being used
    Point<3U> bctr = e_ptr_->BaryCenter(); // of element with property value of 1
    for ( auto i{0U}; i<e_ptr_->Nodes(); ++i ) {
          Point<3U> offset = e_ptr_->N(i)->Coordinate() - bctr;
          // computing the nodal value from the nodes offset from the barycentre where the property value is one,
          // assuming a unit property gradient (1,1,1)
          const Point<3U> grad(1.,1.,1.);
          double n_value = 1. + dotProduct( offset, grad );
          // storing the new value
          e_ptr_->N(i)->Store( key, makeScalar(ANY,n_value) );
      }
 }




/// using 'create_Prism_Hexa_VSet' dataset, tests fem operations against provided data in text file
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
    // centroid of linear prism in parametric space
    Point<3U> rst( 1./3.,1./3., 0. );
    Point<3U> xyz = e_ptr_->RstToXYZ( rst );
    
    // barycenter of element 1 in physical space (double precision on Clang
    _test( xyz == Point<3U>(1.3333333333333335,1.3333333333333335,1.5) );
    return xyz == Point<3U>(1.3333333333333335,1.3333333333333335,1.5);
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
    vector<double> xyz{ 1./3., 1./3., 0.5 }; // offset by 0.5 from barycentre along Z-axis
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
    e_ptr_->PropertyValueAt( tkey, xyz, test_tensor );
    _test( ts == test_tensor );
    passed_all_tests = ( ts == test_tensor );
    
    // TODO: test for spatially variable values (x,y,z gradients or other)

    return passed_all_tests;
 }
  
  
  
    /// returns the value of any property interpolated to the element's center of gravity
bool FiniteElementPolicy_Test::Test_PropertyValueAtBaryCenter()
 {
    bool passed_all_tests{ true };

    const csmp::Index skey  = model_ptr_->Database().StorageKey("nodal variable 1");
    const csmp::Index skey2 = model_ptr_->Database().StorageKey("nodal variable 2");
    const csmp::Index vkey  = model_ptr_->Database().StorageKey("nodal vector 1");
    const csmp::Index tkey  = model_ptr_->Database().StorageKey("nodal tensor 1");

    ScalarVariable     sc(ANY,2.), sc_test;
    VectorVariable<3U> vc(ANY,ANY,ANY,1.,2.,3.), vc_test;
    TensorVariable<3U> ts(ANY,ANY,ANY,1.,0.,0.,0.,2.,0.,0.,0.,3.), ts_test;
    model_ptr_->InputPropertyValue( "nodal variable 1", sc );
    model_ptr_->InputPropertyValue( "nodal vector 1", vc );
    model_ptr_->InputPropertyValue( "nodal tensor 1", ts );
    
    // barycentre
    const vector<double> bctr{ 1./3., 1./3., 0. }; // of linear prism element
    
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
    TensorVariable<3U> test_tensor;
    e_ptr_->PropertyValueAtBaryCenter( tkey, test_tensor );
    //e_ptr_->Read( tkey, tensor );
    _test( ts == test_tensor );
    passed_all_tests = ( ts == test_tensor );
    
    // test for gradient field
    assert( e_ptr_->FE_Type() == ISOPARAMETRIC_LINEAR_PRISM );
    vector<double> xyz{ 1.33333, 1.33333, 1.5 }; // barycentre of prism element 1

    const double at_value{ e_ptr_->PropertyValueAt( skey2, xyz ) },
                 bc_value{ e_ptr_->PropertyValueAtBaryCenter( skey2 ) };

    _test( approximatelyEqual( bc_value, 1. ) );
    _test( approximatelyEqual( at_value, bc_value ) );
    passed_all_tests = approximatelyEqual( bc_value, 1. );
    passed_all_tests = approximatelyEqual( at_value, bc_value );

    return passed_all_tests;
 }
  
  
  
    /// returns the value of any property interpolated to the integration point of interest
bool FiniteElementPolicy_Test::Test_PropertyValueAtIntegrationPoint()
 {
    bool passed_all_tests{ true };

    const csmp::Index skey  = model_ptr_->Database().StorageKey("nodal variable 1");
    const csmp::Index skey2 = model_ptr_->Database().StorageKey("nodal variable 2");
    const csmp::Index vkey  = model_ptr_->Database().StorageKey("nodal vector 1");
    const csmp::Index tkey  = model_ptr_->Database().StorageKey("nodal tensor 1");

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
    TensorVariable<3U> test_tensor;
    e_ptr_->PropertyValueAtIntegrationPoint( tkey, 3U, test_tensor );
    _test( ts == test_tensor );
    passed_all_tests = ( ts == test_tensor );
    
    // same test for scalar gradient field
    e_ptr_->CoordinateMatrix();
    for ( uint32_t ip{0U}; ip<e_ptr_->IntegrationPoints(); ++ip )
      {
         vector<double>  xyz;
         e_ptr_->FE()->IntegrationPoint( ip, xyz );
         const double at_value = e_ptr_->PropertyValueAt( skey2, xyz );
         const double ip_value = e_ptr_->PropertyValueAtIntegrationPoint( skey2, ip );
         const double epsilon{ 1.0e-7 };
         _test( approximatelyEqual( at_value, ip_value, epsilon ) );
         passed_all_tests = approximatelyEqual( at_value, ip_value, epsilon );
      }

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
    // barycentres of the first 3 faces of the element will be tested
    bool passed_all_tests{ true };
    
    // face 0
    uint32_t face_id{ 0U };
    set<Node<3U>*> fnodes = e_ptr_->CornerNodesOfFace( face_id );
    // computing face barycentre
    Point<3U> fbctr; fbctr = 0.;
    for ( const auto nit : fnodes ) {
         fbctr[0] += nit->x();
         fbctr[1] += nit->y();
         fbctr[2] += nit->z();
      }
    fbctr /= static_cast<double>(fnodes.size());
    // TEST
    Point<3U> pt = e_ptr_->FaceBaryCenter( face_id );
    _test( pt == fbctr );
    passed_all_tests = ( pt == fbctr );

    // face 1
    face_id = 1U;
    fnodes = e_ptr_->CornerNodesOfFace( face_id );
    // computing face barycentre
    fbctr = 0.;
    for ( const auto nit : fnodes ) {
         fbctr[0] += nit->x();
         fbctr[1] += nit->y();
         fbctr[2] += nit->z();
      }
    fbctr /= static_cast<double>(fnodes.size());
    // TEST
    pt = e_ptr_->FaceBaryCenter( face_id );
    _test( pt == fbctr );
    passed_all_tests = ( pt == fbctr );

    // face 2
    face_id = 2U;
    fnodes = e_ptr_->CornerNodesOfFace( face_id );
    // computing face barycentre
    fbctr = 0.;
    for ( const auto nit : fnodes ) {
         fbctr[0] += nit->x();
         fbctr[1] += nit->y();
         fbctr[2] += nit->z();
      }
    fbctr /= static_cast<double>(fnodes.size());
    // TEST
    pt = e_ptr_->FaceBaryCenter( face_id );
    _test( pt == fbctr );
    passed_all_tests = ( pt == fbctr );

    return passed_all_tests;
 }



} // csmp
