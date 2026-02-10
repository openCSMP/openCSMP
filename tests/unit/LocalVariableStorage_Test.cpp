#include "LocalVariableStorage_Test.h"

#include "LocalVariableStorage.h"
#include "Model.h"
#include "Region.h"
#include "Element.h"
#include "LinearTriangle.h"
#include "LinearTetrahedron.h"
#include "vsetMakers.h"

using namespace std;

namespace csmp{

/**
 @fn  void LocalVariableStorage_Test::runTest()

 @brief Executes the test operation.
 @todo (1-T) Test integration point access (both simplex & fv)

 @author  P. Lang
 @date  9/24/2012
 */
template<uint32_t dim>
void LocalVariableStorage_Test::runTest()
  {
    // initialize storage with 3 of each (scalars, vectors, tensors)
    LinearTriangle triangle;
    Element<dim>   e(&triangle);

    // LOCAL VARIABLES ONLY
    
    // one var of each type, array has length 20
    const size_t scalars1(1), vectors1(1), tensors1(1), arrays1(1), arrayLength(20),flaggedArrays1(1),flaggedArrayLength(40);
    const size_t totalDataDepth1( scalars1 + vectors1*dim + tensors1*dim*dim + arrays1*arrayLength + flaggedArrays1*flaggedArrayLength );
    const size_t totalFlagDepth1( scalars1 + vectors1*dim + tensors1*dim + arrays1 + flaggedArrays1*flaggedArrayLength );
    LocalVariables lv1( scalars1, vectors1, tensors1, arrays1, arrayLength, flaggedArrays1, flaggedArrayLength, totalDataDepth1, totalFlagDepth1 );
    e.ResizePropertyStorage(lv1);

    
    // some working indices
    Index ks1( SCALAR, ELEMENT, 0, 1, 1, 0, 0, lv1 );
    Index kv1( VECTOR, ELEMENT, 0, dim, dim, 1, 1, lv1 );
    Index kt1( TENSOR, ELEMENT, 0, dim*dim, dim, dim+1, dim+1, lv1 );
    Index ka1( ARRAY,  ELEMENT, 0, arrayLength, 1, dim+dim*dim+1, dim+dim+1, lv1 );
    Index kfa1( FLAGGEDARRAY,  ELEMENT, 0, flaggedArrayLength, flaggedArrayLength, dim+dim*dim+1+arrayLength, dim+dim+1+1, lv1 );


    // some working variables
    ScalarVariable s0;
    const ScalarVariable s1( ROBIN, 1. );
    VectorVariable<dim> v0;
    const VectorVariable<dim> v1( ANY, 2. );
    TensorVariable<dim> t0;
    const TensorVariable<dim> t1( NEUMANN, 3. );
    ArrayVariable a0( arrayLength );
    const ArrayVariable a1( arrayLength, 4., DIRICH );
    FlaggedArrayVariable fa0( flaggedArrayLength );
    const FlaggedArrayVariable fa1( flaggedArrayLength, 5., DIRICH );

    // STORE/READ

    // scalars
    e.Store( ks1, s1 );
    _test( e.Read(ks1) == s1() );
    e.Read( ks1, s0 );
    _test( s0 == s1 );

    // vectors
    e.Store( kv1, v1 );
    e.Read( kv1, v0 );
    _test( v0 == v1 );

    // tensors
    e.Store( kt1, t1 );
    e.Read( kt1, t0 );
    _test( t0 == t1 );

    // arrays
    e.Store( ka1, a1 );
    e.Read( ka1, a0 );
    _test( a0 == a1 );

    // flagged arrays
    e.Store( kfa1, fa1 );
    e.Read( kfa1, fa0 );
    _test( fa0 == fa1 );

    // STATUS
  
    _test( e.Status(ks1) == s1.Flag() );
    _test( e.Status(ka1) == a1.Flag() );
    for( uint32_t d(0); d < flaggedArrayLength; ++d )
        _test( e.Status(kfa1,d) == fa1.Flag(d) );
    for( uint32_t d(0); d < dim; ++d )
      {
        _test( e.Status( kv1, d ) == v1.Flag(d) );
        _test( e.Status( kt1, d ) == t1.Flag(d) );
      }

    e.Status(ks1, INIT_GUESS );
    e.Status(ka1, INIT_COND);
    for( uint32_t d(0); d < flaggedArrayLength; ++d )
        e.Status(kfa1,d, FIELD_DATA);
    for( uint32_t d(0); d < dim; ++d )
      {
      e.Status( kv1, d, CONSTANT_FLUX );
      e.Status( kt1, d, PERIODIC );
      }

    _test( e.Status(ks1) == INIT_GUESS );
    _test( e.Status(ka1) == INIT_COND );
    for( uint32_t d(0); d < flaggedArrayLength; ++d )
        _test( e.Status(kfa1,d) == FIELD_DATA );
    for( uint32_t d(0); d < dim; ++d )
      {
      _test( e.Status( kv1, d ) == CONSTANT_FLUX );
      _test( e.Status( kt1, d ) == PERIODIC );
      }

      
    // reset
    e.Store( ks1, s1 );
    e.Store( kv1, v1 );
    e.Store( kt1, t1 );
    e.Store( ka1, a1 );
    e.Store( kfa1, fa1 );


    // WITHIN
//      _test( e.IsWithinRange( kt1, 2., 4. ) );
    _test( e.IsWithinRange( kt1, 0., 4. ) ); // SKM FIX - since every element is tested, commented out version fails on off-diagonal zeros
    _test( !e.IsWithinRange( kt1, 9., 10. ) );
    _test( !e.IsWithinRange( kt1, 1., 2. ) );

    // DELETE
    e.DeleteProperty(kt1);
    ka1 = Index( ARRAY, ELEMENT, 0, arrayLength, 1, dim+1, dim+1, lv1 );
    e.Read( ka1, a0 );
    _test( a0 == a1 );

    // ADD
    e.AddProperty(kt1);
    ka1 = Index( ARRAY, ELEMENT, 0, arrayLength, 1, dim+dim*dim+1, dim+dim+1, lv1 );
    e.Read( ka1, a0 );
    _test( a0 == a1 );


// FORMERLY COMMENTED OUT?

    // TODO: adding Properties
    /*
    Index ks2, kv3, kt3;
    e.AddProperty(ks3);
    e.AddProperty(kv3);
    e.AddProperty(kt3);

    // storing
    e.Store( ks3, s1 );
    e.Store( kv3, v1 );
    e.Store( kt3, t1 );

    // reading scalars
    e.Read( ks1, s0 );
    _test( s0 == s2 );
    e.Read( ks2, s0 );
    _test( s0 == s3 );
    e.Read( ks3, s0 );
    _test( s0 == s1 );

    // reading vectors
    e.Read( kv1, v0 );
    _test( v0 == v2 );
    e.Read( kv2, v0 );
    _test( v0 == v3 );
    e.Read( kv3, v0 );
    _test( v0 == v1 );

    // reading tensors
    e.Read( kt1, t0 );
    _test( t0 == t2 );
    e.Read( kt2, t0 );
    _test( t0 == t3 );
    e.Read( kt3, t0 );
    _test( t0 == t1 );
    */
    
    // empty
    _test( !e.EmptyLVS() );

// TODO: Is there a problem with  ResizePropertyStorage() the following code fails for v0,t0 (2D), and t0 (3D) in DEBUG mode
/*
    // does copy construction of element preserve variable storage
    Element<dim>  e2(e);
    
if ( dim == 3 ) {
  cout <<"\n3D model";
  }

    // reading scalars
    e2.Read( ks1, s0 );
    _test( s0 == s1 );

    // reading vectors
    e2.Read( kv1, v0 );
    _test( v0 == v1 );

    // reading tensors
    e2.Read( kt1, t0 );
    _test( t0 == t1 );
*/

/* TODO: following code does not compile, constructor arguments are missing

    // initialize storage with 3 of each (scalars, vectors, tensors) and 2 arrays with 14 and 6
    Element<dim> e3;
    e3.ResizePropertyStorage( 3, 3, 3, 2, 20 );

    Index ka1_( ARRAY, ELEMENT, 0, 14, 0 ); // offset to data start / size / offset to flag
    Index ka2_( ARRAY, ELEMENT, 14, 6, 1 );

    ArrayVariable av1( 14, 10., PLAIN );
    ArrayVariable av2( 6, 11., DIRICH );

    e3.Store( ka1_, av1 );
    e3.Store( ka2_, av2 );

    ArrayVariable av0;
    e3.Read( ka1_, av0 );
    _test( av0 == av1 );
    e3.Read( ka2_, av0 );
    _test( av0 == av2 );

    // deleting scalar, vector and tensor
    e.DeleteProperty(ks1);
    e.DeleteProperty(kv1);
    e.DeleteProperty(kt1);

    // reading
    e3.Read( ka1_, av0 );
    _test( av0 == av1 );
    e3.Read( ka2_, av0 );
    _test( av0 == av2 );

    Index ka3_( ARRAY, ELEMENT, 20, 8, 2 );
    e3.AddProperty(ka3_);
    ArrayVariable av3(8);
    e3.Read( ka3_, av0 );
    _test( av0 == av3 );

    Index ka4_( ARRAY, ELEMENT, 28, 16, 3 );
    e3.AddProperty(ka4_);
    av3.Size(16);
    e3.Read( ka4_, av0 );
    _test( av0 == av3 );
    av3 = ArrayVariable( 16, 3. );
    av3.Flag(DIRICH);
    e3.Store( ka4_, av3 );
    av0.Size(16);
    e3.Read( ka4_, av0 );
    _test( av0 == av3 );

    // deleting array
    av1 = ArrayVariable( 14, 6., PLAIN );
    av2 = ArrayVariable( 6, 12., DIRICH );
    av3 = ArrayVariable( 8, 24., PLAIN );
    ArrayVariable av4( 16, 48., DIRICH );
    e3.Store( ka1_, av1 );
    e3.Store( ka2_, av2 );
    e3.Store( ka3_, av3 );
    e3.Store( ka4_, av4 );
    e3.DeleteProperty(ka3_); // now indices have to be updated for offset, so that a1 = a1, a2 = a2, a3 = a4

    // updating arrays
    ka4_ = ka3;
    ka4_.dataDepth = 16;
    ka4_.flagOffset = 2;

    // read & test
    av0.Size(16);
    e3.Read( ka4_, av0 );
    _test( av0 == av4 );
    av0.Size(14);
    e3.Read( ka1_, av0 );
    _test( av0 == av1 );
    av0.Size(6);
    e3.Read( ka2_, av0 );
    _test( av0 == av2 );

    // delete first scalar (first entry in data vector)
    e3.Store( ks2_, s2 );
    e3.DeleteProperty(ks1);
    _test( e3.Read(ks1) == s2() );

    // within
    _test( e3.IsWithinRange( ka4, 47., 50. ) );
    _test( !e3.IsWithinRange( ka4, 41., 45. ) );
    _test( !e3.IsWithinRange( ka2, 9., 10. ) );
    _test( !e3.IsWithinRange( ka2, 15., 18. ) );
    _test( e3.IsWithinRange( ka2, 10., 18. ) );

    // Resize from db
    PropertyDatabase<dim> pdb("CSMP-variables-vsTestLocked.txt");
    e3.ResizePropertyStorage( pdb, ELEMENT );
    _test( e3.LVS().scalars == 1 );
    _test( e3.LVS().vectors == 1 );
    _test( e3.LVS().tensors == 0 );
    _test( e3.LVS().arrays == 2 );
    _test( e3.LVS().arrayLength == 26 );

*/
  } // end runTest






void LocalVariableStorage_Test::run3D()
  {
    // initialize storage with 3 of each (scalars, vectors, tensors)
    const size_t      dim(3);
    LinearTetrahedron tetrahedron;
    Element<3>        e(&tetrahedron);

    // LOCAL VARIABLES ONLY

    // one var of each type, array has length 20
    const size_t scalars1(1), vectors1(1), tensors1(1), arrays1(1), arrayLength(20),flaggedArrays1(1),flaggedArrayLength(40);
    const size_t totalDataDepth1( scalars1 + vectors1*dim + tensors1*dim*dim + arrays1*arrayLength + flaggedArrays1*flaggedArrayLength );
    const size_t totalFlagDepth1( scalars1 + vectors1*dim + tensors1*dim + arrays1 + flaggedArrays1*flaggedArrayLength);
    LocalVariables lv1( scalars1, vectors1, tensors1, arrays1, arrayLength, flaggedArrays1, flaggedArrayLength, totalDataDepth1, totalFlagDepth1 );
    e.ResizePropertyStorage(lv1);

    // some working indices
    const Index ks1( SCALAR, ELEMENT, 0, 1, 1, 0, 0, lv1 );
    const Index kv1( VECTOR, ELEMENT, 0, dim, dim, 1, 1, lv1 );
    const Index kt1( TENSOR, ELEMENT, 0, dim*dim, dim, dim+1, dim+1, lv1 );
    Index ka1( ARRAY,  ELEMENT, 0, arrayLength, 1, dim+dim*dim+1, dim+dim+1, lv1 ); // needs to be assignable, see below
    const Index kfa1( FLAGGEDARRAY,  ELEMENT, 0, flaggedArrayLength, flaggedArrayLength, dim+dim*dim+1+arrayLength, dim+dim+1+1, lv1 );


    // some working variables
    ScalarVariable s0;
    const ScalarVariable s1( ROBIN, 1. );
    VectorVariable<3> v0;
    const VectorVariable<3> v1( ANY, DIRICH, ROBIN, 2., 1., 0. );
    const VectorVariable<3> v2( ANY, DIRICH, ROBIN, 2., 1.1, 0. );
    TensorVariable<3> t0;
    const TensorVariable<3> t1( NEUMANN, DIRICH, NEUMANN, 3., 4., 5., 6., 7., 8., 9., 10., 11. );
    ArrayVariable a0( arrayLength );
    const ArrayVariable a1( arrayLength, 4., DIRICH );
    FlaggedArrayVariable fa0( flaggedArrayLength );
    const FlaggedArrayVariable fa1( flaggedArrayLength, 5., DIRICH );

    // STORE/READ

    // scalars
    e.Store( ks1, s1 );
    _test( e.Read(ks1) == s1() );
    e.Read( ks1, s0 );
    _test( s0 == s1 );

    // vectors
    e.Store( kv1, v1 );
    e.Read( kv1, v0 );
    _test( v0 == v1 );

    // tensors
    e.Store( kt1, t1 );
    e.Read( kt1, t0 );
    _test( t0 == t1 );

    // arrays
    e.Store( ka1, a1 );
    e.Read( ka1, a0 );
    _test( a0 == a1 );

    // flagged arrays
    e.Store( kfa1, fa1 );
    e.Read( kfa1, fa0 );
    _test( fa0 == fa1 );

    // STATUS

    _test( e.Status(ks1) == s1.Flag() );
    _test( e.Status(ka1) == a1.Flag() );
    for( uint32_t d(0); d < flaggedArrayLength; ++d )
        _test( e.Status(kfa1,d) == fa1.Flag(d) );
    for( uint32_t d(0); d < dim; ++d )
      {
        _test( e.Status( kv1, d ) == v1.Flag(d) );
        _test( e.Status( kt1, d ) == t1.Flag(d) );
      }

    e.Status(ks1, INIT_GUESS );
    e.Status(ka1, INIT_COND );
    for( uint32_t d(0); d < flaggedArrayLength; ++d )
        e.Status(kfa1,d, FIELD_DATA);
    for( uint32_t d(0); d < dim; ++d )
      {
      e.Status( kv1, d, CONSTANT_FLUX );
      e.Status( kt1, d, PERIODIC );
      }

    _test( e.Status(ks1) == INIT_GUESS );
    _test( e.Status(ka1) == INIT_COND );
    for( uint32_t d(0); d < flaggedArrayLength; ++d )
        _test( e.Status(kfa1,d) == FIELD_DATA );
    for( uint32_t d(0); d < dim; ++d )
      {
      _test( e.Status( kv1, d ) == CONSTANT_FLUX );
      _test( e.Status( kt1, d ) == PERIODIC );
      }


    // reset
    e.Store( ks1, s1 );
    e.Store( kv1, v1 );
    e.Store( kt1, t1 );
    e.Store( ka1, a1 );
    e.Store( kfa1, fa1 );


    // DELETE
    e.DeleteProperty(kt1);
    ka1 = Index( ARRAY, ELEMENT, 0, arrayLength, 1, dim+1, dim+1, lv1 );
    e.Read( ka1, a0 );
    _test( a0 == a1 );
    e.Read( ks1, s0 );
    _test( s0 == s1 );
    e.Read( kv1, v0 );
    _test( v0 == v1 );

    // ADD
    e.AddProperty(kt1);
    ka1 = Index( ARRAY, ELEMENT, 0, arrayLength, 1, dim+dim*dim+1, dim+dim+1, lv1 );
    e.Read( ka1, a0 );
    _test( a0 == a1 );
    e.Store( kt1, t1 );
    e.Read( ks1, s0 );
    _test( s0 == s1 );
    e.Read( kv1, v0 );
    _test( v0 == v1 );
    e.Read( kt1, t0 );
    _test( t0 == t1 );

    e.Store( kv1, v2 );
    e.Read( kv1, v0 );
    _test( v0 != v1 );

  } // end run3D
  
  
  
  
void LocalVariableStorage_Test::run3D_with_templatized_INDEX()
  {
    // initialize storage with 3 of each (scalars, vectors, tensors)
    const size_t      dim(3);
    LinearTetrahedron tetrahedron;
    Element<3>        e(&tetrahedron);

    // LOCAL VARIABLES ONLY

    // one var of each type, array has length 20
    const size_t scalars1(1), vectors1(1), tensors1(1), arrays1(1), arrayLength(20),flaggedArrays1(1),flaggedArrayLength(40);
    const size_t totalDataDepth1( scalars1 + vectors1*dim + tensors1*dim*dim + arrays1*arrayLength + flaggedArrays1*flaggedArrayLength );
    const size_t totalFlagDepth1( scalars1 + vectors1*dim + tensors1*dim + arrays1 + flaggedArrays1*flaggedArrayLength);
    const LocalVariables lv1( scalars1, vectors1, tensors1, arrays1, arrayLength, flaggedArrays1, flaggedArrayLength, totalDataDepth1, totalFlagDepth1 );
    e.ResizePropertyStorage(lv1);

    // some working indices
    const INDEX<SCALAR,ELEMENT>        ks1( 0, 1, 1, 0, 0, lv1, IntegrationPointVariables() );
    const INDEX<VECTOR,ELEMENT>        kv1( 0, dim, dim, 1, 1, lv1, IntegrationPointVariables() );
    const INDEX<TENSOR,ELEMENT>        kt1( 0, dim*dim, dim, dim+1, dim+1, lv1, IntegrationPointVariables() );
    INDEX<ARRAY,ELEMENT>               ka1( 0, arrayLength, 1, dim+dim*dim+1, dim+dim+1, lv1, IntegrationPointVariables() ); // assignable
    const INDEX<FLAGGEDARRAY,ELEMENT>  kfa1( 0, flaggedArrayLength, flaggedArrayLength,
                                             dim+dim*dim+1+arrayLength, dim+dim+1+1, lv1, IntegrationPointVariables() );

    // some working variables
    ScalarVariable s0;
    const ScalarVariable s1( ROBIN, 1. );
    VectorVariable<3> v0;
    const VectorVariable<3> v1( ANY, DIRICH, ROBIN, 2., 1., 0. );
    const VectorVariable<3> v2( ANY, DIRICH, ROBIN, 2., 1.1, 0. );
    TensorVariable<3> t0;
    const TensorVariable<3> t1( NEUMANN, DIRICH, NEUMANN, 3., 4., 5., 6., 7., 8., 9., 10., 11. );
    ArrayVariable a0( arrayLength );
    const ArrayVariable a1( arrayLength, 4., DIRICH );
    FlaggedArrayVariable fa0( flaggedArrayLength );
    const FlaggedArrayVariable fa1( flaggedArrayLength, 5., DIRICH );

    // STORE/READ

    // scalars
    e.Store( ks1, s1 );
    _test( e.Read(ks1) == s1() );
    e.Read( ks1, s0 );
    _test( s0 == s1 );

    // vectors
    e.Store( kv1, v1 );
    e.Read( kv1, v0 );
    _test( v0 == v1 );

    // tensors
    e.Store( kt1, t1 );
    e.Read( kt1, t0 );
    _test( t0 == t1 );

    // arrays
    e.Store( ka1, a1 );
    e.Read( ka1, a0 );
    _test( a0 == a1 );

    // flagged arrays
    e.Store( kfa1, fa1 );
    e.Read( kfa1, fa0 );
    _test( fa0 == fa1 );

    // STATUS

    _test( e.Status(ks1) == s1.Flag() );
    _test( e.Status(ka1) == a1.Flag() );
    for( uint32_t d(0); d < flaggedArrayLength; ++d )
        _test( e.Status(kfa1,d) == fa1.Flag(d) );
    for( uint32_t d(0); d < dim; ++d )
      {
        _test( e.Status( kv1, d ) == v1.Flag(d) );
        _test( e.Status( kt1, d ) == t1.Flag(d) );
      }

    e.Status(ks1, INIT_GUESS );
    e.Status(ka1, INIT_COND );
    for( uint32_t d(0); d < flaggedArrayLength; ++d )
        e.Status(kfa1,d, FIELD_DATA);
    for( uint32_t d(0); d < dim; ++d )
      {
      e.Status( kv1, d, CONSTANT_FLUX );
      e.Status( kt1, d, PERIODIC );
      }

    _test( e.Status(ks1) == INIT_GUESS );
    _test( e.Status(ka1) == INIT_COND );
    for( uint32_t d(0); d < flaggedArrayLength; ++d )
        _test( e.Status(kfa1,d) == FIELD_DATA );
    for( uint32_t d(0); d < dim; ++d )
      {
      _test( e.Status( kv1, d ) == CONSTANT_FLUX );
      _test( e.Status( kt1, d ) == PERIODIC );
      }


    // reset
    e.Store( ks1, s1 );
    e.Store( kv1, v1 );
    e.Store( kt1, t1 );
    e.Store( ka1, a1 );
    e.Store( kfa1, fa1 );


    // DELETE
    e.DeleteProperty(kt1);
    ka1 = INDEX<ARRAY,ELEMENT>( 0, arrayLength, 1, dim+1, dim+1, lv1, IntegrationPointVariables() );
    e.Read( ka1, a0 );
    _test( a0 == a1 );
    e.Read( ks1, s0 );
    _test( s0 == s1 );
    e.Read( kv1, v0 );
    _test( v0 == v1 );

    // ADD
    e.AddProperty(kt1);
    ka1 = INDEX<ARRAY,ELEMENT>( 0, arrayLength, 1, dim+dim*dim+1, dim+dim+1, lv1, IntegrationPointVariables() );
    e.Read( ka1, a0 );
    _test( a0 == a1 );
    e.Store( kt1, t1 );
    e.Read( ks1, s0 );
    _test( s0 == s1 );
    e.Read( kv1, v0 );
    _test( v0 == v1 );
    e.Read( kt1, t0 );
    _test( t0 == t1 );

    e.Store( kv1, v2 );
    e.Read( kv1, v0 );
    _test( v0 != v1 );

  } // end run3D_with_templatized_INDEX
  
  
  
  
/**
    Test constexpr VARIABLE_TYPE variableType( const Var& )
    
    @attention this method does not even need to be called
 */
void LocalVariableStorage_Test::Test_variableType()
 {
    //constexpr double     db=12.;
    ScalarVariable       sc;
    VectorVariable<1>    vc1;
    VectorVariable<2>    vc2;
    TensorVariable<3>    ts3;
    ArrayVariable        ary5( 5 );
    FlaggedArrayVariable fary5( 4 );
    
    //static_assert( variableType(db)    == SCALAR,       "variableType: double should trigger default: ScalarVariable" );
    static_assert( variableType(sc)    == SCALAR,       "variableType: failed on ScalarVariable" );
    static_assert( variableType(vc1)   == VECTOR,       "variableType: failed on VectorVariable" );
    static_assert( variableType(vc2)   == VECTOR,       "variableType: failed on VectorVariable" );
    static_assert( variableType(ts3)   == TENSOR,       "variableType: failed on TensorVariable" );
    static_assert( variableType(ary5)  == ARRAY,        "variableType: failed on ArrayVariable" );
    static_assert( variableType(fary5) == FLAGGEDARRAY, "variableType: failed on FlaggedArrayVariable" );
 
 } // end Test_variableType
  
  


 /// ReadVector, ReadTensor, ReadArray directly returning csmp variables (skm 12/24)
void LocalVariableStorage_Test::runExtensionsForDeclarativeProgramming()
 {
    VSet<3U> vset;
    create_Pyramid_Hexa_VSet( vset, true ); // skewed model
    Model<3U>  model( vset, "LocalVariableStorage_Test-variables.txt" );
    Region<3U> domain = model.Region("Model");
    
    // testing read operations on element, node, element integration points and FV integration points
    // element
    const csmp::Index ev1 = model.Database().StorageKey("element vector 1");
    const csmp::Index et1 = model.Database().StorageKey("element tensor 1");
    const csmp::Index ea1 = model.Database().StorageKey("element array 1");
    // node
    const csmp::Index nv1 = model.Database().StorageKey("nodal vector 1");
    const csmp::Index nt1 = model.Database().StorageKey("nodal tensor 1");
    const csmp::Index na1 = model.Database().StorageKey("nodal array 1");

    // inputs
    VectorVariable<3U> vc(ANY,PLAIN,ROBIN,1.,2.,3.);
    TensorVariable<3U> ts(ANY,PLAIN,ROBIN,1.,2.,3.,4.,5.,6.,7.,8.,9);
    ArrayVariable ary( 4, 0., ANY ), ary_test;
    for ( uint32_t i{0u}; i<ary.Size(); ++i ) ary(i) = static_cast<double>(i);

    // element & node
    model.InputPropertyValue("element vector 1", makeVector(ANY,PLAIN,ROBIN,1.,2.,3.) );
    model.InputPropertyValue("element tensor 1", makeTensor(ANY,PLAIN,ROBIN,1.,2.,3.,4.,5.,6.,7.,8.,9) );
    model.InputPropertyValue("element array 1", ary );

    model.InputPropertyValue("nodal vector 1", makeVector(ANY,PLAIN,ROBIN,1.,2.,3.) );
    model.InputPropertyValue("nodal tensor 1", makeTensor(ANY,PLAIN,ROBIN,1.,2.,3.,4.,5.,6.,7.,8.,9) );
    model.InputPropertyValue("nodal array 1", ary );
    
    // 1. test of functions reading element and node variables
    for ( const auto& e : domain.CellVector() ) {
         _test( vc == e->ReadVector(ev1) );
         _test( ts == e->ReadTensor(et1) );
         e->Read( ea1, ary_test );
         _test( ary == ary_test );
         vector<double> a1_vec1 = e->ReadArray(ea1), a1_vec2(4);
         // testing return value initialisation
         for ( uint32_t i{0u}; i<ary.Size(); ++i )
              _test( a1_vec1[i] == ary[i] );
         // testing ReadArrayEntry()
         for ( uint32_t i{0u}; i<a1_vec1.size(); ++i )
           a1_vec2[i] = e->ReadArrayEntry( ea1, i );
         // checking consistency 
         _test( a1_vec1 == a1_vec2 );
      }
    for ( const auto& n : domain.NodeVector() ) {
         _test( vc == n->ReadVector(nv1) );
         _test( ts == n->ReadTensor(nt1) );
         n->Read( na1, ary_test );
         _test( ary == ary_test );
         vector<double> a1_vec1 = n->ReadArray(na1), a1_vec2(4);
         for ( uint32_t i{0u}; i<ary.Size(); ++i )
           a1_vec2[i] = n->ReadArrayEntry( na1, i );
         _test( a1_vec1 == a1_vec2 );
         for ( uint32_t i{0u}; i<ary.Size(); ++i )
              _test( ary[i] == a1_vec1[i] ); 
      }

    // 2. test of functions reading variables from element integration points
    // ----------------------------------------------------------------------
    // initialising different array variables for up to 8 integration points (hexahedron)
    vector<VectorVariable<3U>> ip_vector_variables  = { {ANY,ANY,ANY,0.,1.,2.}, {PLAIN,PLAIN,PLAIN,3.,4.,5.},
                                                        {ROBIN,ROBIN,ROBIN,6.,7.,8.}, {DIRICH,DIRICH,DIRICH,9.,10.,11.},
                                                        {INIT_COND,INIT_COND,INIT_COND,12.,13.,14.}, {INIT_GUESS,INIT_GUESS,INIT_GUESS,15.,16.,17.},
                                                        {PERIODIC,PERIODIC,PERIODIC,18.,19.,20.}, {FIELD_DATA,FIELD_DATA,FIELD_DATA,21.,22.,23.} };

    vector<TensorVariable<3U>> ip_tensor_variables  = { {ANY,ANY,ANY,0.,1.,2.,3.,4.,5.,6.,7.,8.},
                                                        {PLAIN,PLAIN,PLAIN,9.,10.,11.,12.,13.,14.,15.,16.,17.},
                                                        {ROBIN,ROBIN,ROBIN,18.,19.,20.,21.,22.,23.,24.,25.,26.},
                                                        {DIRICH,DIRICH,DIRICH,27.,28.,29.,30.,31.,32.,33.,34.,35.},
                                                        {INIT_COND,INIT_COND,INIT_COND,36.,37.,38.,39.,40.,41.,42.,43.,44.},
                                                        {INIT_GUESS,INIT_GUESS,INIT_GUESS,45.,46.,47.,48.,49.,50.,51.,52.,53.},
                                                        {PERIODIC,PERIODIC,PERIODIC,54.,55.,56.,57.,58.,59.,60.,61.,62.},
                                                        {FIELD_DATA,FIELD_DATA,FIELD_DATA,63.,64.,65.,66.,67.,68.,69.,70.,71.} };

    vector<ArrayVariable> ip_array_variables  = { {{0.,1.,2.,3.},ANY}, {{4.,5.,6.,7.},PLAIN},
                                                  {{8.,9.,10.,11.},ROBIN}, {{12.,13.,14.,15.},DIRICH},
                                                  {{16.,17.,18.,19.},INIT_COND}, {{20.,21.,22.,23.},INIT_GUESS},
                                                  {{24.,25.,26.,27.},PERIODIC}, {{28.,29.,30.,31.},FIELD_DATA} };

    // element integration points
    // --------------------------
    const csmp::Index eiv1 = model.Database().StorageKey("eip vector 1");
    const csmp::Index eit1 = model.Database().StorageKey("eip tensor 1");
    const csmp::Index eia1 = model.Database().StorageKey("eip array 1");

    // assigning the values
    for ( auto& e : domain.CellVector() )
      for ( uint32_t i{0u}; i<e->IntegrationPoints(); ++i ) {
           e->Store( i, eiv1, ip_vector_variables[i] );
           e->Store( i, eit1, ip_tensor_variables[i] );
           e->Store( i, eia1, ip_array_variables[i] );
        }
    // testing that they are read correctly
    for ( const auto& e : domain.CellVector() )
      for ( uint32_t i{0u}; i<e->IntegrationPoints(); ++i ) {
           auto vc_vec = e->ReadVector( i, eiv1 );
           _test( vc_vec == ip_vector_variables[i] );
           auto ts_vec = e->ReadTensor( i, eit1 );
           _test( ts_vec == ip_tensor_variables[i] );
           auto ar_vec = e->ReadArray( i, eia1 );
           for ( uint32_t j{0u}; j<ar_vec.size(); ++j )
             _test( ar_vec[j] == ip_array_variables[i][j] );
        }

    // FV facet integration points
    // ---------------------------
    // facet integration points
    const csmp::Index fipv1 = model.Database().StorageKey("fip vector 1");
    const csmp::Index fipt1 = model.Database().StorageKey("fip tensor 1");
    const csmp::Index fipa1 = model.Database().StorageKey("fip array 1");
    
    // adding the necessary 4 extra values to the variable arrays
    // vectors
    ip_vector_variables.push_back( {ROBIN,PLAIN,ANY,24.,25.,26} );
    ip_vector_variables.push_back( {ROBIN,PLAIN,ANY,27.,28.,29} );
    ip_vector_variables.push_back( {ROBIN,PLAIN,ANY,30.,31.,32.} );
    ip_vector_variables.push_back( {ROBIN,PLAIN,ANY,33.,34.,35.} );
    // tensors
    ip_tensor_variables.push_back( {ANY,ANY,ANY,72.,73.,74.,75.,76.,77.,78.,79.,80.} );
    ip_tensor_variables.push_back( {PLAIN,PLAIN,PLAIN,81.,82.,83.,84.,85.,86.,87.,88.,89.} );
    ip_tensor_variables.push_back( {ROBIN,ROBIN,ROBIN,90.,91.,92.,93.,94.,95.,96.,97.,98.} );
    ip_tensor_variables.push_back( {DIRICH,DIRICH,DIRICH,99.,100.,101.,102.,103.,104.,105.,106.,107.} );
    // arrays
    ip_array_variables.push_back( {{32.,33.,34.,34.},ANY} );
    ip_array_variables.push_back( {{35.,36.,37.,38.},PLAIN} );
    ip_array_variables.push_back( {{39.,40.,41.,42.},ROBIN} );
    ip_array_variables.push_back( {{43.,44.,45.,46.},DIRICH} );

    // assigning values to the <=12 facet integration points
    for ( auto& e : domain.CellVector() ) {
          const auto facets(e->Facets());
          const uint32_t fip{0u};
          // looping over the FV facets
          for ( uint32_t i{0u}; i<facets; ++i ) {
               e->Store( i, fip, fipv1, ip_vector_variables[i] );
               e->Store( i, fip, fipt1, ip_tensor_variables[i] );
               e->Store( i, fip, fipa1, ip_array_variables[i] );
            }
      }
    // testing that they are read correctly
    for ( const auto& e : domain.CellVector() ) {
          const auto facets(e->Facets());
          const uint32_t fip{0u};
          // looping over the FV facets
          for ( uint32_t i{0u}; i<facets; ++i ) {
               auto vc_vec = e->ReadVector( i, fip, fipv1 );
               _test( vc_vec == ip_vector_variables[i] );
               auto ts_vec = e->ReadTensor( i, fip, fipt1 );
               _test( ts_vec == ip_tensor_variables[i] );
               auto ar_vec = e->ReadArray( i, fip, fipa1 );
               for ( uint32_t j{0u}; j<ar_vec.size(); ++j )
                 _test( ar_vec[j] == ip_array_variables[i][j] );
            }
       }

 } // end runExtensionsForDeclarativeProgramming




void LocalVariableStorage_Test::runExtensionsForDeclarativeProgramming_INDEX()
 {
    VSet<3U> vset;
    create_Pyramid_Hexa_VSet( vset, true ); // skewed model
    Model<3U>  model( vset, "LocalVariableStorage_Test-variables.txt" );
    Region<3U> domain = model.Region("Model");
    
    // testing read operations on element, node, element integration points and FV integration points
    // element
    const csmp::INDEX<VECTOR,ELEMENT> ev1(model.Database().StorageKey("element vector 1"));
    const csmp::INDEX<TENSOR,ELEMENT> et1(model.Database().StorageKey("element tensor 1"));
    const csmp::INDEX<ARRAY,ELEMENT>  ea1(model.Database().StorageKey("element array 1"));
    // node
    const csmp::INDEX<VECTOR,NODE> nv1(model.Database().StorageKey("nodal vector 1"));
    const csmp::INDEX<TENSOR,NODE> nt1(model.Database().StorageKey("nodal tensor 1"));
    const csmp::INDEX<ARRAY,NODE>  na1(model.Database().StorageKey("nodal array 1"));

    // inputs
    VectorVariable<3U> vc(ANY,PLAIN,ROBIN,1.,2.,3.);
    TensorVariable<3U> ts(ANY,PLAIN,ROBIN,1.,2.,3.,4.,5.,6.,7.,8.,9);
    ArrayVariable ary( 4, 0., ANY ), ary_test;
    for ( uint32_t i{0u}; i<ary.Size(); ++i ) ary(i) = static_cast<double>(i);

    // element & node
    model.InputPropertyValue("element vector 1", makeVector(ANY,PLAIN,ROBIN,1.,2.,3.) );
    model.InputPropertyValue("element tensor 1", makeTensor(ANY,PLAIN,ROBIN,1.,2.,3.,4.,5.,6.,7.,8.,9) );
    model.InputPropertyValue("element array 1", ary );

    model.InputPropertyValue("nodal vector 1", makeVector(ANY,PLAIN,ROBIN,1.,2.,3.) );
    model.InputPropertyValue("nodal tensor 1", makeTensor(ANY,PLAIN,ROBIN,1.,2.,3.,4.,5.,6.,7.,8.,9) );
    model.InputPropertyValue("nodal array 1", ary );
    
    // 1. test of functions reading element and node variables
    for ( const auto& e : domain.CellVector() ) {
         _test( vc == e->ReadVector(ev1) );
         _test( ts == e->ReadTensor(et1) );
         e->Read( ea1, ary_test );
         _test( ary == ary_test );
         vector<double> a1_vec1 = e->ReadArray(ea1), a1_vec2(4);
         // testing return value initialisation
         for ( uint32_t i{0u}; i<ary.Size(); ++i )
              _test( a1_vec1[i] == ary[i] );
         // testing ReadArrayEntry()
         for ( uint32_t i{0u}; i<a1_vec1.size(); ++i )
           a1_vec2[i] = e->ReadArrayEntry( ea1, i );
         // checking consistency 
         _test( a1_vec1 == a1_vec2 );
      }
    for ( const auto& n : domain.NodeVector() ) {
         _test( vc == n->ReadVector(nv1) );
         _test( ts == n->ReadTensor(nt1) );
         n->Read( na1, ary_test );
         _test( ary == ary_test );
         vector<double> a1_vec1 = n->ReadArray(na1), a1_vec2(4);
         for ( uint32_t i{0u}; i<ary.Size(); ++i )
           a1_vec2[i] = n->ReadArrayEntry( na1, i );
         _test( a1_vec1 == a1_vec2 );
         for ( uint32_t i{0u}; i<ary.Size(); ++i )
              _test( ary[i] == a1_vec1[i] ); 
      }

    // 2. test of functions reading variables from element integration points
    // ----------------------------------------------------------------------
    // initialising different array variables for up to 8 integration points (hexahedron)
    vector<VectorVariable<3U>> ip_vector_variables  = { {ANY,ANY,ANY,0.,1.,2.}, {PLAIN,PLAIN,PLAIN,3.,4.,5.},
                                                        {ROBIN,ROBIN,ROBIN,6.,7.,8.}, {DIRICH,DIRICH,DIRICH,9.,10.,11.},
                                                        {INIT_COND,INIT_COND,INIT_COND,12.,13.,14.}, {INIT_GUESS,INIT_GUESS,INIT_GUESS,15.,16.,17.},
                                                        {PERIODIC,PERIODIC,PERIODIC,18.,19.,20.}, {FIELD_DATA,FIELD_DATA,FIELD_DATA,21.,22.,23.} };

    vector<TensorVariable<3U>> ip_tensor_variables  = { {ANY,ANY,ANY,0.,1.,2.,3.,4.,5.,6.,7.,8.},
                                                        {PLAIN,PLAIN,PLAIN,9.,10.,11.,12.,13.,14.,15.,16.,17.},
                                                        {ROBIN,ROBIN,ROBIN,18.,19.,20.,21.,22.,23.,24.,25.,26.},
                                                        {DIRICH,DIRICH,DIRICH,27.,28.,29.,30.,31.,32.,33.,34.,35.},
                                                        {INIT_COND,INIT_COND,INIT_COND,36.,37.,38.,39.,40.,41.,42.,43.,44.},
                                                        {INIT_GUESS,INIT_GUESS,INIT_GUESS,45.,46.,47.,48.,49.,50.,51.,52.,53.},
                                                        {PERIODIC,PERIODIC,PERIODIC,54.,55.,56.,57.,58.,59.,60.,61.,62.},
                                                        {FIELD_DATA,FIELD_DATA,FIELD_DATA,63.,64.,65.,66.,67.,68.,69.,70.,71.} };

    vector<ArrayVariable> ip_array_variables  = { {{0.,1.,2.,3.},ANY}, {{4.,5.,6.,7.},PLAIN},
                                                  {{8.,9.,10.,11.},ROBIN}, {{12.,13.,14.,15.},DIRICH},
                                                  {{16.,17.,18.,19.},INIT_COND}, {{20.,21.,22.,23.},INIT_GUESS},
                                                  {{24.,25.,26.,27.},PERIODIC}, {{28.,29.,30.,31.},FIELD_DATA} };

    // element integration points
    // --------------------------
    const csmp::INDEX<VECTOR,ELEMENT_INTEGRATION_POINT> eiv1(model.Database().StorageKey("eip vector 1"));
    const csmp::INDEX<TENSOR,ELEMENT_INTEGRATION_POINT> eit1(model.Database().StorageKey("eip tensor 1"));
    const csmp::INDEX<ARRAY,ELEMENT_INTEGRATION_POINT>  eia1(model.Database().StorageKey("eip array 1"));

    // assigning the values
    for ( auto& e : domain.CellVector() )
      for ( uint32_t i{0u}; i<e->IntegrationPoints(); ++i ) {
           e->Store( i, eiv1, ip_vector_variables[i] );
           e->Store( i, eit1, ip_tensor_variables[i] );
           e->Store( i, eia1, ip_array_variables[i] );
        }
    // testing that they are read correctly
    for ( const auto& e : domain.CellVector() )
      for ( uint32_t i{0u}; i<e->IntegrationPoints(); ++i ) {
           auto vc_vec = e->ReadVector( i, eiv1 );
           _test( vc_vec == ip_vector_variables[i] );
           auto ts_vec = e->ReadTensor( i, eit1 );
           _test( ts_vec == ip_tensor_variables[i] );
           auto ar_vec = e->ReadArray( i, eia1 );
           for ( uint32_t j{0u}; j<ar_vec.size(); ++j )
             _test( ar_vec[j] == ip_array_variables[i][j] );
        }

    // FV facet integration points
    // ---------------------------
    // facet integration points
    const csmp::INDEX<VECTOR,FACET_INTEGRATION_POINT> fipv1(model.Database().StorageKey("fip vector 1"));
    const csmp::INDEX<TENSOR,FACET_INTEGRATION_POINT> fipt1(model.Database().StorageKey("fip tensor 1"));
    const csmp::INDEX<ARRAY,FACET_INTEGRATION_POINT>  fipa1(model.Database().StorageKey("fip array 1"));
    
    // adding the necessary 4 extra values to the variable arrays
    // vectors
    ip_vector_variables.push_back( {ROBIN,PLAIN,ANY,24.,25.,26} );
    ip_vector_variables.push_back( {ROBIN,PLAIN,ANY,27.,28.,29} );
    ip_vector_variables.push_back( {ROBIN,PLAIN,ANY,30.,31.,32.} );
    ip_vector_variables.push_back( {ROBIN,PLAIN,ANY,33.,34.,35.} );
    // tensors
    ip_tensor_variables.push_back( {ANY,ANY,ANY,72.,73.,74.,75.,76.,77.,78.,79.,80.} );
    ip_tensor_variables.push_back( {PLAIN,PLAIN,PLAIN,81.,82.,83.,84.,85.,86.,87.,88.,89.} );
    ip_tensor_variables.push_back( {ROBIN,ROBIN,ROBIN,90.,91.,92.,93.,94.,95.,96.,97.,98.} );
    ip_tensor_variables.push_back( {DIRICH,DIRICH,DIRICH,99.,100.,101.,102.,103.,104.,105.,106.,107.} );
    // arrays
    ip_array_variables.push_back( {{32.,33.,34.,34.},ANY} );
    ip_array_variables.push_back( {{35.,36.,37.,38.},PLAIN} );
    ip_array_variables.push_back( {{39.,40.,41.,42.},ROBIN} );
    ip_array_variables.push_back( {{43.,44.,45.,46.},DIRICH} );

    // assigning values to the <=12 facet integration points
    for ( auto& e : domain.CellVector() ) {
          const auto facets(e->Facets());
          const uint32_t fip{0u};
          // looping over the FV facets
          for ( uint32_t i{0u}; i<facets; ++i ) {
               e->Store( i, fip, fipv1, ip_vector_variables[i] );
               e->Store( i, fip, fipt1, ip_tensor_variables[i] );
               e->Store( i, fip, fipa1, ip_array_variables[i] );
            }
      }
    // testing that they are read correctly
    for ( const auto& e : domain.CellVector() ) {
          const auto facets(e->Facets());
          const uint32_t fip{0u};
          // looping over the FV facets
          for ( uint32_t i{0u}; i<facets; ++i ) {
               auto vc_vec = e->ReadVector( i, fip, fipv1 );
               _test( vc_vec == ip_vector_variables[i] );
               auto ts_vec = e->ReadTensor( i, fip, fipt1 );
               _test( ts_vec == ip_tensor_variables[i] );
               auto ar_vec = e->ReadArray( i, fip, fipa1 );
               for ( uint32_t j{0u}; j<ar_vec.size(); ++j )
                 _test( ar_vec[j] == ip_array_variables[i][j] );
            }
       }

 } // end runExtensionsForDeclarativeProgramming
  
  


void LocalVariableStorage_Test::run()
  {
    runTest<1>();
    runTest<2>();
    runTest<3>();
    run3D();
    runExtensionsForDeclarativeProgramming();
    run3D_with_templatized_INDEX();
    runExtensionsForDeclarativeProgramming_INDEX();

    //cout << endl << "Size of LocalVariableStorage: " << sizeof(LocalVariableStorage<3,Element<3> >) << endl;
  }

} // csmp
