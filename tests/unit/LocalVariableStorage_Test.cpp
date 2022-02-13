#include "LocalVariableStorage_Test.h"

#include "LocalVariableStorage.h"
#include "Element.h"
#include "LinearTriangle.h"
#include "LinearTetrahedron.h"

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
    for( size_t d(0); d < flaggedArrayLength; ++d )
        _test( e.Status(kfa1,d) == fa1.Flag(d) );
    for( size_t d(0); d < dim; ++d )
      {
        _test( e.Status( kv1, d ) == v1.Flag(d) );
        _test( e.Status( kt1, d ) == t1.Flag(d) );
      }

    e.Status(ks1, INIT_GUESS );
    e.Status(ka1, INIT_COND);
    for( size_t d(0); d < flaggedArrayLength; ++d )
        e.Status(kfa1,d, FIELD_DATA);
    for( size_t d(0); d < dim; ++d )
      {
      e.Status( kv1, d, CONSTANT_FLUX );
      e.Status( kt1, d, PERIODIC );
      }

    _test( e.Status(ks1) == INIT_GUESS );
    _test( e.Status(ka1) == INIT_COND );
    for( size_t d(0); d < flaggedArrayLength; ++d )
        _test( e.Status(kfa1,d) == FIELD_DATA );
    for( size_t d(0); d < dim; ++d )
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
    for( size_t d(0); d < flaggedArrayLength; ++d )
        _test( e.Status(kfa1,d) == fa1.Flag(d) );
    for( size_t d(0); d < dim; ++d )
      {
        _test( e.Status( kv1, d ) == v1.Flag(d) );
        _test( e.Status( kt1, d ) == t1.Flag(d) );
      }

    e.Status(ks1, INIT_GUESS );
    e.Status(ka1, INIT_COND );
    for( size_t d(0); d < flaggedArrayLength; ++d )
        e.Status(kfa1,d, FIELD_DATA);
    for( size_t d(0); d < dim; ++d )
      {
      e.Status( kv1, d, CONSTANT_FLUX );
      e.Status( kt1, d, PERIODIC );
      }

    _test( e.Status(ks1) == INIT_GUESS );
    _test( e.Status(ka1) == INIT_COND );
    for( size_t d(0); d < flaggedArrayLength; ++d )
        _test( e.Status(kfa1,d) == FIELD_DATA );
    for( size_t d(0); d < dim; ++d )
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
    for( size_t d(0); d < flaggedArrayLength; ++d )
        _test( e.Status(kfa1,d) == fa1.Flag(d) );
    for( size_t d(0); d < dim; ++d )
      {
        _test( e.Status( kv1, d ) == v1.Flag(d) );
        _test( e.Status( kt1, d ) == t1.Flag(d) );
      }

    e.Status(ks1, INIT_GUESS );
    e.Status(ka1, INIT_COND );
    for( size_t d(0); d < flaggedArrayLength; ++d )
        e.Status(kfa1,d, FIELD_DATA);
    for( size_t d(0); d < dim; ++d )
      {
      e.Status( kv1, d, CONSTANT_FLUX );
      e.Status( kt1, d, PERIODIC );
      }

    _test( e.Status(ks1) == INIT_GUESS );
    _test( e.Status(ka1) == INIT_COND );
    for( size_t d(0); d < flaggedArrayLength; ++d )
        _test( e.Status(kfa1,d) == FIELD_DATA );
    for( size_t d(0); d < dim; ++d )
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
  
  
  
  
  
  
  


void LocalVariableStorage_Test::run()
  {
    runTest<1>();
    runTest<2>();
    runTest<3>();
    run3D();
    run3D_with_templatized_INDEX();

    //cout << endl << "Size of LocalVariableStorage: " << sizeof(LocalVariableStorage<3,Element<3> >) << endl;
  }

} // csmp
