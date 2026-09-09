#include "CopyReplaceVisitor_Test.h"

#include "CopyReplaceVisitor.h"
#include "PropertyHandle.h"
#include "CSMP_definitions.h"
#include "Element.h"
#include "Model.h"
#include "Region.h"
#include "vsetMakers.h"

namespace csmp {

  
void CopyReplaceVisitor_Test::run()
{
    VSet<2U> vset2D;
    ModelTopology topo = create_MeshPatchWithLineElements_VSet( vset2D );
    vset2D.RemoveData("element variable"); // not needed here
    bool treat_domains_as_regions{true};
    Model<2U> model2D( topo, vset2D, "CSMP-1phase-variables.txt", treat_domains_as_regions );

    // Test for 3D only, creating auxilliaries
    enum{DIM=2U};
    ScalarVariable scalarValue( ROBIN, 2. );
    VectorVariable<DIM> vectorValue( ROBIN, 2. );
    TensorVariable<DIM> tensorValue( ROBIN, 2. );

    // nodal scalar
    PropertyHandle<DIM> nodeScalar1( model2D, "node scalar 1", SCALAR, NODE );
    PropertyHandle<DIM> nodeScalar2( model2D, "node scalar 2", SCALAR, NODE );
    model2D.InputPropertyValue( "node scalar 1", scalarValue );
    CopyReplaceVisitor<ScalarVariable,DIM> cpyNodeScalar1to2( model2D.Database(), "node scalar 1", "node scalar 2");
    model2D.Accept(cpyNodeScalar1to2);
    TestNodes( model2D, scalarValue, "node scalar 2" );

    // nodal vector
    PropertyHandle<DIM> nodeVector1( model2D, "node vector 1", VECTOR, NODE );
    PropertyHandle<DIM> nodeVector2( model2D, "node vector 2", VECTOR, NODE );
    model2D.InputPropertyValue( "node vector 1", vectorValue );
    CopyReplaceVisitor<VectorVariable<DIM>,DIM> cpyNodeVector1to2( model2D.Database(), "node vector 1", "node vector 2" );
    model2D.Accept(cpyNodeVector1to2);
    TestNodes( model2D, vectorValue, "node vector 2" );

    // nodal tensor
    PropertyHandle<DIM> nodeTensor1( model2D, "node tensor 1", TENSOR, NODE );
    PropertyHandle<DIM> nodeTensor2( model2D, "node tensor 2", TENSOR, NODE );
    model2D.InputPropertyValue( "node tensor 1", tensorValue );
    CopyReplaceVisitor<TensorVariable<DIM>,DIM> cpyNodeTensor1to2( model2D.Database(), "node tensor 1", "node tensor 2" );
    model2D.Accept(cpyNodeTensor1to2);
    TestNodes( model2D, tensorValue, "node tensor 2" );

    // element scalar
    PropertyHandle<DIM> elementScalar1( model2D, "element scalar 1", SCALAR, ELEMENT );
    PropertyHandle<DIM> elementScalar2( model2D, "element scalar 2", SCALAR, ELEMENT );
    model2D.InputPropertyValue( "element scalar 1", scalarValue );
    CopyReplaceVisitor<ScalarVariable,DIM> cpyElementScalar1to2( model2D.Database(), "element scalar 1", "element scalar 2" );
    model2D.Accept(cpyElementScalar1to2);
    TestElements( model2D, scalarValue, "element scalar 2" );

    // element vector
    PropertyHandle<DIM> elementVector1( model2D, "element vector 1", VECTOR, ELEMENT );
    PropertyHandle<DIM> elementVector2( model2D, "element vector 2", VECTOR, ELEMENT );
    model2D.InputPropertyValue( "element vector 1", vectorValue );
    CopyReplaceVisitor<VectorVariable<DIM>,DIM> cpyElementVector1to2( model2D.Database(), "element vector 1", "element vector 2" );
    model2D.Accept(cpyElementVector1to2);
    TestElements( model2D, vectorValue, "element vector 2" );

    // element tensor
    PropertyHandle<DIM> elementTensor1( model2D, "element tensor 1", TENSOR, ELEMENT );
    PropertyHandle<DIM> elementTensor2( model2D, "element tensor 2", TENSOR, ELEMENT );
    model2D.InputPropertyValue( "element tensor 1", tensorValue );
    CopyReplaceVisitor<TensorVariable<DIM>,DIM> cpyElementTensor1to2( model2D.Database(), "element tensor 1", "element tensor 2" );
    model2D.Accept(cpyElementTensor1to2);
    TestElements( model2D, tensorValue, "element tensor 2" );

} // run




template<class V>
void CopyReplaceVisitor_Test::TestNodes( Model<2>& model, const V& value, const char* propertyName )
{
  V cache;
  Index key( model.Database().StorageKey(propertyName) );
  Region<2U>& domain{ model.Region("Model") };
  
  for( auto it( domain.NodesBegin() ); it != domain.NodesEnd(); ++it )
  {
    (*it)->Read( key, cache );
    _test( cache == value );
  }  
}


template<class V>
void CopyReplaceVisitor_Test::TestElements( Model<2>& model, const V& value, const char* propertyName )
{
  V cache;
  Index key( model.Database().StorageKey(propertyName) );
  Region<2U>& domain{ model.Region("Model") };

  for( auto it( domain.CellsBegin() ); it != domain.CellsEnd(); ++it )
  {
    (*it)->Read( key, cache );
    _test( cache == value );
  }  
}

template<class V>
void CopyReplaceVisitor_Test::TestElementIntegrationPoints( Model<2>& model, const V& value, const char* propertyName )
{
  V cache;
  Index key( model.Database().StorageKey(propertyName) );
  Region<2U>& domain{ model.Region("Model") };

  for( auto it( domain.CellsBegin() ); it != domain.CellsEnd(); ++it )
    for( uint32_t ip{0U}; ip < (*it)->IntegrationPoints(); ++ip )
      {
        (*it)->Read( key, cache );
        _test( cache == value );
      }  
}


} // csmp
