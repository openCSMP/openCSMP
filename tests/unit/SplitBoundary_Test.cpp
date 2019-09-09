#include "SplitBoundary_Test.h"
#include "Boundary.h"
#include "VTU_Interface.h"
#include "ANSYS_Model2D.h"
#include "ANSYS_Model3D.h"
#include "EclipseModel.h"
#include "VTK_Interface.h"

using namespace std;

namespace csmp
{

template<size_t dim>
void shiftSplitBoundary( SplitBoundary<dim>& splitboundary, INTERFACE_SIDE side, double64 xShift, double64 yShift, double64 zShift )
{
  const typename vector<InterFace<dim>*>::const_iterator facesEnd( splitboundary.ElementsEnd() );
  if ( dim == 2 )
    for ( typename vector<InterFace<dim>*>::const_iterator it = splitboundary.ElementsBegin(); it != facesEnd; ++it )
    {
      const size_t nodes( (*it)->Nodes() );
      for ( size_t i = 0; i < nodes; ++i )
      {
        if ( (*it)->N( i, INSIDE )->Idx() != (*it)->N( i, OUTSIDE )->Idx() )
        {
          (*it)->N( i, side )->x( (*it)->N( i, side )->x() + xShift );
          (*it)->N( i, side )->y( (*it)->N( i, side )->y() + yShift );
        }
      }
    }
  else if ( dim == 3 )
    for ( typename vector<InterFace<dim>*>::const_iterator it = splitboundary.ElementsBegin(); it != facesEnd; ++it )
    {
      const size_t nodes( (*it)->Nodes() );
      for ( size_t i = 0; i < nodes; ++i )
      {
        if ( (*it)->N( i, INSIDE )->Idx() != (*it)->N( i, OUTSIDE )->Idx() )
        {
          (*it)->N( i, side )->x( (*it)->N( i, side )->x() + xShift );
          (*it)->N( i, side )->y( (*it)->N( i, side )->y() + yShift );
          (*it)->N( i, side )->z( (*it)->N( i, side )->z() + zShift );
        }
      }
    }
}


template<size_t dim>
void shiftSplitBoundary( SplitBoundary<dim>& splitboundary, double64 shift )
{
  VectorVariable<dim> displacementPerpedicularToInterface( ANY, 0.0 );
  Point<dim> displacementToBaryCenter;
  Point<dim> baryCenter;
  size_t null_neighbors( 0 );
  Element<dim>* parentElement( NULL );

  const typename vector<InterFace<dim>*>::const_iterator facesEnd( splitboundary.ElementsEnd() );
  if ( dim == 2 )
    for ( typename vector<InterFace<dim>*>::const_iterator it = splitboundary.ElementsBegin(); it != facesEnd; ++it )
    {
      if ( (*it)->BaseElement() != NULL )
        (*it)->BaseElement()->UnitNormal( displacementPerpedicularToInterface );
      else
        (*it)->UnitNormal( displacementPerpedicularToInterface, INSIDE );
      displacementPerpedicularToInterface *= shift;
      const size_t nodes( (*it)->Nodes() );
      for ( size_t i = 0; i < nodes; ++i )
      {
        if ( (*it)->N( i, INSIDE ) != (*it)->N( i, OUTSIDE ) )
        {
          null_neighbors = 0;
          for ( size_t k = 0; k<(*it)->N( i, OUTSIDE )->Parents(); k++ )
          {
            parentElement = (*it)->N( i, OUTSIDE )->Parent( k );
            for ( size_t j = 0; j<parentElement->Neighbors(); j++ )
            {
              if ( parentElement->Neighbor( j ) == NULL )
                null_neighbors++;
            }
          }
          if ( null_neighbors > 1 )
          {
            baryCenter = (*it)->Parent( OUTSIDE )->BaryCenter();
            displacementToBaryCenter = baryCenter - (*it)->N( i, OUTSIDE )->Coordinate();
            displacementToBaryCenter.NormalizeLengthTo( 1.0 );
            displacementToBaryCenter *= shift;
            (*it)->N( i, OUTSIDE )->x( (*it)->N( i, OUTSIDE )->x() + displacementToBaryCenter[0] );
            (*it)->N( i, OUTSIDE )->y( (*it)->N( i, OUTSIDE )->y() + displacementToBaryCenter[1] );
          }
          (*it)->N( i, OUTSIDE )->x( (*it)->N( i, OUTSIDE )->x() + displacementPerpedicularToInterface[0] );
          (*it)->N( i, OUTSIDE )->y( (*it)->N( i, OUTSIDE )->y() + displacementPerpedicularToInterface[1] );

          null_neighbors = 0;
          for ( size_t k = 0; k<(*it)->N( i, INSIDE )->Parents(); k++ )
          {
            parentElement = (*it)->N( i, INSIDE )->Parent( k );
            for ( size_t j = 0; j<parentElement->Neighbors(); j++ )
            {
              if ( parentElement->Neighbor( j ) == NULL )
                null_neighbors++;
            }
          }
          if ( null_neighbors > 1 )
          {
            baryCenter = (*it)->Parent( INSIDE )->BaryCenter();
            displacementToBaryCenter = baryCenter - (*it)->N( i, INSIDE )->Coordinate();
            displacementToBaryCenter.NormalizeLengthTo( 1.0 );
            displacementToBaryCenter *= shift;
            (*it)->N( i, INSIDE )->x( (*it)->N( i, INSIDE )->x() + displacementToBaryCenter[0] );
            (*it)->N( i, INSIDE )->y( (*it)->N( i, INSIDE )->y() + displacementToBaryCenter[1] );
          }
          (*it)->N( i, INSIDE )->x( (*it)->N( i, INSIDE )->x() - displacementPerpedicularToInterface[0] );
          (*it)->N( i, INSIDE )->y( (*it)->N( i, INSIDE )->y() - displacementPerpedicularToInterface[1] );
        }
      }
    }
  else if ( dim == 3 )
    for ( typename vector<InterFace<dim>*>::const_iterator it = splitboundary.ElementsBegin(); it != facesEnd; ++it )
    {
      const size_t nodes( (*it)->Nodes() );
      for ( size_t i = 0; i < nodes; ++i )
      {
        if ( (*it)->BaseElement() != NULL )
          (*it)->BaseElement()->UnitNormal( displacementPerpedicularToInterface );
        else
          (*it)->UnitNormal( displacementPerpedicularToInterface, INSIDE );
        displacementPerpedicularToInterface *= shift;
        if ( (*it)->N( i, INSIDE ) != (*it)->N( i, OUTSIDE ) )
        {
          null_neighbors = 0;
          for ( size_t k = 0; k<(*it)->N( i, OUTSIDE )->Parents(); k++ )
          {
            parentElement = (*it)->N( i, OUTSIDE )->Parent( k );
            for ( size_t j = 0; j<parentElement->Neighbors(); j++ )
            {
              if ( parentElement->Neighbor( j ) == NULL )
                null_neighbors++;
            }
          }
          if ( null_neighbors > 1 )
          {
            baryCenter = (*it)->Parent( OUTSIDE )->BaryCenter();
            displacementToBaryCenter = baryCenter - (*it)->N( i, OUTSIDE )->Coordinate();
            displacementToBaryCenter.NormalizeLengthTo( 1.0 );
            displacementToBaryCenter *= shift;
            (*it)->N( i, OUTSIDE )->x( (*it)->N( i, OUTSIDE )->x() + displacementToBaryCenter[0] );
            (*it)->N( i, OUTSIDE )->y( (*it)->N( i, OUTSIDE )->y() + displacementToBaryCenter[1] );
            (*it)->N( i, OUTSIDE )->z( (*it)->N( i, OUTSIDE )->z() + displacementToBaryCenter[2] );
          }
          (*it)->N( i, OUTSIDE )->x( (*it)->N( i, OUTSIDE )->x() + displacementPerpedicularToInterface[0] );
          (*it)->N( i, OUTSIDE )->y( (*it)->N( i, OUTSIDE )->y() + displacementPerpedicularToInterface[1] );
          (*it)->N( i, OUTSIDE )->z( (*it)->N( i, OUTSIDE )->z() + displacementPerpedicularToInterface[2] );

          null_neighbors = 0;
          for ( size_t k = 0; k<(*it)->N( i, INSIDE )->Parents(); k++ )
          {
            parentElement = (*it)->N( i, INSIDE )->Parent( k );
            for ( size_t j = 0; j<parentElement->Neighbors(); j++ )
            {
              if ( parentElement->Neighbor( j ) == NULL )
                null_neighbors++;
            }
          }
          if ( null_neighbors > 1 )
          {
            baryCenter = (*it)->Parent( INSIDE )->BaryCenter();
            displacementToBaryCenter = baryCenter - (*it)->N( i, INSIDE )->Coordinate();
            displacementToBaryCenter.NormalizeLengthTo( 1.0 );
            displacementToBaryCenter *= shift;
            (*it)->N( i, INSIDE )->x( (*it)->N( i, INSIDE )->x() + displacementToBaryCenter[0] );
            (*it)->N( i, INSIDE )->y( (*it)->N( i, INSIDE )->y() + displacementToBaryCenter[1] );
            (*it)->N( i, INSIDE )->z( (*it)->N( i, INSIDE )->z() + displacementToBaryCenter[2] );
          }
          (*it)->N( i, INSIDE )->x( (*it)->N( i, INSIDE )->x() - displacementPerpedicularToInterface[0] );
          (*it)->N( i, INSIDE )->y( (*it)->N( i, INSIDE )->y() - displacementPerpedicularToInterface[1] );
          (*it)->N( i, INSIDE )->z( (*it)->N( i, INSIDE )->z() - displacementPerpedicularToInterface[2] );
        }
      }
    }
}


template<size_t dim>
void shiftInterfaceTips( Region<dim>& region, double64 shift )
{
  Point<dim> displacement;
  Point<dim> baryCenter;

  if ( dim == 2 )
    for ( typename vector<Element<dim>*>::const_iterator it = region.PerimeterElementsBegin(); it != region.ElementsEnd(); ++it )
    {
      baryCenter = (*it)->BaryCenter();
      const size_t nodes( (*it)->Nodes() );
      for ( size_t i = 0; i < nodes; ++i )
      {
        displacement = baryCenter - (*it)->N( i )->Coordinate();
        displacement.NormalizeLengthTo( 1.0 );
        displacement *= shift;
        (*it)->N( i )->x( (*it)->N( i )->x() + displacement[0] );
        (*it)->N( i )->y( (*it)->N( i )->y() + displacement[1] );
      }
    }
  else if ( dim == 3 )
    for ( typename vector<Element<dim>*>::const_iterator it = region.PerimeterElementsBegin(); it != region.ElementsEnd(); ++it )
    {
      baryCenter = (*it)->BaryCenter();
      const size_t nodes( (*it)->Nodes() );
      for ( size_t i = 0; i < nodes; ++i )
      {
        displacement = baryCenter - (*it)->N( i )->Coordinate();
        displacement.NormalizeLengthTo( 1.0 );
        displacement *= shift;
        (*it)->N( i )->x( (*it)->N( i )->x() + displacement[0] );
        (*it)->N( i )->y( (*it)->N( i )->y() + displacement[1] );
        (*it)->N( i )->z( (*it)->N( i )->z() + displacement[2] );
      }
    }
}


template<size_t dim>
void shiftRegion( Region<dim>& region, double64 xShift, double64 yShift, double64 zShift )
{
  if ( dim == 2 )
    for ( typename vector<Node<dim>*>::const_iterator it = region.NodesBegin(); it != region.NodesEnd(); ++it )
    {
      (*it)->x( (*it)->x() + xShift );
      (*it)->y( (*it)->y() + yShift );
    }
  else if ( dim == 3 )
    for ( typename vector<Node<dim>*>::const_iterator it = region.NodesBegin(); it != region.NodesEnd(); ++it )
    {
      (*it)->x( (*it)->x() + xShift );
      (*it)->y( (*it)->y() + yShift );
      (*it)->z( (*it)->z() + zShift );
    }
}


template<size_t dim>
void shiftRegionAboveLine( Region<dim>& region, size_t x_or_y_or_z, double64 line_coordinate, double64 shift, double64 eps )
{
  if ( dim == 2 )
  {
    if ( x_or_y_or_z == 0 )
    {
      for ( typename vector<Node<dim>*>::const_iterator it = region.NodesBegin(); it != region.NodesEnd(); ++it )
        if ( (*it)->x() > (line_coordinate + eps) )
          (*it)->x( (*it)->x() + shift );
    }
    else if ( x_or_y_or_z == 1 )
    {
      for ( typename vector<Node<dim>*>::const_iterator it = region.NodesBegin(); it != region.NodesEnd(); ++it )
        if ( (*it)->y() > (line_coordinate + eps) )
          (*it)->y( (*it)->y() + shift );
    }
  }
  else if ( dim == 3 )
  {
    if ( x_or_y_or_z == 0 )
    {
      for ( typename vector<Node<dim>*>::const_iterator it = region.NodesBegin(); it != region.NodesEnd(); ++it )
        if ( (*it)->x() > (line_coordinate + eps) )
          (*it)->x( (*it)->x() + shift );
    }
    else if ( x_or_y_or_z == 1 )
    {
      for ( typename vector<Node<dim>*>::const_iterator it = region.NodesBegin(); it != region.NodesEnd(); ++it )
        if ( (*it)->y() > (line_coordinate + eps) )
          (*it)->y( (*it)->y() + shift );
    }
    else {
      for ( typename vector<Node<dim>*>::const_iterator it = region.NodesBegin(); it != region.NodesEnd(); ++it )
        if ( (*it)->z() > (line_coordinate + eps) )
          (*it)->z( (*it)->z() + shift );
    }
  }
}


template<size_t dim>
void shiftRegionBelowLine( Region<dim>& region, size_t x_or_y_or_z, double64 line_coordinate, double64 shift, double64 eps )
{
  if ( dim == 2 )
  {
    if ( x_or_y_or_z == 0 )
    {
      for ( typename vector<Node<dim>*>::const_iterator it = region.NodesBegin(); it != region.NodesEnd(); ++it )
        if ( (*it)->x() < (line_coordinate + eps) )
          (*it)->x( (*it)->x() + shift );
    }
    else if ( x_or_y_or_z == 1 )
    {
      for ( typename vector<Node<dim>*>::const_iterator it = region.NodesBegin(); it != region.NodesEnd(); ++it )
        if ( (*it)->y() < (line_coordinate + eps) )
          (*it)->y( (*it)->y() + shift );
    }
  }
  else if ( dim == 3 )
  {
    if ( x_or_y_or_z == 0 )
    {
      for ( typename vector<Node<dim>*>::const_iterator it = region.NodesBegin(); it != region.NodesEnd(); ++it )
        if ( (*it)->x() < (line_coordinate + eps) )
          (*it)->x( (*it)->x() + shift );
    }
    else if ( x_or_y_or_z == 1 )
    {
      for ( typename vector<Node<dim>*>::const_iterator it = region.NodesBegin(); it != region.NodesEnd(); ++it )
        if ( (*it)->y() < (line_coordinate + eps) )
          (*it)->y( (*it)->y() + shift );
    }
    else {
      for ( typename vector<Node<dim>*>::const_iterator it = region.NodesBegin(); it != region.NodesEnd(); ++it )
        if ( (*it)->z() < (line_coordinate + eps) )
          (*it)->z( (*it)->z() + shift );
    }
  }
}


template<size_t dim>
void scaleRegionSymmetricOverZero( Region<dim>& region, double64 xScale, double64 yScale, double64 zScale )
{
  if ( dim == 2 )
    for ( typename vector<Node<dim>*>::const_iterator it = region.NodesBegin(); it != region.NodesEnd(); ++it )
    {
      (*it)->x( (*it)->x()*xScale );
      (*it)->y( (*it)->y()*yScale );
    }
  else if ( dim == 3 )
    for ( typename vector<Node<dim>*>::const_iterator it = region.NodesBegin(); it != region.NodesEnd(); ++it )
    {
      (*it)->x( (*it)->x()*xScale );
      (*it)->y( (*it)->y()*yScale );
      (*it)->z( (*it)->z()*zScale );
    }
}


template<size_t dim>
void scaleRegion( Region<dim>& region, double64 xScale, double64 yScale, double64 zScale )
{
  Point<dim> min_point;
  Point<dim> max_point;
  Point<dim> mid_point;

  region.MinMaxCoordinates( min_point, max_point );
  mid_point = (min_point + max_point) / 2.0;

  if ( dim == 2 )
    for ( typename vector<Node<dim>*>::const_iterator it = region.NodesBegin(); it != region.NodesEnd(); ++it )
    {
      (*it)->x( mid_point[0] + ((*it)->x() - mid_point[0])*xScale );
      (*it)->y( mid_point[1] + ((*it)->y() - mid_point[1])*yScale );
    }
  else if ( dim == 3 )
    for ( typename vector<Node<dim>*>::const_iterator it = region.NodesBegin(); it != region.NodesEnd(); ++it )
    {
      (*it)->x( mid_point[0] + ((*it)->x() - mid_point[0])*xScale );
      (*it)->y( mid_point[1] + ((*it)->y() - mid_point[1])*yScale );
      (*it)->z( mid_point[2] + ((*it)->z() - mid_point[2])*zScale );
    }
}


template<size_t dim>
void SplitBoundary_Test::NodeParents( const Region<dim>& region, size_t minParentCount /* = 2 */ )
{
  const typename vector<Node<dim>*>::const_iterator nodesEnd( region.NodesEnd() );
  for ( typename vector<Node<dim>*>::const_iterator node( region.NodesBegin() ); node != nodesEnd; ++node )
    _test( (*node)->Parents() >= minParentCount );
}


template<size_t dim>
void SplitBoundary_Test::ElementNodes( const Region<dim>& region )
{
  const typename vector<Element<dim>*>::const_iterator elementsEnd( region.ElementsEnd() );
  for ( typename vector<Element<dim>*>::const_iterator element( region.ElementsBegin() ); element != elementsEnd; ++element )
    _test( (*element)->NodeConnectorSize() == (*element)->Nodes() );
}


/// could be considered a BoundaryInterface/Element/Node test
template<size_t dim>
void csmp::SplitBoundary_Test::CheckRemovedLowDimParents( Boundary<dim>& boundary )
{
  const typename vector<Node<dim>*>::const_iterator nodesEnd( boundary.NodesEnd() );
  for ( typename vector<Node<dim>*>::const_iterator node( boundary.NodesBegin() ); node != nodesEnd; ++node )
    for ( size_t parent( 0 ); parent < (*node)->Parents(); ++parent )
      _test( !isLowDim( (*node)->Parent( parent ) ) );
}


bool isLowDim( Element<3>* ePtr )
{
  return ePtr->FE()->IsSurfaceElement();
}


bool isLowDim( Element<2>* ePtr )
{
  return ePtr->FE()->IsLineElement();
}


template<size_t dim>
bool SplitBoundary_Test::NoNeighborNull( const csmp::InterFace<dim>& interFace )
{
  for ( size_t n( 0 ); n < interFace.Neighbors(); ++n )
    if ( interFace.Neighbor( n ) == NULL )
      return false;
  return true;
}


/// Test whether the splitted nodes share parent elements, located on different sides of interfaces
template<size_t dim>
void SplitBoundary_Test::TestSplitNodeAssignment( const Model<dim>& model )
{
  set<size_t> outerParentsNodes;
  set<size_t> innerParentsNodes;
  set<size_t> not_duplicated_nodes;
  not_duplicated_nodes.clear();

  // loop over SplitBoundaries
  for ( typename Model<dim>::splitBoundaryConstIterator spbit = model.SplitBoundariesBegin(); spbit != model.SplitBoundariesEnd(); ++spbit )
  {
    // loop over InterFaces
    for ( typename std::vector<InterFace<dim>* >::const_iterator ifit = spbit->second.ElementsBegin(); ifit != spbit->second.ElementsEnd(); ++ifit )
    {
      if ( NoNeighborNull( *(*ifit) ) )
      {
        // interior interfaces

        const Element<dim>* const outerParent( (*ifit)->Parent( OUTSIDE ) );
        for ( size_t n( 0 ); n < outerParent->Nodes(); ++n )
          outerParentsNodes.insert( outerParent->N( n )->Idx() );

        const Element<dim>* const innerParent( (*ifit)->Parent( INSIDE ) );
        for ( size_t n( 0 ); n < innerParent->Nodes(); ++n )
          innerParentsNodes.insert( innerParent->N( n )->Idx() );

        // if these test fails, the boundary was not split correctly
        for ( set<size_t>::const_iterator it( outerParentsNodes.begin() ); it != outerParentsNodes.end(); ++it )
          _test( innerParentsNodes.find( *it ) == innerParentsNodes.end() );

        for ( set<size_t>::const_iterator it( innerParentsNodes.begin() ); it != innerParentsNodes.end(); ++it )
          _test( outerParentsNodes.find( *it ) == outerParentsNodes.end() );

        outerParentsNodes.clear();
        innerParentsNodes.clear();

        for ( size_t n( 0 ); n < (*ifit)->Nodes(); ++n )
        {

          _test( (*ifit)->N( n, INSIDE )->Idx() != (*ifit)->N( n, OUTSIDE )->Idx() );

          set<size_t> parents;
          for ( size_t p( 0 ); p < (*ifit)->N( n, OUTSIDE )->Parents(); ++p )
            parents.insert( (*ifit)->N( n, OUTSIDE )->Parent( p )->Idx() );

          for ( size_t p( 0 ); p < (*ifit)->N( n, INSIDE )->Parents(); ++p )
            _test( parents.find( (*ifit)->N( n, INSIDE )->Parent( p )->Idx() ) == parents.end() );
          parents.clear();

          for ( size_t p( 0 ); p < (*ifit)->N( n, INSIDE )->Parents(); ++p )
            parents.insert( (*ifit)->N( n, INSIDE )->Parent( p )->Idx() );
          for ( size_t p( 0 ); p < (*ifit)->N( n, OUTSIDE )->Parents(); ++p )
            _test( parents.find( (*ifit)->N( n, OUTSIDE )->Parent( p )->Idx() ) == parents.end() );
          parents.clear();

        }

      }
      else {
        // perimeter interfaces
        for ( size_t i = 0; i<(*ifit)->Nodes(); ++i )
          if ( (*ifit)->N( i, INSIDE )->Idx() == (*ifit)->N( i, OUTSIDE )->Idx() )
            not_duplicated_nodes.insert( (*ifit)->N( i )->Idx() );
      }
    }

    // check whether the not duplicated nodes are sorted correctly
    set<size_t>::const_iterator nidx( not_duplicated_nodes.begin() );
    if ( not_duplicated_nodes.size() == spbit->second.PerimeterNodes() )
      for ( typename std::vector<Node<dim>* >::const_iterator nit = spbit->second.PerimeterNodesBegin(); nit != spbit->second.NodesEnd(); ++nit, ++nidx )
        _test( (*nidx) == (*nit)->Idx() );
  }
}


template<size_t dim>
void SplitBoundary_Test::TestUnitNormals( Model<dim>& model, const std::string& test_name )
{
  // -----------------------------
  // SplitBoundary Normals
  // -----------------------------

  std::vector<double64> zero_vec( dim, 0.0 );

  VectorVariable<dim> nrml_in( zero_vec );
  VectorVariable<dim> nrml_out( zero_vec );
  VectorVariable<dim> nrml_sum( zero_vec );
  VectorVariable<dim> nrml_zero( zero_vec );

  model.InputPropertyValue( "element vector", nrml_zero );

  Index normal_idx( model.Database().StorageKey( "element vector" ) );

  // loop over SplitBoundaries
  for ( typename Model<dim>::splitBoundaryConstIterator spbit = model.SplitBoundariesBegin(); spbit != model.SplitBoundariesEnd(); ++spbit )
  {
    // loop over InterFaces
    for ( typename std::vector<InterFace<dim>* >::const_iterator ifit = spbit->second.ElementsBegin(); ifit != spbit->second.ElementsEnd(); ++ifit )
    {
      (*ifit)->UnitNormal( nrml_in, INSIDE );
      (*ifit)->UnitNormal( nrml_out, OUTSIDE );

      nrml_sum = nrml_in + nrml_out;

      if(dim < 3) _test( nrml_sum == nrml_zero );
      (*ifit)->Parent( INSIDE )->Store( normal_idx, nrml_in );
      (*ifit)->Parent( OUTSIDE )->Store( normal_idx, nrml_out );
    }
  }

  return;
}


template<size_t dim>
void SplitBoundary_Test::VisualiseSplitBoundaries( Model<dim>& model, const std::string& test_name )
{
  // visualization
  VTU_Interface<dim> vtu( model );
  vtu.OmitZeroInFileName( false );
  list<string> outputProps;
  outputProps.push_back( "nodal id" );
  outputProps.push_back( "nodal variable" );
  outputProps.push_back( "element variable" );

  // -----------------------------
  // SplitBoundary Normals
  // -----------------------------

  const ScalarVariable matrixValue( ANY, 0.0 );
  ScalarVariable       interfaceValue( ANY, 0.0 );
  ScalarVariable       interfaceValueWrite( ANY, 0.0 );
  ScalarVariable       interfaceValueRead( ANY, 0.0 );

  model.InputPropertyValue( "nodal id", matrixValue );
  model.InputPropertyValue( "nodal variable", matrixValue );
  model.InputPropertyValue( "element variable", matrixValue );

  Index node_idx( model.Database().StorageKey( "nodal id" ) );
  Index node_prop_idx( model.Database().StorageKey( "nodal variable" ) );
  Index element_prop_idx( model.Database().StorageKey( "element variable" ) );

  // loop over SplitBoundaries
  for ( typename Model<dim>::splitBoundaryConstIterator spbit = model.SplitBoundariesBegin(); spbit != model.SplitBoundariesEnd(); ++spbit )
  {
    interfaceValue += 1.0;

    // loop over InterFaces
    for ( typename std::vector<InterFace<dim>* >::const_iterator ifit = spbit->second.ElementsBegin(); ifit != spbit->second.ElementsEnd(); ++ifit )
    {
      (*ifit)->Parent( INSIDE )->Read( element_prop_idx, interfaceValueRead );
      interfaceValueWrite = +1;
      (*ifit)->Parent( INSIDE )->Store( element_prop_idx, interfaceValueWrite );
      (*ifit)->Parent( OUTSIDE )->Read( element_prop_idx, interfaceValueRead );
      interfaceValueWrite = -1;
      (*ifit)->Parent( OUTSIDE )->Store( element_prop_idx, interfaceValueWrite );

      for ( size_t i = 0; i<(*ifit)->Nodes(); ++i )
      {
        if ( (*ifit)->N( i, OUTSIDE )->Idx() != (*ifit)->N( i, INSIDE )->Idx() )
        {
          interfaceValueWrite = -1.0;
          (*ifit)->N( i, INSIDE )->Store( node_prop_idx, interfaceValueWrite );
          interfaceValueWrite = (*ifit)->N( i, INSIDE )->Idx();
          (*ifit)->N( i, INSIDE )->Store( node_idx, interfaceValueWrite );

          interfaceValueWrite = 1.0;
          (*ifit)->N( i, OUTSIDE )->Store( node_prop_idx, interfaceValueWrite );
          interfaceValueWrite = (*ifit)->N( i, OUTSIDE )->Idx();
          (*ifit)->N( i, OUTSIDE )->Store( node_idx, interfaceValueWrite );
        }
      }
    }
  }

  if ( verbose_ ) vtu.OutputDataToVTU( test_name.c_str(), outputProps, "Model", 0.0 );

  for ( typename Model<dim>::regionIterator it = model.UniqueRegionsBegin(); it != model.UniqueRegionsEnd(); it++ )
    vtu.OutputDataToVTU( (*it).first.c_str(), outputProps, (*it).first.c_str(), 0.0 );
  return;
}


template<size_t dim>
void SplitBoundary_Test::PullApartSplitboundaries( Model<dim>& model, std::vector<std::string>& interfaces, double64 displacement )
{
  // pull apart SplitBoundaries
  std::string split_boundary_name;

  size_t split_boundaries( model.SplitBoundaries() );
  for ( size_t i = 0; i<split_boundaries; i++ )
  {
    split_boundary_name = "SPLITBOUNDARY_" + interfaces[i];
    shiftSplitBoundary( model.SplitBoundary( split_boundary_name.c_str() ), displacement );
    shiftInterfaceTips( model.Region( interfaces[i].c_str() ), displacement );
  }
}


template<size_t dim>
void SplitBoundary_Test::PullApartSplitboundaries( Model<dim>& model, double64 displacement )
{
  //pPull apart SplitBoundaries
  for ( typename Model<dim>::splitBoundaryIterator spbit = model.SplitBoundariesBegin(); spbit != model.SplitBoundariesEnd(); ++spbit )
    shiftSplitBoundary( (*spbit).second, displacement );
}


template<size_t dim>
void SplitBoundary_Test::LoadModel( const std::string& model_name )
{
  // Model initialization
  const std::string variables_file("SplitBoundary_Test-variables.txt");
  Model<dim>* model = NULL;
  if ( dim == 2U )
    model = dynamic_cast<Model<dim>*>(new ANSYS_Model2D( model_name.c_str(), model_name.c_str(), variables_file.c_str(), false, true, true, true ));
  else if ( dim == 3U )
    model = dynamic_cast<Model<dim>*>(new ANSYS_Model3D( model_name.c_str(), model_name.c_str(), variables_file.c_str(), false, true, true, true ));

  // visualization
  VTU_Interface<dim> vtu( *model );
  vtu.OmitZeroInFileName( false );
  list<string> outputProps;
  outputProps.push_back( "element variable" );
  ScalarVariable regionValue( ANY, 0.0 );

  if ( verbose_ ) cout << "\n\n\nSplitBoundary_Test::PrepareModel: the following interface / interface(s) sets will be considered:\n\n";
  model->InputPropertyValue( "element variable", regionValue );
  for ( typename Model<dim>::regionIterator it = model->RegionsBegin(); it != model->RegionsEnd(); it++ )
  {
    regionValue += 1.0;
    (*it).second.InputPropertyValue( "element variable", regionValue );
    vtu.OutputDataToVTU( model_name.c_str(), outputProps, (*it).first.c_str(), 0.0 );
  }
  if ( verbose_ ) vtu.OutputDataToVTU( model_name.c_str(), outputProps, "Model", 0.0 );

  model->OutputToBinaryFile( model_name.c_str() );
}


template<size_t dim>
void SplitBoundary_Test::LoadModel( const std::string& model_name,
                                    std::vector<std::string>& regions )
{
  // Model initialization
  const std::string variables_file("SplitBoundary_Test-variables.txt");
  Model<dim>* model = NULL;
  if ( dim == 2U )
    model = dynamic_cast<Model<dim>*>(new ANSYS_Model2D( model_name.c_str(), model_name.c_str(), variables_file.c_str(), false, true, true, true ));
  else if ( dim == 3U )
    model = dynamic_cast<Model<dim>*>(new ANSYS_Model3D( model_name.c_str(), model_name.c_str(), variables_file.c_str(), false, true, true, true ));

  // visualization
  VTU_Interface<dim> vtu( *model );
  vtu.OmitZeroInFileName( false );
  list<string> outputProps;
  outputProps.push_back( "element variable" );
  ScalarVariable regionValue( ANY, 0.0 );

  string  spliboundary_regions_file( model_name );

  // load regions from the input file
  inputFromFile( std::string( spliboundary_regions_file + "-disconnected-interface-regions.txt" ).c_str(), regions );

  if ( verbose_ ) cout << "\n\n\nSplitBoundary_Test::PrepareModel: the following interface / interface(s) sets will be considered:\n\n";
  model->InputPropertyValue( "element variable", regionValue );
  for ( std::vector<string>::const_iterator it = regions.begin(); it != regions.end(); it++ )
  {
    regionValue += 1.0;
    model->Region( (*it).c_str() ).InputPropertyValue( "element variable", regionValue );
    vtu.OutputDataToVTU( model_name.c_str(), outputProps, (*it).c_str(), 0.0 );
  }
  if ( verbose_ ) vtu.OutputDataToVTU( model_name.c_str(), outputProps, "Model", 0.0 );

  model->OutputToBinaryFile( model_name.c_str() );
}


template<size_t dim>
void SplitBoundary_Test::LoadContiguousModel( const std::string& model_name,
                                              std::vector<std::string>& interfaces )
{
  // 0. Model initialization
  const string variables_file("SplitBoundary_Test-variables.txt");
  Model<dim>* model = NULL;
  if ( dim == 2U )
    model = dynamic_cast<Model<dim>*>(new ANSYS_Model2D( model_name.c_str(), model_name.c_str(), variables_file.c_str(), false, true, true, true ));
  else if ( dim == 3U )
    model = dynamic_cast<Model<dim>*>(new ANSYS_Model3D( model_name.c_str(), model_name.c_str(), variables_file.c_str(), false, true, true, true ));

  string  spliboundary_regions_file( model_name );

  // 1. InterFace sets
  std::set<std::string> interface_basic_sets;
  std::set<std::string> interface_sets;
  inputFromFile( std::string( spliboundary_regions_file + "-disconnected-interface-regions.txt" ).c_str(), interface_basic_sets );

  // 2. splitting input regions if they are discontigouos
  bool discontiguous_regions( false );
  for ( set<string>::const_iterator it = interface_basic_sets.begin(); it != interface_basic_sets.end(); ++it )
    if ( !model->IsContiguous( (*it).c_str() ) ) {
      if ( verbose_ ) cerr << "\n\tSplitBoundary_Test::PrepareModel: discovered discontiguous region: " << (*it);
      discontiguous_regions = true;
    }
  if ( discontiguous_regions ) {
    set<string>  original_region_names;
    for ( typename RegionInterface<dim, Region>::regionIterator it = model->UniqueRegionsBegin(); it != model->UniqueRegionsEnd(); ++it )
      if ( (*it).first != "Model" )
        original_region_names.insert( (*it).first.c_str() );

    // partitioning regions without revisiting new partitions that can inserted into region map
    for ( set<string>::const_iterator it = original_region_names.begin(); it != original_region_names.end(); ++it )
      model->PartitionRegionIntoContiguousSubRegions( (*it).c_str() );
  }

  // 2. preparing low dimensional regions for making SplitBoundaries around
  etablishContiguosRegionsList( *model, interface_basic_sets, interface_sets );

  model->MergeRegions( interface_sets, "interfaces" );

  for ( auto name : interface_sets )
    model->RemoveRegion( name.c_str(), false );

  interfaces.clear();
  interfaces.push_back( "interfaces" );
  outputToFile( std::string( spliboundary_regions_file + "-connected-interface-regions.txt" ).c_str(), interfaces );

  if ( verbose_ ) cout << "\n\n\nSplitBoundary_Test::PrepareModel: the following interface(s) / interface sets will be considered:\n\n";

  // 3. creating SplitBoundaries
  for ( std::vector<string>::const_iterator it = interfaces.begin(); it != interfaces.end(); it++ )
    model->InsertSplitBoundary( (*it).c_str() );

  // 4. creating lower-dimensional stand-alone meshes from SplitBoundary objects, and 
  //    insert them into a new sub-region (simply named by 'SPLITBOUNDARY_SURFACE')
  model->InsertRegionFromSplitBoundaries();

  // 5. Writing it 
  model->OutputToBinaryFile( model_name.c_str() );
}


/// Input name of regions to split
void SplitBoundary_Test::inputFromFile( const char* file_name,
                                        std::set<string>& interface_basic_set )
{
  assert( file_name != NULL );

  ifstream ifs( file_name );
  if ( !ifs.is_open() )
    throw csmp::Exception( ERROR, "SplitBoundary_Test::inputFromFile",
                           file_name,
                           "file specifying regions is missing" );
  if ( !interface_basic_set.empty() )
    interface_basic_set.clear();

  // reading header line printing it to screen and swallowing empty line thereafter
  char text[256];
  ifs.getline( text, 256 );
  if ( verbose_ ) cout << "\nSplitBoundary_Test::inputFromFile: file header: " << text << endl;
  ifs.getline( text, 256 );

  int n_interfaces( 0 );
  ifs >> n_interfaces;
  assert( n_interfaces > 0 );
  assert( n_interfaces < 10000 );

  string  interface_name;
  for ( int n = 0; n<n_interfaces; n++ ) {
    ifs >> interface_name;
    if ( !interface_name.empty() )
      interface_basic_set.insert( interface_name );
    else
      throw csmp::Exception( ERROR, "SplitBoundary_Test::inputFromFile:", "encountered empty region name." );
  }

  ifs.close();

  if ( verbose_ ) cout << "\nSplitBoundary_Test::inputFromFile: region names stored in '" << file_name << "' read successfully." << endl;

} // end inputFromFile


/// Input name of regions to split
template<size_t dim>
void SplitBoundary_Test::etablishContiguosRegionsList( Model<dim>& model,
                                                       const std::set<string>& interface_basic_set,
                                                       std::set<string>& interface_sets )
{
  if ( verbose_ ) cout << "\nSplitBoundary_Test::etablishContiguosRegionsList:" << endl;

  if ( !interface_sets.empty() )
    interface_sets.clear();
  bool first_call( true );
  for ( typename map<string, csmp::Region<dim> >::const_iterator
        git = model.UniqueRegionsBegin(); git != model.UniqueRegionsEnd(); git++ )
  {
    for ( set<string>::const_iterator it = interface_basic_set.begin(); it != interface_basic_set.end(); ++it )
      if ( (*git).first.find( *it ) != string::npos )
      {
        interface_sets.insert( (*git).first );
        if ( first_call ) {
          cout << "\n\tregion(s) incorporated into the input list: ";
          first_call = false;
        }
        cout << (*git).first << " ";
        break;
      }
  }
  if ( verbose_ ) cout << "\n\n";

} // end inputFromFile


/// Input name of regions to split
void SplitBoundary_Test::inputFromFile( const char* file_name,
                                        std::vector<string>& interfaces )
{
  assert( file_name != NULL );

  ifstream ifs( file_name );
  if ( !ifs.is_open() )
    throw csmp::Exception( ERROR, "SplitBoundary_Test::inputFromFile",
                           file_name,
                           "file specifying regions is missing" );
  if ( !interfaces.empty() )
    interfaces.clear();

  // reading header line printing it to screen and swallowing empty line thereafter
  char text[256];
  ifs.getline( text, 256 );
  if ( verbose_ ) cout << "\nSplitBoundary_Test::inputFromFile: file header: " << text << endl;
  ifs.getline( text, 256 );

  int n_interfaces( 0 );
  ifs >> n_interfaces;
  assert( n_interfaces > 0 );
  assert( n_interfaces < 10000 );

  string  interface_name;
  for ( int n = 0; n<n_interfaces; n++ ) {
    ifs >> interface_name;
    if ( !interface_name.empty() )
      interfaces.push_back( interface_name );
    else
      throw csmp::Exception( ERROR, "SplitBoundary_Test::inputFromFile:", "encountered empty region name." );
  }

  ifs.close();

  if ( verbose_ ) cout << "\nSplitBoundary_Test::inputFromFile: region names stored in '" << file_name << "' read successfully." << endl;

} // end inputFromFile


/// Write contiguous regions
void SplitBoundary_Test::outputToFile( const char* file_name,
                                       const std::vector<string>& interfaces )
{
  assert( file_name != NULL );

  ofstream ofs( file_name );
  if ( !ofs.is_open() )
    throw csmp::Exception( ERROR, "SplitBoundary_Test::outputToFile",
                           file_name, "file specifying permeability-model input variables is missing" );

  ofs << "'" << file_name << "' interface regions to be included.\n\n";

  ofs << interfaces.size() << " ";
  for ( std::vector<string>::const_iterator it = interfaces.begin(); it != interfaces.end(); ++it )
    ofs << (*it) << " ";

  ofs << "\n";
  ofs.close();

  if ( verbose_ ) cout << "\nSplitBoundary_Test::outputToFile: '" << file_name << "' written successfully." << endl;

} // end outputToFile


/// Write interface regions
void SplitBoundary_Test::outputToFile( const char* file_name,
                                       const std::set<string>& interfaces )
{
  assert( file_name != NULL );

  ofstream ofs( file_name );
  if ( !ofs.is_open() )
    throw csmp::Exception( ERROR, "SplitBoundary_Test::outputToFile",
                           file_name, "file specifying permeability-model input variables is missing" );

  ofs << "'" << file_name << "' interface regions to be included.\n\n";

  ofs << interfaces.size() << " ";
  for ( std::set<string>::const_iterator it = interfaces.begin(); it != interfaces.end(); ++it )
    ofs << (*it) << " ";

  ofs << "\n";
  ofs.close();

  if ( verbose_ ) cout << "\nSplitBoundary_Test::outputToFile: '" << file_name << "' written successfully." << endl;

} // end outputToFile


/// TESTS
/// SPLITBOUNDARY BETWEEN REGIONS
template<size_t dim>
void SplitBoundary_Test::test_splitboundary_between_regions( const std::string& model_name )
{
  std::ostringstream ostr;
  std::string dimension( "" );
  ostr << dim;
  dimension += ostr.str();
  dimension += "D";

  if ( verbose_ ) std::cerr << "\nStart " << dimension << " SplitBoundary Test: SplitBoundary between Regions\n";

  std::string test_name( "SPLITBOUNDARY_TEST_BETWEEN_REGIONS_" );
  test_name += dimension;
  test_name += "_";
  test_name += model_name;

  // load Model
  const std::string variables_file( "SplitBoundary_Test-variables.txt" );

  Model<dim>* modelIN = NULL;

  if ( dim == 2U )
    modelIN = dynamic_cast<Model<dim>*>(new ANSYS_Model2D( model_name.c_str(), model_name.c_str(), variables_file.c_str(), false, true, true, true ));
  else if ( dim == 3U )
    // SKM FIX: irregular = true
    modelIN = dynamic_cast<Model<dim>*>(new ANSYS_Model3D( model_name.c_str(), model_name.c_str(), variables_file.c_str(), true, true, true, true ));

  // validating the model
  if ( verbose_ ) cout << "\nSplitBoundary_Test<" << dim << ">::test_splitboundary_between_regions: model contains the regions:";
  for ( auto it = modelIN->RegionsBegin(); it != modelIN->RegionsEnd(); it++ )
    if ( verbose_ ) cout << "\n\tboundary: " << (*it).first;
  if ( verbose_ ) cout << endl;
  if ( verbose_ ) cout << "\nSplitBoundary_Test<" << dim << ">::test_splitboundary_between_regions: model contains the boundaries:";
  for ( auto it = modelIN->BoundariesBegin(); it != modelIN->BoundariesEnd(); it++ )
    if ( verbose_ ) cout << "\n\tboundary: " << (*it).first;
  if ( verbose_ ) cout << endl;

  // create SplitBoundaries
  std::vector<std::string> regions;
  regions.reserve( modelIN->UniqueRegions() );

  const pair<int32, int32>  model_dim = modelIN->Region( "Model" ).SpatialDimensions();
  for ( typename std::map<std::string, csmp::Region<dim> >::iterator
        it = modelIN->UniqueRegionsBegin(); it != modelIN->UniqueRegionsEnd(); ++it ) {
    const pair<int32, int32>  sub_dim = (*it).second.SpatialDimensions();
    if ( sub_dim.second == model_dim.second ) // check whether the highest dimension of the region is equal to the highest dimension of the model
      regions.push_back( (*it).second.Name() );
  }

  // search the existing unique sub-regions
  set<pair<string, string>>	discovered;
  deque<string>	current_regions;
  vector<pair<string, string>>  region_final_pairs;
  for ( size_t i = 0U; i < regions.size(); i++ ) {
    string root = regions[i];
    // starting at the first region
    current_regions.push_back( root );
    while ( !current_regions.empty() ) {
      std::string current_region( *current_regions.begin() );
      set<string> neighbors;
      // for all neighbor sub-regions of the current region
      const csmp::Region<dim>&  gref1( modelIN->Region( current_region ) );
      for ( size_t j = 0; j < regions.size(); j++ ) {
        if ( current_region.compare( regions[j] ) == 0 ) continue;
        const csmp::Region<dim>&  gref2( modelIN->Region( regions[j] ) );
        const size_t  shared_nodes( sharedNodes( gref1, gref2 ) );
        if ( shared_nodes > 1 )
          neighbors.insert( regions[j] );
      }

      for ( auto neighbour_region : neighbors ) {
        // if this neighbor is new one        
        auto new_region = discovered.insert( make_pair( current_region, neighbour_region ) );
        if ( new_region.second ) {
          discovered.insert( make_pair( neighbour_region, current_region ) );
          current_regions.push_back( neighbour_region );
          region_final_pairs.push_back( make_pair( current_region, neighbour_region ) );
        }
      }
      // removing the sub-region from the discovered (but not yet explored) deque
      current_regions.pop_front();
    }
  }

  // do this sequentially according to neighbours, otherwise boundaries are not assgiend properly
  for ( auto it : region_final_pairs ) // for each of the boundary patches discovered, a uniquely named SplitBoundary object is created
    modelIN->InsertSplitBoundary( it.first, it.second );
  
  // create lower-dimensional stand-alone meshes from SplitBoundary objects, and 
  // insert them into a new sub-region (simply named by 'SPLITBOUNDARY_SURFACE')
  modelIN->InsertRegionFromSplitBoundaries();

  // read Model from Binary
  modelIN->OutputToBinaryFile( model_name.c_str() );
  Model<dim> model( model_name.c_str() );
  cout << "\nNodes: " << model.Mesh().Nodes() << "\n";
  cout << "\nNode Groups: " << model.Mesh().NodeGroups() << "\n";
  cout << "\nElements: " << model.Mesh().Elements() << "\n";
  cout << "\nElement Groups: " << model.Mesh().ElementGroups() << "\n";
  cout << "\nFaces: " << model.Mesh().Faces() << "\n";
  cout << "\nFace Groups: " << model.Mesh().FaceGroups() << "\n";
  cout << "\nInterfaces: " << model.Mesh().InterFaces() << "\n";
  cout << "\nInterface Groups: " << model.Mesh().InterFaceGroups() << "\n";

  // tests
  TestSplitNodeAssignment( model );
  TestUnitNormals( model, test_name.c_str() );
  double64 displacement( 0.001 );
  PullApartSplitboundaries( model, displacement );
  VisualiseSplitBoundaries( model, test_name );

  return;
}


/// SPLITBOUNDARY AROUND REGIONS
template<size_t dim>
void SplitBoundary_Test::test_splitboundary_around_regions( const std::string& model_name )
{
  std::ostringstream ostr;
  std::string dimension( "" );
  ostr << dim;
  dimension += ostr.str();
  dimension += "D";

  if ( verbose_ ) std::cerr << "\nStart " << dimension << " SplitBoundary Test: SplitBoundary around Regions\n";
    
  // load Model
  std::vector<std::string> interfaces;
  LoadContiguousModel<dim>( model_name, interfaces );

  // read from binary
  string spliboundary_regions_file( model_name );
  spliboundary_regions_file += "-connected-interface-regions.txt";
  inputFromFile( spliboundary_regions_file.c_str(), interfaces );

  Model<dim> model_out( model_name.c_str() );
  cout << "\nNodes: " << model_out.Mesh().Nodes() << "\n";
  cout << "\nNode Groups: " << model_out.Mesh().NodeGroups() << "\n";
  cout << "\nElements: " << model_out.Mesh().Elements() << "\n";
  cout << "\nElement Groups: " << model_out.Mesh().ElementGroups() << "\n";
  cout << "\nFaces: " << model_out.Mesh().Faces() << "\n";
  cout << "\nFace Groups: " << model_out.Mesh().FaceGroups() << "\n";
  cout << "\nInterfaces: " << model_out.Mesh().InterFaces() << "\n";
  cout << "\nInterface Groups: " << model_out.Mesh().InterFaceGroups() << "\n";
  
  // visualization
  std::string test_name( "SPLITBOUNDARY_TEST_AROUND_REGIONS_" );
  test_name += dimension;
  test_name += "_";
  test_name += model_name;
  VisualiseSplitBoundaries( model_out, test_name );
  
  if ( verbose_ ) std::cerr << "\nFinish " << dimension << " SplitBoundary Test: SplitBoundary around Regions\n";
}


template<size_t dim>
void SplitBoundary_Test::detect_and_create_splitboundaries( const std::string& model_name )
{
  const bool verbose( false );
  const string variables_file( "SplitBoundary_Test-variables.txt" );

  // 1. convert ansys model into CSMP model
  Model<dim>* model = NULL;
  if ( dim == 2U )
    model = dynamic_cast<Model<dim>*>(new ANSYS_Model2D( model_name.c_str(), model_name.c_str(), variables_file.c_str(), false, true, true, true ));
  else if ( dim == 3U )
    model = dynamic_cast<Model<dim>*>(new ANSYS_Model3D( model_name.c_str(), model_name.c_str(), variables_file.c_str(), false, true, true, true ));

  // 2. build CSMP SplitBoundary
  // ---------------------------------------------------------------------------------------
  model->DetectAndCreateSplitBoundaries();

  // 3. see whether the split boundary survives being writting to and recovered from file
  // ------------------------------------------------------------------------------------------
  model->OutputToBinaryFile( model_name.c_str() );

  csmp::Model<dim> model_out( model_name.c_str() );
  cout << "\nNodes: " << model_out.Mesh().Nodes() << "\n";
  cout << "\nNode Groups: " << model_out.Mesh().NodeGroups() << "\n";
  cout << "\nElements: " << model_out.Mesh().Elements() << "\n";
  cout << "\nElement Groups: " << model_out.Mesh().ElementGroups() << "\n";
  cout << "\nFaces: " << model_out.Mesh().Faces() << "\n";
  cout << "\nFace Groups: " << model_out.Mesh().FaceGroups() << "\n";
  cout << "\nInterfaces: " << model_out.Mesh().InterFaces() << "\n";
  cout << "\nInterface Groups: " << model_out.Mesh().InterFaceGroups() << "\n";

  // 4. visualising
  std::string test_name( "DETECTED_SPLITBOUNDARY_TEST_FROM_" );
  test_name += model_name;
  VisualiseSplitBoundaries( model_out, test_name );

  return;
}

template<size_t dim>
void SplitBoundary_Test::detect_and_create_splitboundaries_from_constructor( const std::string& model_name )
{
  const bool verbose( false );
  const string variables_file( "SplitBoundary_Test-variables.txt" );
  
  // 1. convert ansys model into CSMP model
  Model<dim>* model = NULL;
  const bool create_splitboundaries( true );
  if ( dim == 2U )
    model = dynamic_cast<Model<dim>*>(new ANSYS_Model2D( model_name.c_str(), model_name.c_str(), variables_file.c_str(), false, true, true, true, create_splitboundaries ));
  else if ( dim == 3U )
    model = dynamic_cast<Model<dim>*>(new ANSYS_Model3D( model_name.c_str(), model_name.c_str(), variables_file.c_str(), false, true, true, true, create_splitboundaries ));
  
  // 2. create lower - dimensional stand - alone meshes from SplitBoundary objects, and
  //    insert them into a new sub-region (simply named by 'SPLITBOUNDARY_SURFACE')
  model->InsertRegionFromSplitBoundaries();

  // 3. see whether the split boundary survives being writting to and recovered from file
  // ------------------------------------------------------------------------------------------
  model->OutputToBinaryFile( model_name.c_str() );

  csmp::Model<dim> model_out( model_name.c_str() );
  cout << "\nNodes: " << model_out.Mesh().Nodes() << "\n";
  cout << "\nNode Groups: " << model_out.Mesh().NodeGroups() << "\n";
  cout << "\nElements: " << model_out.Mesh().Elements() << "\n";
  cout << "\nElement Groups: " << model_out.Mesh().ElementGroups() << "\n";
  cout << "\nFaces: " << model_out.Mesh().Faces() << "\n";
  cout << "\nFace Groups: " << model_out.Mesh().FaceGroups() << "\n";
  cout << "\nInterfaces: " << model_out.Mesh().InterFaces() << "\n";
  cout << "\nInterface Groups: " << model_out.Mesh().InterFaceGroups() << "\n";

  // 4. visualising
  std::string test_name( "CREATED_SPLITBOUNDARY_TEST_FROM_" );
  test_name += model_name;
  double64 displacement( 0.001 );
  PullApartSplitboundaries( model_out, displacement );
  VisualiseSplitBoundaries( model_out, test_name );   

  return;
}

void SplitBoundary_Test::run()
{
  // test splitboundary between 2D regions
  test_splitboundary_between_regions<2U>( "BoxHalfs2D" );
  test_splitboundary_between_regions<2U>( "ThreeZones2D" );

  // test splitboundary between 3D regions
  test_splitboundary_between_regions<3U>( "BoxHalfs3D" );
  test_splitboundary_between_regions<3U>( "ThreeZones3D" );

  // test splitboundary around interfaces
  // JC: working on the QC process which is requried for the following models
  //test_splitboundary_around_regions<2U>( "UnitSquareFracs_xline" );
  //test_splitboundary_around_regions<2U>( "UnitSquareFracs_yline" );
  //test_splitboundary_around_regions<2U>( "UnitSquareFracs_orthogonal" );
  //test_splitboundary_around_regions<2U>( "UnitSquareFracs_irregular" );

  // test splitboundary from constructor of ansys model
  detect_and_create_splitboundaries_from_constructor<2U>( "Jura-slope1" );
  detect_and_create_splitboundaries_from_constructor<3U>( "Dyke_Split" );

  // test splitboundary for complex ansys models
  test_splitboundary_between_regions<3U>( "lamination" );
  test_splitboundary_between_regions<2U>( "kueper_one_interface" );

  return;
}

} // csmp
