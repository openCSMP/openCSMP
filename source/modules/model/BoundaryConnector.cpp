#include "BoundaryConnector.h"

#include "Model.h"
#include "Element.h"
#include "InterFace.h"
#include "Face.h"
#include "Node.h"
#include "ErrorHandler.h"
#include "Region.h"
#include "Boundary.h"
#include "PL_Utilities.h"
#include "SplitBoundary.h"
#include "variableOperations.h"


using namespace std;

namespace csmp
{


/** Establishes Face<dim> parents from equivalent Element<dim> and initializes Face nodes

If provided with an eligible Element type (surface for 3D, line for 1D), this method
finds the higher dimensional parent element for a given Face. It uses the Face::Assign
member which initializes the internal node container.

@attention Assure that member indices are updated before this method is called.
@attention This method will probably not work if face is at a split boundary, due to duplicate nodes

@section ordering Parent Ordering 

If only one parent element attached to face nodes is found(i.e. the face is located at a model
boundary), this becomes the inner element, thus the face normal is pointing outward of it.

If two parent elements are found, the element with the lower index is made the inner
parent element of the face, hence unit normal points toward the element with the higher index.

@param  face The csmp::Face for which parents are to be found and nodes to be established
@param  element The csmp::Element which equals the face in FEM type and connectivity
@return Number of parents found, zero if none.
*/
template<size_t dim>
size_t BoundaryConnector<dim>::FaceParentsFromElementEquivalent( Face<dim>* face, Element<dim>* const element,
                                                                 std::map<size_t,int>& elementUnitNormalOrientation )
  {
    ErrorHandler& errorHandler( ErrorHandler::Instance() );
    // element fem type equals face fem type
    if( face->FE_Type() != element->FE_Type() )
      errorHandler.notice( csmp::ERROR, "BoundaryConnector<dim>::FaceParentsFromElementEquivalent",
                                        "Finite element type of face and element are not equivalent" );
    // face=element node count
    const size_t elementNodes( element->Nodes() );
    // potential parent element idx container, implemented as map because:
    // automatically orders ids, so the first entry is the inner one (see convention)
    // further it disallows duplicate entries, thus counteracting potential bug.
    // we keep track of how many nodes share same parent element, if number is equal
    // to total number of face nodes the parent element is a valid one
    map<size_t,pair<size_t,Element<dim>*> > potentialParents;
    //  ^^^         ^^^^^   ^^^
    //  idx         count   ptr
    // looping over nodes adding potential parents, but not current (face equivalent)element itself
    for( size_t i = 0; i < elementNodes; ++i )
      {
        const size_t nodeParentCount( element->N(i)->Parents() );
        for( size_t j = 0; j < nodeParentCount; ++j )
          if( element->N(i)->Parent(j)->Idx() != element->Idx() )
            AddPotentialFaceParentElement( element->N(i)->Parent(j), potentialParents );                
      } // face nodes
    // let's see how many parents we found
    std::pair<Element<dim>*,Element<dim>*> parents( make_pair<Element<dim>*,Element<dim>*>( NULL, NULL ) );
    const size_t parentCount( EligibleFaceParentElements( elementNodes, parents, potentialParents ) );
    // we expect at least one(boundary) and max two(interior) and assign them to the face
    if( parentCount == 1 )
      {
        // we need to provide the face node vector
        vector<Node<dim>*> faceNodes( elementNodes, NULL );
        for( size_t i = 0; i < elementNodes; ++i )
          faceNodes[i] = element->N(i);
          // SKM FIX
        face->Assign( parents.first, static_cast<Element<dim>*>(nullptr) );
      }
    else if( parentCount == 2 )
      {
        OuterAndInnerParent( element, parents, elementUnitNormalOrientation );
        face->Assign( parents.first, parents.second );
      }
    else
      errorHandler.notice( csmp::FATAL_ERROR, "BoundaryConnector<dim>::FaceParentsFromElementEquivalent",
                                              "Invalid number of parent elements found" );
    // done
    return parentCount;
  } // FaceParentsFromElementEquivalent



/// adds element with count = 1 if not in map yet, otherwise increases count for element by 1
template<size_t dim>
size_t BoundaryConnector<dim>::AddPotentialFaceParentElement( Element<dim>* potentialParent, 
                                                              map<size_t,std::pair<size_t,Element<dim>*> >& potentialParents )
  {
    if( potentialParent->Placement() != ELEMENT )
      return 0;
    size_t potentialElementRequests( 1 );
    // check whether element already registered as potential parent
    typename map<size_t,pair<size_t,Element<dim>*> >::iterator potentialParentEntry( potentialParents.find( potentialParent->Idx()  ) );
    // if element is new count = 1, else increase count by one
    if( potentialParentEntry == potentialParents.end() )
      potentialParents.insert( make_pair( potentialParent->Idx(), make_pair( potentialElementRequests, potentialParent ) ) );
    else
      {
        potentialElementRequests = potentialParentEntry->second.first;
        potentialParentEntry->second.first = ++potentialElementRequests;
      }
    return potentialElementRequests;
  } // AddPotentialParentElement



/** searches the map for possible parents with occurrences equal to face nodes and inserts to pair 

The inner element(pair.first) will be the one with the lower index.
Throws if more than two elements are found.

@return The number of eligible parents (0, 1 or 2)
*/
template<size_t dim>
size_t BoundaryConnector<dim>::EligibleFaceParentElements( size_t faceNodes, 
                                                           pair<Element<dim>*,Element<dim>*>& faceParentElements,
                                                           const map<size_t,std::pair<size_t,Element<dim>*> >& potentiaFaceParentElements )
  {
//    ErrorHandler& errorHandler( ErrorHandler::Instance() );
    size_t parentCount( 0 );
    const typename map<size_t,pair<size_t,Element<dim>*> >::const_iterator potentialParentsEnd( potentiaFaceParentElements.end() );
    for( typename map<size_t,pair<size_t,Element<dim>*> >::const_iterator it = potentiaFaceParentElements.begin(); it != potentialParentsEnd; ++it )
      {
        if( it->second.first == faceNodes )
          {
            if( parentCount == 0 )
              faceParentElements.first = it->second.second;
            else if( parentCount == 1 )
              faceParentElements.second = it->second.second;
            else
              throw csmp::Exception( csmp::FATAL_ERROR, "BoundaryConnector<dim>::EligibleFaceParentElements", 
                                                        "More than two eligible parent elements found" );
            ++parentCount;
          }
      } // potentiaFaceParentElements
    return parentCount;
  }


/** Returns true and sets provided node vector if a face of the element lies on the given boundary

@attention For equidimensional elements only!
*/
template<size_t dim>
bool BoundaryConnector<dim>::BoundaryFaceNodes( const Element<dim>& element, BOX_BOUNDARY boxBoundary, vector<Node<dim>*>& boundaryFaceNodes, size_t& bFace )
  {
    bFace = 999;
    bool hasBoundaryFace(true);
    boundaryFaceNodes.clear();
    const size_t faces( element.FE()->Faces() );
    for( size_t face(0); face < faces; ++face )
      {
        const size_t nodesPerFace( element.FE()->NodesPerFace(face) );
        vector<size_t> nodesOfFace( nodesPerFace, 0 );
        hasBoundaryFace = true;
        element.FE()->NodesOfFace( face, nodesOfFace );
        for( size_t node(0); node < nodesPerFace; ++node )
          if( !belongsToSide( boxBoundary, element.N( nodesOfFace[node] )->AtBoundary() ) )
            { 
              hasBoundaryFace = false; 
              break; 
            }
        if( hasBoundaryFace )
          {
            bFace = face;
            for( size_t node = 0; node < nodesPerFace; ++node )
                boundaryFaceNodes.push_back( element.N( nodesOfFace[node] ) );
            break;
          }
      }
    return hasBoundaryFace;
  }


/// Auxiliary
template<size_t dim>
size_t connectingFace( Element<dim>* const element, vector<size_t>& fnids, const set<size_t>& boundaryFaceNodeIds )
  {
    const set<size_t>::const_iterator boundaryFaceNodeIdsEnd( boundaryFaceNodeIds.end() );
    size_t face(9999);
    for( size_t i = 0; i < element->Faces(); ++i )
      {
        element->FE()->NodesOfFace( i, fnids );
        const size_t faceNodes( fnids.size() );
        assert( faceNodes == boundaryFaceNodeIds.size() );
        for( size_t j = 0; j < faceNodes; ++j )
          if( boundaryFaceNodeIds.find( element->N( fnids[j] )->Idx() ) == boundaryFaceNodeIdsEnd )
            break;
          else if( j == faceNodes-1 )
            face = i;
        if( face == i )
          break;
      } // parent faces
    return face;
  }


double pointRelativeToPlane( const Point<1>&, const Point<1>&, const vector<double64>& )
  { throw; return -1.; }

double pointRelativeToPlane( const Point<2>& pointOnPlane, const Point<2>& pointToCheck,
                             const vector<double64>& planeUnitNormal )
  { return planeUnitNormal[0]*(pointOnPlane[0]-pointToCheck[0])+planeUnitNormal[1]*(pointOnPlane[1]-pointToCheck[1]); }

double pointRelativeToPlane( const Point<3>& pointOnPlane, const Point<3>& pointToCheck,
                             const vector<double64>& planeUnitNormal )
  { return planeUnitNormal[0]*(pointOnPlane[0]-pointToCheck[0])+planeUnitNormal[1]*(pointOnPlane[1]-pointToCheck[1])+planeUnitNormal[2]*(pointOnPlane[2]-pointToCheck[2]); }


/// Sets inner parent (wrt to outward pointing normal of element) as first and outer as second
template<size_t dim>
void BoundaryConnector<dim>::OuterAndInnerParent( Element<dim>* element, std::pair<Element<dim>*,Element<dim>*>& parents,
                                                  std::map<size_t,int>& elementUnitNormalOrientation )
  {
    if( !parents.second || !parents.first )
      return;
    // 0. get the face UN
    vector<double64> boundaryFaceUN(3,0.);
    element->CoordinateMatrix();
    element->FE()->UnitNormal( boundaryFaceUN ); 
    // 1. node ids of boundary face
    const size_t elementNodes( element->Nodes() );
    set<size_t> boundaryFaceNodeIds;
    for( size_t i = 0; i < elementNodes; ++i )
      boundaryFaceNodeIds.insert( element->N(i)->Idx() ); 
    const set<size_t>::const_iterator boundaryFaceNodeIdsEnd( boundaryFaceNodeIds.end() );
    // 2. non plane element node
    Node<dim> *parentsFirstNonFaceNode(NULL);
    for ( size_t i = 0; i < parents.first->Nodes(); ++i )
      if( boundaryFaceNodeIds.find( parents.first->N(i)->Idx() ) == boundaryFaceNodeIdsEnd )
        { parentsFirstNonFaceNode = parents.first->N(i); break; }
    assert( parentsFirstNonFaceNode ); 
    // 4. check relative location
    const Point<dim> facePoint( element->N(0)->Coordinate() );
    const Point<dim> parentsFirstNonFacePoint( parentsFirstNonFaceNode->Coordinate() );
    assert( elementUnitNormalOrientation.find( element->Idx() ) != elementUnitNormalOrientation.end() );
    if ( (static_cast<double>( elementUnitNormalOrientation.find( element->Idx() )->second ) *
         pointRelativeToPlane( facePoint, parentsFirstNonFacePoint, boundaryFaceUN ) ) < 0. )
      {
        Element<dim>* cachePtr(NULL);
        cachePtr = parents.first;
        parents.first = parents.second;
        parents.second = cachePtr;
      }
  }



template<size_t dim>
void elementUnitNormalVector( VectorVariable<dim>& nv, Element<dim>& element) 
  {
    vector<double> N( dim, 0. );
    element.CoordinateMatrix();
    element.FE()->UnitNormal( N );
    for( size_t d(0); d < dim; ++d )
      nv(d) = N[d];
  }

/**

Starting at one lower dimensional element, we establish the orientation of the neighbor elements unit normal.
This allows us to have a consistent outer parent assignment later on for the case where we wish to build a boundary
within a region based on a lower dimensional one. Due to the limitations that i.e. line elements can only have two neighbors,
which is of course not representative in intesecting structures, we have a non generic approach here. Hence:

@attention This will not work properly if a region turns more than 90 degrees!

@param [out] elementUnitNormalOrientation The map will be filled with the element idx and either OUTER or INNER (see enum) depending on the
                                          orientation of the elements unit normal wrt to the reference element

*/
template<size_t dim>
void BoundaryConnector<dim>::GlobalOrientation( const Region<dim>& region, std::map<size_t,int>& elementUnitNormalOrientation )
  {
    elementUnitNormalOrientation.clear();
    VectorVariable<dim> RNV( PLAIN, 0. ), NV( PLAIN, 0. );
    // first element UN orientation is arbitrary reference orientation
    const typename vector<Element<dim>*>::const_iterator elementsEnd( region.ElementsEnd() );
    typename vector<Element<dim>*>::const_iterator element( region.ElementsBegin() );
    elementUnitNormalOrientation.insert( make_pair( (*element)->Idx(), OUTER ) );
    elementUnitNormalVector( RNV, *(*element) );
    // establish UN orientation of region elements wrt RNV
    for ( ++element; element != elementsEnd; ++element )
      {
        elementUnitNormalVector( NV, *(*element) );
        if( dotProduct(NV,RNV) < 0. )
          elementUnitNormalOrientation.insert( make_pair( (*element)->Idx(), INNER ) );
        else
          elementUnitNormalOrientation.insert( make_pair( (*element)->Idx(), OUTER ) );
      } 
    assert( elementUnitNormalOrientation.size() == region.Elements() );
    return;
  }



/// returns whether initial unit normal of interface produces a positive scalar product shared node to face BC
template<size_t dim>
bool csmp::BoundaryConnector<dim>::OnOutside( InterFace<dim>& interFace, Face<dim>& face,  Node<dim>& sharedNode )
  {
    const Point<dim> bcFace( face.BaryCenter() ), coNode( sharedNode.Coordinate() );
    VectorVariable<dim> nodeToFace( PLAIN, 0. ), interFaceUnitNormal( PLAIN, 0. );
    interFace.UnitNormal( interFaceUnitNormal, INSIDE);
    for( size_t d(0); d < dim; ++d )
      nodeToFace(d) = coNode[d]-bcFace[d];
    return( dotProduct(nodeToFace,interFaceUnitNormal) < 0. );
  }


/// returns whether initial unit normal of interface produces a positive scalar product shared node to interface BC
template<size_t dim>
bool csmp::BoundaryConnector<dim>::OnOutside( InterFace<dim>& interFace, InterFace<dim>& interFaceInQuestion,  Node<dim>& sharedNode )
  {
    const Point<dim> bcInterFace( interFaceInQuestion.BaryCenter() ), coNode( sharedNode.Coordinate() );
    VectorVariable<dim> nodeToInterFace( PLAIN, 0. ), interFaceUnitNormal( PLAIN, 0. );
    interFace.UnitNormal( interFaceUnitNormal, INSIDE);
    for( size_t d(0); d < dim; ++d )
      nodeToInterFace(d) = coNode[d]-bcInterFace[d];
    return( dotProduct(nodeToInterFace,interFaceUnitNormal) < 0. );
  }

/// returns whether initial unit normal of interface produces a positive scalar product shared node to element BC
template<size_t dim>
bool csmp::BoundaryConnector<dim>::OnOutside( InterFace<dim>& interFace, Element<dim>& element,  Node<dim>& sharedNode )
  {
    const Point<dim> bcElement( element.BaryCenter() ), coNode( sharedNode.Coordinate() );
    VectorVariable<dim> nodeToElement( PLAIN, 0. ), interFaceUnitNormal( PLAIN, 0. );
    interFace.UnitNormal( interFaceUnitNormal, INSIDE);
    for( size_t d(0); d < dim; ++d )
      nodeToElement(d) = coNode[d]-bcElement[d];
    return( dotProduct(nodeToElement,interFaceUnitNormal) < 0. );
  }


template<size_t dim>
bool BoundaryConnector<dim>::BoundaryNode( Node<dim>& node, Model<dim>& model )
  {
    for( typename Model<dim>::boundaryConstIterator bit = model.BoundariesBegin(); bit != model.BoundariesEnd(); ++bit )
      if( std::binary_search( bit->second.NodesBegin(), bit->second.NodesEnd(), &node ) )
        return true;
    return false;
  }

template<size_t dim>
bool BoundaryConnector<dim>::SplitBoundaryNode( Node<dim>& node, Model<dim>& model )
  {
    for( typename Model<dim>::splitBoundaryConstIterator bit = model.SplitBoundariesBegin(); bit != model.SplitBoundariesEnd(); ++bit )
      {
        if( std::binary_search( bit->second.NodesBegin(), bit->second.NodesEnd(), &node ) )
          return true;
      }
    return false;
  }

/// removes elements: relies on unique element numbering
template<size_t dim>
void BoundaryConnector<dim>::RemoveElementNeighborConnectivity( Element<dim>& innerElement, Element<dim>& outerElement )
  {
    Element<dim>* const nullPointer(NULL);
    for( size_t i = 0; i < innerElement.Neighbors(); ++i )
      if( innerElement.Neighbor(i) )
        if( innerElement.Neighbor(i) == &outerElement )
          innerElement.Assign( i, nullPointer );
    for( size_t i = 0; i < outerElement.Neighbors(); ++i )
      if( outerElement.Neighbor(i) )
        if( outerElement.Neighbor(i) == &innerElement )
          outerElement.Assign( i, nullPointer );
  }



// explicits
template class BoundaryConnector<3U>;
template class BoundaryConnector<2U>;
template class BoundaryConnector<1U>;


} // csmp
