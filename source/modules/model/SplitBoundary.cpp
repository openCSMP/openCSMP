#include "SplitBoundary.h"
#include "BoundaryConnector.h"

#include "Node.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"

#include "Model.h"
#include "Region.h"
#include "Boundary.h"

#include "FiniteVolumeStencilManager.h"

#include "ErrorHandler.h"
#include "CSMP_highLevelUtilities.h"
#include "PL_Utilities.h"
#include "FEM_Data.h"

#include "Visitor.h"

//#define SPLITBOUNDARY_DEBUG

using namespace std;


namespace csmp {



template<size_t dim>
SplitBoundary<dim>::SplitBoundary( std::string splitboundaryname, const PropertyDatabase<dim>& pref )
 : ModelSubDomain<dim,InterFace>(splitboundaryname,pref)
 {
    this->ResizePropertyStorage( this->pref_.LocalVariablesAt(Placement()) );
 }


template<size_t dim>
SplitBoundary<dim>::SplitBoundary( const SplitBoundary& ed )
 : ModelSubDomain<dim,InterFace>(ed)
 {
 }


template<size_t dim>
SplitBoundary<dim>::SplitBoundary( SplitBoundary&& ed )
 : ModelSubDomain<dim,InterFace>(ed)
 {
 }


template<size_t dim>
SplitBoundary<dim>& SplitBoundary<dim>::operator=( const SplitBoundary<dim>& ed )
 {
    if ( &ed != this ) {
         *this = ed;
      }
    return *this;
 }


/** "All InterFaces" re-constructor of boundary from all faces in the model
*/
template<size_t dim>
SplitBoundary<dim>::SplitBoundary( const PropertyDatabase<dim>& pref,
                                   MeshManager<dim>& mesh,
                                   const SubDomainInfo& info )
 : ModelSubDomain<dim,InterFace>(info.name,pref)
 {
    // building the interface vector
    // -----------------------------
    this->elmt_vec_.reserve( info.interior_elmts.size() + info.perimeter_elmts.size() );
    // assigning pointers to the interior interfaces
    for ( auto it=info.interior_elmts.begin(); it!=info.interior_elmts.end(); ++it )
      this->elmt_vec_.push_back( &(*next(mesh.InterFacesBegin(),(*it))) );
   
    // assigning pointers to the perimeter interfaces
    for ( auto it=info.perimeter_elmts.begin(); it!=info.perimeter_elmts.end(); ++it )
      this->elmt_vec_.push_back( &(*next(mesh.InterFacesBegin(),(*it))) );
   
    // sorting the subvectors for future searching
    const auto perimeterInterFacesBegin( next(this->elmt_vec_.begin(), info.interior_elmts.size()) );
    sort( this->elmt_vec_.begin(), perimeterInterFacesBegin );
    sort( perimeterInterFacesBegin, this->elmt_vec_.end() );

    // building the vector of vectors of those faces of the interfaces that lie on the subdomain perimeter
    // ---------------------------------------------------------------------------------------------------
    this->bd_face_vec_.reserve( info.perimeter_faces.size() );
    for ( auto it=info.perimeter_faces.begin(); it!=info.perimeter_faces.end(); ++it ) {
          const size_t perimeter_faces((*it).size());
          std::vector<ONE_BYTE_NUMBER> face_vec;
          face_vec.reserve(perimeter_faces);
          for ( auto fit=(*it).begin(); fit!=(*it).end(); ++fit )
            face_vec.push_back( static_cast<ONE_BYTE_NUMBER>( (*fit) ) );
      }

    // building the node vector
    // ------------------------
    this->node_vec_.reserve( info.interior_nodes.size() + info.perimeter_nodes.size() );
    // assigning pointers to the interior nodes
    for ( auto it=info.interior_nodes.begin(); it!=info.interior_nodes.end(); ++it )
      this->node_vec_.push_back( &(*next( mesh.NodesBegin(),(*it))) );
   
    // assigning pointers to the perimeter nodes
    for ( auto it=info.perimeter_nodes.begin(); it!=info.perimeter_nodes.end(); ++it )
      this->node_vec_.push_back( &(*next( mesh.NodesBegin(),(*it))) );
   
    // sorting the subvectors for future searching
    const auto perimeterNodesBegin( next(this->node_vec_.begin(), info.interior_nodes.size()) );
    sort( this->node_vec_.begin(), perimeterNodesBegin );
    sort( perimeterNodesBegin, this->node_vec_.end() );
   
    // allocating the storage for subdomain properties
    // -----------------------------------------------
    this->ResizePropertyStorage( pref.LocalVariablesAt(SPLIT_BOUNDARY) );
 }
 
 
 
/** complete construction of boundary using a master boundary that must contain all faces used for the construction
*/
template<size_t dim>
SplitBoundary<dim>::SplitBoundary( const PropertyDatabase<dim>& dbase,
                                   const SplitBoundary<dim>&  all_interfaces,
                                   const SubDomainInfo& info )
 : ModelSubDomain<dim,InterFace>(dbase,all_interfaces,info)
 {
    // property storage is handled by model subdomain
 }



/// Does not delet interfaces, Delete() has to be called for this
template<size_t dim>
SplitBoundary<dim>::~SplitBoundary()
 {
 }



template<size_t dim>
IntegrationPointVariables  SplitBoundary<dim>::InterFaceIntegrationPointVariables() const
  { return this->pref_.IntegrationPointVariablesAt(INTER_FACE); }

template<size_t dim>
LocalVariables  SplitBoundary<dim>::InterFaceVariables() const
  { return this->pref_.LocalVariablesAt(INTER_FACE); }





// TODO: check what the impact of this method is
template<size_t dim>
void  SplitBoundary<dim>::DetachElementsFromNeighbors() const
{
  for ( typename vector<csmp::InterFace<dim>*>::const_iterator it( this->ElementsBegin() ); it != this->ElementsEnd(); ++it )
    (*it)->Detach();
}

// LOCAL VARIABLE STORAGE INTERFACE

template<size_t dim>
bool SplitBoundary<dim>::ValidVariable( const char* variableName ) const
  {
    const PLACEMENT p( this->pref_.Placement(variableName) );
    if( p == NODE || p == INTER_FACE || p == SPLIT_BOUNDARY )
      return true;
    return false;
  }

// VISITORS INTERFACE

/// visitation of a split boundary
template<size_t dim>
void SplitBoundary<dim>::Accept( Visitor<dim>& v )
  {
    if ( v.ApplicationLevel() == MODEL || v.ApplicationLevel() == SPLIT_BOUNDARY )
      v.Visit(this);

    switch( v.ApplicationTarget() ) {
        case MODEL:
          throw csmp::Exception( ERROR, "Region<dim>::Accept",
                                "ApplicationTarget MODEL; Visitor should have never arrived at this SplitBoundary");
          break;
        case SPLIT_BOUNDARY:
          return;
          // element, face and interface are treated the same
        case INTER_FACE:
          for ( typename vector<InterFace<dim>*>::iterator
               it=this->ElementsBegin(); it!=this->ElementsEnd(); it++ )
            (*it)->Accept( v );
          return;
        case NODE:
          for ( typename vector<csmp::Node<dim>*>::iterator
               nd_it=this->NodesBegin(); nd_it!=this->NodesEnd(); nd_it++ )
            (*nd_it)->Accept( v );
          return;
        default:
          throw csmp::Exception( ERROR, "SplitBoundary<dim>::Accept",
                                "ApplicationTarget was not resolved; nothing was done");
    }
  } // end Accept


/// Returns the position of and adjacent region relative to the boundary. Relies on element Idx
template<size_t dim>
INTERFACE_SIDE SplitBoundary<dim>::RegionLocation( const Region<dim>& region )
  {
    assert( !this->elmt_vec_.empty() );
    const typename vector<Element<dim>*>::const_iterator regionElementsEnd( region.ElementsEnd() );
    const typename vector<InterFace<dim>*>::const_iterator sbElementsEnd( this->ElementsEnd() );
    for( typename vector<InterFace<dim>*>::const_iterator ifit( this->ElementsBegin() ); ifit != sbElementsEnd; ++ifit )
      {
        const size_t innerParentIdx( (*ifit)->Parent( INSIDE )->Idx() ),
                     outerParentIdx( (*ifit)->Parent( OUTSIDE )->Idx() );
        for( typename vector<Element<dim>*>::const_iterator eit( region.ElementsBegin() ); eit != regionElementsEnd; ++eit )
          {
            if( (*eit)->Idx() == innerParentIdx )
              return INSIDE;
            if( (*eit)->Idx() == outerParentIdx )
              return OUTSIDE;
          } // region elements
      } // split boundary interfaces
    throw csmp::Exception( ERROR, "SplitBoundary<dim>::RegionLocation", "Region seems not to be adjacent to split boundary!" );
    // shouldn't get here
    return OUTSIDE;
  }




// -----------------------------------------------
// Binary input/output
// -----------------------------------------------



/**
 @fn  void Boundary<dim>::Out( FILE* fp ) const

 @brief Outs the boundary to binary fp.
 @attention Uses current index numbering.

 We use the following order her:
 0. boundary domain data
 1. face count
 2. flag
 3. face fem type
 4. inner, outer parent element id of face
    (in the case of a null parent, i.e. at model perimeter, we also store the face id of the inner element)
 5. fem face data
    -scalar data count, property name & scalar FEM_Data
    -vector data count, property name & vector FEM_Data
    -tensor data count, property name & tensor FEM_Data
    -array data count, property name & array FEM_Data
    -flagged array data count, property name & flagged array FEM_Data

 @author  P. Lang
 @date  9/29/2012

 @tparam  dim Dimension
 @param [in,out]  fp  If non-null, the file pointer to the binary output file

 @todo (2-C) Return values not used properly.
 */
template<size_t dim>
bool SplitBoundary<dim>::Out( FILE* fp ) const
  {
  // split-boundary variables
  domainVariablesOut( fp, *this, this->pref_ ); /// @todo (3-D) Use FEM_Data instead?

  // number of interfaces
  size_t bytes( sizeof(size_t) );
  const size_t interfaceCount( this->Elements() );
  fwrite( (void*) &interfaceCount, bytes, 1, fp );

  // fem type of interfaces
  CSMP_FEM_TYPE interfaceType;
  bytes = sizeof(CSMP_FEM_TYPE);
  InterFace<dim>* interface(NULL);
  for ( size_t f(0); f < interfaceCount; ++f )
    {
        interface = this->elmt_vec_[f];
        if( !interface ) // we check in this loop only for nullptrs
            return false;
        interfaceType = interface->FE()->ElementType();
        fwrite( (void*) &interfaceType, bytes, 1, fp );
    }

  // interface parents
  size_t idx(NULL_IDX);
  bytes = sizeof(size_t);
  
  // higher-dimensional elements will be present on the inside and the outside of the interface
  // because interfaces can only be created from internal model boundaries
  for ( size_t f(0); f < interfaceCount; ++f )
    {
        interface = this->elmt_vec_[f];
      
        assert( interface->Parent(INSIDE) != nullptr );
        assert( interface->Parent(OUTSIDE) != nullptr );

        idx = interface->Parent(INSIDE)->Idx();
        fwrite( &idx, bytes, 1, fp );
        const size_t innerParentFaceId( interface->ParentFaceID( INSIDE ) );
        fwrite( &innerParentFaceId, bytes, 1, fp );

        idx = interface->Parent(OUTSIDE)->Idx();
        fwrite( &idx, bytes, 1, fp );
        const size_t outerParentFaceId( interface->ParentFaceID( OUTSIDE ) );
        fwrite( &outerParentFaceId, bytes, 1, fp );

        if( interface->BaseElement() )
        {
            idx = interface->BaseElement()->Idx();
            fwrite( &idx, bytes, 1, fp );
        }
        else
        {
            idx = NULL_IDX;
            fwrite( &idx, bytes, 1, fp );
        }
        const size_t nodesCount( interface->Nodes() );
        fwrite( &nodesCount, bytes, 1, fp );
        for( size_t fn(0); fn < nodesCount; ++fn )
        {
            const size_t localInnerNodeIdx( interface->ParentNodeNumber( fn, INSIDE ) );
            fwrite( &localInnerNodeIdx, bytes, 1, fp );

            const size_t localOuterNodeIdx( interface->ParentNodeNumber( fn, OUTSIDE ) );
            fwrite( &localOuterNodeIdx, bytes, 1, fp );
        }
    }

  // interface variable count: scalar, vector, tensor, array, flagged array
  this->RenumberElements();
  /*
  Out<ScalarVariable>( fp, FACE, SCALAR );
  Out<VectorVariable<dim> >( fp, FACE, VECTOR );
  Out<TensorVariable<dim> >( fp, FACE, TENSOR );
  Out<ArrayVariable>( fp, FACE, ARRAY );
  Out<FlaggedArrayVariable>( fp, FACE, FLAGGEDARRAY );
  */
  return true;
  }




/**
 @fn  bool Boundary<dim>::In( const FiniteElementManager& femManager, FILE* fp ) const

 @brief INS Boundary from the given binary file pointer.

 see Out( FILE* fp )

 @author  P. Lang
 @date  9/29/2012

 @param [in,out]  fp  If non-null, the fp.

 @return  true if it succeeds, false if it fails.

 @todo (2-C) Return values not used properly.
 */
template<size_t dim>
bool SplitBoundary<dim>::In( MeshManager<dim>& meshManager,
                             const FiniteElementManager& femManager,
                             const Region<dim>& modelRegion,
                             FILE* fp )
  {
    // splitboundary variables
    domainVariablesIn( fp, *this, this->pref_ );
    // number of faces
    size_t bytes( sizeof(size_t) );
    size_t interfaceCount(0);
    fread( (void*) &interfaceCount, bytes, 1, fp );

    // fem type of interfaces
    bytes = sizeof(CSMP_FEM_TYPE);
    vector<CSMP_FEM_TYPE> interfaceTypes(interfaceCount);
    for ( size_t f(0); f < interfaceCount; ++f )
    fread( (void*) &interfaceTypes[f], bytes, 1, fp );

    // interface parents
    bytes = sizeof(size_t);
    std::vector<std::vector<size_t> >                     interfaceParents( interfaceCount );
    std::vector<std::vector<std::pair<size_t, size_t> > > interfaceParentNodes( interfaceCount );
    size_t interfaceNodes(1);
    size_t localInnerNodeIdx(1);
    size_t localOuterNodeIdx(1);
    for ( size_t f(0); f < interfaceCount; ++f )
    {
        interfaceParents[f].resize( 5, NULL_IDX );

        fread( (void*) &interfaceParents[f][0], bytes, 1, fp );
        fread( (void*) &interfaceParents[f][1], bytes, 1, fp );

        fread( (void*) &interfaceParents[f][2], bytes, 1, fp );
        fread( (void*) &interfaceParents[f][3], bytes, 1, fp );

        fread( (void*) &interfaceParents[f][4], bytes, 1, fp );

        fread( (void*) &interfaceNodes, bytes, 1, fp );
        interfaceParentNodes[f].resize( interfaceNodes, make_pair( NULL_IDX, NULL_IDX ) );
        for( size_t fn(0); fn < interfaceNodes; ++fn )
        {
            fread( (void*) &localInnerNodeIdx, bytes, 1, fp );
            fread( (void*) &localOuterNodeIdx, bytes, 1, fp );
            interfaceParentNodes[f][fn] = make_pair( localInnerNodeIdx, localOuterNodeIdx );
        }
    }

    std::map<size_t,Element<dim>*> elementIdPtr;
    const typename vector<Element<dim>*>::const_iterator eitEnd( modelRegion.ElementsEnd() );
    for( typename vector<Element<dim>*>::const_iterator eit( modelRegion.ElementsBegin() ); eit != eitEnd; ++eit  )
    elementIdPtr[ (*eit)->Idx() ] = (*eit);

    // creating splitboundary
    CreateFrom( meshManager, femManager, elementIdPtr, interfaceTypes, interfaceParents, interfaceParentNodes );

    // fem data
    // interface variable count: scalar, vector, tensor, array, flagged array
    /*
    this->UpdateMemberIndexes();
    In<ScalarVariable>( fp, FACE, SCALAR );
    In<VectorVariable<dim> >( fp, FACE, VECTOR );
    In<TensorVariable<dim> >( fp, FACE, TENSOR );
    In<ArrayVariable>( fp, FACE, ARRAY );
    In<FlaggedArrayVariable>( fp, FACE, FLAGGEDARRAY );
    */

    return true;
  }





/**
 @brief Creates splitboundary from input vectors (used in binary IO)

 @author  P. Lang
 @author  R. Manasipov
 @date  01/08/2014

 @param  Dimension
 @param femManager  Manager for finite elements.
 @param interfaceTypes   List of types of the interfaces.
 @param interfaceParents The interface parents.
 @param interfaceParentNodes Local indexes of interface nodes in parent elements

 @return  Failed/Succeeded.
 */
template<size_t dim>
bool SplitBoundary<dim>::CreateFrom( MeshManager<dim>&                                  meshManager,
                                     const FiniteElementManager&                        femManager,
                                     const std::map<size_t,csmp::Element<dim>*>&        elementIdPtr,
                                     const std::vector<CSMP_FEM_TYPE>&                  interfaceTypes,
                                     const std::vector<std::vector<size_t> >&           interfaceParents,
                                     const std::vector<std::vector<std::pair<size_t,size_t> > >&  interfaceParentNodes )
  {

      // LVS
      const LocalVariables lvsInterFaces( InterFaceVariables() );
      const IntegrationPointVariables lvsIntegrationPoints( InterFaceIntegrationPointVariables() );

      // running pointer to fem type of new interfaces as acquired from region element
      FiniteElement* femPtr(NULL);
      Element<dim>* innerParentPtr( NULL );
      Element<dim>* outerParentPtr( NULL );
      Element<dim>* baseElementPtr ( NULL );

      // preparing container for a max of total region element count
      this->elmt_vec_.reserve( interfaceTypes.size() );

      // looping over regions elements, assuring that it's an eligible interface type, creating new interface with variable storage,
      // establishing connectivity and inserting into boundary element container
      for( size_t f(0); f < interfaceTypes.size(); ++f )
        {
            femPtr = femManager.E( interfaceTypes[f] );

            // which is used to create the new interface using the variables prepared above (NULL is FV Stencil)
            //interfacePtr = new InterFace<dim>( femPtr, NULL, lvsInterFaces, lvsIntegrationPoints );
            InterFace<dim>* interfaceObj = meshManager.PushBackIfUnique( InterFace<dim>( femPtr, NULL, lvsInterFaces, lvsIntegrationPoints ) );

            // inner parent element
            if( interfaceParents[f][0] != NULL_IDX )
            {
                innerParentPtr = elementIdPtr.find( interfaceParents[f][0] )->second; // dangerous
                interfaceObj->Assign( innerParentPtr, interfaceParents[f][1], INSIDE );
            }
            else
                throw csmp::Exception( ERROR, "SplitBoundary<dim>::CreateFrom", "Inner interface parent cannot be NULL" );

            // outer parent element
            if( interfaceParents[f][2] != NULL_IDX )
            {
                outerParentPtr = elementIdPtr.find( interfaceParents[f][2] )->second; // dangerous
                interfaceObj->Assign( outerParentPtr, interfaceParents[f][3], OUTSIDE );
            }
            else
                throw csmp::Exception( ERROR, "SplitBoundary<dim>::CreateFrom", "Outer interface parent cannot be NULL" );

            // assign base element if it exist
            if( interfaceParents[f][4] != NULL_IDX )
            {
                baseElementPtr = elementIdPtr.find( interfaceParents[f][4] )->second; // dangerous
                interfaceObj->Assign( baseElementPtr );
            }

            // assign corresponding parent nodes
            const size_t interfaceNodes( interfaceParentNodes[f].size() );
            for( size_t fn(0); fn < interfaceNodes; ++fn )
                interfaceObj->Assign( fn, interfaceParentNodes[f][fn] );

            // push back into face container
            this->elmt_vec_.emplace_back( interfaceObj );

        } // region elements

      // free
      vector<InterFace<dim>*>( this->elmt_vec_ ).swap( this->elmt_vec_ );

      // establishing splitboundary essentials
      Initialize( true /* update neighbor connectivity*/, false /* do not update indexes */ );

      // done
      return true;

  }











// ------------------------------------------------------------------
// Building blocks
// -------------------------------------------------------------------


/**
Dispatched initialize method to establish node vector, perimeter entities,
and entity sorting
*/
template<size_t dim>
void SplitBoundary<dim>::Initialize( bool updateNeighborConnectivity, bool updateIndexes )
  {
    // establishing splitboundary node container
    this->CreateNodePointerVector();

    // identify entities on splitboundary perimeter
    if( updateNeighborConnectivity )
      this->EstablishNeighborConnectivity();

    this->IdentifyPerimeter();

    // initialize indices
    if(updateIndexes)
        this->UpdateMemberIndexes();
  }


template<size_t dim>
void SplitBoundary<dim>::ConnectFiniteVolumeStencils( const FiniteVolumeStencilManager<dim>& fvm_mgr )
 {
    for ( typename vector<InterFace<dim>*>::iterator it=this->elmt_vec_.begin(); it!=this->elmt_vec_.end(); it++ )
        (*it)->Assign( fvm_mgr.Stencil( (*it)->FE_Type() ) );

 }


template<size_t dim>
void SplitBoundary<dim>::CreateNodePointerVector()
 {
    assert( !this->elmt_vec_.empty() );

    if ( !this->node_vec_.empty() )
        this->node_vec_.clear();

    // creating the node index vector
    set<csmp::Node<dim>*>  nodes_set;
    for ( typename vector<InterFace<dim>*>::const_iterator it = this->elmt_vec_.begin(); it!=this->elmt_vec_.end(); it++ )
        for ( typename vector<Node<dim>*>::size_type i=0U; i<(*it)->Nodes(); i++ )
        {
            nodes_set.insert( (*it)->N(i, INSIDE ) );
            nodes_set.insert( (*it)->N(i, OUTSIDE ) );
            if ( (*it)->hasBase() )
                nodes_set.insert( (*it)->N(i, MIDDLE ) );
        }

    this->node_vec_.assign( nodes_set.begin(), nodes_set.end() );
 }




/**
    Detects of how many spatial dimensions interface types are contained in model.
    It returns a pair: first value gives number of different spatial dimensions contained,
    second value returns the highest spatial dimension contained.

    @author SKM 1/11/2013

*/
template<size_t dim>
pair<int32,int32>  SplitBoundary<dim>::InterFaceSpatialDimensions() const
 {
    return this->SpatialDimensions();

 } // end ElementSpatialDimensions
















/// Creates a split boundary from a boundary. Requires unique indices.
template<size_t dim>
bool  SplitBoundary<dim>::CreateFrom( Model<dim>& model,
                                      Boundary<dim>& boundary )
  {
    //LVS
    const LocalVariables lvsInterFace( InterFaceVariables() );
    const IntegrationPointVariables lvsIntegrationPoint( InterFaceIntegrationPointVariables() );

    // allocating SubDomain element container
    this->elmt_vec_.clear();
    this->elmt_vec_.reserve( boundary.Elements() );

    // check whether all faces have inner AND outer parents.
    FiniteElement*  femPtr( NULL );
    Element<dim>*   innerElement( NULL );
    Element<dim>*   outerElement( NULL );

    const typename vector<Face<dim>*>::const_iterator facesEnd( boundary.ElementsEnd() );
    for( typename vector<Face<dim>*>::const_iterator fit( boundary.ElementsBegin() ); fit != facesEnd; ++fit )
      {
        // getting parents
        innerElement = (*fit)->Parent(INSIDE);
        outerElement = (*fit)->Parent(OUTSIDE);

        // the InterFace that is being build from the current face
        femPtr = model.FE_Manager().E( (*fit)->FE_Type() );
        InterFace<dim>* interfaceObj = model.Mesh().PushBackIfUnique( InterFace<dim>( femPtr, NULL, lvsInterFace, lvsIntegrationPoint ) );

        interfaceObj->Assign( innerElement, outerElement );
        interfaceObj->Idx( (*fit)->Idx() );

        this->elmt_vec_.emplace_back( interfaceObj );

        BoundaryConnector<dim>::RemoveElementNeighborConnectivity( *innerElement, *outerElement );

        // resetting
        innerElement = NULL;
        outerElement = NULL;

      } // interfaces from faces

    // free excessive allocated capacity
    vector<InterFace<dim>*>( this->elmt_vec_ ).swap( this->elmt_vec_ );

    // initialize splitboundary essentials
    Initialize( true /* update neighbor connectivity*/, false /* do  not update member indexes */);

    // split nodes and assign to corresponding elements
    Split( model, boundary );

    // establishing splitboundary essentials
    Initialize( false /* update neighbor connectivity*/, false /* do not update indexes */ );
    
// TODO: model region needs to be rebuild after this operation, since we have additional nodes now?

    return true;
  }
















template<size_t dim>
void SplitBoundary<dim>::Split( Model<dim>& model,
                                Boundary<dim>& boundary )
{
    // Update Connectivity ( Replace all old nodes with duplicated ones )

#ifdef SPLITBOUNDARY_DEBUG
    clock_t start = clock();
#endif

    // Collect nodes that are going to be duplicated:

    std::set<Node<dim>*> nodesToDuplicate;

    // we don't want to find nodes shared with the base itself
    std::vector<Node<dim>*> boundaryNodeCache( boundary.NodeVector() );
    boundary.NodeVector().clear();

    // duplicate perimeter nodes of boundary only if they are located on already existing boundary or splitboundary
    for( typename vector<Face<dim>*>::const_iterator fit( boundary.PerimeterElementsBegin() ); fit != boundary.ElementsEnd(); ++fit )
      for( size_t n(0); n < (*fit)->Nodes(); ++n )
        {
          if( std::binary_search( nodesToDuplicate.begin(), nodesToDuplicate.end(), (*fit)->N(n) ) )
            continue;
          if( BoundaryConnector<dim>::BoundaryNode( *(*fit)->N(n), model ) )
            { nodesToDuplicate.insert( (*fit)->N(n) ); continue; }
          if( BoundaryConnector<dim>::SplitBoundaryNode( *(*fit)->N(n), model ) )
            { nodesToDuplicate.insert( (*fit)->N(n) ); continue; }
        }

    // restore boundary nodes
    boundary.NodeVector().swap( boundaryNodeCache );

    // duplicate all interior nodes of boundary
    for( typename vector<Node<dim>*>::const_iterator nit( boundary.NodesBegin() ); nit != boundary.PerimeterNodesBegin(); ++nit )
        nodesToDuplicate.insert( (*nit) );


    // UpdateConnectivity:
    std::set<Node<dim>*> updatedNodes;
    Node<dim>*           originalNode  ( NULL );
    size_t               np(0);

    for( typename std::vector<InterFace<dim>*>::const_iterator ifit( this->ElementsBegin() ); ifit != this->ElementsEnd(); ++ifit )
    {
      for( size_t ifn(0); ifn < (*ifit)->Nodes(); ++ifn )
      {
        originalNode  = (*ifit)->N( ifn, INSIDE );

        if( updatedNodes.find( originalNode ) != updatedNodes.end() )
            continue;
        else
            updatedNodes.insert( originalNode );

        bool duplicate( false );
        if( nodesToDuplicate.find( originalNode ) != nodesToDuplicate.end() )
            duplicate = true;

        if( duplicate )
          {
            /// TODO: WARNING: this operation inserts Node into Mesh Manager container without checking whether such Node already exist or not.
            /// Therefore this function should be called with caution and only if all previos steps approves it.
            Node<dim>* duplicatedNode = model.Mesh().PushBack( Node<dim>( *originalNode ) );
            updatedNodes.insert( duplicatedNode );

            // referencing index and checking if node connections already updated
            const size_t originalNodeIndex( originalNode->Idx() );

            // replacing all element nodes
            const Region<dim>&  mref(model.Region("Model"));
            for( typename vector<Element<dim>*>::const_iterator eit( mref.ElementsBegin() ); eit != mref.ElementsEnd(); ++eit )
              for( size_t en(0); en < (*eit)->Nodes(); ++en )
                if( (*eit)->N(en)->Idx() == originalNodeIndex )
                  if( BoundaryConnector<dim>::OnOutside( *(*ifit), *(*eit), *duplicatedNode ) )
                    (*eit)->Assign( en, duplicatedNode );

            // rm all inner elements from duplicated node parent list
            np = 0;
            while( np < duplicatedNode->Parents() )
              {
                if( !BoundaryConnector<dim>::OnOutside( *(*ifit), *duplicatedNode->Parent(np), *duplicatedNode ) )
                  { duplicatedNode->Unassign( duplicatedNode->Parent(np) ); np = 0; }
                else
                  ++np;
              }

            // rm all outer elements from original counterpart node parent list
            np = 0;
            while( np < originalNode->Parents() )
              {
                if( BoundaryConnector<dim>::OnOutside( *(*ifit), *originalNode->Parent(np), *originalNode ) )
                  { originalNode->Unassign( originalNode->Parent(np) ); np = 0; }
                else
                  ++np;
              }

          } // original-duplicated node combo loop

        }// check if node need's to be duplicated
    }

  #ifdef SPLITBOUNDARY_DEBUG

      clock_t end = clock();

      // ------------------------------------------------------------------------------------------
      // Output Elapsed Time
      // ------------------------------------------------------------------------------------------

      unsigned long millisec ( (end - start) * 1000 / CLOCKS_PER_SEC);

      cerr<<"\n"<<"Splitboundary<dim>::Split: "<<name_.first<<"_"<<name_.second<<endl;
      cerr<<"\nElapsed Time = "<<millisec<<" ms ("<<(double)(millisec)/1000.<<" sec; "<<(double)(millisec)/60000.<<" min; "<<(double)(millisec)/3600000.<<" hours)"<<endl;
      cerr.flush();

  #endif

    return;
}








// CALCULATIONS



template<size_t dim>
double64  SplitBoundary<dim>::Perimeter() const
 {
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );
     if ( dim != 3U )
       csmp_error.notice( ERROR, "SplitBoundary<dim>::Perimeter",
                                 "result would not be meaningful" );

     double64        perimeter_length(0.);
     vector<size_t>  fnids;
     size_t          n(0U);
     for ( typename vector<InterFace<dim>*>::const_iterator
           it=this->PerimeterElementsBegin(); it!=this->ElementsEnd(); it++, n++ )
       for ( size_t i=0U; i<this->PerimeterFaces(n); i++ ) {
            (*it)->FE()->NodesOfFace( this->PerimeterFace(n,i), fnids );
            perimeter_length += ((*it)->N(fnids[1])->Coordinate() -
                                 (*it)->N(fnids[0])->Coordinate()).Length();
         }

     return perimeter_length;
 }



template<size_t dim>
double64  SplitBoundary<dim>::Area() const
 {
     double64  integrated_area(0.);
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );

     if ( dim == 3U ) {
         for ( typename vector<InterFace<dim>*>::const_iterator
               it=this->elmt_vec_.begin(); it!=this->elmt_vec_.end(); it++ )
           if ( (*it)->FE()->IsSurfaceElement() )
             integrated_area += (*it)->Volume();
       }
     else if ( dim == 2U ) {
         for ( typename vector<InterFace<dim>*>::const_iterator
               it=this->elmt_vec_.begin(); it!=this->elmt_vec_.end(); it++ )
           if ( (*it)->FE()->IsLineElement() )
             integrated_area += (*it)->Volume();
       }
     else
     csmp_error.notice( ERROR, "SplitBoundary<dim>::Area", "not defined in 1D");

     return integrated_area;
 }


template<size_t dim>
double64 SplitBoundary<dim>::SurfaceIntegral( const PropertyDatabase<dim>& p, const char* property ) const
 {
     csmp::Index prop_key = p.StorageKey(property);

     if ( prop_key.place == ELEMENT_INTEGRATION_POINT or prop_key.place == REGION ) {
          throw csmp::Exception( ERROR, "SplitBoundary<dim>::SurfaceIntegral",
                            property, "placed on IntegrationPoint or Region cannot be assigned on boundary");
          return std::numeric_limits<double64>::quiet_NaN();
       }
      if ( prop_key.type == TENSOR ) {
          throw csmp::Exception( ERROR, "SplitBoundary<dim>::SurfaceIntegral",
                            property, "is a tensor property; this method does not know how to integrate it");
          return std::numeric_limits<double64>::quiet_NaN();
       }

     double64  property_integral(0.);

     // 1. if the property is a scalar
     if ( prop_key.type == SCALAR ) {
           if ( prop_key.place == ELEMENT ) {
                throw csmp::Exception( ERROR, "SplitBoundary<dim>::SurfaceIntegral",
                                  property, "is an Element property; this method does not know how to integrate it");
             }
           else if ( prop_key.place == FACE or prop_key.place == INTER_FACE ) {
                for ( typename vector<InterFace<dim>*>::const_iterator
                      it=this->elmt_vec_.begin(); it!=this->elmt_vec_.end(); it++ )
                  property_integral += (*it)->Volume() * (*it)->Read( prop_key );
             }
           else if ( prop_key.place == NODE ) { // for nodes on first side of interface
                ScalarVariable  sc;
                for ( typename vector<InterFace<dim>*>::const_iterator
                      it=this->elmt_vec_.begin(); it!=this->elmt_vec_.end(); it++ ) {
                     (*it)->PropertyValueAtBaryCenter( prop_key, sc );
                     property_integral += (*it)->Volume() * sc.Value();
                  }
             }
           else {
                throw csmp::Exception( FATAL_ERROR, "SplitBoundary<dim>::SurfaceIntegral",
                                                    "Property placement not recognized");
             }
       }
     // 2. if the property is a vector
     if ( prop_key.type == VECTOR ) {
           // the average of the values projected onto the normal are being used.
           if ( prop_key.place == FACE or prop_key.place == INTER_FACE ) {
                VectorVariable<dim>  unrml, vc;
                for ( typename vector<InterFace<dim>*>::const_iterator
                      it=this->elmt_vec_.begin(); it!=this->elmt_vec_.end(); it++ ) {
                     (*it)->UnitNormal( unrml );
                     (*it)->Read( prop_key, vc );
                     property_integral += unrml & vc;
                  }
             }
           else if ( prop_key.place == NODE ) { // for nodes on first side of interface
                VectorVariable<dim>  unrml, vc;
                for ( typename vector<InterFace<dim>*>::const_iterator
                      it=this->elmt_vec_.begin(); it!=this->elmt_vec_.end(); it++ ) {
                     (*it)->UnitNormal( unrml );
                     (*it)->PropertyValueAtBaryCenter( prop_key, vc );
                     property_integral += unrml & vc;
                  }
             }
           else {
                throw csmp::Exception( FATAL_ERROR, "SplitBoundary<dim>::SurfaceIntegral",
                                                    "Property placement not recognized");
             }
       }

    return property_integral;

 } // end SurfaceIntegral





/**
     Property assignment to nodes on either side of the interface.
*/
template<size_t dim>
template<class Var>
void SplitBoundary<dim>::InputNodePropertyValue( const char* input_prop, const Var& new_value, SUBDOMAIN_PART part, INTERFACE_SIDE innerOuter )
  {
    Index ipKey( this->pref_.StorageKey(input_prop) );

    if( ipKey.place != NODE )
      throw csmp::Exception(  ERROR, "SplitBoundary<dim>::InputNodePropertyValue:", "This method applies to node properties only!" );

    if ( part == COMPLETE ) {
        const typename vector<InterFace<dim>*>::const_iterator ifEnd( this->ElementsEnd() );
        for( typename vector<InterFace<dim>*>::const_iterator ifit( this->ElementsBegin() ); ifit != ifEnd; ++ifit )
          for( size_t n(0); n < (*ifit)->Nodes(); ++n )
            (*ifit)->N( n, innerOuter )->Store( ipKey, new_value );
      }
    else if ( part == INTERIOR ) {
      const typename vector<InterFace<dim>*>::const_iterator ifEnd( this->PerimeterElementsBegin() );
      for( typename vector<InterFace<dim>*>::const_iterator ifit( this->ElementsBegin() ); ifit != ifEnd; ++ifit )
        for( size_t n(0); n < (*ifit)->Nodes(); ++n )
          (*ifit)->N( n, innerOuter )->Store( ipKey, new_value );
      }
    else if ( part == PERIMETER ) {
      const typename vector<InterFace<dim>*>::const_iterator ifEnd( this->ElementsEnd() );
      for( typename vector<InterFace<dim>*>::const_iterator ifit( this->PerimeterElementsBegin() ); ifit != ifEnd; ++ifit )
        for( size_t n(0); n < (*ifit)->Nodes(); ++n )
          (*ifit)->N( n, innerOuter )->Store( ipKey, new_value );
      }
    else
    cerr << "\nSplitBoundary<dim>::InputNodePropertyValue: subdomain part not recognized; nothing was done.\n";

  } // InputPropertyValue

template void SplitBoundary<1>::InputNodePropertyValue( const char*, const ScalarVariable&, SUBDOMAIN_PART, INTERFACE_SIDE );
template void SplitBoundary<2>::InputNodePropertyValue( const char*, const ScalarVariable&, SUBDOMAIN_PART, INTERFACE_SIDE );
template void SplitBoundary<3>::InputNodePropertyValue( const char*, const ScalarVariable&, SUBDOMAIN_PART, INTERFACE_SIDE );

template void SplitBoundary<1>::InputNodePropertyValue( const char*, const VectorVariable<1>&, SUBDOMAIN_PART, INTERFACE_SIDE );
template void SplitBoundary<2>::InputNodePropertyValue( const char*, const VectorVariable<2>&, SUBDOMAIN_PART, INTERFACE_SIDE );
template void SplitBoundary<3>::InputNodePropertyValue( const char*, const VectorVariable<3>&, SUBDOMAIN_PART, INTERFACE_SIDE );

template void SplitBoundary<1>::InputNodePropertyValue( const char*, const TensorVariable<1>&, SUBDOMAIN_PART, INTERFACE_SIDE );
template void SplitBoundary<2>::InputNodePropertyValue( const char*, const TensorVariable<2>&, SUBDOMAIN_PART, INTERFACE_SIDE );
template void SplitBoundary<3>::InputNodePropertyValue( const char*, const TensorVariable<3>&, SUBDOMAIN_PART, INTERFACE_SIDE );







// SCREEN OUTPUT

template<size_t dim>
void SplitBoundary<dim>::Out() const
 {
    cout <<"\nSplitBoundary<dim>::Out(): ";
    cout <<" member elements: interior="<< this->InteriorElements();
    cout <<", boundary="<< this->elmt_vec_.size()-this->InteriorElements() <<": "<< endl;

    for ( typename vector<InterFace<dim>*>::const_iterator
          it=this->elmt_vec_.begin(); it!=this->elmt_vec_.end(); it++ ) {
         if ( (*it) == NULL )
           throw csmp::Exception( ERROR, "SplitBoundary<dim>::Out",
                                 "member element pointer not initialised");
      }

    cout <<"\n\n edge elements and their edge faces (current local numbering): "<< endl;
    vector<vector<ONE_BYTE_NUMBER> >::const_iterator  bit(this->bd_face_vec_.begin());
    for ( size_t i=this->InteriorElements(); i<this->elmt_vec_.size(); i++, bit++ ) {
         cout <<"\nelement "<< i <<": edge face numbers: ";
         for ( vector<ONE_BYTE_NUMBER>::const_iterator
               ft=(*bit).begin(); ft!=(*bit).end(); ft++ ) cout << (*ft) <<" ";
      }

    cout <<"\n\n edge nodes: "<< this->node_vec_.size() - this->first_bd_node_ <<" (current local numbering):"<< endl;
    for ( size_t i=this->first_bd_node_; i<this->node_vec_.size(); i++ ) {
         if ( this->node_vec_[i] == NULL )
           throw csmp::Exception( ERROR, "SplitBoundary<dim>::Out", "member node pointer not initialised.");
         else cout << this->node_vec_[i]->Idx() <<" ";
      }

    cout << endl;
 }












template class SplitBoundary<1>;
template class SplitBoundary<2>;
template class SplitBoundary<3>;

} // end csmp

