#include "PolygonCell.h"

#include "MeshingTools.h"
#include "PolygonGrid.h"
#include "CSMP_ElementSpecifications.h"

#include "ErrorHandler.h"

namespace csmp {

// POLYGON

PolygonCell::PolygonCell()
:grid_(nullptr)
{
}

/// custom constructor
PolygonCell::PolygonCell( PolygonGridManager& pgm )
: grid_( &pgm )
{
   assert( grid_ != nullptr );
}



// TODO: what happens to pgm pointer
PolygonCell::PolygonCell( const PolygonCell&  poly )
: elements_   ( poly.elements_ ),
  faces_      ( poly.faces_ ),
  custom_face_order_( poly.custom_face_order_ ),
  num_faces_  ( poly.num_faces_ ),
  face_nodes_ ( poly.face_nodes_ ),
  face_nodes_in_custom_order_( poly.face_nodes_in_custom_order_ ),
  num_nodes_  ( poly.num_nodes_ ),
  extra_nodes_( poly.extra_nodes_ ),
  nodes_      ( poly.nodes_ ),
  nodes_in_custom_order_( poly.nodes_in_custom_order_ ),
  custom_node_order_( poly.custom_node_order_ ),
  // do we want this to be the same pointer?
  grid_       ( poly.grid_ )
{
   assert( grid_ != nullptr );
}


PolygonCell& PolygonCell::operator=( const PolygonCell&  cell )
{
    if ( &cell != this )
    {
        grid_  = ( cell.grid_ == nullptr )  ? nullptr : cell.grid_;

        elements_           = cell.elements_;
        num_faces_          = cell.num_faces_;
        faces_              = cell.faces_;
        custom_face_order_  = cell.custom_face_order_;
        face_nodes_         = cell.face_nodes_;
        face_nodes_in_custom_order_ = cell.face_nodes_in_custom_order_;

        num_nodes_          = cell.num_nodes_;
        extra_nodes_        = cell.extra_nodes_;
        nodes_              = cell.nodes_;
        custom_node_order_  = cell.custom_node_order_;
        nodes_in_custom_order_ = cell.nodes_in_custom_order_;
      }
    return *this;
}


PolygonCell::~PolygonCell()
{
}

// NODES


size_t PolygonCell
::GetNumNodes() const
{
    return num_nodes_;
}


size_t PolygonCell
::GetNumExtraNodes() const
{
    return extra_nodes_.size();
}


size_t PolygonCell
::GetNodeGlobalIdOriginalOrder( size_t nid ) const
{
    if( nid < num_nodes_ )
        return nodes_[ nid ]->GetIdx();
    return extra_nodes_[ nid - num_nodes_ ]->GetIdx();
}


size_t PolygonCell
::GetNodeGlobalIdCustomOrder( size_t nid ) const
{
    if( nid < num_nodes_ )
        return nodes_[ custom_node_order_[nid] ]->GetIdx();
    return extra_nodes_[ nid - num_nodes_ ]->GetIdx();
}


const csmp::Point<3U>& PolygonCell
::GetPointOriginalOrder( size_t nid ) const
{
    if( nid < num_nodes_ )
        return nodes_[ nid ]->GetPoint();
    return extra_nodes_[ nid - num_nodes_ ]->GetPoint();
}


const csmp::Point<3U>& PolygonCell
::GetPointCustomOrder( size_t nid ) const
{
    if( nid < num_nodes_ )
        return nodes_[ custom_node_order_[nid] ]->GetPoint();
    return extra_nodes_[ nid - num_nodes_ ]->GetPoint();
}


const csmp::Point<3U>& PolygonCell
::GetExtraPoint( size_t nid ) const
{
    return extra_nodes_[nid]->GetPoint();
}


GridNode* PolygonCell
::GetNodeOriginalOrder( size_t nid )
{
    if( nid < num_nodes_ )
        return nodes_[ nid ];
    return extra_nodes_[ nid - num_nodes_ ];
}



GridNode* PolygonCell
::GetNodeCustomOrder( size_t nid )
{
    if( nid < num_nodes_ )
        return nodes_[ custom_node_order_[nid] ];
    return extra_nodes_[ nid - num_nodes_ ];
}



GridNode* PolygonCell
::GetExtraNode( size_t nid )
{
    return extra_nodes_[nid];
}



void PolygonCell::AddExtraNode( const csmp::Point<3U>& pt )
{
    extra_nodes_.reserve( extra_nodes_.size() + 1 );
    extra_nodes_.push_back( grid_->AddNode( pt ) );
std::cerr <<"\ncalled PolygonCell::AddExtraNode:";
}


/// cell centroid ( if exist: by convention it's a first node of extra nodes arrays )
size_t PolygonCell
::GetCellCentroidGlobalId() const
{
    if( !extra_nodes_.empty() )
        return extra_nodes_[0]->GetIdx();
    return std::numeric_limits<uint32_t>::max();
}


const csmp::Point<3U>& PolygonCell
::GetCellCentroidPoint() const
{
    csmp::ErrorHandler& csmp_error( ErrorHandler::Instance() );
    if( !extra_nodes_.empty() )
        return extra_nodes_[0]->GetPoint();
    csmp_error.notice( csmp::ERROR,
                       "GetCellCentroidPoint()",
                       "Centroid node doesn't exist!!!");
    return nodes_[0]->GetPoint();
}



GridNode* PolygonCell
::GetCellCentroidNode()
{
    if( !extra_nodes_.empty() )
        return extra_nodes_[0];
    return NULL;
}



void PolygonCell
::InitializeNodeOrder( size_t num_nodes )
{
    num_nodes_ = num_nodes;
    custom_node_order_.clear();
    custom_node_order_.resize( num_nodes_ );
    nodes_in_custom_order_.resize( num_nodes_ );
    for( size_t nid = 0 ; nid < num_nodes_; ++nid )
    {
        custom_node_order_[nid] = nid;
        nodes_in_custom_order_[nid] = nodes_[nid];
    }
}



void PolygonCell
::AssignNodeOrder( size_t cnid, size_t onid )
{
    custom_node_order_.resize( num_nodes_ );
    nodes_in_custom_order_.resize( num_nodes_ );
    custom_node_order_[cnid] = onid;
    nodes_in_custom_order_[cnid] = nodes_[onid];
}

// FACES

/// returns how many faces make up the polygon cell
size_t PolygonCell::GetNumFaces() const
{
    return num_faces_;
}



/// returns 2 values for each spatial dimension (x-1, x+1, y-1, y+1, z-1, z+1) ?
size_t PolygonCell::GetNumSubFaces( size_t fid ) const
{
    assert( faces_[fid].size() <= 6 );
    return faces_[fid].size();
}



/**
    returns the number of surface nodes, but I left the name this way as well as FaceType and so on 
    because I wanted to have an access to the final instances which are Faces. 
    So if there is even an undivided face the access happens through the subface id, which is  0 for a single Face
   
    returns the number of nodes that make up the specific face

    /// faces[i][j] is a 6 x 5 matrix of faces defined by face-type and corresponding node pointers as entries (type, value pairs)
    std::vector<std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode<dim>*> > > > faces_;

*/
size_t PolygonCell::GetNumFaceNodes( size_t fid, size_t sfid ) const
{
  assert( fid < faces_.size() );
  assert( !faces_[fid].empty() );
  assert( sfid < faces_[fid].size() );
    return faces_[fid][sfid].second.size();
}


size_t PolygonCell
::GetFaceNodeGlobalId( size_t fid, size_t sfid, size_t nid ) const
{
    return faces_[fid][sfid].second[nid]->GetIdx();
}


const csmp::Point<3U>& PolygonCell
::GetFacePoint( size_t fid, size_t sfid, size_t nid ) const
{
    return faces_[fid][sfid].second[nid]->GetPoint();
}


GridNode* PolygonCell
::GetFaceNode( size_t fid, size_t sfid, size_t nid )
{
    return faces_[fid][sfid].second[nid];
}


const csmp::CSMP_FEM_TYPE& PolygonCell
::GetFaceType( size_t fid, size_t sfid ) const
{
    return faces_[fid][sfid].first;
}


size_t PolygonCell
::GetFaceDim( size_t fid, size_t sfid ) const
{
    if( CSMP_ElementSpecifications::SurfaceElement( GetFaceType( fid, sfid ) ) )
        return 2U;
    if( CSMP_ElementSpecifications::LineElement( GetFaceType( fid, sfid ) ) )
        return 1U;
    return 0U;
}


size_t PolygonCell
::GetNumPolygonFaceNodes( size_t fid ) const
{
    return face_nodes_[fid].size();
}


size_t PolygonCell::GetPolygonFaceNodeLocalId( size_t fid, size_t nid ) const
{
    return face_nodes_[fid][nid];
}


size_t PolygonCell::GetPolygonFaceNodeGlobalId( size_t fid, size_t nid ) const
{
    return nodes_[face_nodes_[fid][nid]]->GetIdx();
}


const csmp::Point<3U>& PolygonCell::GetPolygonFacePoint( size_t fid, size_t nid ) const
{
    return nodes_[face_nodes_[fid][nid]]->GetPoint();
}


GridNode* const PolygonCell::GetPolygonFaceNode( size_t fid, size_t nid ) const
{
    return nodes_[face_nodes_[fid][nid]];
}


void PolygonCell::AssignPolygonFaceNode( size_t fid, size_t nid, size_t cell_nid )
{
    face_nodes_[fid][nid] = cell_nid;
}


void PolygonCell::AddPolygonFaceNode( size_t fid, size_t cell_nid )
{
   face_nodes_[fid].reserve( face_nodes_[fid].size() + 1 );
   face_nodes_[fid].push_back( cell_nid );
}


void PolygonCell::ClearPolygonFaceNodes( size_t fid )
{
    face_nodes_[fid].clear();
}


void PolygonCell::ResizePolygonFaceNodes( size_t fid, size_t size )
{
    face_nodes_[fid].resize( size );
}


void PolygonCell::InitializeFaceOrder( size_t num_faces )
{
    num_faces_ = num_faces;
    custom_face_order_.clear();
    face_nodes_in_custom_order_.clear();
    custom_face_order_.resize( num_faces_ );
    face_nodes_in_custom_order_.resize( num_faces_ );
    face_nodes_.resize( num_faces_ );
    faces_.resize( num_faces_ );
    for( size_t fid = 0 ; fid < num_faces_; ++fid )
    {
        custom_face_order_[ fid ] = fid;
        const size_t num_face_nodes( face_nodes_[fid].size() );
        for ( size_t nid = 0 ; nid < num_face_nodes; ++nid ) {
             face_nodes_in_custom_order_[fid].reserve( face_nodes_in_custom_order_[fid].size() + 1 );
             face_nodes_in_custom_order_[fid].push_back( GetCustomNodeId( face_nodes_[fid][nid] ) );
          }
    }
}



void PolygonCell::AssignFaceOrder( size_t cfid , size_t ofid )
{
    custom_face_order_.resize( num_faces_ );
    custom_face_order_[cfid] = ofid;
    face_nodes_in_custom_order_.resize( num_faces_ );
    face_nodes_in_custom_order_[cfid].clear();
    face_nodes_.resize( num_faces_ );
    faces_.resize( num_faces_ );
    const size_t num_face_nodes( face_nodes_[ofid].size() );
    for( size_t nid = 0 ; nid < num_face_nodes; ++nid ) {
         face_nodes_in_custom_order_[cfid].reserve( face_nodes_in_custom_order_[cfid].size() + 1 );
         face_nodes_in_custom_order_[cfid].push_back( GetCustomNodeId( face_nodes_[ofid][nid] ) );
      }
}



void PolygonCell::AddFace( size_t fid, const std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> >& face )
{
    faces_[fid].reserve( faces_[fid].size() + 1 );
    faces_[fid].push_back( face );
}


void PolygonCell::EraseFace( size_t fid )
{
    faces_.erase( faces_.begin() + fid );
}


void PolygonCell::EraseFace( size_t fid, size_t position )
{
    faces_[fid].erase( faces_[fid].begin() + position );
}



void PolygonCell::ClearFaces( size_t fid )
{
    faces_[fid].clear();
}



void PolygonCell::ClearFaces( )
{
    faces_.clear();
}





// ELEMENTS


size_t PolygonCell
::GetNumElements() const
{
    return elements_.size();
}


size_t PolygonCell
::GetNumElementNodes( size_t eid ) const
{
    return elements_[eid].second.size();
}


size_t PolygonCell
::GetElementNodeGlobalId( size_t eid, size_t nid ) const
{
    return elements_[eid].second[nid]->GetIdx();
}


const csmp::Point<3U>& PolygonCell
::GetElementPoint( size_t eid, size_t nid ) const
{
    return elements_[eid].second[nid]->GetPoint();
}


GridNode* PolygonCell
::GetElementNode( size_t eid, size_t nid )
{
    return elements_[eid].second[nid];
}


const csmp::CSMP_FEM_TYPE& PolygonCell
::GetElementType( size_t eid ) const
{
    return elements_[eid].first;
}



size_t PolygonCell
::GetElementDim( size_t eid ) const
{
    if( CSMP_ElementSpecifications::VolumeElement( elements_[eid].first ) )
        return 3U;
    if( CSMP_ElementSpecifications::SurfaceElement( elements_[eid].first ) )
        return 2U;
    if( CSMP_ElementSpecifications::LineElement( elements_[eid].first ) )
        return 1U;
    return 0U;
}



size_t PolygonCell
::GetElementDim( ) const
{
    size_t cell_dim(0);
    const size_t num_elmts( elements_.size() );
    for( size_t eid = 0; eid < num_elmts; ++eid )
        cell_dim = std::max( cell_dim, GetElementDim( eid ) );
    return cell_dim;
}



void PolygonCell::AddCell( const std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> >& cell  )
{
    elements_.reserve( elements_.size() + 1 );
    elements_.push_back( cell );
}



void PolygonCell::EraseCell( size_t position )
{
    elements_.erase( elements_.begin() + position );
}



void PolygonCell::ClearCells( )
{
    elements_.clear();
}


// INDEXING

size_t PolygonCell
::GetOriginalNodeId( size_t custom_nid ) const
{
    return custom_node_order_[ custom_nid ];
}



size_t PolygonCell
::GetCustomNodeId( size_t original_nid ) const
{
    ErrorHandler& csmp_error( ErrorHandler::Instance() );

    for( size_t nid = 0; nid<num_nodes_; ++nid )
        if( custom_node_order_[nid] == original_nid )
            return nid;
    csmp_error.notice( ERROR,
                      "CustomNodeId()",
                      "Original node id is out of range!!!");
    return -1;
}



size_t PolygonCell::GetOriginalFaceId( size_t custom_fid ) const
{
    return custom_face_order_[ custom_fid ];
}



size_t PolygonCell
::GetCustomFaceId( size_t original_fid ) const
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    const size_t num_faces( 6U );
    for( size_t fid = 0; fid<num_faces; ++fid )
        if( custom_face_order_[fid] == original_fid )
            return fid;
    csmp_error.notice(csmp::ERROR,
                      "CustomFaceId()",
                      "Original face id is out of range!!!");
    return -1;
}




void PolygonCell::Out() const
 {
    std::cout<<"\nPolygonCell::Out: ";
    /// grid manager
    //PolygonGridManager<dim>* grid_;
    std::cout<<"\n\tparent grid manager 'PolygonGridManager' is not output...\n";

    std::cout <<"\nPolygonCell::Out:\n";
    std::cout <<"\tnodes: "<< num_nodes_ <<", extra nodes: "<< extra_nodes_.size();
    std::cout <<", custom node order:\n\t";
    for ( auto i{0}; i<custom_node_order_.size(); i++ ) std::cout << custom_node_order_[i] <<" ";

    std::cout <<"\n\tnodes:\n";
    for ( typename std::vector<GridNode*>::const_iterator it=nodes_.begin(); it!=nodes_.end(); ++it )
      std::cout << *(*it) <<" ";

    std::cout <<"\n\tnodes in custom order:\n";
    for ( typename std::vector<GridNode*>::const_iterator it=nodes_in_custom_order_.begin(); it!=nodes_in_custom_order_.end(); ++it )
      std::cout << *(*it) <<" ";

    std::cout <<"\n\textra nodes:\n";
    for ( typename std::vector<GridNode*>::const_iterator it=extra_nodes_.begin(); it!=extra_nodes_.end(); ++it )
      std::cout << *(*it) <<" ";
   
    /// cells
    std::cout <<"\n\telements:\n";
    //std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode<dim>*> > > elements_;
    for ( auto it=elements_.begin(); it!=elements_.end(); it++ ) {
         std::cout <<"\n"<< parseFiniteElementType((*it).first) <<": ";
         for ( typename std::vector<GridNode*>::const_iterator nit=(*it).second.begin();
              nit!=(*it).second.end(); nit++ )
           std::cout << *(*nit) <<" ";
      }

    /// faces
    std::cout <<"\n\tfaces ("<< num_faces_ << ") and their nodes:\n";
    for ( auto i{0}; i<face_nodes_.size(); i++ ) {
         std::cout <<"\tface "<< i+1 <<": ";
         for ( size_t j=0U; j<face_nodes_[i].size(); j++ )
           std::cout << face_nodes_[i][j] <<" ";
         std::cout <<"\n";
      }

    std::cout <<"\n\tfaces ("<< num_faces_ << ") and their nodes in custom order:\n";
    for ( auto i{0}; i<face_nodes_in_custom_order_.size(); i++ ) {
         std::cout <<"\tface "<< i+1 <<": ";
         for ( size_t j=0U; j<face_nodes_in_custom_order_[i].size(); j++ )
           std::cout << face_nodes_in_custom_order_[i][j] <<" ";
         std::cout <<"\n";
      }

    std::cout <<"\n\tface order :\n";
    for ( auto i{0}; i<custom_face_order_.size(); i++ )
      std::cout << custom_face_order_[i] <<" ";
    std::cout <<"\n";

    /// faces
    // std::vector<std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode<dim>*> > > > faces_;
    std::cout <<"\n\tfaces:\n";
    for ( auto fit=faces_.begin(); fit!=faces_.end(); fit++ )
      for ( auto it=(*fit).begin(); it!=(*fit).end(); it++ ) {
           std::cout <<"\n"<< parseFiniteElementType((*it).first) <<": ";
           for ( typename std::vector<GridNode*>::const_iterator nit=(*it).second.begin();
                nit!=(*it).second.end(); nit++ )
             std::cout << *(*nit) <<" ";
        }

    std::cout <<"\n\n";

 } // end Out

} // end namespace csmp
