#include "PolygonGrid.h"

#include "ErrorHandler.h"

namespace csmp {


// GRID NODE

template<size_t dim>
GridNode<dim>::GridNode( size_t index, const Point<dim>& pt )
: idx_( index ),
  pt_( pt )
 {
 }


template<size_t dim>
GridNode<dim>::GridNode()
: idx_( UINT_MAX ),
 pt_ ( 0.0 )
{
}

template<size_t dim>
GridNode<dim>::GridNode( const GridNode<dim>& gn )
:idx_( gn.idx_ ),
 pt_ ( gn.pt_  )
{
}

template<size_t dim>
GridNode<dim>& GridNode<dim>::operator=( const GridNode<dim>& gn )
{
    if( &gn != this )
    {
        idx_ = gn.idx_;
        pt_  = gn.pt_;
    }
    return *this;
}

template<size_t dim>
GridNode<dim>::~GridNode()
{
}

template<size_t dim>
bool GridNode<dim>::operator==( const GridNode<dim>& gn )
{
    if( gn.pt_ == gn.pt_ )
        return true;
    return false;
}

template<size_t dim>
size_t GridNode<dim>::GetIdx() const
{
    return idx_;
}

template<size_t dim>
bool GridNode<dim>::IsNullIdx() const
{
   // SKM FIX
    return ( idx_ == UINT_MAX );
}

template<size_t dim>
const csmp::Point<dim>& GridNode<dim>::GetPoint() const
{
    return pt_;
}

template<size_t dim>
csmp::Point<dim>& GridNode<dim>::GetPoint()
{
    return pt_;
}

template<size_t dim>
void GridNode<dim>::AssignNullIdx()
{
   // SKM FIX: assignment of a negative number to an unsigned int lead to undefined results
    idx_ = UINT_MAX;
}

template<size_t dim>
void GridNode<dim>::AssignIdx( size_t idx )
{
    idx_ = idx;
}

template<size_t dim>
void GridNode<dim>::AssignPoint( const csmp::Point<dim>& pt )
{
    pt_ = pt;
}

template class GridNode<1U>;
template class GridNode<2U>;
template class GridNode<3U>;


template<size_t dim>
std::ostream&  operator<<( std::ostream& stream, const GridNode<dim>& gn )
 {
     stream <<" idx: "<< gn.GetIdx() <<", xyz: ";
     for ( size_t i=0U; i<dim; i++ ) std::cout << gn.GetPoint()[i] <<" ";
     return stream;
 }

template std::ostream&  operator<<( std::ostream& stream, const GridNode<1U>& );
template std::ostream&  operator<<( std::ostream& stream, const GridNode<2U>& );
template std::ostream&  operator<<( std::ostream& stream, const GridNode<3U>& );
 


// GRID FACE

template<size_t dim>
GridFace<dim>::GridFace()
:type_(csmp::UNKNOWN),
 idx_(UINT_MAX)
{
}

template<size_t dim>
GridFace<dim>::GridFace( const GridFace<dim>& gf )
:type_ ( gf.type_  ),
 nodes_( gf.nodes_ ),
 idx_( gf.idx_ )
{
}

template<size_t dim>
GridFace<dim>& GridFace<dim>::operator=( const GridFace<dim>& gf )
{
    if( this != &gf )
    {
        type_   =  gf.type_;
        nodes_  =  gf.nodes_;
        idx_    =  gf.idx_;
    }
    return *this;
}

template<size_t dim>
GridFace<dim>::~GridFace()
{
}

template<size_t dim>
bool GridFace<dim>::operator==( const GridFace<dim>& gf )
{
    if( gf.type_ == type_ )
    {
        const size_t gf_num_nodes( gf.GetNumNodes() );
        const size_t num_nodes( GetNumNodes() );
        if( gf_num_nodes != num_nodes )
            return false;

        std::set<csmp::Point<dim> > gf_pts;
        for( size_t nid=0;nid<gf_num_nodes; ++nid)
            gf_pts.insert( gf.GetNode( nid )->GetPoint() );
        std::set<csmp::Point<dim> > pts;
        for( size_t nid=0;nid<num_nodes; ++nid)
            pts.insert( GetNode( nid )->GetPoint() );
        typename std::set<csmp::Point<dim> >::const_iterator gf_it = gf_pts.begin();
        typename std::set<csmp::Point<dim> >::const_iterator it = pts.begin();
        typename std::set<csmp::Point<dim> >::const_iterator itEnd( pts.end() );
        for( ; it != itEnd; ++it, ++gf_it )
            if( *it != *gf_it )
                return false;
        return true;
    }
    return false;
}

template<size_t dim>
const CSMP_FEM_TYPE& GridFace<dim>::GetType() const
{
    return type_;
}

template<size_t dim>
size_t GridFace<dim>::GetNumNodes() const
{
    return nodes_.size();
}

template<size_t dim>
GridNode<dim>* GridFace<dim>::GetNode( size_t nid )
{
    return nodes_[ nid ];
}

template<size_t dim>
GridNode<dim>* const GridFace<dim>::GetNode( size_t nid ) const
{
    return nodes_[ nid ];
}

template<size_t dim>
void GridFace<dim>::AssignType( const csmp::CSMP_FEM_TYPE& type )
{
    type_ = type;
}

template<size_t dim>
void GridFace<dim>::AssignIdx( size_t idx )
{
    idx_ = idx;
}

template<size_t dim>
void GridFace<dim>::AddNode( GridNode<dim>* gn )
{
    nodes_.reserve( nodes_.size() + 1 );
    nodes_.push_back( gn );
}

template class GridFace<1U>;
template class GridFace<2U>;
template class GridFace<3U>;


// GRID ELEMENT

template<size_t dim>
GridElement<dim>::GridElement()
:type_(csmp::UNKNOWN),
 idx_(UINT_MAX)
{
}

template<size_t dim>
GridElement<dim>::GridElement( const GridElement& ge )
:type_ ( ge.type_  ),
 nodes_( ge.nodes_ ),
 idx_  ( ge.idx_   )
{
}

template<size_t dim>
GridElement<dim>& GridElement<dim>::operator=( const GridElement& ge )
{
    if( this != &ge )
    {
        type_   =  ge.type_;
        nodes_  =  ge.nodes_;
        idx_    =  ge.idx_;
    }
    return *this;
}

template<size_t dim>
GridElement<dim>::~GridElement()
{
}

template<size_t dim>
bool GridElement<dim>::operator==( const GridElement<dim>& ge )
{
    if( ge.type_ == type_ )
    {
        const size_t gf_num_nodes( ge.GetNumNodes() );
        const size_t num_nodes( GetNumNodes() );
        if( gf_num_nodes != num_nodes )
            return false;

        std::set<csmp::Point<dim> > gf_pts;
        for( size_t nid=0;nid<gf_num_nodes; ++nid)
            gf_pts.insert( ge.GetNode( nid )->GetPoint() );
        std::set<csmp::Point<dim> > pts;
        for( size_t nid=0;nid<num_nodes; ++nid)
            pts.insert( GetNode( nid )->GetPoint() );
        typename std::set<csmp::Point<dim> >::const_iterator gf_it = gf_pts.begin();
        typename std::set<csmp::Point<dim> >::const_iterator it = pts.begin();
        typename std::set<csmp::Point<dim> >::const_iterator itEnd( pts.end() );
        for( ; it != itEnd; ++it, ++gf_it )
            if( *it != *gf_it )
                return false;
        return true;
    }
    return false;
}

template<size_t dim>
const csmp::CSMP_FEM_TYPE& GridElement<dim>::GetType() const
{
    return type_;
}

template<size_t dim>
size_t GridElement<dim>::GetNumNodes() const
{
    return nodes_.size();
}

template<size_t dim>
GridNode<dim>* GridElement<dim>::GetNode( size_t nid )
{
    return nodes_[ nid ];
}

template<size_t dim>
GridNode<dim>* const GridElement<dim>::GetNode( size_t nid ) const
{
    return nodes_[ nid ];
}

template<size_t dim>
void GridElement<dim>::AssignType( const csmp::CSMP_FEM_TYPE& type )
{
    type_ = type;
}

template<size_t dim>
void GridElement<dim>::AssignIdx( size_t idx )
{
    idx_ = idx;
}

template<size_t dim>
void GridElement<dim>::AddNode( GridNode<dim>* gn )
{
    nodes_.reserve( nodes_.size() + 1 );
    nodes_.push_back( gn );
}

template class GridElement<1U>;
template class GridElement<2U>;
template class GridElement<3U>;


// POLYGON GRID

template<size_t dim>
PolygonGrid<dim>::PolygonGrid()
{

}

template<size_t dim>
PolygonGrid<dim>::PolygonGrid( const PolygonGrid<dim>& pg )
:nodes_     ( pg.nodes_     ),
 faces_     ( pg.faces_     ),
 elements_  ( pg.elements_  )
{
}

template<size_t dim>
PolygonGrid<dim>& PolygonGrid<dim>::operator=( const PolygonGrid<dim>& pg )
{
    if( &pg != this )
    {
        nodes_      = pg.nodes_;
        faces_      = pg.faces_;
        elements_   = pg.elements_;
    }
    return *this;
}

template<size_t dim>
PolygonGrid<dim>::~PolygonGrid()
{
    /// delete nodes
    const size_t num_nodes( nodes_.size() );
    for(size_t i =0; i<num_nodes; ++i)
        delete nodes_[i];
    /// delete faces
    const size_t num_faces( faces_.size() );
    for(size_t i =0; i<num_faces; ++i)
        delete faces_[i];
    /// delete elements
    const size_t num_elements( elements_.size() );
    for(size_t i =0; i<num_elements; ++i)
        delete elements_[i];
}

template<size_t dim>
GridNode<dim>* PolygonGrid<dim>::AddNode( csmp::GridNode<dim>& gn )
{
    typename std::deque<GridNode<dim>*>::iterator nit = find_if( nodes_.begin(), nodes_.end(), [&]( csmp::GridNode<dim>* n )->bool { return ( gn == *n ); } );
    if( nit != nodes_.end() )
        return (*nit);
    gn.AssignIdx( nodes_.size() );
    nodes_.push_back( new GridNode<dim>(gn) );
    return nodes_.back();
}

template<size_t dim>
GridFace<dim>* PolygonGrid<dim>::AddFace( csmp::GridFace<dim>& gf )
{
    typename std::deque<GridFace<dim>*>::iterator nit = find_if( faces_.begin(), faces_.end(), [&]( csmp::GridFace<dim>* f )->bool { return ( gf == *f ); } );
    if( nit != faces_.end() )
        return (*nit);
    gf.AssignIdx( faces_.size() );
    faces_.push_back( new GridFace<dim>(gf) );
    return faces_.back();
}

template<size_t dim>
GridElement<dim>* PolygonGrid<dim>::AddElement( csmp::GridElement<dim>& ge )
{
    typename std::deque<GridElement<dim>*>::iterator nit = find_if( elements_.begin(), elements_.end(), [&]( csmp::GridElement<dim>* e )->bool { return ( ge == *e ); } );
    if( nit != elements_.end() )
        return (*nit);
    ge.AssignIdx( elements_.size() );
    elements_.push_back( new GridElement<dim>(ge) );
    return elements_.back();
}

template class PolygonGrid<1U>;
template class PolygonGrid<2U>;
template class PolygonGrid<3U>;


// POLYGON GRID MANAGER

template<size_t dim>
PolygonGridManager<dim>::PolygonGridManager()
{
}


template<size_t dim>
PolygonGridManager<dim>::PolygonGridManager( const PolygonGridManager<dim>& pgm )
: points_    ( pgm.points_     ),
  grid_nodes_( pgm.grid_nodes_ )
{
}


template<size_t dim>
PolygonGridManager<dim>& PolygonGridManager<dim>::operator=( const PolygonGridManager<dim>& pgm )
{
    if( &pgm != this )
      {
          points_     = pgm.points_;
          grid_nodes_ = pgm.grid_nodes_;
      }
    return *this;
}


template<size_t dim>
PolygonGridManager<dim>::~PolygonGridManager()
{
    Clear();
}


/** 
    deleting the grid nodes that the pointers in the grid_nodes vector point to.
*/
template<size_t dim>
void PolygonGridManager<dim>::Clear()
{
    /// delete nodes
    const size_t num_nodes( grid_nodes_.size() );
    for(size_t i=0; i<num_nodes; ++i)
      {
          if( grid_nodes_[i] != NULL )
            {
                delete grid_nodes_[i];
                grid_nodes_[i] = NULL;
            }
      }
    /// clear points
    grid_nodes_.clear();
    points_.clear();
}


template<size_t dim>
size_t PolygonGridManager<dim>::GetNumNodes() const
{
    return points_.size();
}

template<size_t dim>
size_t PolygonGridManager<dim>::GetNodeId( const csmp::Point<dim>& pt ) const
{
    return points_.at( pt );
}

template<size_t dim>
const csmp::Point<dim>& PolygonGridManager<dim>::GetPoint( size_t nid ) const
{
    return grid_nodes_[nid]->GetPoint();
}

template<size_t dim>
csmp::GridNode<dim>* PolygonGridManager<dim>::GetNode( size_t nid )
{
    return grid_nodes_[nid];
}

template<size_t dim>
csmp::GridNode<dim>* PolygonGridManager<dim>::GetNode( const csmp::Point<dim>& pt )
{
    return grid_nodes_[points_.at( pt )];
}



/**
    Methods dynamically allocates nodes on the heap,
    storing the pointers to them in the grid_nodes_ array, 
    but also returning them its caller.
    
    @todo SKM: when vector is resized all the stored pointers will be invalidated during rereconstruction; avoid!
*/
template<size_t dim>
csmp::GridNode<dim>* PolygonGridManager<dim>::AddNode( const csmp::Point<dim>& pt )
{
    const size_t nid( points_.size() );
    std::pair<typename std::map<csmp::Point<dim>,size_t>::iterator,bool>
            pit = points_.insert( std::make_pair( pt, nid ) );
  
    if( pit.second == true )
      {
//          assert( grid_nodes_.capacity() > grid_nodes_.size() );
          grid_nodes_.push_back( new csmp::GridNode<dim>(nid,pt) );
      }
    return grid_nodes_[ (*pit.first).second ];
}

/* SKM REPLACED

template<size_t dim>
csmp::GridNode<dim>* PolygonGridManager<dim>::AddNode( const csmp::Point<dim>& pt )
{
    const size_t nid( points_.size() );
    std::pair<typename std::map<csmp::Point<dim>,size_t>::iterator,bool>
            pit = points_.insert( std::make_pair( pt, nid ) );
  
    if( pit.second == true )
      {
          // TODO: Why is there no custom constructor for GridNode that could be called with new() ?
          csmp::GridNode<dim> gn;
          gn.AssignIdx( nid );
          gn.AssignPoint( pt );
          //grid_nodes_.reserve( grid_nodes_.size() + 1 );
          assert( grid_nodes_.capacity() > grid_nodes_.size() );
          grid_nodes_.push_back( new csmp::GridNode<dim>(nid) );
      }
    return grid_nodes_[ (*pit.first).second ];
}
*/



template<size_t dim>
void PolygonGridManager<dim>::Out(std::ostream& os) const
 {
    os<<"\nPolygonGridManager<"<< dim <<">::Out: points: "<< points_.size() <<", grid nodes: "<< grid_nodes_.size() <<"\n";
    os<<"\tpoints and their indices:";
   
    for ( auto it=points_.begin(); it!=points_.end(); it++ )
      os << (*it).first <<", "<< (*it).second;
    os<<"\n";

    os <<"\n\tgrid nodes:\n";
    for ( auto it=grid_nodes_.begin(); it!=grid_nodes_.end(); ++it )
      os << *(*it) <<" ";
    os<<std::endl;
 } // end Out




template class PolygonGridManager<1U>;
template class PolygonGridManager<2U>;
template class PolygonGridManager<3U>;

} // end namespace csmp









