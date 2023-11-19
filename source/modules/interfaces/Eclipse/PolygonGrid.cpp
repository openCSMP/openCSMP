#include "PolygonGrid.h"

#include "ErrorHandler.h"

using namespace std;

namespace csmp {

// GRID NODE

GridNode::GridNode( size_t index, const Point<3U>& pt )
: Point(pt), idx_( index )
 {
 }


GridNode::GridNode( const GridNode& gn )
: Point(*this), idx_( gn.idx_ )
{
}


GridNode& GridNode::operator=( const GridNode& gn )
{
    if ( &gn != this ) {
         idx_ = gn.idx_;
      }
    return *this;
}


bool GridNode::operator==( const GridNode& gn )
{
    if( this->Point<3>::operator==(gn) == true ) return true;
    return false;
}


size_t GridNode::GetIdx() const
{
    return idx_;
}


bool GridNode::IsNullIdx() const
{
   // SKM FIX
    return ( idx_ == UINT_MAX );
}


const csmp::Point<3U>& GridNode::GetPoint() const
{
    return *this;
}


csmp::Point<3U>& GridNode::GetPoint()
{
    return *this;
}


void GridNode::AssignNullIdx()
{
   // SKM FIX: assignment of a negative number to an unsigned int lead to undefined results
    idx_ = UINT_MAX;
}


void GridNode::AssignIdx( size_t idx )
{
    idx_ = idx;
}


void GridNode::AssignPoint( const csmp::Point<3U>& pt )
{
    static_cast<Point<3U> >(*this) = pt;
}






ostream&  operator<<( ostream& stream, const GridNode& gn )
 {
     stream <<" idx: "<< gn.GetIdx() <<", xyz: ";
     for ( auto i{0U}; i<3U; i++ ) cout << gn[i] <<" ";
     return stream;
 }
 





// GRID FACE

GridFace::GridFace()
:type_(csmp::UNKNOWN),
 idx_(UINT_MAX)
{
}



GridFace::GridFace( const GridFace& gf )
:type_ ( gf.type_  ),
 nodes_( gf.nodes_ ),
 idx_( gf.idx_ )
{
}


GridFace& GridFace::operator=( const GridFace& gf )
{
    if( this != &gf )
    {
        type_   =  gf.type_;
        nodes_  =  gf.nodes_;
        idx_    =  gf.idx_;
    }
    return *this;
}



GridFace::~GridFace()
{
}



bool GridFace::operator==( const GridFace& gf )
{
    if( gf.type_ == type_ )
    {
        const size_t gf_num_nodes( gf.GetNumNodes() );
        const size_t num_nodes( GetNumNodes() );
        if( gf_num_nodes != num_nodes )
            return false;

        set<csmp::Point<3U> > gf_pts;
        for( size_t nid=0;nid<gf_num_nodes; ++nid)
            gf_pts.insert( gf.GetNode( nid )->GetPoint() );
        set<csmp::Point<3U> > pts;
        for( size_t nid=0;nid<num_nodes; ++nid)
            pts.insert( GetNode( nid )->GetPoint() );
        typename set<csmp::Point<3U> >::const_iterator gf_it = gf_pts.begin();
        typename set<csmp::Point<3U> >::const_iterator it = pts.begin();
        typename set<csmp::Point<3U> >::const_iterator itEnd( pts.end() );
        for( ; it != itEnd; ++it, ++gf_it )
            if( *it != *gf_it )
                return false;
        return true;
    }
    return false;
}


const CSMP_FEM_TYPE& GridFace::GetType() const
{
    return type_;
}


size_t GridFace::GetNumNodes() const
{
    return nodes_.size();
}


GridNode* GridFace::GetNode( size_t nid )
{
    return nodes_[ nid ];
}


GridNode* const GridFace::GetNode( size_t nid ) const
{
    return nodes_[ nid ];
}


void GridFace::AssignType( const csmp::CSMP_FEM_TYPE& type )
{
    type_ = type;
}


void GridFace::AssignIdx( size_t idx )
{
    idx_ = idx;
}

void GridFace::AddNode( GridNode* gn )
{
    nodes_.reserve( nodes_.size() + 1 );
    nodes_.push_back( gn );
}










// GRID ELEMENT

GridElement::GridElement()
:type_(csmp::UNKNOWN),
 idx_(UINT_MAX)
{
}


GridElement::GridElement( const GridElement& ge )
:type_ ( ge.type_  ),
 nodes_( ge.nodes_ ),
 idx_  ( ge.idx_   )
{
}


GridElement& GridElement::operator=( const GridElement& ge )
{
    if( this != &ge )
    {
        type_   =  ge.type_;
        nodes_  =  ge.nodes_;
        idx_    =  ge.idx_;
    }
    return *this;
}


GridElement::~GridElement()
{
}


bool GridElement::operator==( const GridElement& ge )
{
    if( ge.type_ == type_ )
    {
        const size_t gf_num_nodes( ge.GetNumNodes() );
        const size_t num_nodes( GetNumNodes() );
        if( gf_num_nodes != num_nodes )
            return false;

        set<csmp::Point<3U> > gf_pts;
        for( size_t nid=0;nid<gf_num_nodes; ++nid)
            gf_pts.insert( ge.GetNode( nid )->GetPoint() );
        set<csmp::Point<3U> > pts;
        for( size_t nid=0;nid<num_nodes; ++nid)
            pts.insert( GetNode( nid )->GetPoint() );
        typename set<csmp::Point<3U> >::const_iterator gf_it = gf_pts.begin();
        typename set<csmp::Point<3U> >::const_iterator it = pts.begin();
        typename set<csmp::Point<3U> >::const_iterator itEnd( pts.end() );
        for( ; it != itEnd; ++it, ++gf_it )
            if( *it != *gf_it )
                return false;
        return true;
    }
    return false;
}


const csmp::CSMP_FEM_TYPE& GridElement::GetType() const
{
    return type_;
}


size_t GridElement::GetNumNodes() const
{
    return nodes_.size();
}


GridNode* GridElement::GetNode( size_t nid )
{
    return nodes_[ nid ];
}


GridNode* const GridElement::GetNode( size_t nid ) const
{
    return nodes_[ nid ];
}


void GridElement::AssignType( const csmp::CSMP_FEM_TYPE& type )
{
    type_ = type;
}


void GridElement::AssignIdx( size_t idx )
{
    idx_ = idx;
}


void GridElement::AddNode( GridNode* gn )
{
    nodes_.reserve( nodes_.size() + 1 );
    nodes_.push_back( gn );
}






// POLYGON GRID

PolygonGrid::PolygonGrid()
{
}




PolygonGrid::PolygonGrid( const PolygonGrid& pg )
:nodes_     ( pg.nodes_     ),
 faces_     ( pg.faces_     ),
 elements_  ( pg.elements_  )
{
}


PolygonGrid& PolygonGrid::operator=( const PolygonGrid& pg )
{
    if( &pg != this )
    {
        nodes_      = pg.nodes_;
        faces_      = pg.faces_;
        elements_   = pg.elements_;
    }
    return *this;
}


PolygonGrid::~PolygonGrid()
{
    /// delete nodes
    const size_t num_nodes( nodes_.size() );
    for(auto i =0; i<num_nodes; ++i)
        delete nodes_[i];
    /// delete faces
    const size_t num_faces( faces_.size() );
    for(auto i =0; i<num_faces; ++i)
        delete faces_[i];
    /// delete elements
    const size_t num_elements( elements_.size() );
    for(auto i =0; i<num_elements; ++i)
        delete elements_[i];
}


GridNode* PolygonGrid::AddNode( GridNode& gn )
{
    typename deque<GridNode*>::iterator nit = find_if( nodes_.begin(), nodes_.end(), [&]( GridNode* n )->bool { return ( gn == *n ); } );
    if( nit != nodes_.end() )
        return (*nit);
    gn.AssignIdx( nodes_.size() );
    nodes_.push_back( new GridNode(gn) );
    return nodes_.back();
}


GridFace* PolygonGrid::AddFace( GridFace& gf )
{
    typename deque<GridFace*>::iterator nit = find_if( faces_.begin(), faces_.end(), [&]( GridFace* f )->bool { return ( gf == *f ); } );
    if( nit != faces_.end() )
        return (*nit);
    gf.AssignIdx( faces_.size() );
    faces_.push_back( new GridFace(gf) );
    return faces_.back();
}


GridElement* PolygonGrid::AddElement( GridElement& ge )
{
    typename deque<GridElement*>::iterator nit = find_if( elements_.begin(), elements_.end(), [&]( GridElement* e )->bool { return ( ge == *e ); } );
    if( nit != elements_.end() )
        return (*nit);
    ge.AssignIdx( elements_.size() );
    elements_.push_back( new GridElement(ge) );
    return elements_.back();
}







// POLYGON GRID MANAGER


PolygonGridManager::PolygonGridManager()
{
}



PolygonGridManager::PolygonGridManager( const PolygonGridManager& pgm )
: points_    ( pgm.points_     ),
  grid_nodes_( pgm.grid_nodes_ )
{
}



PolygonGridManager& PolygonGridManager::operator=( const PolygonGridManager& pgm )
{
    if( &pgm != this )
      {
          points_     = pgm.points_;
          grid_nodes_ = pgm.grid_nodes_;
      }
    return *this;
}



PolygonGridManager::~PolygonGridManager()
{
   Clear();
}


/** 
    deleting the grid nodes that the pointers in the grid_nodes vector point to.
*/
void PolygonGridManager::Clear()
{
    /// delete nodes
    const size_t num_nodes( grid_nodes_.size() );
    for(size_t i{0U}; i<num_nodes; ++i)
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


size_t PolygonGridManager::GetNumNodes() const
{
    return grid_nodes_.size();
}


const csmp::Point<3U>& PolygonGridManager::GetPoint( size_t nid ) const
{
    return grid_nodes_[nid]->GetPoint();
}


GridNode* PolygonGridManager::GetNode( size_t nid )
{
    return grid_nodes_[nid];
}


GridNode* PolygonGridManager::GetNode( const csmp::Point<3U>& pt )
{
    return grid_nodes_[points_.at( pt )];
}



/**
    Methods dynamically allocates nodes on the heap,
    storing the pointers to them in the grid_nodes_ array, 
    but also returning them to its caller.
    
    @todo SKM: when vector is resized all the stored pointers will be invalidated during rereconstruction; avoid!
*/

GridNode* PolygonGridManager::AddNode( const csmp::Point<3U>& pt )
{
    const size_t nid( points_.size() );
    pair<typename map<csmp::Point<3U>,size_t>::iterator,bool>
            pit = points_.insert( make_pair( pt, nid ) );
  
    if( pit.second == true )
      {
//          assert( grid_nodes_.capacity() > grid_nodes_.size() );
          grid_nodes_.push_back( new GridNode(nid,pt) );
      }
    return grid_nodes_[ (*pit.first).second ];
}

/* SKM REPLACED

template<uint32_t dim>
GridNode* PolygonGridManager::AddNode( const csmp::Point<dim>& pt )
{
    const size_t nid( points_.size() );
    pair<typename map<csmp::Point<dim>,size_t>::iterator,bool>
            pit = points_.insert( make_pair( pt, nid ) );
  
    if( pit.second == true )
      {
          // TODO: Why is there no custom constructor for GridNode that could be called with new() ?
          GridNode gn;
          gn.AssignIdx( nid );
          gn.AssignPoint( pt );
          //grid_nodes_.reserve( grid_nodes_.size() + 1 );
          assert( grid_nodes_.capacity() > grid_nodes_.size() );
          grid_nodes_.push_back( new GridNode(nid) );
      }
    return grid_nodes_[ (*pit.first).second ];
}
*/




void PolygonGridManager::Out() const
 {
    cout<<"\nPolygonGridManager::Out: points: "<< points_.size() <<", grid nodes: "<< grid_nodes_.size() <<"\n";
    cout<<"\tpoints and their indices:";
   
    for ( auto it=points_.begin(); it!=points_.end(); it++ )
      cout << (*it).first <<", "<< (*it).second;
    cout<<"\n";

    cout <<"\n\tgrid nodes:\n";
    for ( typename vector<GridNode*>::const_iterator it=grid_nodes_.begin(); it!=grid_nodes_.end(); ++it )
      cout << *(*it) <<" ";
    cout<<"\n";

 } // end Out

} // end namespace csmp
