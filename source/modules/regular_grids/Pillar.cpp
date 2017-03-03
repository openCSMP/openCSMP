#include "Pillar.h"
#include "ErrorHandler.h"

namespace csmp {

// PILLARS

Pillar::Pillar()
:a_(0.0),
 b_(0.0)
{
}

Pillar::Pillar( const csmp::Point<3U>& pStart, const csmp::Point<3U>& pEnd )
{
    AssignEnds( pStart, pEnd );
}

Pillar::Pillar( const Pillar& pr )
: a_(pr.a_),
  b_(pr.b_),
  points_(pr.points_),
  attached_cells_(pr.attached_cells_)
{
}

Pillar& Pillar::operator=( const Pillar& pr )
{
    if( this != &pr )
    {
        a_              = pr.a_;
        b_              = pr.b_;
        points_         = pr.points_;
        attached_cells_ = pr.attached_cells_;
    }
    return *this;
}

Pillar::~Pillar()
{
    points_.clear();
    attached_cells_.clear();
}

// ASSIGNMENTS & ACCESS

void Pillar::AssignEnds( const csmp::Point<3U>& pStart, const csmp::Point<3U>& pEnd )
{
    a_ = pStart;
    b_ = ( pEnd - pStart );
}

void Pillar::AssignEnds()
{
    const size_t num_points( points_.size() );
    a_ = points_[0];
    b_ = points_[num_points-1];
}

void Pillar::AssignAttachedCellId( size_t i, size_t j, size_t offset )
{
    attached_cells_.insert( std::make_pair( std::make_pair(i,j), offset ) );
}

size_t Pillar::GetNumPoints() const
{
    return points_.size();
}

/**
    Copy of point is needed because it is going to be modified.
*/
void Pillar::AddPoint( csmp::Point<3U> zpt )
{
   ResolveXYcoords( zpt );
   points_.push_back( zpt );
}

const csmp::Point<3U>& Pillar::GetPoint( size_t zid ) const
{
    return points_[ zid ];
}

// CELL IDS

size_t Pillar
::NodeIndexIncrementI( size_t nid ) const
{
    return ( nid%2 );
}

size_t Pillar
::NodeIndexIncrementJ( size_t nid ) const
{
    return ( (nid/2)%2 );
}

size_t Pillar
::NodeIndexIncrementK( size_t nid ) const
{
    return ( nid/4 );
}

size_t Pillar::GetNumAttachedCells() const
{
    return attached_cells_.size();
}

const csmp::Point<3U>& Pillar::GetPoint( size_t xid, size_t yid, size_t zid, size_t nid ) const
{
    const size_t pzid( zid*2 + NodeIndexIncrementK( nid ) );
    const size_t offset( attached_cells_.at( std::make_pair(xid,yid) ) );
    return points_[ pzid*attached_cells_.size() + offset ];
}

void Pillar::AssignPoint( size_t xid, size_t yid, size_t zid, size_t nid, const csmp::Point<3U>& pt )
{
    const size_t pzid( zid*2 + NodeIndexIncrementK( nid ) );
    const size_t offset( attached_cells_.at( std::make_pair(xid,yid) ) );
    points_[ pzid*attached_cells_.size() + offset ] = pt;
}

// RESOLVE MISSING COORDINATES

/**
    Calculates the XY coordinates of the supplied point that is supposed to lie on the subvertical Pillar
    from its original ones. Then these new coordinates are assigned.

    r = a + t*b;
    x: r[0] = a[0] + t*b[0]
    y: r[1] = a[1] + t*b[1]
    z: r[2] = a[2] + t*b[2]
*/
void Pillar::ResolveXYcoords( csmp::Point<3U>& r )
{
    // not safe: if( b_ != 0.0 ) return;
    if ( std::fabs(b_[0]) <= std::numeric_limits<double>::epsilon() and
         std::fabs(b_[1]) <= std::numeric_limits<double>::epsilon() and
         std::fabs(b_[2]) <= std::numeric_limits<double>::epsilon() )
      return;

    const double t = (r[2U] - a_[2U]) / b_[2U];
    r[0U] = a_[0U] + t*b_[0U];
    r[1U] = a_[1U] + t*b_[1U];
}

/// simultaneously resolves the XY coordinates for all pillar beads for which a Z value has already been assigned
void Pillar::ResolveXYcoords()
{
    const size_t num_points( points_.size() );
    for( size_t pid = 0; pid<num_points; ++pid )
        ResolveXYcoords( points_[pid] );
}

csmp::Point<3U> Pillar::ResolveXYcoords( double z )
{
    /// calculate x and y coordinates based on z
    csmp::Point<3U> pt(0.,0.,z);
    ResolveXYcoords( pt );
    return pt;
}




void Pillar::Out(std::ostream& os) const
 {
    os <<"\nPillar::Out:\n";
    os <<"\torigin: "<< a_ <<", direction: "<< b_;
    os <<"\n\tstored points ("<< points_.size() <<"):\n\t";
    for ( auto it = points_.begin(); it!=points_.end(); ++it )
      os << (*it) <<" ";

    os <<"\n\n\tindices of attached cells ("<< attached_cells_.size() <<"):\n\t";
    for ( auto it = attached_cells_.begin(); it!=attached_cells_.end(); ++it )
      os <<"\t"<< (*it).first.first <<","<< (*it).first.second <<": "<< (*it).second <<"\n";
   
    os <<"\n";
 }





} // end namespace csmp









