#include "Pillar.h"
#include "ErrorHandler.h"

namespace csmp {

// PILLARS

Pillar::Pillar()
{
}

Pillar::Pillar( const csmp::Point<3U>& pStart, const csmp::Point<3U>& pEnd )
  : start_(pStart), end_(pEnd)
{
}

Pillar::~Pillar()
{
}
  
  void Pillar::Reserve(size_t zpoints)
  {
    points_.reserve(zpoints);
  }

size_t Pillar::GetNumPoints() const
{
    return points_.size();
}

void Pillar::AddZCoord( double z )
{
   points_.push_back( z );
}
  
void Pillar::SetZCoords(const std::vector<double> &zcoords)
{
  std::vector<double> newzcoords(zcoords.begin(), zcoords.end());
  points_ = std::move(newzcoords);
}
  
double Pillar::GetZCoord( size_t zid ) const
{
  return points_[zid];
}
  
size_t Pillar::FindPoint(double z) const
{
auto it = std::lower_bound(points_.begin(), points_.end(), z);
if (it != points_.end() && *it == z) {
    return std::distance(points_.begin(), it);
}
else {
    throw csmp::Exception(ERROR, "Pillar::FindPoint", "Could not find zcoord in pillar");
}
}

void Pillar::SetFirstNodeNum(size_t node_num)
{
firstNodeNum_ = node_num;
}
  
size_t Pillar::FirstNodeNum() const
{
return firstNodeNum_;
}


csmp::Point<3U> Pillar::GetPoint( size_t zid ) const
{
  Point<3U> p;
  assert( zid < points_.size() );
  double z = points_[zid];
  double t = (z - start_[2]) / (end_[2] - start_[2]);
  p[0] = start_[0] * (1.0 - t) + end_[0] * t;
  p[1] = start_[1] * (1.0 - t) + end_[1] * t;
  p[2] = z;
  return p;
}

void Pillar::Out() const
 {
    std::cout <<"\nPillar::Out:\n";
    std::cout <<"\tstart: "<< start_ <<", end: "<< end_;
    std::cout <<"\n\tstored points ("<< points_.size() <<"):\n\t";
   for ( auto z : points_ )
      std::cout << z <<" ";
    std::cout <<"\n";
 }

} // end namespace csmp

