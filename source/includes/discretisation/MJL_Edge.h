#ifndef MJL_EDGE_H
#define MJL_EDGE_H

#include "MJL_Point.h"

/// Lazlo p. 88
namespace mjl {

class Edge {
  public:
    Edge();
    Edge( const mjl::Point&, const mjl::Point& );
    Edge( const mjl::Edge& );
    ~Edge() {};
    Edge&  operator=( const mjl::Edge& m );
    void       Set( const mjl::Point&, const mjl::Point& );
    void       Set( double x1, double y1, double x2, double y2 );
    Edge&      Rot();
    Edge&      Flip();
    void       Move( double dx, double dy );
    void       Scale( double xfac, double yfac );
    void       NormalizeTo( double n );
    bool       IfCloserThanMerge( double tolerance, mjl::Point& );

    bool       operator==( const mjl::Edge& ) const;
    bool       operator!=( const mjl::Edge& ) const;
    bool       operator<( const mjl::Edge& )  const;

    mjl::Point Point( double t ) const;  // gives xy-coordinates of parametric value t along line
    double     Y( double x ) const;      // gives y-coordinate for an x-coordinate along infinite line
    mjl::Point Origin() const;
    mjl::Point Destination() const;
    void       MidPoint( mjl::Point& ) const;
    mjl::Point MidPoint() const;
    bool       IsVertical() const;
    double     Slope() const;
    double     DistanceFrom( const mjl::Point& ) const;
    double     Length() const;
    
    LOCATION      Classify( const mjl::Point& ) const;
    INTERSECTION  Intersect( const Edge&, double& ) const;
    INTERSECTION  Cross( const Edge&, double& ) const;
    INTERSECTION  CrossingPoint( double tolerance, const Edge&, mjl::Point& ) const;
    INTERSECTION  SKM_CrossingPoint( double tolerance, const Edge&, mjl::Point& ) const;
    EDGE_CLASSIFICATION  Type( const mjl::Point& ) const;
    double  AngleTo( const Edge& ) const;
    double  DotProduct( const mjl::Point&, const mjl::Point& ) const;

    void Out() const;
    
    mjl::Point org_, dest_;
};


bool                aimsAt( const Edge&, const Edge&, LOCATION, INTERSECTION ); // MJL, p. 156
double              signedAngle( const Point&, const Edge& ); // MJL, p. 120
INTERSECTION        crossingPoint( double tolerance, const Edge&, const Edge&, Point& ); // MJL, p. 162
EDGE_CLASSIFICATION edgeType( const mjl::Point&, const Edge& ); // MJL, p. 119
    

// inlined methods

inline Edge::Edge()
  :  org_(0.,0.), dest_(0.,0.)
 {
 }


inline Edge& Edge::operator=( const Edge& m )
 {
    if ( &m == this ) return *this;
    org_  = m.org_;
    dest_ = m.dest_;
    return *this;
 }


inline  bool Edge::operator==( const Edge& m ) const
 {
    return ((org_==m.org_) && (dest_==m.dest_));
 }

inline  bool Edge::operator!=( const Edge& m ) const
 {
    return ((org_!=m.org_) && (dest_!=m.dest_));
 }


inline  bool Edge::operator<( const Edge& m ) const
 {
    return ((org_<m.org_) || (org_==m.org_ && (dest_<m.dest_)) );
 }



inline Edge::Edge( const mjl::Point& org, const mjl::Point& dest )
  :  org_(org), dest_(dest)
 {
 }


inline Edge::Edge( const Edge& m )
   : org_(m.org_), dest_(m.dest_)
 {
 }


inline void  Edge::Set( const mjl::Point& org, const mjl::Point& dest )
 {
    org_  = org;
    dest_ = dest;
 }


inline void  Edge::Set( double x1, double y1, double x2, double y2 )
 {
    org_(0)  = x1;
    org_(1)  = y1;
    dest_(0) = x2;
    dest_(1) = y2;
 }



inline double Edge::DotProduct( const mjl::Point& p, const mjl::Point& q ) const
 {
    return( p[0] * q[0] + p[1] * q[1] );
 }



/// flip the direction of an edge
inline Edge&  Edge::Flip()
 {
    return Rot().Rot();
 }



inline  mjl::LOCATION  Edge::Classify( const mjl::Point& p ) const
 {
    return p.Classify( org_, dest_ );
 }    



inline mjl::Point  Edge::Point( double t ) const
 {
    mjl::Point tmp( dest_ - org_ );
    tmp *= t;
    return  mjl::Point( org_ + tmp );
 }    



inline mjl::Point  Edge::Origin() const { return org_; }



inline mjl::Point  Edge::Destination() const { return dest_; }



inline void Edge::MidPoint( mjl::Point& mp ) const 
 { 
    mp.Set( (org_[0]+dest_[0])/2., (org_[1]+dest_[1])/2. ); 
 }



inline mjl::Point  Edge::MidPoint() const 
 { 
    return mjl::Point( (org_[0]+dest_[0])/2., (org_[1]+dest_[1])/2. ); 
 }


/**

@return SKEW_CROSS if the line segments intersect (in constrast with the infinite lines)
COLLINEAR, PARALLEL
SKEW_NO_CROSS as appropriate
t = x value of intersection if intersection occurs
*/
inline  INTERSECTION  Edge::Cross( const Edge& e, double& t ) const
 {
    double s;
    INTERSECTION crossType = e.Intersect( *this, s );
    if ((crossType == COLLINEAR || crossType == PARALLEL)) return crossType;
    if ((s < 0.) || (s > 1.)) return SKEW_NO_CROSS;
    Intersect( e, t );
    if ((0. <= t) && (t <= 1.)) return SKEW_CROSS;
    return SKEW_NO_CROSS;
 }


inline  bool  Edge::IsVertical() const
 {
    return( org_[0] == dest_[0] );
 }


inline double Edge::Slope() const
 {
    if ( org_[0] != dest_[0] ) 
      return (dest_[1] - org_[1]) / (dest_[0] - org_[0]);
      
    return DBL_MAX;  
 }


inline void Edge::NormalizeTo( double n )
 {
     mjl::Point pn( dest_ - org_ );
     pn = (n / pn.Length()) * pn;
     org_.Set( 0., 0. );
     dest_ = pn; 
 }


/** for a given x value Y() returns the corresponding y value
on the line defined by the Edge*/
inline double Edge::Y( double x ) const
 {
    return Slope() * (x - org_[0]) + org_[1];
 }




inline  void  Edge::Move( double dx, double dy )
 {
    org_.Move( dx, dy );
    dest_.Move( dx, dy );
 }
 
 
inline  void  Edge::Scale( double xfac, double yfac )
 {
    org_.Scale( xfac, yfac );
    dest_.Scale( xfac, yfac );
 }




inline  double  Edge::Length() const
 {
     double dx = dest_[0] - org_[0];
     double dy = dest_[1] - org_[1];
     
     return std::hypot( dx, dy );
 }


}

#endif






















