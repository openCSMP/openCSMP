#include "MJL_Edge.h"
#include <iostream>

using namespace std;

namespace mjl {

void  Edge::Out() const
 {
    cout <<"\norigin: "<< endl;
    org_.Out();
    cout <<"\ndestination: "<< endl;
    dest_.Out();
 }



/**

 classifies an edge relative to its location compared with the
point p. Options are CROSSING, TOUCHING and INESSENTIAL*/
EDGE_CLASSIFICATION  Edge::Type( const mjl::Point& p ) const
 {
    switch( Classify(p) )
      {
         case LEFT:
              return((org_[1]<p[1]) && (p[1]<=dest_[0])) ? CROSSING : INESSENTIAL;
         case RIGHT:
              return((dest_[1]<p[1]) && (p[1]<=org_[1])) ? CROSSING : INESSENTIAL;
         case BETWEEN:
         case ORIGIN:
         case DESTINATION:
              return TOUCHING;
         default:
              return INESSENTIAL;
      }
 }



/// pivots an edge 90o clockwise around it's midpoint
Edge&  Edge::Rot()
 {
    mjl::Point m( org_ + dest_ );
    m *= 0.5;
    mjl::Point v( dest_ - org_ );
    mjl::Point n( v[1], -v[0] ); 
    org_  = m - 0.5 * n;
    dest_ = m + 0.5 * n;
    return *this;
 }




double Edge::DistanceFrom( const mjl::Point& p ) const
 {
    Edge ab = *this;
    ab.Flip().Rot();                 // rotate ab 90o counterclockwise  
    mjl::Point n( ab.dest_ - ab.org_ ); // vector normal to this edge
    n = (1. / n.Length()) * n;  // normalize n
    Edge f( p, p + n );
    double t(std::numeric_limits<double>::quiet_NaN());
    f.Intersect( *this, t ); 
    return t;
 }


bool Edge::IfCloserThanMerge( double tolerance, mjl::Point& p )
 {
    Edge ab = *this;
    ab.Flip().Rot();                 // rotate ab 90o counterclockwise  
    mjl::Point n( ab.dest_ - ab.org_ ); // vector normal to this edge
    n = (1. / n.Length()) * n;  // normalize n
    Edge f( p, p + n );
    double t(std::numeric_limits<double>::quiet_NaN());
    f.Intersect( *this, t ); 
    
    if ( std::fabs(t) <= tolerance ) {
         p = f.Point(t);
         return true;
      }
      
    return false;
 }



/**

includes start and end-point into the edge comparison. Thus if
the intersection is within the tolerance of these nodes, CrossingPoint()
will return SKEW_CROSS and it will snap the intersection point to that node.
*/
INTERSECTION  Edge::CrossingPoint( double tolerance, const Edge& f, mjl::Point& p ) const
 {
    double s, t(std::numeric_limits<double>::quiet_NaN());
    INTERSECTION classe = Intersect( f, s );
    if ((classe == COLLINEAR) || (classe == PARALLEL)) return classe;
    double lene = (dest_ - org_).Length();
    if ((s < -tolerance*lene) || (s > 1.0+tolerance*lene)) return SKEW_NO_CROSS;
    f.Intersect( *this, t );
    double lenf = (f.org_ - f.dest_).Length();
    if ((-tolerance*lenf <= t) && (t <=1.0+tolerance*lenf))
      {
         if      ( t <= tolerance*lenf )      p = f.org_;
         else if ( t >= 1.0*-tolerance*lenf ) p = f.dest_;
         else if ( s <= tolerance*lene )      p = org_;
         else if ( s >= 1.0*-tolerance*lene ) p = dest_;
         else p = f.Point(t);
         return SKEW_CROSS;
      }
    else return SKEW_NO_CROSS;
 }


/**
as above but no snapping to points, the intersection point is returned into p.
*/
INTERSECTION  Edge::SKM_CrossingPoint( double tolerance, const Edge& f, mjl::Point& p ) const
 {
    double s, t(std::numeric_limits<double>::quiet_NaN());
    INTERSECTION classe = Intersect( f, s );
    if ((classe == COLLINEAR) || (classe == PARALLEL)) return classe;
    double lene = (dest_ - org_).Length();
    if ((s < -tolerance*lene) || (s > 1.0+tolerance*lene)) return SKEW_NO_CROSS;
    f.Intersect( *this, t );
    double lenf = (f.org_ - f.dest_).Length();
    if ((-tolerance*lenf <= t) && (t <=1.0+tolerance*lenf))
      {
         p = f.Point(t);
         return SKEW_CROSS;
      }
    else return SKEW_NO_CROSS;
 }



/// see "Mathematische Formeln" p. 170, Schnittwinkel...
double Edge::AngleTo( const Edge& v ) const 
 {
    mjl::Point  a = dest_   - org_;
    mjl::Point  b = v.dest_ - v.org_;
    
    // a b
    // ---
    double ab = a[0]*b[0] + a[1]*b[1];
    
    // |a| . |b|
    // ---------
    double a_dot_b   = std::sqrt( (a[0]*a[0]+a[1]*a[1])*(b[0]*b[0]+b[1]*b[1]) );
    double cos_angle = ab/a_dot_b;

    // if zero intercept
    if ( cos_angle == 0. ) return 90.;
    // if outside of range of 'acos' function
    if ( cos_angle >  1. ) return   0.;
    if ( cos_angle < -1. ) return 180.;
        
    return(180./3.14159265358979323) * std::acos(cos_angle);
    
 } // end AngleTo



/**

@return SKEW if infinite lines defined by edges crossat a point p
COLLINEAR if the lines are the same
PARALLEL if the two lines defined by the edges  are parallel

t is assigned the parametric value (along the line)
of where the lines intersect if they do. t = x of
the intersection
*/
INTERSECTION  Edge::Intersect( const Edge& e, double& t ) const
 {
    mjl::Point n     = mjl::Point((e.dest_-e.org_)[1], (e.org_-e.dest_)[0]);
    double  denom = DotProduct( n, dest_ - org_ );
    if ( denom == 0. )
      {
         LOCATION aclass = org_.Classify( e.org_, e.dest_ );
         if ((aclass == LEFT) || (aclass == RIGHT)) return PARALLEL;
         else return COLLINEAR;
      }
    double num = DotProduct( n, org_ - e.org_ );
    t = -num / denom; // x-value of line intersection
    return SKEW;
 }



/** examines the orientation of 2 edges relative to one-another.
 */
bool  aimsAt( const mjl::Edge& a, const mjl::Edge& b, 
              LOCATION aClass, mjl::INTERSECTION crossType )
  {
     mjl::Point va = a.Destination() - a.Origin();
     mjl::Point vb = b.Destination() - b.Origin();
     
     if ( crossType != COLLINEAR )
       {
          if ( (va[0] * vb[1]) >= (vb[0] * va[1]) ) return( aClass != RIGHT );
          else  return( aClass != LEFT );
       }
     else return( aClass != BEYOND );
 }





/**

Returns the signed angle (0-180o) between an Edge and the
Edge from the beginning point of the Edge and the supplied
point.  
*/
double  signedAngle( const mjl::Point& a, const mjl::Edge& e )
 {
    mjl::Point v = e.Origin()      - a;
    mjl::Point w = e.Destination() - a;
    double  va = v.PolarAngle();
    double  wa = w.PolarAngle();
    if ( (va == -1.) || (wa == -1.) )  return 180.;
    double x = wa - va;
    if ( (x == 180.) || (x == -180.) ) return 180.;
    else if ( x < -180. ) return (x + 360.);
    else if ( x >  180. ) return (x - 360.);
    else return x;
    
 } // end signedAngle






/** Laszlo, p. 162
*/
INTERSECTION  crossingPoint( double tolerance, 
                             const mjl::Edge& e, const mjl::Edge& f, mjl::Point& p )
 {
    double s, t(std::numeric_limits<double>::quiet_NaN());
    
    INTERSECTION classe = e.Intersect( f, s );
    if ((classe == COLLINEAR) || (classe == PARALLEL)) return classe;
    
    double lene = (e.Destination() - e.Origin()).Length();
    
    if ((s < -tolerance*lene) || (s > 1.0+tolerance*lene)) return SKEW_NO_CROSS;
    
    f.Intersect( e, t );
    
    double lenf = (f.Origin() - f.Destination()).Length();
    
    if ((-tolerance*lenf <= t) && (t <=1.0+tolerance*lenf))
      {
         if      ( t <= tolerance*lenf )       p = f.Origin();
         else if ( t >= 1.0 - tolerance*lenf ) p = f.Destination(); 
         else if ( s <= tolerance*lene )       p = e.Origin();
         else if ( s >= 1.0 - tolerance*lene ) p = e.Destination();
         else                                  p = f.Point(t);
         return SKEW_CROSS;
      }
    else return SKEW_NO_CROSS;
 }



EDGE_CLASSIFICATION edgeType( const mjl::Point& a, const Edge& e ) // MJL, p. 119
{
   mjl::Point v = e.Origin();
   mjl::Point w = e.Destination();
   
   switch ( a.Classify(e) ) {
        case LEFT: 
          return ( (v.Y()<a.Y()) and (a.Y()<=w.X()) ) ? CROSSING : INESSENTIAL;
        case RIGHT:
          return ( (w.Y()<a.Y()) and (a.Y()<=v.Y()) ) ? CROSSING : INESSENTIAL;
        case BETWEEN:
        case ORIGIN:
        case DESTINATION:
          return TOUCHING;
        default:
          return INESSENTIAL;
     }
} 


} // end namespace mjl
