#include "MJL_functions.h"

using namespace std;

namespace mjl {

/// Laszlo, p. 87
bool pointInConvexPolygon( const mjl::Point& s, Polygon& p )
 {
    if ( p.Size() == 1U ) return (s == p.Point());
    if ( p.Size() == 2U ) {
         LOCATION c = s.Classify( p.Edge() );
         return( (c==BETWEEN) || (c==ORIGIN) || (c==DESTINATION) );
      }

    Vertex*  org = p.V();
    
    for ( unsigned int i=0U; i<p.Size(); i++, p.Advance(CLOCKWISE) )
      if ( s.Classify( p.Edge() ) == LEFT ) {
           p.SetV(org);
           return false;
        }
        
    return true;

 } // end pointInConvexPolygon





/// Ray shooting method for arbitrary polygons
POINT_CLASSIFICATION pointInPolygonRS( const mjl::Point a, const Polygon& p ) // MJL, p. 118
 {
    int parity(0);
    
    for ( unsigned int i=0; i<p.Size(); i++, p.Advance(CLOCKWISE) ) {
         Edge e = p.Edge();
         switch ( edgeType(a,e) ) {
              case TOUCHING: return BOUNDARY;
              case CROSSING: 
                parity = 1 - parity;
              default:
                cout <<"\nmjl::pointInPolygonRS: unhandled case INESSENTIAL in switch statement: "<< edgeType(a,e) << endl;
           }
      }
    
    return ( parity ? INSIDE : OUTSIDE );
      
 } // end pointInPolygonRS
 
 
 
 
/// signed angle method for arbitrary polygons
POINT_CLASSIFICATION pointInPolygonSA( const mjl::Point a, const Polygon& p ) // MJL, p. 120ff
 {
    double total(0.);
    
    for ( unsigned int i=0; i<p.Size(); i++, p.Advance(CLOCKWISE) ) {
         Edge e = p.Edge();
         double alpha = signedAngle( a, e );
         if ( alpha <= (180. + EPSILON1) and alpha >= (180. - EPSILON1) ) return BOUNDARY;
         total += alpha;
      }
    
    return ( (total < -180.) ? INSIDE : OUTSIDE );
    
 } // end pointInPolygonSA 



bool firstContainsSecondRS( const Polygon& first, const Polygon& second ) // SKM
 {
    for ( unsigned int i=0; i<second.Size(); i++, second.Advance(CLOCKWISE) ) {
         mjl::Point a = second.Point();
         if ( pointInPolygonRS( a, first ) != mjl::INSIDE ) return false; 
      }
      
    return true;
 }



/// Lazlo, p. 102 
INTERSECTION  lineTriangle3DIntersect( const Edge3D& e, const Triangle3D& p, double& t )
 {
     INTERSECTION  aclass = e.Intersect( p, t );
     
     if ( aclass == PARALLEL or aclass == COLLINEAR ) return aclass;
     
     int h, v;
     
     // choosing the projection plane (SKM fix of original code which gave
     // degenerate solutions in some cases) 
     if ( fabs(p.n().DotProduct(Point3D(0.,0.,1.))) > EPSILON1 ) { 
          h = 0; v = 1;
       }
     else if ( fabs(p.n().DotProduct(Point3D(1.,0.,0.))) > EPSILON1 ) {
          h = 1; v = 2;
       }
     else {
          h = 2; v = 0;
       }
 
     Polygon*  pp = project( p, h, v ); // pp->Out();
 
     mjl::Point3D  q(e.Point(t));
     
     bool answer = pointInConvexPolygon( mjl::Point( q[h], q[v] ), *pp ); 
     
     delete pp;
     
     return ( answer ? SKEW_CROSS : SKEW_NO_CROSS );
 
 } // end lineTriangle3DIntersect




/**

assumes that the projection of the triangle into the h,v plane
is non-degenerate, i.e. it is a triangle
@warning watch out! - the returned triangle pointer must be deleted
by the user of the function*/
Polygon*  project( const Triangle3D& p, int h, int v )
 {
    // project vertices of triangle p into the plane defined by h,v
    mjl::Point  pts[3] = { mjl::Point( p[0][h], p[0][v] ),
                           mjl::Point( p[1][h], p[1][v] ),
                           mjl::Point( p[2][h], p[2][v] ) };
                               
    // insert first two projected vertices into polygon
    Polygon*  pp = new Polygon();  
    pp->Insert( pts[0] );
    pp->Insert( pts[1] );
    
    // insert last projected vertex into polygon
    if ( pts[2].Classify( pts[0], pts[1] ) == LEFT )
    pp->Advance( CLOCKWISE );
    pp->Insert( pts[2] );
    
    return pp;
    
 } // end project






/// builds star-shaped polygon around vector which must contain its first point
Polygon* starShapedPolygon( const std::vector<mjl::Point>& pts ) // MJL, p. 110
 {
    if ( pts.empty() ) {
         cout <<"\nstarShapedPolygon: supplied point-list is empty; nothing was done."<< endl;
         return nullptr;
      }

    std::vector<mjl::Point>::const_iterator  it=pts.begin();
    
    Polygon* p = new Polygon();
    p->Insert( (*it++) );
    
    Vertex* origin = p->V();
    mjl::Point originPt = origin->Point();
    
    while ( it != pts.end() ) {
         p->SetV(origin);
         p->Advance(CLOCKWISE);
         while ( polarCmp( originPt, (*it), (*p->V()) ) < 0 )
           p->Advance(CLOCKWISE);
         p->Advance(COUNTER_CLOCKWISE);
         p->Insert( (*it) );
         it++;
      }
    
    return p;    
 
 } // end starShapedPolygon





mjl::Point somePoint;

Polygon* insertionHull( const std::vector<mjl::Point>& s ) // MJL, p. 114
 {
    extern mjl::Point somePoint;
    Polygon* p = new Polygon();
    p->Insert( s[0] );
    
    for ( unsigned int i=1U; i<s.size(); i++ ) {
         if ( pointInConvexPolygon( s[i], *p ) ) continue;
         somePoint = s[i];
         p->LeastVertex( closestToPolygonCmp );
         supportingLine( s[i], p, mjl::LEFT );
         Vertex *l = p->V();
         supportingLine( s[i], p, mjl::RIGHT );
         delete p->Split( l );
         p->Insert( s[i] );
      } 
      
    return p;
 
 } // end insertionHull




/// Laszlo, p. 158
Polygon*  convexPolygonIntersect( Polygon& P, Polygon& Q )
 {
    Polygon*               R(nullptr);
    mjl::Point             iPnt, startPnt;
    POLYGON_CLASSIFICATION inflag(UNKNOWN);
    int                    phase(1);
    size_t                 maxItns = 2 * (P.Size() + Q.Size());
    
    for ( int i=1; (i<=static_cast<int>(maxItns)) or (phase==2); i++ ) {
         Edge p = P.Edge();
         Edge q = Q.Edge();
         LOCATION pclass = p.dest_.Classify(q);
         LOCATION qclass = q.dest_.Classify(p);
         INTERSECTION crossType = crossingPoint( EPSILON1, p, q, iPnt );
         if ( crossType == SKEW_CROSS ) {
              if ( phase == 1 ) {
                   phase = 2;
                   R = new Polygon();
                   R->Insert( iPnt );
                   startPnt = iPnt;
                }
              else if ( iPnt != R->Point() ) {
                   if ( iPnt != startPnt ) R->Insert( iPnt );
                   else return R;
                }
              if      ( pclass == RIGHT ) inflag = P_IS_INSIDE;
              else if ( qclass == RIGHT ) inflag = Q_IS_INSIDE;
              else inflag = UNKNOWN;
           }
         else if ( crossType == COLLINEAR and pclass != BEHIND and qclass != BEHIND )
           inflag = UNKNOWN;
         bool pAIMSq = aimsAt( p, q, pclass, crossType );
         bool qAIMSp = aimsAt( q, p, qclass, crossType );
         if ( pAIMSq and qAIMSp ) {
              if ( inflag == Q_IS_INSIDE or (inflag == UNKNOWN and pclass == LEFT) )
                advance( P, *R, false );
              else
                advance( Q, *R, false );
           }
         else if ( pAIMSq ) advance( P, *R, inflag == P_IS_INSIDE );
         else if ( qAIMSp ) advance( Q, *R, inflag == Q_IS_INSIDE );
         else {
              if ( inflag == Q_IS_INSIDE or (inflag == UNKNOWN and pclass == LEFT) )
                advance( P, *R, false );
              else
                advance( Q, *R, false );
           }
      } // end for
      
    if ( pointInConvexPolygon( P.Point(), Q ) ) return new Polygon( P );
    else if ( pointInConvexPolygon( Q.Point(), P ) ) return new Polygon( Q );
    return new Polygon();
 
 } // end convexPolygonIntersect   

 
 
 
 

enum { UPPER, LOWER };
 
 
/** If a and b intersect eachother their union is return as potentially non-convex
polygon into res. 
 */
Polygon*  merge( Polygon& L, Polygon& R )
 {
    Vertex *l1(0), *r1(0), *l2(0), *r2(0);
    Vertex* vl = L.LeastVertex( rightToLeftCmp );
    Vertex* vr = R.LeastVertex( leftToRightCmp );
    bridge( L, R, l1, r1, UPPER );
    L.SetV(vl);
    R.SetV(vr);
    bridge( L, R, l2, r2, LOWER );
    L.SetV(l1);
    L.Split(r1);
    R.SetV(r2);
    delete R.Split(l2);
    return &R;
 }
 
 
 void bridge( Polygon& L, Polygon& R, Vertex* vl, Vertex* vr, int type )
  {
     LOCATION sides[2] = { LEFT, RIGHT };
     int      indx = (type == UPPER) ? 0 : 1;
     do {
          vl = L.V();
          vr = R.V();
          supportingLine( L.Point(), &R, sides[indx] );
          supportingLine( R.Point(), &L, sides[1-indx] );
       } 
     while ( vl != L.V() or vr != R.V() );
     
  } // end bridge   

 
 
void supportingLine( const mjl::Point& s, Polygon* p, LOCATION side )
 {
    ORIENTATION  rotation = (side==LEFT) ? CLOCKWISE : COUNTER_CLOCKWISE;
    Vertex*      a = p->V();
    Vertex*      b = p->Neighbor( rotation );
    LOCATION     c = b->Classify( s, *a );
    
    while( c == side or c == BEYOND or c == BETWEEN ) {
         p->Advance( rotation );
         a = p->V();
         b = p->Neighbor( rotation );
         c = b->Classify( s, *a );
      }
 }




bool clipPolygonToEdge( const Polygon& s, const Edge& e, Polygon& result ) // MJL, p. 127
 {
    Polygon*    p = &result;
    mjl::Point  crossingPt;
    
    for ( int i=0; i<static_cast<int>(s.Size()); s.Advance(CLOCKWISE), i++ ) {
         mjl::Point org  = s.Point();
         mjl::Point dest = s.Cw()->Point();
         bool orgIsInside  = ( org.Classify(e)  != mjl::LEFT );
         bool destIsInside = ( dest.Classify(e) != mjl::LEFT );
         if ( orgIsInside != destIsInside ) {
              double t;
              e.Intersect( s.Edge(), t );
              crossingPt = e.Point( t );
           }
         // case 1
         if ( orgIsInside and destIsInside ) p->Insert( dest );
         // case 2
         else if ( orgIsInside and !destIsInside ) {
              if ( org != crossingPt ) p->Insert( crossingPt );
           }
         // case 3
         else if ( !orgIsInside and !destIsInside )
           ;
         // case 4
         else {
              p->Insert( crossingPt );
              if ( dest != crossingPt ) p->Insert( dest );
           } 
      }
    
    return ( p->Size() > 0u );  
      
 } // end 
 
 

/**

the supplied edge is clipped where it extends beyond the polygon, p
if e does not intersect p, the method returns false
variation of clipPolygonToEdge(), MJL, p. 127

@warning works only for convex polygons which yield exactly two
crossing points.
*/
bool clipEdgeToPolygon( const Polygon& s, Edge& e )
 {
    mjl::Point crossingPt[2u];
    int   crosses(0);
    
    for ( int i=0; i<static_cast<int>(s.Size()); s.Advance(CLOCKWISE), i++ ) {
         mjl::Point org  = s.Point();
         mjl::Point dest = s.Cw()->Point();
         bool  orgIsInside  = ( org.Classify(e)  != mjl::LEFT );
         bool  destIsInside = ( dest.Classify(e) != mjl::LEFT );
         if ( orgIsInside != destIsInside ) {
              double t;
              e.Intersect( s.Edge(), t );
              crossingPt[crosses++] = e.Point( t );
           }
         if ( crosses >= 2 ) {
              e.Set( crossingPt[0], crossingPt[1] );
              return true;
           }
      }
    
    return false;
      
 } // end clipEdgeToPolygon



} // end namespace csp
















