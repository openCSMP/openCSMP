#ifndef MJL_POLYGON_H
#define MJL_POLYGON_H

#include "MJL_Vertex.h"
#include "MJL_Edge.h"

#include <list>

namespace mjl {

class Polygon {
  public: // points // should be unique and ordered clockwise
    Polygon( std::list<mjl::Point>::const_iterator, 
             std::list<mjl::Point>::const_iterator );  
    explicit Polygon( mjl::Vertex* );
    Polygon();
    Polygon( const Polygon& );
    virtual ~Polygon();
    Polygon&      operator=( const Polygon& p ); 
    mjl::Vertex*  SetV( mjl::Vertex* );
    mjl::Vertex*  V() const;
    size_t        Size() const;
    bool          Empty() const;
    mjl::Point    Point() const;
    mjl::Edge     Edge() const;
    mjl::Vertex*  Cw() const;
    mjl::Vertex*  Ccw() const;
    mjl::Vertex*  Neighbor( ORIENTATION );
    mjl::Vertex*  Advance( ORIENTATION );
    const mjl::Vertex* Advance( ORIENTATION ) const;
    mjl::Vertex*  Insert( const mjl::Point& );
    void          Erase();
    void          Remove();
    void          Revert(); ///< rebuilds itself in reverse order
    Polygon*      Split( mjl::Vertex* );
                  // MJL, p.87
    mjl::Vertex*  LeastVertex( int (*cmp )( const mjl::Point& a, const mjl::Point& b ) );
    mjl::Point    CenterOfGravity() const;
    void          CenterOfGravity( double&, double& ) const;
    double        Perimeter() const;
    bool          CheckAngles( double tolerated_bound ) const;
    mjl::Point    InsidePointSA(); ///< uses signed angle > and returns edge midpoint 
    void          Scale( double factor );
    void          Move( double dx, double dy );
    void          BoundingRectangle( mjl::Point& cnr_min, mjl::Point& cnr_max ) const;

    bool          IsPositivelyOriented() const;
    void          OrientPositively();

    void          Out() const;
    void          Out( const char* ) const;
    void          OutputCoordinatesTo( std::list<mjl::Point>& ) const;

  private:
    mutable mjl::Vertex* v_;
    size_t  size_;
    
    void Resize(); 
};      

// nonmember functions
int leftToRightCmp( const mjl::Point& a, const mjl::Point& b );
int rightToLeftCmp( const mjl::Point& a, const mjl::Point& b );
int closestToPolygonCmp( const mjl::Point& a, const mjl::Point& b );



// member functions

inline Polygon::Polygon()
 : v_(0), size_(0u)
  {
  }

  
inline Polygon::Polygon( const Polygon& p )
 {
    *this = p;
 }



inline Polygon::Polygon( Vertex* v )
 : v_(v)
  {
     Resize();
  }


inline mjl::Vertex* Polygon::V() const
 {
    return v_;
 }
 

inline size_t Polygon::Size() const
 { return size_; }
 
 
inline bool  Polygon::Empty() const 
 { return (size_==0u); }

 
inline mjl::Point Polygon::Point() const
 {
    return v_->Point();
 } 
 

inline Edge  Polygon::Edge() const
 {
    return mjl::Edge( Point(), v_->Cw()->Point() );
 }
 
 
/// @return vertex successor
inline mjl::Vertex* Polygon::Cw() const
 {
    return v_->Cw();
 }


/// @return predecessor of vertex
inline mjl::Vertex* Polygon::Ccw() const
 {
    return v_->Ccw();
 }
 
 
inline mjl::Vertex* Polygon::Neighbor( ORIENTATION rotation )
 {
    return v_->Neighbor( rotation );
 }
 
 
inline mjl::Vertex* Polygon::Advance( ORIENTATION rotation )
 {
    return v_ = v_->Neighbor( rotation );
 }

inline const mjl::Vertex* Polygon::Advance( ORIENTATION rotation ) const
 {
    return v_ = v_->Neighbor( rotation );
 }



inline mjl::Vertex* Polygon::SetV( Vertex* v )
 {
    return v_ = v;
 }



inline mjl::Vertex* Polygon::Insert( const mjl::Point& p )
 {
    if ( size_++ == 0U ) v_ = new mjl::Vertex(p);
    else v_ = v_->Insert( new mjl::Vertex(p) );
    
    return v_;
 }


inline void Polygon::Remove()
 {
    mjl::Vertex* v = v_;
    v_ = (--size_ == 0U) ? 0 : v_->Ccw();
    delete v->Remove(); 
 }
 
 

inline Polygon*  Polygon::Split( mjl::Vertex* b )
 {
    mjl::Vertex* bp = v_->Split(b);
    Resize();
    return new Polygon(bp);
 }
 



// old functionality recaptured


inline mjl::Point Polygon::CenterOfGravity() const
 {
     double x, y;
     CenterOfGravity( x, y );
     assert( !isnan(x) );
     assert( !isnan(y) );
     return mjl::Point( x, y );
 }





// nonmember functions

inline int leftToRightCmp( const mjl::Point& a, const mjl::Point& b )
 {
    if ( a < b ) return -1;
    if ( a > b ) return  1;
    return 0;
 }

inline int rightToLeftCmp( const mjl::Point& a, const mjl::Point& b )
 {
    return leftToRightCmp( b, a );
 }


inline int closestToPolygonCmp( const mjl::Point& a, const mjl::Point& b )
 {
     extern mjl::Point somePoint;
     double distA = (somePoint - a).Length();
     double distB = (somePoint - b).Length();
     if ( distA < distB )     return -1;
     else if ( distA > distB ) return 1;
     return 0;
 }




} 

#endif














