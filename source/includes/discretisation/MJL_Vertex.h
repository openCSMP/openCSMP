#ifndef MJL_VERTEX_H
#define MJL_VERTEX_H

#include "MJL_Point.h"
#include "MJL_Node.h"

namespace mjl {

enum ORIENTATION { CLOCKWISE, COUNTER_CLOCKWISE };

class Vertex : public Node, public mjl::Point {
  public:
    Vertex( double, double );
    explicit    Vertex( const mjl::Point& );
    Vertex*     Cw() const;
    Vertex*     Ccw() const;
    Vertex*     Neighbor( ORIENTATION ) const;
    mjl::Point  Point() const;
    Vertex*     Insert( Vertex* );
    Vertex*     Remove();
    void        Splice( Vertex* );
    Vertex*     Split( Vertex* );
    bool        IsConvex() const; // MJL, p. 198
};

bool adjacent( const Vertex*, const Vertex* );


inline bool adjacent( const Vertex* v, const Vertex* w )
 {
    return ( w==v->Cw() or w==v->Ccw() );
 }


// methods

inline Vertex::Vertex( double x, double y )
   : mjl::Point(x,y)
 {
 }


inline Vertex::Vertex( const mjl::Point& p )
   : mjl::Point(p)
 {
 }


inline Vertex* Vertex::Cw() const
 {
    return static_cast<Vertex*>(next_);
 }
 
 
inline Vertex* Vertex::Ccw() const
 {
    return static_cast<Vertex*>(prev_);
 }
 
 
inline Vertex* Vertex::Neighbor( ORIENTATION rotation ) const 
 {
    return ( (rotation == CLOCKWISE) ? Cw() : Ccw() );
 }
 
 
inline mjl::Point Vertex::Point() const
 {
    return static_cast<mjl::Point>(*this);
 }
 
 
inline Vertex* Vertex::Insert( Vertex* v )
 {
    return static_cast<Vertex*>(Node::Insert(v));
 }
 
 
inline Vertex*  Vertex::Remove()
 {
    return static_cast<Vertex*>(Node::Remove());
 }
 
 
inline void  Vertex::Splice( Vertex* b )
 {
    Node::Splice(b);
 }
 
 

} // end namespace mjl


#endif


