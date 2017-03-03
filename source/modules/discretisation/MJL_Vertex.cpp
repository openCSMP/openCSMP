#include "MJL_Vertex.h"

namespace mjl {

Vertex*  Vertex::Split( Vertex* b )
 {
    Vertex* bp = b->Ccw()->Insert( new Vertex( b->Point()) );
    Insert( new Vertex(mjl::Point()) );
    Splice( bp );
    
    return bp;
 }


bool  Vertex::IsConvex() const
 {
    Vertex*  u = Ccw();
    Vertex*  w = Cw();
    LOCATION c = w->Classify( *u, *w );
    return ( (c == BEYOND) or (c == RIGHT) );
 }

} // end namespace mjl
