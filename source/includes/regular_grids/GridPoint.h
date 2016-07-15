#ifndef CSMP_GRID_POINT_H
#define CSMP_GRID_POINT_H

#include "CSMP_definitions.h"

namespace csmp {

class GridPoint {
  public:
    GridPoint() : i(0), j(0), data(1) { data[0]=std::numeric_limits<double64>::quiet_NaN(); };
    GridPoint( size_t in, size_t jn, std::vector<double64>& d )
      : i(in), j(jn), data(d) {};
    GridPoint( size_t in, size_t jn, size_t n )
      : i(in), j(jn), data(n) {};
    GridPoint( size_t in, size_t jn, double64 d )
      : i(in), j(jn), data(1) { data[0]=d; };
    ~GridPoint() {};
    GridPoint( const GridPoint& gp ) { *this=gp; };
    GridPoint& operator=( const GridPoint& gp )
      {
         if ( &gp != this ) {
              i=gp.i;
              j=gp.j;
              data=gp.data;
           }
         return *this;
      }
    size_t   I()                  const { return i; };
    size_t   J()                  const { return j; };
    double64 Value( size_t i )      const { return data[i]; };
    void   I( size_t i_inp )            { i=i_inp; };
    void   J( size_t j_inp )            { j=j_inp; };
    double64 operator[]( size_t i )  const { return data[i]; };
    size_t Size()               const { return data.size(); };
    void   Value( size_t i, double64 val )
      {
         if ( i < data.size() ) data[i] = val;
         else
           {
              data.reserve( i+1 );
              for ( size_t j=data.size()-1; j<=i; j++ ) data.push_back( val );
           }
      };
      
  private:
    size_t   i, j;
    std::vector<double64> data;

};

} // csmp
// ------------------------------------------------------------------

#endif /* GridPoint.h */


