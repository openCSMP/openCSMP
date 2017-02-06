#ifndef ELEMENT_GRID_H
#define ELEMENT_GRID_H

#include "CSMP_definitions.h"
#include "FiniteDifferenceGrid.h"

namespace csmp {

/// mapping from triangular finite elements to regular grid (copyright (c) SKM 1996)
class ElementGrid {
  public:
    ElementGrid()  {};
    ElementGrid( const ElementGrid& );
    ~ElementGrid() {};
    ElementGrid& operator=( const ElementGrid& );
    void AddPoint( int32 i, int32 j, double64 val );
    std::map<std::pair<int32,int32>,double64>::iterator  Begin();
    std::map<std::pair<int32,int32>,double64>::iterator  End();
    double64 Value( int32 i, int32 j ) const;
    double64 GridAverage() const;
    double64 Sum() const;
    void  MinMax( double64& dmin, double64& dmax ) const;
    void  Out(std::ostream& os) const;
    void  TransferDataToGrid( FiniteDifferenceGrid& ) const;
    bool  Empty() const;
    int32  Size() const;

  private:
    std::map<std::pair<int32,int32>,double64>  grid;
};
  
typedef std::map<std::pair<int32,int32>,double64>::iterator  EGridIterator;
typedef std::map<std::pair<int32,int32>,double64>::const_iterator  cEGridIterator;

  

inline  EGridIterator  ElementGrid::Begin()
 {
    return grid.begin();
 }
 
 
inline  EGridIterator  ElementGrid::End()
 {
    return grid.end();
 }




inline  void ElementGrid::AddPoint( int32 i, int32 j, double64 val )
 {
    grid[ std::make_pair(i,j) ] = val;
 }




inline double64 ElementGrid::Value( int32 i, int32 j ) const 
 {
    assert( grid.count( std::make_pair(i,j) ) > 0 );
    return (*grid.find( std::make_pair(i,j) )).second;
 }




inline  ElementGrid& ElementGrid::operator=( const ElementGrid& g )
 {
    if ( &g != this ) grid = g.grid;
    return *this;
 }




inline  ElementGrid::ElementGrid( const ElementGrid& g )
 : grid(g.grid)
 {
 }


inline  double64 ElementGrid::GridAverage() const
 {
    double64 avg = 0.0F;
    long     n   = 0;
    
    for ( cEGridIterator cit=grid.begin(); cit!=grid.end(); cit++, n++ )
      {
          avg += (*cit).second;
      }
    return avg / static_cast<double64>(n);  
 }



inline  double64 ElementGrid::Sum() const
 {
    double64 sum = 0.0F;
    
    for ( cEGridIterator cit=grid.begin(); cit!=grid.end(); cit++ ) sum += (*cit).second;

    return sum;  
 }



inline  void ElementGrid::MinMax( double64& dmin, double64& dmax ) const
 {
    cEGridIterator  cit(grid.begin());
    dmin = dmax = (*cit).second;
    
    while ( cit!=grid.end() )
      {
         if ( (*cit).second > dmax ) dmax = (*cit).second;
         if ( (*cit).second < dmin ) dmin = (*cit).second;
         cit++;
      }
 }

inline  void ElementGrid::TransferDataToGrid( FiniteDifferenceGrid& g ) const
 { 
    for ( cEGridIterator cit=grid.begin(); cit!=grid.end(); cit++ )
      //  i=vertical=g.row   j=horizontal=g.col 
      g( (*cit).first.first, (*cit).first.second ) = (*cit).second;
 }


inline bool  ElementGrid::Empty() const
 {
    return grid.empty();
 }
 
 
inline int32  ElementGrid::Size() const
 {
    return static_cast<int32>(grid.size());
 }

} // csmp
  
#endif






