#include "ElementGrid.h"

using namespace std;

namespace csmp {

<<<<<<< HEAD

EGridIterator  ElementGrid::Begin()
 {
    return grid.begin();
 }
 
 
EGridIterator  ElementGrid::End()
 {
    return grid.end();
 }




void ElementGrid::AddPoint( int32 i, int32 j, double64 val )
 {
    grid[ std::make_pair(i,j) ] = val;
 }




double64 ElementGrid::Value( int32 i, int32 j ) const
 {
    assert( grid.count( std::make_pair(i,j) ) > 0 );
    return (*grid.find( std::make_pair(i,j) )).second;
 }




ElementGrid& ElementGrid::operator=( const ElementGrid& g )
 {
    if ( &g != this ) grid = g.grid;
    return *this;
 }




ElementGrid::ElementGrid( const ElementGrid& g )
 : grid(g.grid)
 {
 }


double64 ElementGrid::GridAverage() const
 {
    double64 avg = 0.0F;
    long     n   = 0;
    
    for ( cEGridIterator cit=grid.begin(); cit!=grid.end(); cit++, n++ )
      {
          avg += (*cit).second;
      }
    return avg / static_cast<double64>(n);  
 }



double64 ElementGrid::Sum() const
 {
    double64 sum = 0.0F;
    
    for ( cEGridIterator cit=grid.begin(); cit!=grid.end(); cit++ ) sum += (*cit).second;

    return sum;  
 }



void ElementGrid::MinMax( double64& dmin, double64& dmax ) const
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

void ElementGrid::TransferDataToGrid( FiniteDifferenceGrid& g ) const
 { 
    for ( cEGridIterator cit=grid.begin(); cit!=grid.end(); cit++ )
      //  i=vertical=g.row   j=horizontal=g.col 
      g( (*cit).first.first, (*cit).first.second ) = (*cit).second;
 }


bool  ElementGrid::Empty() const
 {
    return grid.empty();
 }
 
 
int32  ElementGrid::Size() const
 {
    return static_cast<int32>(grid.size());
 }


void ElementGrid::Out() const
=======
void ElementGrid::Out(std::ostream& os) const
>>>>>>> b71117cb444b1539e747fa5855ab341062d3c4b3
 {
    pair<int32,int32>  index;
    
    os <<"\nElementGrid::Out(): printing grid...\n" << endl;
    for ( cEGridIterator it=grid.begin(); it!=grid.end(); it++ )
      {
          index = (*it).first;
          os <<"ElementGrid["<< index.first <<","<< index.second <<"] ";
          os << (*it).second << endl;
      }
 }

} // end namespace csmp
