// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "ElementGrid.h"

using namespace std;

namespace csmp {


EGridIterator  ElementGrid::Begin()
 {
    return grid.begin();
 }
 
 
EGridIterator  ElementGrid::End()
 {
    return grid.end();
 }




void ElementGrid::AddPoint( int32_t i, int32_t j, double val )
 {
    grid[ std::make_pair(i,j) ] = val;
 }




double ElementGrid::Value( int32_t i, int32_t j ) const
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


double ElementGrid::GridAverage() const
 {
    double avg = 0.0F;
    long     n   = 0;
    
    for ( cEGridIterator cit=grid.begin(); cit!=grid.end(); cit++, n++ )
      {
          avg += (*cit).second;
      }
    return avg / static_cast<double>(n);  
 }



double ElementGrid::Sum() const
 {
    double sum = 0.0F;
    
    for ( cEGridIterator cit=grid.begin(); cit!=grid.end(); cit++ ) sum += (*cit).second;

    return sum;  
 }



void ElementGrid::MinMax( double& dmin, double& dmax ) const
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
 
 
int32_t  ElementGrid::Size() const
 {
    return static_cast<int32_t>(grid.size());
 }


void ElementGrid::Out() const
 {
    pair<int32_t,int32_t>  index;
    
    cout <<"\nElementGrid::Out(): printing grid...\n" << endl;
    for ( cEGridIterator it=grid.begin(); it!=grid.end(); it++ )
      {
          index = (*it).first;
          cout <<"ElementGrid["<< index.first <<","<< index.second <<"] ";
          cout << (*it).second << endl;
      }
 }

} // end namespace csmp
