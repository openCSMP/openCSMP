#ifndef ELEMENT_GRID_H
#define ELEMENT_GRID_H

#include "CSMP_definitions.h"
#include "FiniteDifferenceGrid.h"

namespace csmp {

/**
    @brief mapping from triangular finite elements to regular grid 
    
    @copyright (c) SKM (1996)
*/
class ElementGrid {
  public:
    ElementGrid()  {};
    ElementGrid( const ElementGrid& );
    ~ElementGrid() {};
    ElementGrid& operator=( const ElementGrid& );
    void AddPoint( int32_t i, int32_t j, double val );
    std::map<std::pair<int32_t,int32_t>,double>::iterator  Begin();
    std::map<std::pair<int32_t,int32_t>,double>::iterator  End();
    double Value( int32_t i, int32_t j ) const;
    double GridAverage() const;
    double Sum() const;
    void     MinMax( double& dmin, double& dmax ) const;
    void     TransferDataToGrid( FiniteDifferenceGrid& ) const;
    bool     Empty() const;
    int32_t    Size() const;

    void     Out() const;

  private:
    std::map<std::pair<int32_t,int32_t>,double>  grid;
};
  
typedef std::map<std::pair<int32_t,int32_t>,double>::iterator  EGridIterator;
typedef std::map<std::pair<int32_t,int32_t>,double>::const_iterator  cEGridIterator;

} // csmp
  
#endif






