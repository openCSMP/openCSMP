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
    void AddPoint( int32 i, int32 j, double64 val );
    std::map<std::pair<int32,int32>,double64>::iterator  Begin();
    std::map<std::pair<int32,int32>,double64>::iterator  End();
    double64 Value( int32 i, int32 j ) const;
    double64 GridAverage() const;
    double64 Sum() const;
<<<<<<< HEAD
    void     MinMax( double64& dmin, double64& dmax ) const;
    void     TransferDataToGrid( FiniteDifferenceGrid& ) const;
    bool     Empty() const;
    int32    Size() const;

    void     Out() const;
=======
    void  MinMax( double64& dmin, double64& dmax ) const;
    void Out() const { Out(std::cout); }
    void  Out(std::ostream& os) const;
    void  TransferDataToGrid( FiniteDifferenceGrid& ) const;
    bool  Empty() const;
    int32  Size() const;
>>>>>>> b71117cb444b1539e747fa5855ab341062d3c4b3

  private:
    std::map<std::pair<int32,int32>,double64>  grid;
};
  
typedef std::map<std::pair<int32,int32>,double64>::iterator  EGridIterator;
typedef std::map<std::pair<int32,int32>,double64>::const_iterator  cEGridIterator;

} // csmp
  
#endif






