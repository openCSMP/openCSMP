#ifndef FEM_TO_GRID_VISITOR_H
#define FEM_TO_GRID_VISITOR_H

#include "CSMP_definitions.h"
#include "Visitor.h"
#include "ElementGrid.h"
#include "FiniteElement.h"

namespace csmp {

template<size_t> class PropertyDatabase;
class FiniteDifferenceGrid;
template<size_t> class DenseMatrix;
template<size_t> class Element;

// works only in 2D, requires contiguous, unique element numbering
template<size_t dim>
class FemToGridVisitor : public Visitor<dim> {
  public:
    FemToGridVisitor( const PropertyDatabase<dim>&,
                      FiniteDifferenceGrid&,
                      const char* var, 
                      size_t elements );

    virtual ~FemToGridVisitor();
    
    virtual void Visit( Element<dim>* );   
    virtual void Visit( Region<dim>* );
  
    void        OutputProperty( const char* prop );
    const char* OutputProperty() const;
    
    bool OverWrite() const;
    void OverWrite( bool yes_or_no );
    void ChangeOutputProperty( const char* var );
    
    FemToGridVisitor& operator=( const FemToGridVisitor& );

  private:
    const PropertyDatabase<dim>&  pref;
    FiniteDifferenceGrid&         grid;
    std::vector<bool>             visited;
    std::vector<ElementGrid>      egrids;
    cEGridIterator                eit, it_end;
    std::vector<double64>         NF;
    std::vector<double64>         xy;
    std::vector<ScalarVariable >  P;
    ScalarVariable                sc;
    DenseMatrix<DM_MIN>           XY, NN;
    csmp::Index                   prop_key; 
    double64                      time_increment;   
    bool                          overwrite;

    bool IsInsideTriangle( double64 x, double64 y, bool update=true );
    bool IsInsideQuadrilateral( double64 x, double64 y );
    
    void MinMaxCoordinates( double64& min_x, double64& max_x, 
                            double64& min_y, double64& max_y );
                            
    void InitializeElementGrid( size_t id, CSMP_FEM_TYPE );
     
};
/**
@class FemToGridVisitor FemToGridVisitor "main_library/FemToGridVisitor.h"

@author S.K. Matthaei
@date 1997

class 'FemToGridVisitor' uses linear interpolation on a triangle
to calculate values of the input "var" on the supplied
regular grid.
*/


} // csmp    


#endif

















