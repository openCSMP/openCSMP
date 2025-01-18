#ifndef FEM_TO_GRID_VISITOR_H
#define FEM_TO_GRID_VISITOR_H

#include "CSMP_definitions.h"
#include "Visitor.h"
#include "ElementGrid.h"
#include "FiniteElement.h"

namespace csmp {

template<uint32_t> class PropertyDatabase;
class FiniteDifferenceGrid;
template<uint32_t> class DenseMatrix;
template<uint32_t> class Element;

// works only in 2D, requires contiguous, unique element numbering
template<uint32_t dim>
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
    std::vector<double>         NF;
    std::vector<double>         xy;
    std::vector<ScalarVariable >  P;
    ScalarVariable                sc;
    DenseMatrix<DM_MIN>           XY, NN;
    csmp::Index                   prop_key; 
    double                      time_increment;   
    bool                          overwrite;

    bool IsInsideTriangle( double x, double y, bool update=true );
    bool IsInsideQuadrilateral( double x, double y );
    
    void MinMaxCoordinates( double& min_x, double& max_x, 
                            double& min_y, double& max_y );
                            
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

















