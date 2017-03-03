#ifndef FEM_FROM_GRID_VISITOR_H
#define FEM_FROM_GRID_VISITOR_H

#include "CSMP_definitions.h"
#include "Visitor.h"
#include "ScalarVariable.h"
#include "ElementGrid.h"
#include "FiniteElement.h"

namespace csmp {

template<size_t> class PropertyDatabase;
class FiniteDifferenceGrid;
template<size_t> class DenseMatrix;
template<size_t> class Element;

template<size_t dim>
class FemFromGridVisitor : public Visitor<dim> {
  public:
    FemFromGridVisitor( const PropertyDatabase<dim>& p, 
                        const FiniteDifferenceGrid& g, 
                        const char* var, size_t elements );

    virtual ~FemFromGridVisitor();
    
    virtual void Visit( Element<dim>* );  // this includes constraint points as target 

  private:
    const PropertyDatabase<dim>&          pref;
    const FiniteDifferenceGrid&  grid;
    std::vector<ElementGrid>         egrids;
    DenseMatrix<DM_MIN>           XY, NN;
    ScalarVariable               val;
    csmp::Index                      key;

    bool IsInsideTriangle( double64 x, double64 y, bool update=true );
    bool IsInsideQuadrilateral( double64 x, double64 y );
    void MinMaxCoordinates( double64& min_x, double64& max_x, double64& min_y, double64& max_y );
    void InitializeElementGrid( size_t id, CSMP_FEM_TYPE fe_type );
};

/**
@class FemFromGridVisitor FemFromGridVisitor "main_library/FemFromGridVisitor.h"

@author S.K. Matthaei
@date 1998

class 'FemFromGridVisitor' uses arithmetic means of values
from supplied grid to assign these to variables placed on
the finite ELEMENTS.
*/

} // csmp


#endif







