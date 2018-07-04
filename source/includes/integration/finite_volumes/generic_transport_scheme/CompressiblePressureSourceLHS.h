#ifndef COMPRESSIBLE_PRESSURE_SOURCE_LHS_H
#define COMPRESSIBLE_PRESSURE_SOURCE_LHS_H

#include "MatrixOperator.h"


namespace csmp {

template<size_t dim>
class CompressiblePressureSourceLHS : public MatrixOperator<dim> {
  public:
    CompressiblePressureSourceLHS(Model<dim>& model, const char*, const char*, const char*, const char*);
    CompressiblePressureSourceLHS( const CompressiblePressureSourceLHS&);
  
    virtual void AccumulateFiniteVolume( Node<dim>&, SparseMatrix& ) const;
    virtual void AccumulateStencil( Element<dim>&, SparseMatrix& ) const;

  private:
    const csmp::INDEX<SCALAR,ELEMENT> key_PHI; // porosity
    const csmp::INDEX<SCALAR,NODE> key_CT;
    const csmp::INDEX<SCALAR,NODE> key_PF0;
    const csmp::INDEX<SCALAR,NODE> key_PF1;
};

}

#endif /* COMPRESSIBLE_PRESSURE_SOURCE_LHS_H */
