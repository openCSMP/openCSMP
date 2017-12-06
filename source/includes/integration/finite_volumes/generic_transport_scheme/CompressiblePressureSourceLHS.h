#ifndef COMPRESSIBLE_PRESSURE_SOURCE_LHS_H
#define COMPRESSIBLE_PRESSURE_SOURCE_LHS_H

#include "MatrixOperator.h"


namespace csmp {

template<size_t dim>
class CompressiblePressureSourceLHS : public MatrixOperator<dim> {
  public:
    CompressiblePressureSourceLHS(Model<dim>& model, const char*, const char*, const char*, const char*);
    CompressiblePressureSourceLHS( const CompressiblePressureSourceLHS&);
  
    virtual void AccumulateFiniteVolume( Node<dim>*, SparseMatrix& );
    virtual void AccumulateStencil( Element<dim>*, SparseMatrix& );

  private:
    const csmp::Index phi_key_, CT_key_, pf0_key_, pf1_key_;
};

}

#endif /* COMPRESSIBLE_PRESSURE_SOURCE_LHS_H */
