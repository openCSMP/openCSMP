#ifndef FLUX_LHS_H
#define FLUX_LHS_H

#include "MatrixOperator.h"


namespace csmp {

template<size_t dim>
class FluxLHS : public MatrixOperator<dim> {
  public:
    FluxLHS(Model<dim>& model, const char*);
    FluxLHS( const FluxLHS&);
  
    virtual void AccumulateFiniteVolume( Node<dim>&, SparseMatrix& ) const;
    virtual void AccumulateStencil( Element<dim>&, SparseMatrix& ) const;

  private:
    const csmp::INDEX<SCALAR,FACET_INTEGRATION_POINT> ff_key_;
};

}

#endif /* FLUX_LHS_H */
