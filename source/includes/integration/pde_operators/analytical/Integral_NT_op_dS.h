#ifndef INTEGRAL_NT_OP_N_DS_H
#define INTEGRAL_NT_OP_N_DS_H

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"

namespace csmp {
/**
@author S.K. Matthaei
@author S. Geiger
@date 2001 */

/// works only if there are faces
template<typename fT,size_t dim>
class Integral_NT_op_dS : public MathOperatorRHS<dim> {
  public:
    Integral_NT_op_dS( const PropertyDatabase<dim>& p, const char* oper, const char* test );
    
    virtual void GetOperands( const SIMPLEX& );

    virtual void ComputeContribution( const SIMPLEX& );
  
  private:
    ScalarVariable          sc;
    VectorVariable<dim>     vc, vc2, pvc;
    std::vector<double>   fn;
    double                face_value, edge;
    bool                    faces;
    SG_BOUNDARY             at_boundary;
};


/* copyright (c) 2001 by Dr. Stephan K. Matthaei & Sebastian Geiger */

} // csmp

#endif
















