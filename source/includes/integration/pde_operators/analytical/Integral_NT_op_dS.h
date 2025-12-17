#ifndef INTEGRAL_NT_OP_N_DS_H
#define INTEGRAL_NT_OP_N_DS_H

#include "MathOperatorRHS.h"

namespace csmp {

template<uint32_t> class Face;

/**
@author S.K. Matthaei
@author S. Geiger
@date 2001 */

/// works only if there are faces
template<uint32_t dim>
class Integral_NT_op_dS : public MathOperatorRHS<dim,Face> {
  public:
    Integral_NT_op_dS( const PropertyDatabase<dim>&, const char* oper, const char* test );
    
    void GetOperands( const Face<dim>& ) override final;

    void ComputeContribution( const Face<dim>& ) override final;
  
  private:
    ScalarVariable        sc;
    VectorVariable<dim>   vc, vc2, pvc;
    std::vector<double>   fn;
    double                face_value, edge;
    bool                  faces;
    SG_BOUNDARY           at_boundary;
};


/* copyright (c) 2001 by Dr. Stephan K. Matthaei & Sebastian Geiger */

} // csmp

#endif
















