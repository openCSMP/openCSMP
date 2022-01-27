#ifndef NT_OP_H
#define NT_OP_H

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"
#include "Operand.h"

namespace csmp {

/**

@brief Assigns point loads, sources etc. to nodes.

@note More or less equivalent to PointSource.

@author S.K. Matthai
@author S. Roberts
@date 1999 

*/
template<size_t dim,class SIMPLEX=Element<dim> >
class NT_op : public MathOperatorRHS<dim> {

  public:

    NT_op( const PropertyDatabase<dim>& p, const char* oper, const char* test );
    virtual void GetOperands( const SIMPLEX& e );
    virtual void ComputeContribution( const SIMPLEX& e );

  private:

    std::vector<ScalarVariable >  M_;
};

// copyright (c) 2000 by Dr. Stephan K. Matthaei & Stephen G. Roberts

} // csmp

#endif
















