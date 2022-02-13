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
template<uint32_t dim,class CELL=Element<dim> >
class NT_op : public MathOperatorRHS<dim> {

  public:

    NT_op( const PropertyDatabase<dim>& p, const char* oper, const char* test );
    virtual void GetOperands( const CELL& e );
    virtual void ComputeContribution( const CELL& e );

  private:

    std::vector<ScalarVariable >  M_;
};

// copyright (c) 2000 by Dr. Stephan K. Matthaei & Stephen G. Roberts

} // csmp

#endif
















