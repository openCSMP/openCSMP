#ifndef NT_OP_H
#define NT_OP_H

#include "MathOperatorRHS.h"

namespace csmp {

template<uint32_t> class Element;

/**

@brief Assigns point loads, sources etc. to nodes.

@note More or less equivalent to PointSource.

@author S.K. Matthai
@author S. Roberts
@date 1999 

*/
template<uint32_t dim, template<uint32_t> class CELL=Element>
class NT_op : public MathOperatorRHS<dim,CELL> {

  public:
    NT_op( const PropertyDatabase<dim>&, const char* oper, const char* test );
    virtual ~NT_op() {}
    
    virtual void GetOperands( const CELL<dim>& );
    virtual void ComputeContribution( const CELL<dim>& );

  private:
    std::vector<ScalarVariable >  M_;
};

// copyright (c) 2000 by Dr. Stephan K. Matthaei & Stephen G. Roberts

} // csmp

#endif
















