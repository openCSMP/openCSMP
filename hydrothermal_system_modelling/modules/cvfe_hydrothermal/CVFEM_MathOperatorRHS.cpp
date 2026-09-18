// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "CVFEM_MathOperatorRHS.h"
#include "Exception.h"
#include "Element.h"

using namespace std;

namespace csmp {

/** default destructor */
template<uint32_t dim, template<uint32_t> class CELL>
CVFEM_MathOperatorRHS<dim,CELL>::~CVFEM_MathOperatorRHS() {}

/** custom constructor */
template<uint32_t dim, template<uint32_t> class CELL>
CVFEM_MathOperatorRHS<dim,CELL>::CVFEM_MathOperatorRHS( const PropertyDatabase<dim>& pref,
                                                        const char* test )
  : MathOperatorRHS<dim,CELL>(pref,test)
 {

 }

/** custom constructor */
template<uint32_t dim, template<uint32_t> class CELL>
CVFEM_MathOperatorRHS<dim,CELL>::CVFEM_MathOperatorRHS( const PropertyDatabase<dim>& pref,
                                                        const char* oper,
                                                        const char* test )
  : MathOperatorRHS<dim,CELL>(pref,oper,test)
 {
 }

/** virtual function for CVFEM_Visitor */
template<uint32_t dim, template<uint32_t> class CELL>
void CVFEM_MathOperatorRHS<dim,CELL>::GetOperandsCVFEM( const CELL<dim>& e, csmp::Index upwind_var_key )
 {

	 throw csmp::Exception( ERROR, "CVFEM_MathOperatorRHS<dim,CELL>::GetOperandsCVFEM", 
                           " not specifically defined for this operator" );

}

/** access function for CVFEM_Visitor to the vector entries*/
template<uint32_t dim, template<uint32_t> class CELL>
std::vector<double> CVFEM_MathOperatorRHS<dim,CELL>::GetContribution()
 {
   return MathOperatorRHS<dim,CELL>::RHS;
 }

template class CVFEM_MathOperatorRHS<1U>;
template class CVFEM_MathOperatorRHS<2U>;
template class CVFEM_MathOperatorRHS<3U>;

template class CVFEM_MathOperatorRHS<1U,Face>;
template class CVFEM_MathOperatorRHS<2U,Face>;
template class CVFEM_MathOperatorRHS<3U,Face>;

} // csmp
