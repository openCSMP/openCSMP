#include "CVFEM_MathOperatorLHS.h"
#include "Exception.h"
#include "Face.h"

using namespace std;

namespace csmp {

/** default constructor */
template<uint32_t dim, template<uint32_t> class CELL>
CVFEM_MathOperatorLHS<dim,CELL>::CVFEM_MathOperatorLHS() {}

/** default destructor */
template<uint32_t dim, template<uint32_t> class CELL>
CVFEM_MathOperatorLHS<dim,CELL>::~CVFEM_MathOperatorLHS() {}

/** custom constructor */
template<uint32_t dim, template<uint32_t> class CELL>
CVFEM_MathOperatorLHS<dim,CELL>::CVFEM_MathOperatorLHS( const PropertyDatabase<dim>& pref,
                                                        const char* basic,
                                                        const char* test )
  : MathOperatorLHS<dim,CELL>(pref,basic,test)
 {
 }

/** custom constructor */
template<uint32_t dim, template<uint32_t> class CELL>
CVFEM_MathOperatorLHS<dim,CELL>::CVFEM_MathOperatorLHS( const PropertyDatabase<dim>& pref,
                                                        const char* oper,
                                                        const char* basic,
                                                        const char* test )
  : MathOperatorLHS<dim,CELL>(pref,oper,basic,test)
 {

 }

/** virtual function for CVFEM_Visitor */
template<uint32_t dim, template<uint32_t> class CELL>
void CVFEM_MathOperatorLHS<dim,CELL>::GetOperandsCVFEM( const CELL<dim>&, csmp::Index upwind_var_key )
 {

	 throw csmp::Exception( ERROR, "CVFEM_MathOperatorLHS<dim>::GetOperandsCVFEM", 
                           " not specifically defined for this operator" );

 }

/** access function for CVFEM_Visitor to the matrix entries*/
template<uint32_t dim, template<uint32_t> class CELL>
DenseMatrix<DM_MIN> CVFEM_MathOperatorLHS<dim,CELL>::GetContribution()
 {
   return MathOperatorLHS<dim,CELL>::LHS;
 }

template class CVFEM_MathOperatorLHS<1U>;
template class CVFEM_MathOperatorLHS<2U>;
template class CVFEM_MathOperatorLHS<3U>;

template class CVFEM_MathOperatorLHS<1U,Face>;
template class CVFEM_MathOperatorLHS<2U,Face>;
template class CVFEM_MathOperatorLHS<3U,Face>;

} // csmp
