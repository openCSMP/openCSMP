#include "CVFEM_MathOperatorLHS.h"
#include "Exception.h"

using namespace std;

namespace csmp {

/** default constructor */
template<uint32_t dim>
CVFEM_MathOperatorLHS<dim>::CVFEM_MathOperatorLHS() {}

/** default destructor */
template<uint32_t dim>
CVFEM_MathOperatorLHS<dim>::~CVFEM_MathOperatorLHS() {}

/** custom constructor */
template<uint32_t dim>
CVFEM_MathOperatorLHS<dim>::CVFEM_MathOperatorLHS( const PropertyDatabase<dim>& pref, 
                                                              const char* basic, 
                                                              const char* test )
  : MathOperatorLHS<dim>(pref,basic,test)
 {

 }

/** custom constructor */
template<uint32_t dim>
CVFEM_MathOperatorLHS<dim>::CVFEM_MathOperatorLHS( const PropertyDatabase<dim>& pref, 
                                                              const char* oper,
                                                              const char* basic, 
                                                              const char* test )
  : MathOperatorLHS<dim>(pref,oper,basic,test)
 {

 }

/** virtual function for CVFEM_Visitor */
template<uint32_t dim>
void CVFEM_MathOperatorLHS<dim>::GetOperandsCVFEM( Element<dim>& e, csmp::Index upwind_var_key )
 {

	 throw csmp::Exception( ERROR, "CVFEM_MathOperatorLHS<dim>::GetOperandsCVFEM", 
                           " not specifically defined for this operator" );

 }

/** access function for CVFEM_Visitor to the matrix entries*/
template<uint32_t dim>
DenseMatrix<DM_MIN> CVFEM_MathOperatorLHS<dim>::GetContribution( )
 {
   return MathOperatorLHS<dim>::LHS;
 }

template class CVFEM_MathOperatorLHS<1U>;
template class CVFEM_MathOperatorLHS<2U>;
template class CVFEM_MathOperatorLHS<3U>;

} // csmp
