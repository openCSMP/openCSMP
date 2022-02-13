#include "CVFEM_MathOperatorRHS.h"
#include "Exception.h"

using namespace std;

namespace csmp {

/** default constructor */
template<uint32_t dim>
CVFEM_MathOperatorRHS<dim>::CVFEM_MathOperatorRHS() {}

/** default destructor */
template<uint32_t dim>
CVFEM_MathOperatorRHS<dim>::~CVFEM_MathOperatorRHS() {}

/** custom constructor */
template<uint32_t dim>
CVFEM_MathOperatorRHS<dim>::CVFEM_MathOperatorRHS( const PropertyDatabase<dim>& pref, 
                                                              const char* test )
  : MathOperatorRHS<dim>(pref,test)
 {

 }

/** custom constructor */
template<uint32_t dim>
CVFEM_MathOperatorRHS<dim>::CVFEM_MathOperatorRHS( const PropertyDatabase<dim>& pref, 
                                                              const char* oper,
                                                              const char* test )
  : MathOperatorRHS<dim>(pref,oper,test)
 {
 }

/** virtual function for CVFEM_Visitor */
template<uint32_t dim>
void CVFEM_MathOperatorRHS<dim>::GetOperandsCVFEM( Element<dim>& e, csmp::Index upwind_var_key )
 {

	 throw csmp::Exception( ERROR, "CVFEM_MathOperatorRHS<dim>::GetOperandsCVFEM", 
                           " not specifically defined for this operator" );

}

/** access function for CVFEM_Visitor to the vector entries*/
template<uint32_t dim>
std::vector<double> CVFEM_MathOperatorRHS<dim>::GetContribution( )
 {
   return MathOperatorRHS<dim>::RHS;
 }

template class CVFEM_MathOperatorRHS<1U>;
template class CVFEM_MathOperatorRHS<2U>;
template class CVFEM_MathOperatorRHS<3U>;

} // csp
