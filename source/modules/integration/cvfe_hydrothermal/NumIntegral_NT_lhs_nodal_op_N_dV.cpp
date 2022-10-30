#include "NumIntegral_NT_lhs_nodal_op_N_dV.h"
#include "Element.h"
#include "Face.h"
#include "Exception.h"

using namespace std;

namespace csmp {

/** custom constructor */
template<uint32_t dim, template<uint32_t> class CELL>
NumIntegral_NT_lhs_nodal_op_N_dV<dim,CELL>::NumIntegral_NT_lhs_nodal_op_N_dV( const PropertyDatabase<dim>& pref,
                                                                              const char* oper,
                                                                              const char* basic,
                                                                              const char* test )
  : MathOperatorLHS<dim,CELL>(pref,oper,basic,test)
 {
    MathOperatorLHS<dim,CELL>::Name("NumIntegral_NT_lhs_nodal_op_N_dV", oper, basic, test );
    
    if ( MathOperatorLHS<dim,CELL>::MaterialOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_NT_lhs_nodal_op_N_dV<dim>::(constructor)", 
                      oper, "Operand must be a scalar property." );

    if ( MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() != NODE )
      throw csmp::Exception( ERROR, "NumIntegral_NT_lhs_nodal_op_N_dV<dim>::(constructor)", 
                             oper, "operand has only been implemented for nodal variables." );

    if ( MathOperatorLHS<dim,CELL>::TestOperandPlacement() != NODE ||
         MathOperatorLHS<dim,CELL>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_NT_lhs_nodal_op_N_dV<dim>::(constructor)", 
                      test, "Dependent variable must be a scalar property placed on the nodes." );
                   
    if ( MathOperatorLHS<dim,CELL>::BasicOperandPlacement() != NODE ||
         MathOperatorLHS<dim,CELL>::BasicOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_NT_lhs_nodal_op_N_dV<dim>::(constructor)", 
                      basic, "Weighting variable must be a scalar property placed on the nodes." );

}





template<uint32_t dim, template<uint32_t> class CELL>
NumIntegral_NT_lhs_nodal_op_N_dV<dim,CELL>::~NumIntegral_NT_lhs_nodal_op_N_dV()
 {
 }




/** read in nodal vector of material operand */
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_NT_lhs_nodal_op_N_dV<dim,CELL>::GetOperands( const CELL<dim>& e )
   { 

      e.NodePropertyVector( MathOperatorLHS<dim,CELL>::MaterialOperandKey(), VAR );

   } // end GetOperands






/** Calculate contributions to the capacitance matrix */
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_NT_lhs_nodal_op_N_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
 {
    MathOperatorLHS<dim,CELL>::LHS.Resize(e.Nodes(),e.Nodes());
    MathOperatorLHS<dim,CELL>::LHS.Zero();

	volume = e.Volume();
	 
	 for ( uint32_t j{0U}; j<e.Nodes(); j++ )
	   MathOperatorLHS<dim,CELL>::LHS(j,j) = (VAR[j]()*volume) / static_cast<double>(e.Nodes());

} // end ComputeContribution


template class NumIntegral_NT_lhs_nodal_op_N_dV<1U>;
template class NumIntegral_NT_lhs_nodal_op_N_dV<2U>;
template class NumIntegral_NT_lhs_nodal_op_N_dV<3U>;

template class NumIntegral_NT_lhs_nodal_op_N_dV<1U,Face>;
template class NumIntegral_NT_lhs_nodal_op_N_dV<2U,Face>;
template class NumIntegral_NT_lhs_nodal_op_N_dV<3U,Face>;

} // csmp
