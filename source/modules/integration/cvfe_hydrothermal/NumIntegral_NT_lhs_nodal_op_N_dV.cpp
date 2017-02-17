#include "NumIntegral_NT_lhs_nodal_op_N_dV.h"
#include "Element.h"

using namespace std;

namespace csmp {

/** default constructor */
template<size_t dim, class SIMPLEX>
NumIntegral_NT_lhs_nodal_op_N_dV<dim,SIMPLEX>::NumIntegral_NT_lhs_nodal_op_N_dV() {}

/** default destructor */
template<size_t dim, class SIMPLEX>
NumIntegral_NT_lhs_nodal_op_N_dV<dim,SIMPLEX>::~NumIntegral_NT_lhs_nodal_op_N_dV() {}

/** custom constructor */
template<size_t dim, class SIMPLEX>
NumIntegral_NT_lhs_nodal_op_N_dV<dim,SIMPLEX>::NumIntegral_NT_lhs_nodal_op_N_dV( const PropertyDatabase<dim>& pref, 
                                                              const char* oper,
                                                              const char* basic, 
                                                              const char* test )
  : MathOperatorLHS<dim>(pref,oper,basic,test)
 {
    MathOperatorLHS<dim>::Name("NumIntegral_NT_lhs_nodal_op_N_dV", oper, basic, test );
    
    if ( MathOperatorLHS<dim>::MaterialOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_NT_lhs_nodal_op_N_dV<dim>::(constructor)", 
                      oper, "Operand must be a scalar property." );

    if ( MathOperatorLHS<dim>::MaterialOperandPlacement() != NODE )
      throw csmp::Exception( ERROR, "NumIntegral_NT_lhs_nodal_op_N_dV<dim>::(constructor)", 
                             oper, "operand has only been implemented for nodal variables." );

    if ( MathOperatorLHS<dim>::TestOperandPlacement() != NODE || 
         MathOperatorLHS<dim>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_NT_lhs_nodal_op_N_dV<dim>::(constructor)", 
                      test, "Dependent variable must be a scalar property placed on the nodes." );
                   
    if ( MathOperatorLHS<dim>::BasicOperandPlacement() != NODE || 
         MathOperatorLHS<dim>::BasicOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_NT_lhs_nodal_op_N_dV<dim>::(constructor)", 
                      basic, "Weighting variable must be a scalar property placed on the nodes." );

}

/** read in nodal vector of material operand */
template<size_t dim, class SIMPLEX>
void NumIntegral_NT_lhs_nodal_op_N_dV<dim,SIMPLEX>::GetOperands( SIMPLEX& e )
   { 

      e.NodePropertyVector( MathOperatorLHS<dim>::MaterialOperandKey(), VAR );

   } // end GetOperands


/** Calculate contributions to the capacitance matrix */
template<size_t dim, class SIMPLEX>
void NumIntegral_NT_lhs_nodal_op_N_dV<dim,SIMPLEX>::ComputeContribution( SIMPLEX& e )
 {
    MathOperatorLHS<dim>::LHS.Resize(e.Nodes(),e.Nodes());
    MathOperatorLHS<dim>::LHS.Zero();

	volume = e.Volume();
	 
	 for ( size_t j=0; j<e.Nodes(); j++ )
	   MathOperatorLHS<dim>::LHS(j,j) = (VAR[j]()*volume) / static_cast<double64>(e.Nodes());

} // end ComputeContribution


template class NumIntegral_NT_lhs_nodal_op_N_dV<1U,Element<1U> >;
template class NumIntegral_NT_lhs_nodal_op_N_dV<2U,Element<2U> >;
template class NumIntegral_NT_lhs_nodal_op_N_dV<3U,Element<3U> >;

} // csmp
