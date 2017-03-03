#include "NumIntegral_NT_rhs_nodal_op_N_dV.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"
#include "Exception.h"

using namespace std;

namespace csmp {

/** default constructor */
template<size_t dim, class SIMPLEX>
NumIntegral_NT_rhs_nodal_op_N_dV<dim,SIMPLEX>::NumIntegral_NT_rhs_nodal_op_N_dV() {}

/** default destructor */
template<size_t dim, class SIMPLEX>
NumIntegral_NT_rhs_nodal_op_N_dV<dim,SIMPLEX>::~NumIntegral_NT_rhs_nodal_op_N_dV() {}

/** custom constructor */
template<size_t dim, class SIMPLEX>
NumIntegral_NT_rhs_nodal_op_N_dV<dim,SIMPLEX>::NumIntegral_NT_rhs_nodal_op_N_dV( const PropertyDatabase<dim>& pref, 
                                                                                 const char* oper,
																				 const char* test )
  : MathOperatorRHS<dim>(pref,oper,test)
 {
    MathOperatorRHS<dim>::Name("NumIntegral_NT_rhs_nodal_op_N_dV", oper, test );

    if ( MathOperatorRHS<dim>::MaterialOperandType() != SCALAR )
      throw csmp::Exception( CSMP_ERROR, "NumIntegral_NT_rhs_nodal_op_N_dV<dim>::(constructor)", 
                      oper, "Operand must be a scalar property." );

    if ( MathOperatorRHS<dim>::MaterialOperandPlacement() != NODE )
      throw csmp::Exception( CSMP_ERROR, "NumIntegral_NT_rhs_nodal_op_N_dV<dim>::(constructor)", 
                             oper, "operand has only been implemented for nodal variables." );

    if ( MathOperatorRHS<dim>::TestOperandPlacement() != NODE || 
         MathOperatorRHS<dim>::TestOperandType() != SCALAR )
      throw csmp::Exception( CSMP_ERROR, "NumIntegral_NT_rhs_nodal_op_N_dV<dim>::(constructor)", 
                      test, "Dependent variable must be a scalar property placed on the nodes." );


}

/** read in nodal vector of material operand */
template<size_t dim, class SIMPLEX>
void NumIntegral_NT_rhs_nodal_op_N_dV<dim,SIMPLEX>::GetOperands( SIMPLEX& e )
   { 

      e.NodePropertyVector( MathOperatorRHS<dim>::MaterialOperandKey(), VAR );

   } // end GetOperands


/** Calculate contributions to the capacitance matrix */
template<size_t dim, class SIMPLEX>
void NumIntegral_NT_rhs_nodal_op_N_dV<dim,SIMPLEX>::ComputeContribution( SIMPLEX& e )
 {

   MathOperatorRHS<dim>::RHS.resize(e.Nodes());
   fill( MathOperatorRHS<dim>::RHS.begin(), MathOperatorRHS<dim>::RHS.end(), 0.0 );
    
   volume = e.Volume();

   for ( size_t j=0; j<e.Nodes(); j++ )
      MathOperatorRHS<dim>::RHS[j] = (VAR[j]()*volume) / static_cast<double64>(e.Nodes());
   
} // end ComputeContribution

template class NumIntegral_NT_rhs_nodal_op_N_dV<1U,Element<1U> >;
template class NumIntegral_NT_rhs_nodal_op_N_dV<2U,Element<2U> >;
template class NumIntegral_NT_rhs_nodal_op_N_dV<3U,Element<3U> >;

} // csmp
