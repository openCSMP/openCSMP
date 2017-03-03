#include "NumIntegral_DNT_rhsop_DN_dV.h"
#include "PropertyDatabase.h"
#include "Element.h"
#include "Face.h"
#include "Exception.h"

using namespace std;

namespace csmp {

template<size_t dim,class SIMPLEX>
void NumIntegral_DNT_rhsop_DN_dV<dim,SIMPLEX>::IgnoreOperand( bool ignore )
  { ignore_operand = ignore; }


template<size_t dim,class SIMPLEX>
NumIntegral_DNT_rhsop_DN_dV<dim,SIMPLEX>::NumIntegral_DNT_rhsop_DN_dV( const PropertyDatabase<dim>& pref,
                                                                  const char* integral_multiplier,
                                                                  const char* test )
                          
  : MathOperatorRHS<dim>(pref,integral_multiplier,test),
    mult_key(pref.StorageKey(integral_multiplier)),    
    DN(dim,3),
    DNT(3,dim),
    OPMAT(dim,1),
    ignore_operand(true)
 {
    MathOperatorRHS<dim>::Name("NumIntegral_DNT_rhsop_DN_dV", integral_multiplier, test );
    
    if ( mult_key.type != SCALAR || mult_key.place != ELEMENT )
    throw csmp::Exception( ERROR, "NumIntegral_DNT_rhsop_DN_dV<dim>::(constructor)", 
                   integral_multiplier, " must be a scalar element property." );
    
    if ( MathOperatorRHS<dim>::TestOperandPlacement() != NODE || 
         MathOperatorRHS<dim>::TestOperandType() != SCALAR )
    throw csmp::Exception( ERROR, "NumIntegral_DNT_rhsop_DN_dV<dim>::(constructor)", 
                   test, "Dependent variable must be a scalar property placed on the nodes." );
 }




template<size_t dim,class SIMPLEX>
NumIntegral_DNT_rhsop_DN_dV<dim,SIMPLEX>::NumIntegral_DNT_rhsop_DN_dV( const PropertyDatabase<dim>& pref,
                                                              const char* integral_multiplier,
                                                              const char* oper,     
                                                              const char* test )
                          
  : MathOperatorRHS<dim>(pref,oper,test),
    mult_key(pref.StorageKey(integral_multiplier)),
    DN(dim,3),
    DNT(3,dim),
    OPMAT(dim,1),
    ignore_operand(false)
 {
    MathOperatorRHS<dim>::Name("NumIntegral_DNT_rhsop_DN_dV", oper, test );
    
    if ( mult_key.type != SCALAR || mult_key.place != ELEMENT )
    throw csmp::Exception( ERROR, "NumIntegral_DNT_rhsop_DN_dV<dim>::(constructor)", 
                   integral_multiplier, " must be a scalar element property." );
    
    if ( MathOperatorRHS<dim>::MaterialOperandPlacement() == NODE ||
         MathOperatorRHS<dim>::MaterialOperandType() != SCALAR )
    throw csmp::Exception( ERROR, "NumIntegral_DNT_rhsop_DN_dV<dim>::(constructor)", 
                   oper, "Dependent variable must be a scalar property." );

    if ( MathOperatorRHS<dim>::TestOperandPlacement() != NODE || 
         MathOperatorRHS<dim>::TestOperandType() != SCALAR )
    throw csmp::Exception( ERROR, "NumIntegral_DNT_rhsop_DN_dV<dim>::(constructor)", 
                   test, "Dependent variable must be a scalar property placed on the nodes." );
 }





template<size_t dim,class SIMPLEX>
void NumIntegral_DNT_rhsop_DN_dV<dim,SIMPLEX>::GetOperands( SIMPLEX& e )
{
   // 0. reading the integral multiplier, e.g., permeability or so
   e.Read( mult_key, multiplier ); 

   // 1. reading nodal property the div^2 of which is integrated
   if ( MathOperatorRHS<dim>::MaterialOperandPlacement() == ELEMENT ) 
     e.Read( MathOperatorRHS<dim>::MaterialOperandKey(), eoperand );
   else 
     e.NodePropertyVector( MathOperatorRHS<dim>::MaterialOperandKey(), noperand );

} // end GetOperands





/**
 
@section arguments Input Arguments 

A reference to the finite-element from which the contribution is 
computed.  
*/
template<size_t dim,class SIMPLEX>
void NumIntegral_DNT_rhsop_DN_dV<dim,SIMPLEX>::ComputeContribution( SIMPLEX& e )
{
    MathOperatorRHS<dim>::RHS.resize(e.Nodes());
    fill( MathOperatorRHS<dim>::RHS.begin(), MathOperatorRHS<dim>::RHS.end(), 0. );
    
    double64 detJ;
    
    if ( ignore_operand ) {
        OPMAT.Resize( e.Nodes(), 1 );
        OPMAT = 1.;
    
	    for ( size_t i=0; i<e.FE()->IntegrationPoints(); i++ ) {
	         detJ = e.dN_AtIntegrationPoint( DN, i );
	         DN.Transposed( DNT );

	         // multiplying BT . B 
	         DNT *= DN;
	             
	         // collapse matrix into righthand vector
	         DNT *= OPMAT;

	         // assembling contribution to right-hand vector
	         for ( size_t j=0; j<e.Nodes(); j++ ) 
	           // multiplying with determinant and weights
	           MathOperatorRHS<dim>::RHS[j] += DNT(j,0) * e.WeightAtIntegrationPoint(i) * detJ * multiplier(); 
	      }
      }
    // the Operand is considered  
    else {
        OPMAT.Resize( e.Nodes(), 1 );
    
	    for ( size_t i=0; i<e.FE()->IntegrationPoints(); i++ ) {
	         detJ = e.dN_AtIntegrationPoint( DN, i );
	         DN.Transposed( DNT );
	         
	         // multiplying BT . B 
	         DNT *= DN;
	             
	         if ( MathOperatorRHS<dim>::MaterialOperandPlacement() == ELEMENT ) OPMAT = eoperand();
	         else
             for ( size_t j=0; j<e.Nodes(); j++ ) OPMAT(j,0) = noperand[j]();
               
	         // collapse matrix into righthand vector format
	         DNT *= OPMAT;

	         // assembling contribution to right-hand vector
	         for ( size_t j=0; j<e.Nodes(); j++ ) 
	           // multiplying with determinant and weights
	           MathOperatorRHS<dim>::RHS[j] += 
	              DNT(j,0) * e.WeightAtIntegrationPoint(i) * detJ * multiplier(); 
	      }
      }

} // end ComputeContribution



// cout <<"\nNumIntegral_DNT_rhsop_DN_dV<dim>::ComputeContribution: Element "<< e.Idx() <<":"<< endl; 
// nicePrint( RHS );

template class NumIntegral_DNT_rhsop_DN_dV<1U,Element<1U> >;
template class NumIntegral_DNT_rhsop_DN_dV<2U,Element<2U> >;
template class NumIntegral_DNT_rhsop_DN_dV<3U,Element<3U> >;

template class NumIntegral_DNT_rhsop_DN_dV<1U,Face<1U> >;
template class NumIntegral_DNT_rhsop_DN_dV<2U,Face<2U> >;
template class NumIntegral_DNT_rhsop_DN_dV<3U,Face<3U> >;

} // csmp
