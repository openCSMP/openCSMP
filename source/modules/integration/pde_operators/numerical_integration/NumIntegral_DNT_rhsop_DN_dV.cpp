#include "NumIntegral_DNT_rhsop_DN_dV.h"
#include "PropertyDatabase.h"
#include "Element.h"
#include "Face.h"
#include "Exception.h"

using namespace std;

namespace csmp {

template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_DNT_rhsop_DN_dV<dim,CELL>::IgnoreOperand( bool ignore )
  { ignore_operand = ignore; }


template<uint32_t dim, template<uint32_t> class CELL>
NumIntegral_DNT_rhsop_DN_dV<dim,CELL>::NumIntegral_DNT_rhsop_DN_dV( const PropertyDatabase<dim>& pref,
                                                                    const char* integral_multiplier,
                                                                    const char* test )
                          
  : MathOperatorRHS<dim,CELL>(pref,integral_multiplier,test),
    mult_key(pref.StorageKey(integral_multiplier)),    
    DN(dim,3),
    DNT(3,dim),
    OPMAT(dim,1),
    ignore_operand(true)
 {
    MathOperatorRHS<dim,CELL>::Name("NumIntegral_DNT_rhsop_DN_dV", integral_multiplier, test );
    
    if ( mult_key.type != SCALAR || mult_key.place != ELEMENT )
    throw csmp::Exception( ERROR, "NumIntegral_DNT_rhsop_DN_dV<dim>::(constructor)", 
                   integral_multiplier, " must be a scalar element property." );
    
    if ( MathOperatorRHS<dim,CELL>::TestOperandPlacement() != NODE || 
         MathOperatorRHS<dim,CELL>::TestOperandType() != SCALAR )
    throw csmp::Exception( ERROR, "NumIntegral_DNT_rhsop_DN_dV<dim>::(constructor)", 
                   test, "Dependent variable must be a scalar property placed on the nodes." );
 }




template<uint32_t dim, template<uint32_t> class CELL>
NumIntegral_DNT_rhsop_DN_dV<dim,CELL>::NumIntegral_DNT_rhsop_DN_dV( const PropertyDatabase<dim>& pref,
                                                                    const char* integral_multiplier,
                                                                    const char* oper,
                                                                    const char* test )
                          
  : MathOperatorRHS<dim,CELL>(pref,oper,test),
    mult_key(pref.StorageKey(integral_multiplier)),
    DN(dim,3),
    DNT(3,dim),
    OPMAT(dim,1),
    ignore_operand(false)
 {
    MathOperatorRHS<dim,CELL>::Name("NumIntegral_DNT_rhsop_DN_dV", oper, test );
    
    if ( mult_key.type != SCALAR || mult_key.place != ELEMENT )
    throw csmp::Exception( ERROR, "NumIntegral_DNT_rhsop_DN_dV<dim>::(constructor)", 
                   integral_multiplier, " must be a scalar element property." );
    
    if ( MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() == NODE ||
         MathOperatorRHS<dim,CELL>::MaterialOperandType() != SCALAR )
    throw csmp::Exception( ERROR, "NumIntegral_DNT_rhsop_DN_dV<dim>::(constructor)", 
                   oper, "Dependent variable must be a scalar property." );

    if ( MathOperatorRHS<dim,CELL>::TestOperandPlacement() != NODE || 
         MathOperatorRHS<dim,CELL>::TestOperandType() != SCALAR )
    throw csmp::Exception( ERROR, "NumIntegral_DNT_rhsop_DN_dV<dim>::(constructor)", 
                   test, "Dependent variable must be a scalar property placed on the nodes." );
 }





template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_DNT_rhsop_DN_dV<dim,CELL>::GetOperands( const CELL<dim>& e )
{
    // this integral is only for numerically integrated isoparametric finite elements
    assert( e.FE()->Isoparametric() == true );

   // 0. reading the integral multiplier, e.g., permeability or so
   e.Read( mult_key, multiplier ); 

   // 1. reading nodal property the div^2 of which is integrated
   if ( MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT ) 
     e.Read( MathOperatorRHS<dim,CELL>::MaterialOperandKey(), eoperand );
   else 
     e.NodePropertyVector( MathOperatorRHS<dim,CELL>::MaterialOperandKey(), noperand );

} // end GetOperands





/**
 
@section arguments Input Arguments 

A reference to the finite-element from which the contribution is 
computed.  
*/
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_DNT_rhsop_DN_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
{
    MathOperatorRHS<dim,CELL>::RHS.resize(e.Nodes());
    fill( MathOperatorRHS<dim,CELL>::RHS.begin(), MathOperatorRHS<dim,CELL>::RHS.end(), 0. );
    
    double detJ;
    
    if ( ignore_operand ) {
        OPMAT.Resize( e.Nodes(), 1 );
        OPMAT = 1.;
    
	    for ( auto i{0U}; i<e.FE()->IntegrationPoints(); i++ ) {
	         detJ = e.dN_AtIntegrationPoint( DN, i );
	         DN.Transposed( DNT );

	         // multiplying BT . B 
	         DNT *= DN;
	             
	         // collapse matrix into righthand vector
	         DNT *= OPMAT;

	         // assembling contribution to right-hand vector
	         for ( auto j{0U}; j<e.Nodes(); j++ )
	           // multiplying with determinant and weights
	           MathOperatorRHS<dim,CELL>::RHS[j] += DNT(j,0) * e.WeightAtIntegrationPoint(i) * detJ * multiplier(); 
	      }
      }
    // the Operand is considered  
    else {
        OPMAT.Resize( e.Nodes(), 1 );
    
	    for ( auto i{0U}; i<e.FE()->IntegrationPoints(); i++ ) {
	         detJ = e.dN_AtIntegrationPoint( DN, i );
	         DN.Transposed( DNT );
	         
	         // multiplying BT . B 
	         DNT *= DN;
	             
	         if ( MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT ) OPMAT = eoperand();
	         else
             for ( auto j{0U}; j<e.Nodes(); j++ ) OPMAT(j,0) = noperand[j]();
               
	         // collapse matrix into righthand vector format
	         DNT *= OPMAT;

	         // assembling contribution to right-hand vector
	         for ( auto j{0U}; j<e.Nodes(); j++ ) 
	           // multiplying with determinant and weights
	           MathOperatorRHS<dim,CELL>::RHS[j] += 
	              DNT(j,0) * e.WeightAtIntegrationPoint(i) * detJ * multiplier(); 
	      }
      }

} // end ComputeContribution



// cout <<"\nNumIntegral_DNT_rhsop_DN_dV<dim>::ComputeContribution: Element "<< e.Idx() <<":"<< endl; 
// nicePrint( RHS );

template class NumIntegral_DNT_rhsop_DN_dV<1U,Element>;
template class NumIntegral_DNT_rhsop_DN_dV<2U,Element>;
template class NumIntegral_DNT_rhsop_DN_dV<3U,Element>;

template class NumIntegral_DNT_rhsop_DN_dV<1U,Face>;
template class NumIntegral_DNT_rhsop_DN_dV<2U,Face>;
template class NumIntegral_DNT_rhsop_DN_dV<3U,Face>;

} // csmp
