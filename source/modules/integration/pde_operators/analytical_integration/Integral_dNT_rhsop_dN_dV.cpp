#include "Integral_dNT_rhsop_dN_dV.h"
#include "PropertyDatabase.h"
#include "Exception.h"
#include "Node.h"
#include "Element.h"
#include "Face.h"

using namespace std;

namespace csmp {

template<size_t dim,class SIMPLEX>
Integral_dNT_rhsop_dN_dV<dim,SIMPLEX>::Integral_dNT_rhsop_dN_dV( const PropertyDatabase<dim>& pref,
                                                          const char*             oper, 
                                                          const char*             test,
                                                          const char*             grad_var ) 
  : MathOperatorRHS<dim>(pref,oper,test),
    grad_key(pref.StorageKey( grad_var )),
    DN(2,3), DNT(3,2), VAR(3,1)
{
    MathOperatorRHS<dim>::Name("Integral_dNT_rhsop_dN_dV", oper, test );

    // testing the Operands 
    if ( MathOperatorRHS<dim>::MaterialOperandPlacement() != ELEMENT and MathOperatorRHS<dim>::MaterialOperandPlacement() != REGION )
    throw csmp::Exception( ERROR, "Integral_dNT_rhsop_dN_dV::(constructor)", 
                    oper, "Operand must be placed on the element or group.");

    if ( MathOperatorRHS<dim>::TestOperandPlacement() != NODE || MathOperatorRHS<dim>::TestOperandType() != SCALAR )
    throw csmp::Exception( ERROR, "Integral_dNT_rhsop_dN_dV::(constructor)", 
                    test, "Testfunction (dependent) variable must be a scalar property placed on the nodes.");

    if ( grad_key.place != NODE || grad_key.type != SCALAR )
    throw csmp::Exception( ERROR, "Integral_dNT_rhsop_dN_dV::(constructor)", 
                    grad_var, "Testfunction (dependent) variable must be a scalar property placed on the nodes.");
}




template<size_t dim,class SIMPLEX>
void Integral_dNT_rhsop_dN_dV<dim,SIMPLEX>::GetOperands( SIMPLEX& e )
{
    // only if the property is an element property  something is done here
    if ( MathOperatorRHS<dim>::MaterialOperandPlacement() == ELEMENT or MathOperatorRHS<dim>::MaterialOperandPlacement() == REGION ) 
      {
         MathOperatorRHS<dim>::MTRL[0].Resize(dim,dim);
         MathOperatorRHS<dim>::MTRL[0].Zero();
      
         if ( MathOperatorRHS<dim>::MaterialOperandType() == SCALAR ) {
              MathOperatorRHS<dim>::MTRL[0].AssignToDiagonal( dim, e.Read( MathOperatorRHS<dim>::MaterialOperandKey() ) ); 
           }
         if ( MathOperatorRHS<dim>::MaterialOperandType() == VECTOR ) {
              VectorVariable<dim>  vc;
              e.Read( MathOperatorRHS<dim>::MaterialOperandKey(), vc );
              MathOperatorRHS<dim>::MTRL[0].AssignToDiagonal( vc ); 
           }
         if ( MathOperatorRHS<dim>::MaterialOperandType() == TENSOR ) {
              TensorVariable<dim>  ts;
              e.Read( MathOperatorRHS<dim>::MaterialOperandKey(), ts );
              MathOperatorRHS<dim>::MTRL[0] = ts;
           }
      }
    else throw csmp::Exception( FATAL_ERROR, "Integral_dNT_rhsop_dN_dV::GetOperands", 
                               "The current finite element is analytically integrated and therefore has no integration points",
                                        "Hence, nodal properties cannot be integrated. Use numerically integrated elements in stead.");
   VAR.Resize(e.Nodes(),1);   
   // the gradient variable
   for ( size_t i=0U; i<e.Nodes(); i++ )
     VAR(i,0) = e.N(i)->Read( grad_key );   
          
} // end GetOperands






template<size_t dim,class SIMPLEX>
void Integral_dNT_rhsop_dN_dV<dim,SIMPLEX>::ComputeContribution( SIMPLEX& e )
 {
    e.dN( DN );
    // transpose the shape function derivative matrix
    DNT.Resize(e.Nodes(),dim); 
    DN.Transposed( DNT );
    
    // calculate the element contribution to RHS (lumping into vector format)
    DNT *= MathOperatorRHS<dim>::MTRL[0];
    DNT *= DN;
    DNT *= VAR;
    
    // assigning element contribution & integrating the matrix
    double64 volume = e.Volume();
    
    MathOperatorRHS<dim>::RHS.resize(e.Nodes());
    for ( size_t i=0; i<e.Nodes(); i++ )
      MathOperatorRHS<dim>::RHS[i] = DNT(i,0) * volume;

} // end ComputeContribution



template class Integral_dNT_rhsop_dN_dV<1U,Element<1U> >;
template class Integral_dNT_rhsop_dN_dV<2U,Element<2U> >;
template class Integral_dNT_rhsop_dN_dV<3U,Element<3U> >;

template class Integral_dNT_rhsop_dN_dV<1U,Face<1U> >;
template class Integral_dNT_rhsop_dN_dV<2U,Face<2U> >;
template class Integral_dNT_rhsop_dN_dV<3U,Face<3U> >;

} // csmp
