#include "Upwind_Integral_dNT_op_dN_dV.h"
#include "Exception.h"
#include "PropertyDatabase.h"
#include "Element.h"
#include "Face.h"

using namespace std;

namespace csmp {

/**
 
@section arguments Input Arguments 
The upwinding is based on the basic operand (which will be the pressure)
variable.
*/
template<size_t dim,class SIMPLEX>
Upwind_Integral_dNT_op_dN_dV<dim,SIMPLEX>::Upwind_Integral_dNT_op_dN_dV( const PropertyDatabase<dim>& pref,
                                                            				const char* oper, 
                                                            				const char* basic, 
                                                            				const char* test,
                                                            				const char* upwind,
                                                            				const char* trigger,
                                                            				const double64 prefactor)
 
  : MathOperatorLHS<dim>(pref,oper,basic,test),
    DN(2,3),
    DNT(3,2),
    uvar_(pref.StorageKey(upwind)),
    tvar_(pref.StorageKey(trigger)),
    prefactor_(prefactor)
{
    string name = "Upwind_Integral_dNT_op_dN_dV ";
    name += upwind;
    char * cname = new char[name.length()+1];
    strcpy(cname, name.c_str());
    MathOperatorLHS<dim>::Name(cname, oper, basic, test);

    // testing the Operands 
    if ( MathOperatorLHS<dim>::MaterialOperandPlacement() != ELEMENT )
    throw csmp::Exception( ERROR, "Upwind_Integral_dNT_op_dN_dV::(constructor)", 
                    oper, "Operand must be placed on the element.");

    if ( MathOperatorLHS<dim>::BasicOperandPlacement() != NODE || MathOperatorLHS<dim>::BasicOperandType() != SCALAR )
    throw csmp::Exception( ERROR, "Upwind_Integral_dNT_op_dN_dV::(constructor)", 
                    test, "Basic (dependent) variable must be a scalar property placed on the nodes.");

    if ( MathOperatorLHS<dim>::TestOperandPlacement() != NODE || MathOperatorLHS<dim>::TestOperandType() != SCALAR )
    throw csmp::Exception( ERROR, "Upwind_Integral_dNT_op_dN_dV::(constructor)", 
                    test, "Testfunction (dependent) variable must be a scalar property placed on the nodes.");
    
    if (uvar_.place != NODE || uvar_.type != SCALAR)
    throw csmp::Exception( ERROR, "Upwind_Integral_dNT_op_dN_dV::(constructor)", 
                    upwind, "Upwind variable must be a scalar property placed on the nodes.");
                    
    if (tvar_.place != NODE || uvar_.type != SCALAR)
    throw csmp::Exception( ERROR, "Upwind_Integral_dNT_op_dN_dV::(constructor)", 
                    trigger, "Upwind variable must be a scalar property placed on the nodes.");
                    
}


/**
 
GetOperands fills a vector of material property matrices with the 
required values for later computation. These matrices always have 
the dimensions spatial-dimension^2 and they will hold either scalar,
vector or tensor properties, depending on what kind of property the 
Operand is.  

When the property is an element property, it will be put into the
first vector entry MTRL[0]. */
template<size_t dim,class SIMPLEX>
void Upwind_Integral_dNT_op_dN_dV<dim,SIMPLEX>::GetOperands( SIMPLEX& e )
 {
 	MathOperatorLHS<dim>::GetOperands(e);
      
  if (uvar_.place == NODE && uvar_.type == SCALAR && tvar_.place == NODE && tvar_.type == SCALAR) {
    e.NodePropertyVector( uvar_, el_uvar);
    e.NodePropertyVector( tvar_, el_tvar);
  }
  else {
    throw csmp::Exception( FATAL_ERROR, "Upwind_Integral_dNT_op_dN_dV<dim>::GetOperands", 
                                      "Only nodal properties allowed" );
  }
    
 } // end GetOperands


template<size_t dim,class SIMPLEX>
void Upwind_Integral_dNT_op_dN_dV<dim,SIMPLEX>::ComputeContribution( SIMPLEX& e )
 {
    e.dN( DN );
    // transpose the shape function derivative matrix
    DNT.Resize(e.Nodes(),dim); 
    DN.Transposed( DNT );
    
    // calculate the element contribution to LHS (E_OP is the Basic Operand)
    DNT *= MathOperatorLHS<dim>::MTRL[0];
    DNT *= DN;
    
    // calculate upwinding coefficients and multiply them with operand matrix   
    for (size_t i = 0; i < e.Nodes(); ++i) {
    	for (size_t j = 0; j < e.Nodes(); ++j) {
    		if (i != j) {
    			const double64 decision = DNT(i, j)*(el_tvar[j]() - el_tvar[i]());
    			if      (decision > 0) DNT(i, j) *= el_uvar[i]();
    			else if (decision < 0) DNT(i, j) *= el_uvar[j]();
    			else                   DNT(i, j) *= 0.5*(el_uvar[i]() + el_uvar[j]());
    		}
    		else                     DNT(i, j) = static_cast<double64>(0.0);
    	}
    }
    
    for (size_t i = 0; i < e.Nodes(); ++i) {
    	DNT(i, i) = -DNT.RowSum(i);
    }
    
    // assigning element contribution
    MathOperatorLHS<dim>::LHS.Resize(e.Nodes(),e.Nodes());
    MathOperatorLHS<dim>::LHS = DNT;
    // analytical integration over area / volume for linear triangle 
    // and tetrahedron elements, respectively
    MathOperatorLHS<dim>::LHS *= e.Volume() * prefactor_;

 //cout <<"\nElement: "<< e.Idx();
 //MathOperatorLHS<dim>::LHS.Out();

} // end ComputeContribution




template class Upwind_Integral_dNT_op_dN_dV<2U,Element<2U> >;
template class Upwind_Integral_dNT_op_dN_dV<3U,Element<3U> >;

template class Upwind_Integral_dNT_op_dN_dV<2U,Face<2U> >;
template class Upwind_Integral_dNT_op_dN_dV<3U,Face<3U> >;

} // csmp
