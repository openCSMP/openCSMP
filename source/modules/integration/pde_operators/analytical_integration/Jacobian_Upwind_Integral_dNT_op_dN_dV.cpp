#include "Jacobian_Upwind_Integral_dNT_op_dN_dV.h"
#include "PropertyDatabase.h"
#include "Exception.h"
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
Jacobian_Upwind_Integral_dNT_op_dN_dV<dim,SIMPLEX>::Jacobian_Upwind_Integral_dNT_op_dN_dV( const PropertyDatabase<dim>& pref,
                                                            				const char* oper, 
                                                            				const char* basic, 
                                                            				const char* test,
                                                            				const char* test_orig,
                                                            				const char* upwind,
                                                            				const char* d_upwind,
                                                            				const char* trigger,
                                                            				const double delta,
                                                            				const double prefactor)
 
  : MathOperatorLHS<dim>(pref,oper,basic,test),
    DN(2,3),
    DNT(3,2),
    test_orig_(pref.StorageKey(test_orig)),
    upwind_(pref.StorageKey(upwind)),
    d_upwind_(pref.StorageKey(d_upwind)),
    trigger_(pref.StorageKey(trigger)),
    delta_(delta),
    prefactor_(prefactor)
{
    string name = "Jacobian_Upwind_Integral_dNT_op_dN_dV ";
    name += upwind;
    name += ", ";
    name += trigger;
    char* cname = new char[name.length()+1];
    strcpy(cname, name.c_str()); 
    MathOperatorLHS<dim>::Name(cname, oper, basic, test );

    // testing the Operands 
    // * Add tests for new operands
    if ( MathOperatorLHS<dim>::MaterialOperandPlacement() != ELEMENT && MathOperatorLHS<dim>::MaterialOperandPlacement() != REGION )
    throw csmp::Exception( ERROR, "Jacobian_Upwind_Integral_dNT_op_dN_dV::(constructor)", 
                    oper, "Operand must be placed on the element or group.");

    if ( MathOperatorLHS<dim>::BasicOperandPlacement() != NODE || MathOperatorLHS<dim>::BasicOperandType() != SCALAR )
    throw csmp::Exception( ERROR, "Jacobian_Upwind_Integral_dNT_op_dN_dV::(constructor)", 
                    test, "Basic (dependent) variable must be a scalar property placed on the nodes.");

    if ( MathOperatorLHS<dim>::TestOperandPlacement() != NODE || MathOperatorLHS<dim>::TestOperandType() != SCALAR )
    throw csmp::Exception( ERROR, "Jacobian_Upwind_Integral_dNT_op_dN_dV::(constructor)", 
                    test, "Testfunction (dependent) variable must be a scalar property placed on the nodes.");
    
    if (upwind_.place != NODE || upwind_.type != SCALAR)
    throw csmp::Exception( ERROR, "Jacobian_Upwind_Integral_dNT_op_dN_dV::(constructor)", 
                    upwind, "Jacobian_Upwind variable must be a scalar property placed on the nodes.");
                    
    if (trigger_.place != NODE || trigger_.type != SCALAR)
    throw csmp::Exception( ERROR, "Jacobian_Upwind_Integral_dNT_op_dN_dV::(constructor)", 
                    trigger, "Jacobian_Upwind variable must be a scalar property placed on the nodes.");
                    
}


/**
 
GetOperands fills a vector of material property matrices with the 
required values for later computation. These matrices always have 
the dimensions spatial-dimension^2 and they will hold either scalar,
vector or tensor properties, depending on what kind of property the 
Operand is.  

When the property is an element property, it will be put into the
first vector entry MTRL[0]. Else, 
*/
template<size_t dim,class SIMPLEX>
void Jacobian_Upwind_Integral_dNT_op_dN_dV<dim,SIMPLEX>::GetOperands( const SIMPLEX& e )
 {
    // this integral is only for analytically integrated finite elements
    assert( e.FE()->UsesLocalCoordinates() == false );

    MathOperatorLHS<dim>::GetOperands(e);
        
    if (upwind_.place == NODE && upwind_.type == SCALAR && trigger_.place == NODE && trigger_.type == SCALAR) {
      e.NodePropertyVector( upwind_, el_upwind);
      e.NodePropertyVector( trigger_, el_trigger);
      e.NodePropertyVector( test_orig_, el_test_orig);
      e.NodePropertyVector( d_upwind_, el_d_upwind);
    }
    else {
      throw csmp::Exception( FATAL_ERROR, "Jacobian_Upwind_Integral_dNT_op_dN_dV<dim>::GetOperands", 
                                   "Only nodal properties allowed" );
    }
    
 } // end GetOperands


/// @todo compute Jacobian
template<size_t dim,class SIMPLEX>
void Jacobian_Upwind_Integral_dNT_op_dN_dV<dim,SIMPLEX>::ComputeContribution( const SIMPLEX& e )
 {
    e.dN( DN );
    // transpose the shape function derivative matrix
    DNT.Resize(e.Nodes(),dim); 
    DN.Transposed( DNT );
    
    // calculate the element contribution to LHS (E_OP is the Basic Operand)
    DNT *= MathOperatorLHS<dim>::MTRL[0];
    DNT *= DN;
    
    // generate delta vector
    // * can be replaced by STL algorithm
    vector<double>  delta_up(e.Nodes());
    for (size_t i = 0; i < e.Nodes(); ++i) {
      delta_up[i] = el_d_upwind[i]() - el_upwind[i]();
    }
    
    // calculate DN as modified DNT by multiplying with p_k-p_i
    DN.Resize(e.Nodes(), e.Nodes());
    for (size_t i = 0; i < e.Nodes(); ++i) {
      for (size_t k = 0; k < e.Nodes(); ++k) {
        DN(i, k) = DNT(i, k)*(el_test_orig[k]() - el_test_orig[i]());
      }
    }
    
    // assemble solution in LHS
    MathOperatorLHS<dim>::LHS.Resize(e.Nodes(),e.Nodes());
    MathOperatorLHS<dim>::LHS.Fill(0.0);
     
    for (size_t i = 0; i < e.Nodes(); ++i) {
      for (size_t k = 0; k < e.Nodes(); ++k) {
        if (i != k) {
          const double decision = DNT(i, k)*(el_trigger[k]() - el_trigger[i]());
          if (decision > 0) { // take j = i
            MathOperatorLHS<dim>::LHS(i, i) += delta_up[i]*DN(i, k);
          } else if (decision < 0) { // take j = k
            MathOperatorLHS<dim>::LHS(i, k) += delta_up[k]*DN(i, k);
          } else { // take j = i and k
            MathOperatorLHS<dim>::LHS(i, i) += 0.5*delta_up[i]*DN(i, k);
            MathOperatorLHS<dim>::LHS(i, k) += 0.5*delta_up[k]*DN(i, k);
          }
        }
      }
    }
    
    // analytical integration over area / volume for linear triangle 
    // and tetrahedron elements, respectively
    MathOperatorLHS<dim>::LHS *= e.Volume() * prefactor_/delta_;
    
    /*
    double fac_upwind;
    
    // loop over j to get an column in the element tangent stiffness matrix
    for (size_t j = 0; j < e.Nodes(); ++j) {
      for (size_t i = 0; i < e.Nodes(); ++i) {
        for (size_t k = 0; k < e.Nodes(); ++k) {
          if (i != k) {
            fac_upwind = 0.0;
            const double decision = DNT(i, k)*(el_trigger[k]() - el_trigger[i]());
            if (decision > 0) {
              if (i == j) fac_upwind = el_d_upwind[i]() - el_upwind[i]();
            }
    			  else if (decision < 0) {
    			    if (k == j) fac_upwind = el_d_upwind[k]() - el_upwind[k]();
    			  }
    			  else {
    			    if (i == j) fac_upwind = 0.5*(el_d_upwind[i]() - el_upwind[i]());
    			    else if (k == j) fac_upwind = 0.5*(el_d_upwind[k]() - el_upwind[k]());
    			  }
    			  // off diagonal contribution
    			  MathOperatorLHS<dim>::LHS(i, j) += fac_upwind*DNT(i, k)*el_test_orig[k](); 
    			  // diagonal contribution
    			  MathOperatorLHS<dim>::LHS(i, j) -= fac_upwind*DNT(i, k)*el_test_orig[i]();
          } // end if (i != k)
        }
      }
    }
    

    // analytical integration over area / volume for linear triangle 
    // and tetrahedron elements, respectively
    MathOperatorLHS<dim>::LHS *= e.Volume() * prefactor_/delta_;
  */

   // cout <<"\nElement: "<< e.Idx();
   // MathOperatorLHS<dim>::LHS.Out();

} // end ComputeContribution


template class Jacobian_Upwind_Integral_dNT_op_dN_dV<2U,Element<2U> >;
template class Jacobian_Upwind_Integral_dNT_op_dN_dV<3U,Element<3U> >;

template class Jacobian_Upwind_Integral_dNT_op_dN_dV<2U,Face<2U> >;
template class Jacobian_Upwind_Integral_dNT_op_dN_dV<3U,Face<3U> >;

} // csmp
