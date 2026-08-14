#include "Jacobian_Integral_dNT_op_dN_dV.h"
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
template<uint32_t dim, template<uint32_t> class CELL>
Jacobian_Integral_dNT_op_dN_dV<dim,CELL>::Jacobian_Integral_dNT_op_dN_dV( const PropertyDatabase<dim>& pref,
                                                                          const char* oper,
                                                                          const char* basic,
                                                                          const char* test,
                                                                          const char* test_orig,
                                                                          const char* lambda,
                                                                          const char* d_lambda,
                                                                          const double delta,
                                                                          const double prefactor )
 
  : MathOperatorLHS<dim,CELL>(pref,oper,basic,test),
    DN(2,3),
    DNT(3,2),
    test_orig_(pref.StorageKey(test_orig)),
    lambda_(pref.StorageKey(lambda)),
    d_lambda_(pref.StorageKey(d_lambda)),
    delta_(delta),
    prefactor_(prefactor)
{
    string name = "Jacobian_Integral_dNT_op_dN_dV ";
    name += lambda;
    char* cname = new char[name.length()+1];
    strcpy(cname, name.c_str());
    MathOperatorLHS<dim,CELL>::Name(cname, oper, basic, test );

    // testing the Operands 
    // * Add tests for new operands
    if ( MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() != ELEMENT &&
         MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() != REGION )
      throw csmp::Exception( ERROR, "Jacobian_Integral_dNT_op_dN_dV::(constructor)", 
                             oper, "Operand must be placed on the element or group.");

    if ( MathOperatorLHS<dim,CELL>::BasicOperandPlacement() != NODE ||
         MathOperatorLHS<dim,CELL>::BasicOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "Jacobian_Integral_dNT_op_dN_dV::(constructor)", 
                             test, "Basic (dependent) variable must be a scalar property placed on the nodes.");

    if ( MathOperatorLHS<dim,CELL>::TestOperandPlacement() != NODE ||
         MathOperatorLHS<dim,CELL>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "Jacobian_Integral_dNT_op_dN_dV::(constructor)", 
                             test, "Testfunction (dependent) variable must be a scalar property placed on the nodes.");
    
    if (lambda_.place != NODE || lambda_.type != SCALAR)
      throw csmp::Exception( ERROR, "Jacobian_Integral_dNT_op_dN_dV::(constructor)", 
                             lambda, "Jacobian_Upwind variable must be a scalar property placed on the nodes.");
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
template<uint32_t dim, template<uint32_t> class CELL>
void Jacobian_Integral_dNT_op_dN_dV<dim,CELL>::GetOperands( const CELL<dim>& e )
 {
    // this integral is only for analytically integrated finite elements
    assert( e.UsesLocalCoordinates() == false );

    MathOperatorLHS<dim,CELL>::GetOperands(e);
        
    if (lambda_.place == NODE && lambda_.type == SCALAR) {
      e.NodePropertyVector( lambda_, el_lambda);
      e.NodePropertyVector( test_orig_, el_test_orig);
      e.NodePropertyVector( d_lambda_, el_d_lambda);
    }
    else {
      throw csmp::Exception( FATAL_ERROR, "Jacobian_Integral_dNT_op_dN_dV<dim>::GetOperands", 
                                        "Only nodal properties allowed" );
    }
   
 } // end GetOperands




/// @todo (3) Compute Jacobian
template<uint32_t dim, template<uint32_t> class CELL>
void Jacobian_Integral_dNT_op_dN_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
 {
    e.dN( DN );
    // transpose the shape function derivative matrix
    DNT.Resize(e.Nodes(),dim); 
    DN.Transposed( DNT );
    
    // calculate the element contribution to LHS (E_OP is the Basic Operand)
    DNT *= MathOperatorLHS<dim,CELL>::MTRL[0];
    DNT *= DN;
    
    // Assemble solution in DN
    res_.Resize(e.Nodes(),e.Nodes());
    res_.Fill(0.0);
  
    const double global_factor = prefactor_ * e.Volume() / (delta_ * static_cast<double>(e.Nodes()) );
    for ( uint32_t j{0U}; j < e.Nodes(); ++j) {
      const double j_factor = el_d_lambda[j]() - el_lambda[j]();
      for (uint32_t i{0U}; i < e.Nodes(); ++i) {
        for (uint32_t k{0U}; k < e.Nodes(); ++k) {
          res_(i, j) += DNT(i, k) * el_test_orig[k]() * j_factor * global_factor;
        }
      }
    }
    
    // assigning element contribution
    MathOperatorLHS<dim,CELL>::LHS.Resize(e.Nodes(),e.Nodes());
    MathOperatorLHS<dim,CELL>::LHS = res_;

} // end ComputeContribution



template class Jacobian_Integral_dNT_op_dN_dV<2U,Element>;
template class Jacobian_Integral_dNT_op_dN_dV<3U,Element>;

template class Jacobian_Integral_dNT_op_dN_dV<2U,Face>;
template class Jacobian_Integral_dNT_op_dN_dV<3U,Face>;

} // csmp
