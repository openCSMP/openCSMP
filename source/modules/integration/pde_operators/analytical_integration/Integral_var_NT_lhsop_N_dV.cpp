#include "Integral_var_NT_lhsop_N_dV.h"
#include "PropertyDatabase.h"
#include "Exception.h"
#include "Element.h"
#include "Face.h"

using namespace std;

namespace csmp {


template<uint32_t dim, template<uint32_t> class CELL>
Integral_var_NT_lhsop_N_dV<dim,CELL>::Integral_var_NT_lhsop_N_dV( const PropertyDatabase<dim>& pref,
                                                                  const char* oper,
                                                                  const char* basic,
                                                                  const char* test,
                                                                  const char* var,
                                                                  const double prefactor )
  : MathOperatorLHS<dim,CELL>(pref,oper,basic,test),
    var_(pref.Parameter(var)),
    prefactor_(prefactor)
 {
    string name = "Integral_var_NT_lhsop_N_dV ";
    name += var;
    char * cname = new char[name.length()+1];
    strcpy(cname, name.c_str());
    MathOperatorLHS<dim,CELL>::Name(cname, oper, basic, test );
    
        // testing the Operands 
    if ( MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() != ELEMENT &&
         MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() != REGION )
    throw csmp::Exception( ERROR, "Integral_var_NT_lhsop_N_dV::(constructor)", 
                           oper, "Operand must be a property placed on the element or region." );

    if ( MathOperatorLHS<dim,CELL>::TestOperandPlacement() != NODE ||
         MathOperatorLHS<dim,CELL>::TestOperandType() != SCALAR )
    throw csmp::Exception( ERROR, "Integral_var_NT_lhsop_N_dV::(constructor)", 
                   test, "Dependent variable must be a scalar property placed on the nodes." );
                   
    if (var_.key.place != NODE || var_.key.type != SCALAR)
      throw csmp::Exception(ERROR, "Integral_var_NT_lhsop_N_dV::(constructor)", 
                   var, "Additional variable must be a scalar property placed on the nodes.");
 }




/** Reads the Operand values from the elements.
*/
template<uint32_t dim, template<uint32_t> class CELL>
void Integral_var_NT_lhsop_N_dV<dim,CELL>::GetOperands( const CELL<dim>& e )
{
    // this integral is only for analytically integrated finite elements
    assert( e.FE()->UsesLocalCoordinates() == false );

   // reading oper and var
   e.Read(MathOperatorLHS<dim,CELL>::MaterialOperandKey(), op_ );
   e.NodePropertyVector(var_.key, vvar_ );
}




/** Computes the volume (area) integral over the testfunction products
multiplied with the Operand.  
*/
template<uint32_t dim, template<uint32_t> class CELL>
void Integral_var_NT_lhsop_N_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
{
  if ( !MathOperatorLHS<dim,CELL>::LumpedFormulation() ) { // consistent formulation
    if (e.FE_Type() != LINEAR_TRIANGLE)
       throw csmp::Exception(ERROR, 
                      "Integral_var_NT_lhsop_N_dV::ComputeContribution",
                      "Consistent formulation only provided for elements of type LinearTriangle");
     ComputeIntegral(e);
  } else { // lumped formulation
    MathOperatorLHS<dim,CELL>::LHS.Resize(e.Nodes(),e.Nodes());
    MathOperatorLHS<dim,CELL>::LHS.Zero();
  
    DenseMatrix<DM_MIN> mat;
    e.IntegralNN(mat);
    
    for (uint32_t i = 0; i < e.Nodes(); ++i) {
      for (uint32_t k = 0; k < e.Nodes(); ++k) {
        MathOperatorLHS<dim,CELL>::LHS(i, i) += mat(i, k) * vvar_[k]();
      }
    }
    MathOperatorLHS<dim,CELL>::LHS *= op_() * prefactor_;
  }
  
  //MathOperatorLHS<dim>::LHS.Out();
} // end ComputeContribution





template<uint32_t dim, template<uint32_t> class CELL>
void Integral_var_NT_lhsop_N_dV<dim,CELL>::ComputeIntegral( const CELL<dim>& e ) {
  assert(e.Nodes() == 3);
  assert(vvar_.size() == 3);

  MathOperatorLHS<dim,CELL>::LHS.Resize(3,3);

  MathOperatorLHS<dim,CELL>::LHS(0, 0) = 6.0*vvar_[0]() + 2.0*vvar_[1]() + 2.0*vvar_[2]();
  MathOperatorLHS<dim,CELL>::LHS(0, 1) = 2.0*vvar_[0]() + 2.0*vvar_[1]() + vvar_[2]();
  MathOperatorLHS<dim,CELL>::LHS(0, 2) = 2.0*vvar_[0]() + vvar_[1]() + 2.0*vvar_[2]();
  MathOperatorLHS<dim,CELL>::LHS(1, 0) = MathOperatorLHS<dim,CELL>::LHS(0, 1);
  MathOperatorLHS<dim,CELL>::LHS(1, 1) = 2.0*vvar_[0]() + 6.0*vvar_[1]() + 2.0*vvar_[2]();
  MathOperatorLHS<dim,CELL>::LHS(1, 2) = vvar_[0]() + 2.0*vvar_[1]() + 2.0*vvar_[2]();
  MathOperatorLHS<dim,CELL>::LHS(2, 0) = MathOperatorLHS<dim,CELL>::LHS(0, 2);
  MathOperatorLHS<dim,CELL>::LHS(2, 1) = MathOperatorLHS<dim,CELL>::LHS(1, 2);
  MathOperatorLHS<dim,CELL>::LHS(2, 2) = 2.0*vvar_[0]() + 2.0*vvar_[1]() + 6.0*vvar_[2]();
  
  MathOperatorLHS<dim,CELL>::LHS *= e.Volume() * prefactor_ * op_() / 60.0;
}

template class Integral_var_NT_lhsop_N_dV<2U,Element>;
template class Integral_var_NT_lhsop_N_dV<3U,Element>;

template class Integral_var_NT_lhsop_N_dV<2U,Face>;
template class Integral_var_NT_lhsop_N_dV<3U,Face>;

} // csmp

