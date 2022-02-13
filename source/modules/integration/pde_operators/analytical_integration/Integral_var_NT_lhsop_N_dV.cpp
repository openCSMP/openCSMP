#include "Integral_var_NT_lhsop_N_dV.h"
#include "PropertyDatabase.h"
#include "Exception.h"
#include "Element.h"
#include "Face.h"

using namespace std;

namespace csmp {


template<uint32_t dim,class CELL>
Integral_var_NT_lhsop_N_dV<dim,CELL>::Integral_var_NT_lhsop_N_dV( const PropertyDatabase<dim>& pref,
                                                                const char* oper,
                                                                const char* basic, 
                                                                const char* test,
                                                                const char* var,
                                                                const double prefactor )
  : MathOperatorLHS<dim>(pref,oper,basic,test),
    var_(pref.Parameter(var)),
    prefactor_(prefactor)
 {
    string name = "Integral_var_NT_lhsop_N_dV ";
    name += var;
    char * cname = new char[name.length()+1];
    strcpy(cname, name.c_str());
    MathOperatorLHS<dim>::Name(cname, oper, basic, test );
    
        // testing the Operands 
    if ( MathOperatorLHS<dim>::MaterialOperandPlacement() != ELEMENT && MathOperatorLHS<dim>::MaterialOperandPlacement() )
    throw csmp::Exception( ERROR, "Integral_var_NT_lhsop_N_dV::(constructor)", 
                   oper, "Operand must be a property placed on the element or group." );

    if ( MathOperatorLHS<dim>::TestOperandPlacement() != NODE || 
         MathOperatorLHS<dim>::TestOperandType() != SCALAR )
    throw csmp::Exception( ERROR, "Integral_var_NT_lhsop_N_dV::(constructor)", 
                   test, "Dependent variable must be a scalar property placed on the nodes." );
                   
    if (var_.key.place != NODE || var_.key.type != SCALAR)
      throw csmp::Exception(ERROR, "Integral_var_NT_lhsop_N_dV::(constructor)", 
                   var, "Additional variable must be a scalar property placed on the nodes.");
 }


/** Reads the Operand values from the elements.
*/
template<uint32_t dim,class CELL>
void Integral_var_NT_lhsop_N_dV<dim,CELL>::GetOperands( const CELL& e )
{
    // this integral is only for analytically integrated finite elements
    assert( e.FE()->UsesLocalCoordinates() == false );

   // reading oper and var
   e.Read(MathOperatorLHS<dim>::MaterialOperandKey(), op_ );
   e.NodePropertyVector(var_.key, vvar_ );
}


/** Computes the volume (area) integral over the testfunction products
multiplied with the Operand.  
*/
template<uint32_t dim,class CELL>
void Integral_var_NT_lhsop_N_dV<dim,CELL>::ComputeContribution( const CELL& e )
{
  if ( !MathOperatorLHS<dim>::LumpedFormulation() ) { // consistent formulation
    if (e.FE_Type() != LINEAR_TRIANGLE)
       throw csmp::Exception(ERROR, 
                      "Integral_var_NT_lhsop_N_dV::ComputeContribution",
                      "Consistent formulation only provided for elements of type LinearTriangle");
     ComputeIntegral(e);
  } else { // lumped formulation
    MathOperatorLHS<dim>::LHS.Resize(e.Nodes(),e.Nodes());
    MathOperatorLHS<dim>::LHS.Zero();
  
    DenseMatrix< DM_MIN> mat;
    e.IntegralNN(mat);
    
    for (auto i = 0; i < e.Nodes(); ++i) {
      for (size_t k = 0; k < e.Nodes(); ++k) {
        MathOperatorLHS<dim>::LHS(i, i) += mat(i, k) * vvar_[k]();
      }
    }
    MathOperatorLHS<dim>::LHS *= op_() * prefactor_;
  }
  
  //MathOperatorLHS<dim>::LHS.Out();
} // end ComputeContribution





template<uint32_t dim,class CELL>
void Integral_var_NT_lhsop_N_dV<dim,CELL>::ComputeIntegral( const CELL& e ) {
  assert(e.Nodes() == 3);
  assert(vvar_.size() == 3);

  MathOperatorLHS<dim>::LHS.Resize(3,3);

  MathOperatorLHS<dim>::LHS(0, 0) = 6.0*vvar_[0]() + 2.0*vvar_[1]() + 2.0*vvar_[2]();
  MathOperatorLHS<dim>::LHS(0, 1) = 2.0*vvar_[0]() + 2.0*vvar_[1]() + vvar_[2]();
  MathOperatorLHS<dim>::LHS(0, 2) = 2.0*vvar_[0]() + vvar_[1]() + 2.0*vvar_[2]();
  MathOperatorLHS<dim>::LHS(1, 0) = MathOperatorLHS< dim>::LHS(0, 1);
  MathOperatorLHS<dim>::LHS(1, 1) = 2.0*vvar_[0]() + 6.0*vvar_[1]() + 2.0*vvar_[2]();
  MathOperatorLHS<dim>::LHS(1, 2) = vvar_[0]() + 2.0*vvar_[1]() + 2.0*vvar_[2]();
  MathOperatorLHS<dim>::LHS(2, 0) = MathOperatorLHS< dim>::LHS(0, 2);
  MathOperatorLHS<dim>::LHS(2, 1) = MathOperatorLHS< dim>::LHS(1, 2);
  MathOperatorLHS<dim>::LHS(2, 2) = 2.0*vvar_[0]() + 2.0*vvar_[1]() + 6.0*vvar_[2]();
  
  MathOperatorLHS< dim>::LHS *= e.Volume() * prefactor_ * op_() / 60.0;
}

template class Integral_var_NT_lhsop_N_dV<2U,Element<2U> >;
template class Integral_var_NT_lhsop_N_dV<3U,Element<3U> >;

template class Integral_var_NT_lhsop_N_dV<2U,Face<2U> >;
template class Integral_var_NT_lhsop_N_dV<3U,Face<3U> >;

} // csmp

