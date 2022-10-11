#include "CVFEM_Visitor.h"
#include "ErrorHandler.h"
#include "Exception.h"
#include "Model.h"
#include "CVFEM_MathOperatorLHS.h"
#include "CVFEM_MathOperatorRHS.h"

using namespace std;

namespace csmp {

/** custom constructor */
template<uint32_t dim, template<uint32_t> class CELL>
CVFEM_Visitor<dim,CELL>::CVFEM_Visitor( Model<dim>& model, const char* variable )
  : with_operand( false ),
  dt(0.0)
 { 
     
	this->ApplicationLevel(REGION);
  this->ApplicationTarget(ELEMENT);

	variable_key = model.Database().StorageKey(variable);
    
    if ( variable_key.place != NODE || variable_key.type != SCALAR )
      throw Exception( ERROR, "CVFEM_Visitor::(constructor)", 
                      variable, " must be a scalar property placed on the nodes." );

  }



/** custom constructor */
template<uint32_t dim, template<uint32_t> class CELL>
CVFEM_Visitor<dim,CELL>::CVFEM_Visitor( Model<dim>& model, const char* operand, const char* variable )
  : with_operand( true ),
  dt(0.0)
 {
	 this->ApplicationLevel(REGION);
   this->ApplicationTarget(ELEMENT);

	 variable_key = model.Database().StorageKey(variable);
   operand_key = model.Database().StorageKey(operand);

    if ( variable_key.place != NODE || variable_key.type != SCALAR )
      throw Exception( ERROR, "CVFEM_Visitor::(constructor)", 
                      variable, " must be a scalar property placed on the nodes." );

    if ( operand_key.type != SCALAR )
      throw Exception( ERROR, "CVFEM_Visitor::(constructor)", 
                      variable, " must be a scalar property." );

  }


/** default destructor */
template<uint32_t dim, template<uint32_t> class CELL>
CVFEM_Visitor<dim,CELL>::~CVFEM_Visitor()
 {}




/** visit function for elem */
template<uint32_t dim, template<uint32_t> class CELL>
void CVFEM_Visitor<dim,CELL>::Visit( Element<dim>* n )
  {
      GetOperands( *n );
      ComputeContribution( *n );
      WriteOperands( *n );
  }
 
 
 
 
 
template<uint32_t dim, template<uint32_t> class CELL>
void CVFEM_Visitor<dim,CELL>::GetOperands( const CELL<dim>& e )
{
  e.NodePropertyVector( variable_key, variable );
  if (with_operand) {
      if (operand_key.place == NODE)
	  {
         e.NodePropertyVector( operand_key, nodal_operand );
	  }
      else if (operand_key.place == ELEMENT)
	  {
         e.Read( operand_key, element_operand );
	  }
      else
	  {
         throw Exception( ERROR, "CVFEM_Visitor::GetOperands", 
                      " Operand variable must be placed on the node or element." );
	  }
  }

    for (lhs_it = lhs_operators.begin(); lhs_it < lhs_operators.end(); lhs_it++)
     {
      e.NodePropertyVector( lhs_it->basic_operand_key, lhs_it->basic_operands );
      lhs_it->lhs_operator->GetOperands(e);
     }

    for (lhs_it_up = lhs_operators_upwind.begin(); lhs_it_up < lhs_operators_upwind.end(); lhs_it_up++)
     {
      e.NodePropertyVector( lhs_it_up->basic_operand_key, lhs_it_up->basic_operands );
//      e.NodePropertyVector( lhs_it_up->upwind_operand_key, lhs_it_up->upwind_operands );
      lhs_it_up->lhs_operator->GetOperandsCVFEM(e, lhs_it_up->upwind_operand_key);
     }

    for (rhs_it = rhs_operators.begin(); rhs_it < rhs_operators.end(); rhs_it++)
      rhs_it->rhs_operator->GetOperands(e);         

    for (rhs_it_up = rhs_operators_upwind.begin(); rhs_it_up < rhs_operators_upwind.end(); rhs_it_up++)
     {
//      e.NodePropertyVector( rhs_it_up->upwind_operand_key, rhs_it_up->upwind_operands );
      rhs_it_up->rhs_operator->GetOperandsCVFEM(e, rhs_it_up->upwind_operand_key);
//      rhs_it_up->rhs_operator->GetOperands(e);
     }
     

} // end GetOperands







template<uint32_t dim, template<uint32_t> class CELL>
void CVFEM_Visitor<dim,CELL>::ComputeContribution( const CELL<dim>& e )
{

   for (lhs_it = lhs_operators.begin(); lhs_it < lhs_operators.end(); lhs_it++)
     {
       lhs_it->lhs_operator->ComputeContribution(e);
       if (lhs_it->lhs_operator->MultiplyWithTimeIncrement())
          lhs_it->lhs_operator->MultiplyWithTimeFactor( dt );
       LHS = lhs_it->lhs_operator->GetContribution();
       for ( uint32_t m = 0; m < e.Nodes(); m++)
          for ( uint32_t n = 0; n < e.Nodes(); n++)
             {
             LHS(m,n) *= lhs_it->basic_operands[n]();
             if (with_operand && operand_key.place == NODE
                 && nodal_operand[m]() != 0.0)
               LHS(m,n) /= nodal_operand[m]();
             else if (with_operand && operand_key.place == ELEMENT
                     && element_operand() != 0.0)
               LHS(m,n) /= element_operand();                
             variable[m]() -= LHS(m,n);
             }
     }

   for (lhs_it_up = lhs_operators_upwind.begin(); lhs_it_up < lhs_operators_upwind.end(); lhs_it_up++)
     {
       lhs_it_up->lhs_operator->ComputeContribution(e);
       if (lhs_it_up->lhs_operator->MultiplyWithTimeIncrement())
          lhs_it_up->lhs_operator->MultiplyWithTimeFactor( dt );
       LHS = lhs_it_up->lhs_operator->GetContribution();
       for ( uint32_t m = 0; m < e.Nodes(); m++)
          for ( uint32_t n = 0; n < e.Nodes(); n++)
             {
             LHS(m,n) *= lhs_it_up->basic_operands[n]();
             if (with_operand && operand_key.place == NODE
                 && nodal_operand[m]() != 0.0)
               LHS(m,n) /= nodal_operand[m]();
             else if (with_operand && operand_key.place == ELEMENT
                     && element_operand() != 0.0)
               LHS(m,n) /= element_operand();                
             variable[m]() -= LHS(m,n);
             }
     }

    for (rhs_it = rhs_operators.begin(); rhs_it < rhs_operators.end(); rhs_it++)
     {
       rhs_it->rhs_operator->ComputeContribution(e);
       if (rhs_it->rhs_operator->MultiplyWithTimeIncrement())
          rhs_it->rhs_operator->MultiplyWithTimeFactor( dt );
       RHS = rhs_it->rhs_operator->GetContribution();
       for (size_t m = 0; m < e.Nodes(); m++)
         {
          if (with_operand && operand_key.place == NODE
              && nodal_operand[m]() != 0.0)
            RHS[m] /= nodal_operand[m]();
          else if (with_operand && operand_key.place == ELEMENT
                   && element_operand() != 0.0)
            RHS[m] /= element_operand();                
          variable[m]() += RHS[m];
         }
     }

    for (rhs_it_up = rhs_operators_upwind.begin(); rhs_it_up < rhs_operators_upwind.end(); rhs_it_up++)
     {
       rhs_it_up->rhs_operator->ComputeContribution(e);
       if (rhs_it_up->rhs_operator->MultiplyWithTimeIncrement())
          rhs_it_up->rhs_operator->MultiplyWithTimeFactor( dt );
       RHS = rhs_it_up->rhs_operator->GetContribution();
       for (size_t m = 0; m < e.Nodes(); m++)
         {
          if (with_operand && operand_key.place == NODE
              && nodal_operand[m]() != 0.0)
            RHS[m] /= nodal_operand[m]();
          else if (with_operand && operand_key.place == ELEMENT
                   && element_operand() != 0.0)
            RHS[m] /= element_operand();                
          variable[m]() += RHS[m];
         }
     }

} // end ComputeContribution







template<uint32_t dim, template<uint32_t> class CELL>
void CVFEM_Visitor<dim,CELL>::WriteOperands( CELL<dim>& e )
 {

   for (auto i = 0; i < e.Nodes(); i++)
     {
   	  if ( e.N(i)->Status( variable_key ) != DIRICH )
        e.N(i)->Store( variable_key, variable[i] );
     }

 } // end WriteOperands





template<uint32_t dim, template<uint32_t> class CELL>
void CVFEM_Visitor<dim,CELL>::Add( const PropertyDatabase<dim>& pref,
                                   CVFEM_MathOperatorLHS<dim,CELL>* lhs_op )
 {

   lhs_operators.push_back( Operator_LHS(pref,lhs_op) );

 } // Add




template<uint32_t dim, template<uint32_t> class CELL>
void CVFEM_Visitor<dim,CELL>::Add( CVFEM_MathOperatorRHS<dim,CELL>* rhs_op )
 {

   if ( rhs_op->AddLater() )
     rhs_operators.push_back( Operator_RHS(rhs_op) );
   else 
     throw Exception( ERROR, "CVFEM_Visitor::Add( CVFEM_MathOperatorRHS<dim>* rhs_op )", 
                      " RHS operator must be a late accumulate." );
   
 } // Add





template<uint32_t dim, template<uint32_t> class CELL>
void CVFEM_Visitor<dim,CELL>::Add( const PropertyDatabase<dim>& pref,
                                   CVFEM_MathOperatorLHS<dim,CELL>* lhs_op, const char* upwind_variable )
 {

   lhs_operators_upwind.push_back( Operator_LHS_Upwind(pref,lhs_op,upwind_variable) );

 } // Add




template<uint32_t dim, template<uint32_t> class CELL>
void CVFEM_Visitor<dim,CELL>::Add( const PropertyDatabase<dim>& pref,
                                   CVFEM_MathOperatorRHS<dim,CELL>* rhs_op,
                                   const char* upwind_variable )
 {

   if ( rhs_op->AddLater() )
     rhs_operators_upwind.push_back( Operator_RHS_Upwind(pref, rhs_op, upwind_variable) );
   else 
     throw Exception( ERROR, "CVFEM_Visitor::Add( CVFEM_MathOperatorRHS<dim>* rhs_op )", 
                      " RHS operator must be a late accumulate." );
   
 } // Add
 




template<uint32_t dim, template<uint32_t> class CELL>
void CVFEM_Visitor<dim,CELL>::SetTimeIncrement( double time_increment )
 {
   dt = time_increment;
 } // SetTimeIncrement( double time_increment )






template<uint32_t dim, template<uint32_t> class CELL>
CVFEM_Visitor<dim,CELL>::Operator_LHS::Operator_LHS( const PropertyDatabase<dim>& pref,
                                                     CVFEM_MathOperatorLHS<dim,CELL>* lhs_op )
 : lhs_operator( lhs_op )
 {
    basic_operand_key = pref.StorageKey(lhs_op->BasicOperandName().c_str());
 } // Operator_LHS( constructor )






template<uint32_t dim, template<uint32_t> class CELL>
CVFEM_Visitor<dim,CELL>::Operator_RHS::Operator_RHS( CVFEM_MathOperatorRHS<dim,CELL>* rhs_op )
 : rhs_operator( rhs_op )
 {} // Operator_RHS( constructor )




template<uint32_t dim, template<uint32_t> class CELL>
CVFEM_Visitor<dim,CELL>::Operator_LHS_Upwind::Operator_LHS_Upwind( const PropertyDatabase<dim>& pref,
                                                                   CVFEM_MathOperatorLHS<dim,CELL>* lhs_op,
                                                                   const char* upwind_variable )
 : lhs_operator(lhs_op),
   basic_operand_key( pref.StorageKey(lhs_op->BasicOperandName().c_str()) ),
   upwind_operand_key( pref.StorageKey(upwind_variable) )
 {
    //basic_operand_key  = pref.StorageKey(lhs_op->BasicOperandName().c_str());
    //upwind_operand_key = pref.StorageKey(upwind_variable);
 } // Operator_LHS( constructor )






template<uint32_t dim, template<uint32_t> class CELL>
CVFEM_Visitor<dim,CELL>::Operator_RHS_Upwind::Operator_RHS_Upwind( const PropertyDatabase<dim>& pref,
                                                                    CVFEM_MathOperatorRHS<dim,CELL>* rhs_op,
                                                                    const char* upwind_variable  )
 : rhs_operator( rhs_op ), upwind_operand_key( pref.StorageKey(upwind_variable) )
 {
     //upwind_operand_key = pref.StorageKey(upwind_variable);
 } // Operator_RHS( constructor )


template class CVFEM_Visitor<1U,Element>;
template class CVFEM_Visitor<2U,Element>;
template class CVFEM_Visitor<3U,Element>;

} // csmp






















