#ifndef CVFEM_VISITOR_H
#define CVFEM_VISITOR_H

#include "Visitor.h"
#include "Model.h"
#include "CVFEM_MathOperatorLHS.h"
#include "CVFEM_MathOperatorRHS.h"

/*   Changelog
     February 2015, Philipp Weis:
	 - initial port to CSMP++ and strong simplification as compared to the csmp5-version.
*/

namespace csmp {

template<size_t dim>
class CVFEM_Visitor : public Visitor<dim> {
  public:
    CVFEM_Visitor( Model<dim>& model, const char* variable);
    CVFEM_Visitor( Model<dim>& model, const char* operand, const char* variable);

    ~CVFEM_Visitor();
    
    virtual void Visit(Element<dim>* n);   
    
    void Add( const PropertyDatabase<dim>& pref,
              CVFEM_MathOperatorLHS<dim>* lhs_op );
    void Add( CVFEM_MathOperatorRHS<dim>* rhs_op );
    void Add( const PropertyDatabase<dim>& pref,
              CVFEM_MathOperatorLHS<dim>* lhs_op,
              const char* upwind_variable );
    void Add( const PropertyDatabase<dim>& pref,
              CVFEM_MathOperatorRHS<dim>* rhs_op,
              const char* upwind_variable );
    
    void SetTimeIncrement( double64 time_increment );
    
  private:
    struct Operator_LHS{ CVFEM_MathOperatorLHS<dim>* lhs_operator;
                         csmp::Index basic_operand_key;
                         std::vector<ScalarVariable> basic_operands;
                         Operator_LHS( const PropertyDatabase<dim>& pref, CVFEM_MathOperatorLHS<dim>* lhs_op );};
    struct Operator_LHS_Upwind{ CVFEM_MathOperatorLHS<dim>* lhs_operator;
                                csmp::Index basic_operand_key;
                                std::vector<ScalarVariable> basic_operands;
                                csmp::Index upwind_operand_key;
                                Operator_LHS_Upwind( const PropertyDatabase<dim>& pref,
                                              CVFEM_MathOperatorLHS<dim>* lhs_op,
                                              const char* upwind_variable );};
    struct Operator_RHS{ CVFEM_MathOperatorRHS<dim>* rhs_operator;
                         Operator_RHS( CVFEM_MathOperatorRHS<dim>* rhs_op );};
    struct Operator_RHS_Upwind{ CVFEM_MathOperatorRHS<dim>* rhs_operator;
                                csmp::Index upwind_operand_key;
                                Operator_RHS_Upwind( const PropertyDatabase<dim>& pref,
                                              CVFEM_MathOperatorRHS<dim>* rhs_op,
                                              const char* upwind_variable ); };
  
    void GetOperands( Element<dim>& e );
    void ComputeContribution( Element<dim>& e );
    void WriteOperands( Element<dim>& e );

    std::vector<Operator_LHS>  lhs_operators;
    std::vector<Operator_LHS_Upwind>  lhs_operators_upwind;
    std::vector<Operator_RHS>  rhs_operators;
    std::vector<Operator_RHS_Upwind>  rhs_operators_upwind;

    typename std::vector<Operator_LHS>::iterator        lhs_it;
    typename std::vector<Operator_LHS_Upwind>::iterator lhs_it_up;
    typename std::vector<Operator_RHS>::iterator        rhs_it;
    typename std::vector<Operator_RHS_Upwind>::iterator rhs_it_up;

    csmp::Index             variable_key, operand_key;
    std::vector<ScalarVariable> variable, nodal_operand;
    ScalarVariable  element_operand;
    
    DenseMatrix<DM_MIN> LHS;
    std::vector<double64>        RHS;
    
    double64 dt;    
    bool with_operand;

};

    /**
     @class CVFEM_Visitor CVFEM_Visitor.h

     @author Philipp Weis, ETH Zuerich
     @section contact Contact
     philipp.weis@erdw.ethz.ch

     @changes changes Latest Changes                                                                                  
  
     @section motivation Motivation
      Specialized visitor that allows the use of CVFEM_MathOperators.
	  Initially used for consistency checks for development of CVFEM scheme.
	  Currently not in use.
	  Could be used in future as replacement for FV calculations

     @section usage Usage
      To be used within the CVFEM scheme (Weis et al., Geofluids, 2014).

     @code
	 Can assemble PDE_operators in the same way as the PDE_integrator.
          
     @endcode
     
     @section dependencies Dependencies
	 CVFEM_MathOperatorLHS
	 CVFEM_MathOperatorRHS
     
     @section issues Known issues
     
     @section testing Testing
     testing was done in the period before publication in 2014.

  */
 } // csmp

#endif
