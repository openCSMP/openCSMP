//
//  IntegralEquation.h
//  CSMP_unit_tests
//
//  Created by Stephan Matthai on 19/7/20.
//  Copyright © 2020 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_INTEGRAL_EQUATION_H
#define CSMP_INTEGRAL_EQUATION_H

#include "MatrixOperator.h"
#include "VectorOperator.h"
#include "Operation.h"

namespace csmp {

enum ADE_TERM : std::int8_t { ADVECTION,
                FE_DIFFUSION,
                DISPERSION,
                SOURCE_SINK,
                COMPRESSIBILITY,
                SORPTION,
                DECAY };
                
template<size_t> class PropertyDatabase;

/**

To build an equation add terms in the desired sequence into the LHS and RHS.

This class depends on the supply of INDEX keys for the transport problem,
from the USER template parameter.

IntegralEquation (policy)
Operations 
+ Add( operator& )
+ Subtract( operator& )
+ MultiplyWith( operator& )
+ AddPostProcess( domain );
Attributes
+ map<pair<OperatorLHS,enum:Operation> >;
+ map<pair<OperatorRHS,enum:Operation> >;
+ enum OPERATION {};

*/
template<size_t dim, class VARIABLE_SET>
class IntegralEquation : public VARIABLE_SET {
  public:
    /// advection with sources and sinks
    explicit IntegralEquation( const PropertyDatabase<dim>& );
    
    /// builds equation for one of the available options and using the variables supplied in the set
    IntegralEquation( const PropertyDatabase<dim>&, const std::set<csmp::ADE_TERM>& ADE_components );
    
    /// deletes all copies of  MathOperators
    ~IntegralEquation();

    /// access to the INDEX keys for the variables used in this equation
    const VARIABLE_SET&  Notation;
    
    /// operators are added in the user-specified sequence (add, subtract, multiply, specified in operators 
    void Add( const Operation&, MatrixOperator<dim>* );
    void Add( const Operation&, VectorOperator<dim>* );
    
    /// post-processing operations involving the solution variable (these will never have to be modified)
    void AddPostProcess( MatrixOperator<dim>* );
    
    typedef typename std::map<Operation,MatrixOperator<dim>*>::const_iterator MatrixOperatorConstIterator;
    typedef typename std::map<Operation,VectorOperator<dim>*>::const_iterator VectorOperatorConstIterator;
    typedef typename std::map<std::string,MatrixOperator<dim>*>::const_iterator PostProcessingOperatorConstIterator;
    
    MatrixOperatorConstIterator LHS_OperatorsBegin() const;
    MatrixOperatorConstIterator LHS_OperatorsEnd()   const;
    
    VectorOperatorConstIterator RHS_OperatorsBegin() const;
    VectorOperatorConstIterator RHS_OperatorsEnd()   const;
    
    PostProcessingOperatorConstIterator PostProcessingOperatorsBegin() const;
    PostProcessingOperatorConstIterator PostProcessingOperatorsEnd() const;
    
    void Out() const;
  
  private:
    // sequence number-operation type, operator
    std::map<Operation,MatrixOperator<dim>*>   mat_operators_;
    std::map<Operation,VectorOperator<dim>*>   vec_operators_;
    std::map<std::string,MatrixOperator<dim>*> postpro_operators_;
};

} // end csmp

#endif /* CSMP_INTEGRAL_EQUATION_H */
