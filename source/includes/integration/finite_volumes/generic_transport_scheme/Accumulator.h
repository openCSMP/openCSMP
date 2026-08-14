#ifndef CSMP_ACCUMULATOR_H
#define CSMP_ACCUMULATOR_H

#include "MatrixOperator.h"
#include "VectorOperator.h"
#include "SparseMatrix.h"

namespace csmp {

/**
    Accumulator
    Operations 
    + SetupSystem()
    + AccumulateInteriorByFV( operator, operation )
    + AccumulateInteriorByFE( range )
    + LateAccumulateInteriorByFV()
    + LateAccumulateInteriorByFE()
    + AccumulatePerimeterByFV( operator, operation )
    + AccumulatePerimeterByFE( range )
    + LateAccumulatePerimeterByFV()
    + LateAccumulatePerimeterByFE()=
    Attributes
    + halo elements

    @attention the Accumulator uses stencil accumulation for the interior and finite volume accumulation for its perimeter.
    @attention the order in which LHS and RHS are accumulated is defined in Operation
*/
template<uint32_t dim, template<uint32_t> class USER>
class Accumulator {
  public:
    /// initalises sparse matrix and vector, taking into account how many solution variables there are
    // TODO: look at cases with multiple solution variables
    void SetUp(); 

    /// accumuates finite volume / hybride FVM _ FEM transport equation into global system
    void Accumulate( double time_increment=1. );
     
    /// compensate (+) balance at inflow boundaries, (-) balance at outflow boundaries, and any potential divergence of flow at no-flow boundaries              
    void BalanceFlowsThroughTruncatedBoundaryFiniteVolumes( SparseMatrix& lhs, std::vector<double>& rhs,
                                                            double time_increment=1.);
  private:  
     /// accumulates individual finite element - finite volume stencils for interior and halo elements
    void AccumulateByStencil( typename std::vector<Element<dim>*>::const_iterator begin, 
                              typename std::vector<Element<dim>*>::const_iterator end,
                              const MatrixOperator<dim>* const, ///< modified to contain operation information
                              SparseMatrix& lhs ) const;

     /// accumulates individual finite element - finite volume stencils for interior and halo elements
    void AccumulateByStencil( typename std::vector<Element<dim>*>::const_iterator begin, 
                              typename std::vector<Element<dim>*>::const_iterator end,
                              const VectorOperator<dim>* const, ///< modified to contain operation information
                              std::vector<double>& rhs ) const;
    
    /// accumulates individual finite volumes
    void AccumulateByFiniteVolume( typename std::vector<Node<dim>*>::const_iterator begin,
                                   typename std::vector<Node<dim>*>::const_iterator end,
                                   const MatrixOperator<dim>* const, ///< modified to contain operation information
                                   SparseMatrix& lhs ) const;

    void AccumulateByFiniteVolume( typename std::vector<Node<dim>*>::const_iterator begin,
                                   typename std::vector<Node<dim>*>::const_iterator end,
                                   const VectorOperator<dim>* const, ///< modified to contain operation information
                                   std::vector<double>& rhs ) const;

  private:
    /// shorthand for accessing the class that this is a policy of
    USER<dim>* User() { return static_cast<USER<dim>*>(this); }
    USER<dim> const* User() const { return static_cast<const USER<dim>*>(this); }
};

} // end csmp


#endif /* CSMP_MATRIX_ACCUMULATOR_H */
