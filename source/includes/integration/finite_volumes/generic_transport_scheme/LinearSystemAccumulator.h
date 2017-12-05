#ifndef CSMP_MATRIX_ACCUMULATOR_H
#define CSMP_MATRIX_ACCUMULATOR_H

#include "DenseMatrix.h"
#include "SparseMatrix.h"
#include "LinearSolver.h"
#include <memory>

namespace csmp {

  template<size_t dim> class Model;
  template<size_t dim> class Region;
  template<size_t dim> class Element;
  template<size_t dim> class Node;

enum ACCUMULATION_MODE {
    ADD_ACCUMULATE,
    ADD_LATER,
    MULTIPLY_ACCUMULATE
};


template<size_t dim>
class MatrixOperator
{
public:
    ACCUMULATION_MODE AccumulationMode() const { return mode_; }

    virtual void AccumulateStencil( Element<dim>& fe, SparseMatrix& mat ) const = 0;
    virtual void AccumulateFiniteVolume( Node<dim>& fv, SparseMatrix& mat ) const = 0;

    virtual ~MatrixOperator() { }

protected:
    MatrixOperator( ACCUMULATION_MODE mode )
       : mode_(mode)
    {
    }

    ACCUMULATION_MODE mode_;
};


template<size_t dim>
class VectorOperator
{
public:
    ACCUMULATION_MODE AccumulationMode() const { return mode_; }

    virtual void AccumulateStencil( Element<dim>& fe, std::vector<double64>& mat ) const = 0;
    virtual void AccumulateFiniteVolume( Node<dim>& fv, std::vector<double64>& mat ) const = 0;

    virtual ~VectorOperator() { }

protected:
    VectorOperator( ACCUMULATION_MODE mode )
       : mode_(mode)
    {
    }

    ACCUMULATION_MODE mode_;    
};


template<size_t dim>
class LinearSystemAccumulator {
  public:
    LinearSystemAccumulator( Model<dim>& model, const char* region );

    void AddOperatorPerimeter( MatrixOperator<dim>* op );
    void AddOperatorInterior( MatrixOperator<dim>* op );
    void AddOperatorPerimeter( VectorOperator<dim>* op );
    void AddOperatorInterior( VectorOperator<dim>* op );

    void AccumulateByStencil();
    void AccumulateByFiniteVolume();

    void SolveSystem();

  private:
    LinearSystemAccumulator( ) = delete;
    LinearSystemAccumulator( const LinearSystemAccumulator& ) = delete;
    LinearSystemAccumulator( LinearSystemAccumulator&& ) = delete;
    LinearSystemAccumulator& operator=( const LinearSystemAccumulator& ) = delete;

    Model<dim>& model_;
    Region<dim>& gref_;

    size_t fvs_;

    std::unique_ptr<Solver> solver_;

    SparseMatrix LHS_;
    std::vector<double64> RHS_;
    std::vector<double64> RESULT_;

    std::vector<MatrixOperator<dim>*> interior_lhs_;
    std::vector<MatrixOperator<dim>*> perimeter_lhs_;
    
    std::vector<VectorOperator<dim>*> interior_rhs_;
    std::vector<VectorOperator<dim>*> perimeter_rhs_;
};


/**
@class LinearSystemAccumulator LinearSystemAccumulator "integration/finite_volumes/generic_transport_scheme/LinearSystemAccumulator.h"

\brief     Linear system accumulator
\details   Part of the Colleoli transport scheme.
\author    
\version   0a
\date      5/12/2017
\pre       high-level class depending on CSMP++ API
\bug
\warning
\copyright The University of Melbourne

@section motivation Motivation


@section design Design Intent


@section applicability Applicability


@section collaborations Collaborations


@section implementation Implementation


@section examples Application Examples

@code

@endcode

*/

} // end csmp


#endif /* CSMP_MATRIX_ACCUMULATOR_H */
