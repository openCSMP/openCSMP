#ifndef CSMP_MATRIX_ACCUMULATOR_H
#define CSMP_MATRIX_ACCUMULATOR_H

#include "GenericTransportScheme.h"
#include "LinearSolver.h"
#include <memory>

namespace csmp {

template<size_t dim>
class LinearSystemAccumulator {
  public:
    LinearSystemAccumulator( Model<dim>& model, const char* region );

    void AddOperatorPerimeter( MatrixOperator<dim>* op );
    void AddOperatorInterior( MatrixOperator<dim>* op );
    void AddOperatorPerimeter( VectorOperator<dim>* op );
    void AddOperatorInterior( VectorOperator<dim>* op );

    void FinaliseOperators();
  
    void TimeIncrement( double64 dt );

    void AccumulateByStencil( SparseMatrix& lhs, std::vector<double64>& rhs );
    void AccumulateByFiniteVolume( SparseMatrix& lhs, std::vector<double64>& rhs );

  private:
    LinearSystemAccumulator( ) = delete;
    LinearSystemAccumulator( const LinearSystemAccumulator& ) = delete;
    LinearSystemAccumulator( LinearSystemAccumulator&& ) = delete;
    LinearSystemAccumulator& operator=( const LinearSystemAccumulator& ) = delete;
  
    Model<dim>& model_;
    const Region<dim>& gref_;

    size_t fvs_;

    std::vector<MatrixOperator<dim>*> interior_lhs_;
    std::vector<MatrixOperator<dim>*> perimeter_lhs_;
    
    std::vector<VectorOperator<dim>*> interior_rhs_;
    std::vector<VectorOperator<dim>*> perimeter_rhs_;
};


/**
@class LinearSystemAccumulator LinearSystemAccumulator "integration/finite_volumes/generic_transport_scheme/LinearSystemAccumulator.h"

\brief     Linear system accumulator
\details   Part of the Colleoli transport scheme.
\author    Andrew J. Bromage
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
