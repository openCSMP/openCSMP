//
//  LinearSystemAccumulator.cpp
//
//  Created by Andrew J. Bromage on 5/12/2017.
//  Copyright (c) 2017 The University of Melbourne. All rights reserved.
//

#include "LinearSystemAccumulator.h"
#include "Model.h"
#include "Region.h"
#include "Element.h"
#include "Node.h"
#include "MatrixOperator.h"
#include "VectorOperator.h"
#include "LinearSolver.h"
#ifdef CSMP_WITH_SAMG_SOLVER
#include "SAMG_Settings.h"
#endif

using namespace std;

namespace csmp {


template<size_t dim>
LinearSystemAccumulator<dim>::LinearSystemAccumulator( Model<dim>& model, const char* region )
  : model_(model),
    gref_(model_.Region(region))
{
}


template<size_t dim>
void LinearSystemAccumulator<dim>::AddOperatorPerimeter( MatrixOperator<dim>* op )
{
    perimeter_lhs_.push_back(op);
}


template<size_t dim>
void LinearSystemAccumulator<dim>::AddOperatorInterior( MatrixOperator<dim>* op )
{
    interior_lhs_.push_back(op);
}


template<size_t dim>
void LinearSystemAccumulator<dim>::AddOperatorPerimeter( VectorOperator<dim>* op )
{
    perimeter_rhs_.push_back(op);
}


template<size_t dim>
void LinearSystemAccumulator<dim>::AddOperatorInterior( VectorOperator<dim>* op )
{
    interior_rhs_.push_back(op);
}


template<size_t dim>
void LinearSystemAccumulator<dim>::FinaliseOperators()
{
}

  
template<size_t dim>
  void LinearSystemAccumulator<dim>::TimeIncrement(double64 dt)
  {
    for (auto op: interior_lhs_) {
      if (op->MultiplyWithTimeIncrement())
        op->TimeIncrement( dt );
    }
    
    for (auto op: interior_rhs_) {
      if (op->MultiplyWithTimeIncrement())
        op->TimeIncrement( dt );
    }
    
    for (auto op: perimeter_lhs_) {
      if (op->MultiplyWithTimeIncrement())
        op->TimeIncrement( dt );
    }
    
    for (auto op: perimeter_rhs_) {
      if (op->MultiplyWithTimeIncrement())
        op->TimeIncrement( dt );
    }
  }

template<size_t dim>
void LinearSystemAccumulator<dim>::AccumulateByStencil(SparseMatrix& lhs, std::vector<double64>& rhs)
{
  auto interior_elmts_end = gref_.InteriorElementsEnd();
  for (auto eit = gref_.InteriorElementsBegin(); eit != interior_elmts_end; ++eit) {
    for (auto op: interior_lhs_) {
      op->AccumulateStencil( *eit, lhs );
    }

    for (auto op: interior_rhs_) {
      op->AccumulateStencil( *eit, rhs );
    }
  }
  
  auto perimeter_elmts_end = gref_.PerimeterElementsEnd();
  for (auto eit = gref_.PerimeterElementsBegin(); eit != perimeter_elmts_end; ++eit) {
    for (auto op: perimeter_lhs_) {
      op->AccumulateStencil( *eit, lhs );
    }
    
    for (auto op: perimeter_rhs_) {
      op->AccumulateStencil( *eit, rhs );
    }
  }
}


template<size_t dim>
void LinearSystemAccumulator<dim>::AccumulateByFiniteVolume(SparseMatrix& lhs, std::vector<double64>& rhs)
{
  auto interior_nodes_end = gref_.InteriorNodesEnd();
  for (auto nit = gref_.InteriorNodesBegin(); nit != interior_nodes_end; ++nit) {
    for (auto op: interior_lhs_) {
      op->AccumulateFiniteVolume( *nit, lhs );
    }
    
    for (auto op: interior_rhs_) {
      op->AccumulateFiniteVolume( *nit, rhs );
    }
  }
  
  auto perimeter_nodes_end = gref_.PerimeterNodesEnd();
  for (auto nit = gref_.PerimeterNodesBegin(); nit != perimeter_nodes_end; ++nit) {
    for (auto op: perimeter_lhs_) {
      op->AccumulateFiniteVolume( *nit, lhs );
    }
    
    for (auto op: perimeter_rhs_) {
      op->AccumulateFiniteVolume( *nit, rhs );
    }
  }
}

  
template class LinearSystemAccumulator<1U>;
template class LinearSystemAccumulator<2U>;
template class LinearSystemAccumulator<3U>;

} // end csmp

