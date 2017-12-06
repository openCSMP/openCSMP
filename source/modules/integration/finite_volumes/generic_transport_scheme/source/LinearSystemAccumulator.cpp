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
#ifdef CSMP_WITH_SAMG_SOLVER
#include "SAMG_Settings.h"
#include "SAMG_Solver.h"
#else
#include "GaussJordan_Solver.h"
#include "Meschach_Solver.h"
#endif

using namespace std;

namespace csmp {


template<size_t dim>
LinearSystemAccumulator<dim>::LinearSystemAccumulator( Model<dim>& model, const char* region )
  : model_(model),
    gref_(model_.Region(region)),
    fvs_(gref_.Nodes()),
    LHS_(fvs_),
    RHS_(fvs_)
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
    auto matrix_operator_ordering = [](const MatrixOperator<dim>* x, const MatrixOperator<dim>* y) {
        return x->AccumulationMode() < y->AccumulationMode();
    };
    std::sort(interior_lhs_.begin(), interior_lhs_.end(), matrix_operator_ordering);
    std::sort(perimeter_lhs_.begin(), perimeter_lhs_.end(), matrix_operator_ordering);

    auto vector_operator_ordering = [](const VectorOperator<dim>* x, const VectorOperator<dim>* y) {
        return x->AccumulationMode() < y->AccumulationMode();
    };
    std::sort(interior_rhs_.begin(), interior_rhs_.end(), vector_operator_ordering);
    std::sort(perimeter_rhs_.begin(), perimeter_rhs_.end(), vector_operator_ordering);
}


template<size_t dim>
void LinearSystemAccumulator<dim>::AccumulateByStencil()
{
}


template<size_t dim>
void LinearSystemAccumulator<dim>::AccumulateByFiniteVolume()
{
}

  
template<size_t dim>
void LinearSystemAccumulator<dim>::EnsureSolver()
{
  if (!solver_) {
#ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings settings;
    settings.Set_iout1(1);
    settings.Set_iout2(0);
    solver_= std::unique_ptr<Solver>(new SAMG_Solver(&settings));
#else
    solver_= std::unique_ptr<Solver>(new CSMP_DEFAULT_LINEAR_SOLVER);
#endif
  }
}

template<size_t dim>
void LinearSystemAccumulator<dim>::SolveSystem()
{
  EnsureSolver();
  
#ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings* settings = static_cast<SAMG_Settings*>(solver_->GetSolverSettings());
    cout <<"\n\nLinearSystemAccumulator::SolveMatrixEquation: calling solver instance: ";
    cout << settings->GetSolverInstance() << endl;
    cout <<"\tSAMG settings:";
    cout <<"\n\t\tiswit  = " << settings->Get_iswit();
    cout <<"\n\t\titypu  = " << settings->Get_ifirst();
    cout <<"\n\t\tlevelx = " << settings->Get_levelx();
    cout << endl;
#endif
  
  solver_->Solve( LHS_, RHS_, RESULT_ );
}
  

template<size_t dim>
void LinearSystemAccumulator<dim>::WriteResultIntoModel()
{
}


template class LinearSystemAccumulator<1U>;
template class LinearSystemAccumulator<2U>;
template class LinearSystemAccumulator<3U>;

} // end csmp

