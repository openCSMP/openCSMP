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



template class LinearSystemAccumulator<1U>;
template class LinearSystemAccumulator<2U>;
template class LinearSystemAccumulator<3U>;

} // end csmp

