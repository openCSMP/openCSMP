// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "FractureBrooksCoreyParameterVisitor.h"
#include "Model.h"
#include "ScalarVariable.h"
#include "VTU_Interface.h"

namespace csmp{


template<uint32_t dim>
FractureBrooksCoreyParameterVisitor<dim>::FractureBrooksCoreyParameterVisitor( Model<dim>& model, const char* frapTag,
                                                                               const char* brooksCoreyLambdaTag, double meanPoreRadius )
  : Visitor<dim>( MODEL, ELEMENT ), model_( model ), frapKey_( model.Database().StorageKey(frapTag) ),
    brooksCoreyLambdaKey_( model.Database().StorageKey(brooksCoreyLambdaTag) ), meanPoreRadius_(meanPoreRadius)
  {
  }



  template<uint32_t dim>
  void FractureBrooksCoreyParameterVisitor<dim>::Visit( Element<dim>* element )
  {
    const double frap( element->Read(frapKey_) );
    double lambda( element->Read( brooksCoreyLambdaKey_ ) );

    if ( frap >= 2 *meanPoreRadius_ )
      lambda = 0.;

    element->Store( brooksCoreyLambdaKey_, makeScalar( PLAIN, lambda ) );
  }

  template class FractureBrooksCoreyParameterVisitor<1>;
  template class FractureBrooksCoreyParameterVisitor<2>;
  template class FractureBrooksCoreyParameterVisitor<3>;

} //csmp
