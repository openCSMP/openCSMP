// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "FracturePecletNumberVisitor.h"
#include "Model.h"
#include "ScalarVariable.h"

namespace csmp{


template<uint32_t dim>
FracturePecletNumberVisitor<dim>::FracturePecletNumberVisitor( Model<dim>& model, const char* velocityTag,
                                                               const char* pecletNumberTag, double nodalSourceSink )
    : Visitor<dim>( MODEL, ELEMENT ), model_( model ), velKey_( model.Database().StorageKey(velocityTag) ),
      pecletNumberKey_( model.Database().StorageKey(pecletNumberTag) ), sourceSink_(nodalSourceSink)
    {
    }


  template<uint32_t dim>
  void FracturePecletNumberVisitor<dim>::Visit( Element<dim>* element )
    {
      VectorVariable<dim> nodalVelocity;
      element->Read(velKey_, nodalVelocity );
      FracNPe_ = nodalVelocity.Length() / sourceSink_;
      element->Store( pecletNumberKey_, FracNPe_ );
    }

  template class FracturePecletNumberVisitor<1>;
  template class FracturePecletNumberVisitor<2>;
  template class FracturePecletNumberVisitor<3>;

}

