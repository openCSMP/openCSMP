// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "CapillaryNumberVisitor.h"
#include "Model.h"
#include "NodeManifold.h"
#include "ScalarVariable.h"
#include "VectorVariable.h"

namespace csmp{


template<uint32_t dim>
CapillaryNumberVisitor<dim>::CapillaryNumberVisitor( Model<dim>& model, const char* velTag,
                                                     const char* viscTag, const char* nCapTag, double ift )
    : Visitor<dim>( MODEL, ELEMENT ), model_(model), velKey_( model.Database().StorageKey(velTag) ),
      viscKey_( model.Database().StorageKey(viscTag) ), nCapKey_( model.Database().StorageKey(nCapTag) ), ift_(ift)
    {
    }


template<uint32_t dim>
void CapillaryNumberVisitor<dim>::Visit( Element<dim>* element )
  {
    const double visc( element->Read(viscKey_) );
    ScalarVariable capillaryNumber( PLAIN, 0.0 );
    VectorVariable<dim> velocity;
    element->Read( velKey_, velocity );

    // Capillary number (Saffman and Taylor, 1958) = v * mu * del(frap) / ift
    capillaryNumber = visc * velocity.Length() / ift_;

    element->Store( nCapKey_, capillaryNumber );
  }

template class CapillaryNumberVisitor<1>;
template class CapillaryNumberVisitor<2>;
template class CapillaryNumberVisitor<3>;

}

