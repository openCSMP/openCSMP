// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "ConductivityVisitor.h"

#include "CapillaryNumberVisitor.h"
#include "Model.h"
#include "ScalarVariable.h"
#include "VectorVariable.h"
#include "ErrorHandler.h"

namespace csmp {


template<uint32_t dim, template<uint32_t> class CELL>
ConductivityVisitor<dim,CELL>::ConductivityVisitor(Model<dim>& model,
                                              const char* specific_saturated_hydraulic_conductivity, // this is without mult. by density
                                              const char* permeability,
                                              const char* viscosity,
                                              const char* density,
                                              const char* saturated_hydraulic_conductivity,
                                              const char* compressibility,
                                              const char* porosity,
                                              const char* specific_saturated_hydraulic_diffusivity)
    : Visitor<dim>( MODEL, ELEMENT ), model_(model),
      sshcKey_( model.Database().StorageKey(specific_saturated_hydraulic_conductivity) ),
      kKey_   ( model.Database().StorageKey(permeability) ),
      muKey_  ( model.Database().StorageKey(viscosity) ),
      rhoKey_ ( ((density==NULL) ? csmp::Index() : model.Database().StorageKey(density))),
      shcKey_ ( ((saturated_hydraulic_conductivity==NULL) ? csmp::Index() : model.Database().StorageKey(saturated_hydraulic_conductivity))),
      ctKey_  ( ((compressibility==NULL) ? csmp::Index() : model.Database().StorageKey(compressibility))),
      phiKey_ ( ((porosity==NULL) ? csmp::Index() : model.Database().StorageKey(porosity))),
      sshdKey_( ((specific_saturated_hydraulic_diffusivity==NULL) ? csmp::Index() : model.Database().StorageKey(specific_saturated_hydraulic_diffusivity)))
{
    if ( sshcKey_.place != ELEMENT )
        throw csmp::Exception( ERROR, "ConductivityVisitor", specific_saturated_hydraulic_conductivity, " must be an element property." );
    if ( kKey_.place != ELEMENT )
        throw csmp::Exception( ERROR, "ConductivityVisitor", permeability, " must be an element property." );
    if ( muKey_.place != NODE )
        throw csmp::Exception( ERROR, "ConductivityVisitor", viscosity, " must be a node property." );
    if ( rhoKey_ != csmp::Index() && rhoKey_.place != ELEMENT )
        throw csmp::Exception( ERROR, "ConductivityVisitor", density, " must be an element property." );
    if ( shcKey_ != csmp::Index() && shcKey_.place != ELEMENT )
        throw csmp::Exception( ERROR, "ConductivityVisitor", saturated_hydraulic_conductivity, " must be an element property." );
    if ( ctKey_ != csmp::Index()  && ctKey_.place != ELEMENT )
        throw csmp::Exception( ERROR, "ConductivityVisitor", compressibility, " must be an element property." );
    if ( phiKey_ != csmp::Index() && phiKey_.place != ELEMENT )
        throw csmp::Exception( ERROR, "ConductivityVisitor", porosity, " must be an element property." );
    if ( sshdKey_ != csmp::Index() && sshdKey_.place != ELEMENT )
        throw csmp::Exception( ERROR, "ConductivityVisitor", specific_saturated_hydraulic_diffusivity, " must be an element property." );

}




template<uint32_t dim, template<uint32_t> class CELL>
void ConductivityVisitor<dim,CELL>::Visit( CELL<dim>* cell )
{
    // first we calculate sshc
    cell->PropertyValueAtBaryCenter( muKey_, result );
    result = cell->Read( kKey_ ) / result();
    cell->Store( sshcKey_, result );

    // if key exists, calculate sshd
    if (sshdKey_!=csmp::Index()){
        resultd = result;
        resultd /= cell->Read(phiKey_);
        resultd /= cell->Read(ctKey_);
        cell->Store( sshdKey_, resultd);
    }

    // if key exists, calculate shc
    if (shcKey_!=csmp::Index()){
        result=result* cell->Read(rhoKey_);
        cell->Store( shcKey_, result );
    }
}


template class ConductivityVisitor<1,Element>;
template class ConductivityVisitor<2,Element>;
template class ConductivityVisitor<3,Element>;

template class ConductivityVisitor<1,Face>;
template class ConductivityVisitor<2,Face>;
template class ConductivityVisitor<3,Face>;

template class ConductivityVisitor<1,InterFace>;
template class ConductivityVisitor<2,InterFace>;
template class ConductivityVisitor<3,InterFace>;

} // end csmp

