// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "BO_LookupVisitor.h"

#include "Model.h"
#include "UndersaturatedBOCorrelations.h"


namespace csmp {

template<uint32_t dim>
BO_LookupVisitor<dim>::BO_LookupVisitor( const Model<dim>& model,
                                         BO_PropertyCalculatorInterface& bo,
                                         ReservoirWaterPropertyCalculatorInterface& water )
                                        
  : Visitor<dim>( MODEL, NODE ), 
    bo_(bo),
    water_(water)
 {
    pressureFluidKey_ =  model.Database().StorageKey( "fluid pressure" );
    tempFluidKey_ = model.Database().StorageKey("fluid temperature");
    
    viscosityOilKey_ =  model.Database().StorageKey( "viscosity oil" );
    densityOilKey_ = model.Database().StorageKey( "density oil" );
    
    viscosityWaterKey_ = model.Database().StorageKey( "viscosity water" );
    densityWaterKey_ = model.Database().StorageKey( "density water" );
 }
 
 

template<uint32_t dim>
void  BO_LookupVisitor<dim>::Visit( Node<dim>* node )
 {
    node->Read(pressureFluidKey_,presFluid_);
    node->Read(tempFluidKey_, tempFluid_);
    
    bo_.Update(presFluid_(), tempFluid_());
    water_.Update(presFluid_(), tempFluid_());
    
    node->Store(viscosityOilKey_, makeScalar(PLAIN, bo_.getUndersaturatedOilViscosity()));
    node->Store(densityOilKey_, makeScalar(PLAIN, bo_.getUndersaturatedOilDensity()));
    node->Store(compressibilityOilKey_, makeScalar(PLAIN, bo_.getUndersaturatedOilCompressibility()));
    
    node->Store(viscosityWaterKey_, makeScalar(PLAIN, water_.getReservoirWaterViscosity()));
    node->Store(densityWaterKey_, makeScalar(PLAIN, water_.getReservoirWaterDensity()));
    node->Store(compressibilityWaterKey_, makeScalar(PLAIN, water_.getReservoirWaterCompressibility()));
    
 }
  
  
  
template<uint32_t dim>
void  BO_LookupVisitor<dim>::Visit( Model<dim>* m )
 {
    csmp::Region<dim>&  sg(m->Region("Model"));
    for ( auto el_it=sg.CellsBegin(); el_it!=sg.CellsEnd(); el_it++ )
      (*el_it)->Accept( *this );
 }

template class BO_LookupVisitor<1>;
template class BO_LookupVisitor<2>;
template class BO_LookupVisitor<3>;

} //end csmp
