#ifndef CSMP_FLUID_COMPRESSIBILITY_SOURCE_LHS_H
#define CSMP_FLUID_COMPRESSIBILITY_SOURCE_LHS_H

#include "MatrixOperator.h"

namespace csmp {

/** 
  Volumetric source and sink terms that arise as the fluid expands or compresses due to pressure changes.
    
     @note this version uses sector-pore volume and total systems compressibilities of rocktypes that are associated with the elements.
     
     @note fluid pressure is either taken from the node(FV) or interpolated with the finite-element interpolation to sector integration points.
*/
template<size_t dim>
class FluidCompressibilitySourceLHS : public MatrixOperator<dim> {
  public:
    FluidCompressibilitySourceLHS( const csmp::INDEX<SCALAR,SECTOR_INTEGRATION_POINT>& key_SPV,
                                   const csmp::INDEX<SCALAR,ELEMENT>& key_PHI,
                                   const csmp::INDEX<SCALAR,ELEMENT>& key_CT,
                                   const csmp::INDEX<SCALAR,NODE>& key_PF0,
                                   const csmp::INDEX<SCALAR,NODE>& key_PF1,
                                   bool interpolate_pf_to_sector_ip );
    
    virtual ~FluidCompressibilitySourceLHS() {}
  
    virtual void AccumulateFiniteVolume( const Node<dim>&, SparseMatrix& ) const;

    virtual void AccumulateStencil( const Element<dim>&, SparseMatrix& ) const;
    // virtual void AccumulateStencil( const Face<dim>&, SparseMatrix& ) const;
    // virtual void AccumulateStencil( const InterFace<dim>&, SparseMatrix& ) const;

  private:
    const csmp::INDEX<SCALAR,SECTOR_INTEGRATION_POINT>& key_SPV_; ///<  sector pore volume

    const csmp::INDEX<SCALAR,ELEMENT>& key_PHI_; ///<  porosity
    const csmp::INDEX<SCALAR,ELEMENT>& key_CT_;  ///<  total system compressibility

    const csmp::INDEX<SCALAR,NODE>& key_PF0_; ///< previous fluid pressure
    const csmp::INDEX<SCALAR,NODE>& key_PF1_; ///< current fluid pressure
    const bool interpolate_pf_to_sector_ip_;  ///<  sector by sector computation of source term using pressure interpolation
};

} // end csmp

#endif /* CSMP_FLUID_COMPRESSIBILITY_SOURCE_LHS_H */

