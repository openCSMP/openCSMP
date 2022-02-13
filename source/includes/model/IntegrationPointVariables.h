#ifndef INTEGRATION_POINT_VARIABLES_H
#define INTEGRATION_POINT_VARIABLES_H

#include "LocalVariables.h"

namespace csmp {


  /// Data class to store physical variable count at given integration points
  struct IntegrationPointVariables {
    IntegrationPointVariables() 
      : ipvSimplex  (),
        ipvSector   (),
        ipvFacet    ()
    {}

    IntegrationPointVariables( const IntegrationPointVariables& ivs ) 
      : ipvSimplex  (ivs.ipvSimplex),
        ipvSector   (ivs.ipvSector),
        ipvFacet    (ivs.ipvFacet)
    {}

    IntegrationPointVariables( uint32_t scalarsSimplex,
                               uint32_t vectorsSimplex,
                               uint32_t tensorsSimplex,
                               uint32_t arrayCountSimplex,
                               uint32_t arrayLengthSimplex,
                               uint32_t flaggedArrayCountSimplex,
                               uint32_t flaggedArrayLengthSimplex,
                               uint32_t totalDataDepthSimplex,
                               uint32_t totalFlagDepthSimplex,
                               uint32_t scalarsSector,
                               uint32_t vectorsSector,
                               uint32_t tensorsSector,
                               uint32_t arrayCountSector,
                               uint32_t arrayLengthSector,
                               uint32_t flaggedArrayCountSector,
                               uint32_t flaggedArrayLengthSector,
                               uint32_t totalDataDepthSector,
                               uint32_t totalFlagDepthSector,
                               uint32_t scalarsFacet,
                               uint32_t vectorsFacet,
                               uint32_t tensorsFacet,
                               uint32_t arrayCountFacet,
                               uint32_t arrayLengthFacet,
                               uint32_t flaggedArrayCountFacet,
                               uint32_t flaggedArrayLengthFacet,
                               uint32_t totalDataDepthFacet,
                               uint32_t totalFlagDepthFacet
                            )
      :
      ipvSimplex( scalarsSimplex,
                  vectorsSimplex,
                  tensorsSimplex,
                  arrayCountSimplex,
                  arrayLengthSimplex,
                  flaggedArrayCountSimplex,
                  flaggedArrayLengthSimplex,
                  totalDataDepthSimplex,
                  totalFlagDepthSimplex
                ),
      ipvSector ( scalarsSector,
                  vectorsSector,
                  tensorsSector,
                  arrayCountSector,
                  arrayLengthSector,
                  flaggedArrayCountSector,
                  flaggedArrayLengthSector,
                  totalDataDepthSector,
                  totalFlagDepthSector
                ),
      ipvFacet  ( scalarsFacet,
                  vectorsFacet,
                  tensorsFacet,
                  arrayCountFacet,
                  arrayLengthFacet,
                  flaggedArrayCountFacet,
                  flaggedArrayLengthFacet,
                  totalDataDepthFacet,
                  totalFlagDepthFacet
                )
    {}

    IntegrationPointVariables& operator=( const IntegrationPointVariables& ivs )
    { 
      if( this != &ivs ) 
      { 
        ipvSimplex = ivs.ipvSimplex;
        ipvSector  = ivs.ipvSector;
        ipvFacet   = ivs.ipvFacet;
      } 
      return *this; 
    }

    bool Empty() const { return ( ipvSimplex.Empty() && ipvSector.Empty() && ipvFacet.Empty() ); }

    LocalVariables ipvSimplex,
                   ipvSector,
                   ipvFacet;
  };


} // csmp

#endif
