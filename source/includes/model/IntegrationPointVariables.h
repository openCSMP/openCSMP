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

    IntegrationPointVariables( size_t scalarsSimplex,
                               size_t vectorsSimplex,
                               size_t tensorsSimplex,
                               size_t arrayCountSimplex,
                               size_t arrayLengthSimplex,
                               size_t flaggedArrayCountSimplex,
                               size_t flaggedArrayLengthSimplex,
                               size_t totalDataDepthSimplex,
                               size_t totalFlagDepthSimplex,
                               size_t scalarsSector,
                               size_t vectorsSector,
                               size_t tensorsSector,
                               size_t arrayCountSector,
                               size_t arrayLengthSector,
                               size_t flaggedArrayCountSector,
                               size_t flaggedArrayLengthSector,
                               size_t totalDataDepthSector,
                               size_t totalFlagDepthSector,
                               size_t scalarsFacet,
                               size_t vectorsFacet,
                               size_t tensorsFacet,
                               size_t arrayCountFacet,
                               size_t arrayLengthFacet,
                               size_t flaggedArrayCountFacet,
                               size_t flaggedArrayLengthFacet,
                               size_t totalDataDepthFacet,
                               size_t totalFlagDepthFacet
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
