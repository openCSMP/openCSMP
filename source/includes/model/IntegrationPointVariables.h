#ifndef INTEGRATION_POINT_VARIABLES_H
#define INTEGRATION_POINT_VARIABLES_H

#include "LocalVariables.h"

namespace csmp {

  /** Data container storing physical variable counts at given integration points
      see Paluszny et al. (2007, Geofluids)
 */
  struct IntegrationPointVariables {
  
    using int_type = LocalVariables::int_type; ///< unsigned integer type that is big enough to hold 'totalDataDepth'

    IntegrationPointVariables()
      : ipvCell(),   ///< variables placed on the integration points of Element, Face or InterFace
        ipvSector(), ///< variables on finite-volume sector integration points
        ipvFacet()   ///< variables on finite-volume facet integration points
    {}

    IntegrationPointVariables( const IntegrationPointVariables& ivs ) 
      : ipvCell  (ivs.ipvCell),
        ipvSector   (ivs.ipvSector),
        ipvFacet    (ivs.ipvFacet)
    {}

    IntegrationPointVariables( int_type scalarsSimplex,
                               int_type vectorsSimplex,
                               int_type tensorsSimplex,
                               int_type arrayCountSimplex,
                               int_type arrayLengthSimplex,
                               int_type flaggedArrayCountSimplex,
                               int_type flaggedArrayLengthSimplex,
                               int_type totalDataDepthSimplex,
                               int_type totalFlagDepthSimplex,
                               int_type scalarsSector,
                               int_type vectorsSector,
                               int_type tensorsSector,
                               int_type arrayCountSector,
                               int_type arrayLengthSector,
                               int_type flaggedArrayCountSector,
                               int_type flaggedArrayLengthSector,
                               int_type totalDataDepthSector,
                               int_type totalFlagDepthSector,
                               int_type scalarsFacet,
                               int_type vectorsFacet,
                               int_type tensorsFacet,
                               int_type arrayCountFacet,
                               int_type arrayLengthFacet,
                               int_type flaggedArrayCountFacet,
                               int_type flaggedArrayLengthFacet,
                               int_type totalDataDepthFacet,
                               int_type totalFlagDepthFacet
                            )
      :
      ipvCell( scalarsSimplex,
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
        ipvCell = ivs.ipvCell;
        ipvSector  = ivs.ipvSector;
        ipvFacet   = ivs.ipvFacet;
      } 
      return *this; 
    }

    bool Empty() const { return ( ipvCell.Empty() && ipvSector.Empty() && ipvFacet.Empty() ); }

    LocalVariables ipvCell,
                   ipvSector,
                   ipvFacet;
  };


} // csmp

#endif
