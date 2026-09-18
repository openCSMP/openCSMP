// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef INTEGRATION_POINT_COUNT_H
#define INTEGRATION_POINT_COUNT_H


namespace csmp {

  /// Data class that tracks number of element/facet/sector integration points in Element/Face/InterFace
  struct IntegrationPointCount {
    IntegrationPointCount() : simplexIntegrationPoints(0), sectorIntegrationPoints(0), facetIntegrationPoints(0) {}
    IntegrationPointCount( uint32_t siIntegrationPoints, uint32_t seIntegrationPoints, uint32_t faIntegrationPoints )
      : simplexIntegrationPoints(siIntegrationPoints),
        sectorIntegrationPoints(seIntegrationPoints), 
        facetIntegrationPoints(faIntegrationPoints) {}

    IntegrationPointCount& operator = ( const IntegrationPointCount& ipc ) {
        if( this != &ipc )
          {
            simplexIntegrationPoints = ipc.simplexIntegrationPoints;
            sectorIntegrationPoints = ipc.sectorIntegrationPoints;
            facetIntegrationPoints = ipc.facetIntegrationPoints;
          }
        return *this;
      }

    uint32_t simplexIntegrationPoints;
    uint32_t sectorIntegrationPoints;
    uint32_t facetIntegrationPoints;
  };



} // csmp

#endif
