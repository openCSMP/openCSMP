#ifndef INTEGRATION_POINT_COUNT_H
#define INTEGRATION_POINT_COUNT_H


namespace csmp{

  /// Data class that tracks number of element/facet/sector integration points in Element/Face/InterFace
  struct IntegrationPointCount{
    IntegrationPointCount() : simplexIntegrationPoints(0), sectorIntegrationPoints(0), facetIntegrationPoints(0) {}
    IntegrationPointCount( size_t siIntegrationPoints, size_t seIntegrationPoints, size_t faIntegrationPoints ) 
      : simplexIntegrationPoints(siIntegrationPoints), 
        sectorIntegrationPoints(seIntegrationPoints), 
        facetIntegrationPoints(faIntegrationPoints) {}

    IntegrationPointCount& operator = ( const IntegrationPointCount& ipc )
    {
      if( this != &ipc )
        {
          simplexIntegrationPoints = ipc.simplexIntegrationPoints;
          sectorIntegrationPoints = ipc.sectorIntegrationPoints;
          facetIntegrationPoints = ipc.facetIntegrationPoints;
        }
      return *this;
    }

    size_t simplexIntegrationPoints;
    size_t sectorIntegrationPoints;
    size_t facetIntegrationPoints;
  };



} // csmp

#endif