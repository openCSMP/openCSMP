// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

//
//  DFM_PermeabilityAndFluxRatioAnalyzer.h
//  CSMP_DFM_Upscaling
//
//

#ifndef DFM_PERMEABILITY_AND_FLUX_RATIO_ANALYZER_H
#define DFM_PERMEABILITY_AND_FLUX_RATIO_ANALYZER_H

#include "RegionMonitor.h"
#include "StatisticalAnalyzer.h"
#include "Standard_IO_Handler.h"
#include "Model.h"
#include "CSMP_physical_constants.h"

namespace csmp {

/// Flow-based analysis of DFM: measure k_equivalent and qf/qm ratio (also for multilayer rock matrix) / 'excess permeability' of fractured rock
template<uint32_t dim>
class DFM_PermeabilityAndFluxRatioAnalyzer {
  public:
    DFM_PermeabilityAndFluxRatioAnalyzer( const char* model_name, bool verbose=false );
  
    /// equivalent permeability and fracture-matrix flux ratio (qf/qf) output to file '*-equivalent-properties.text'
    void DiagonalTensorAnalysis( const char* model_name, bool with_vtk_output );
  
    /// equivalent permeability and fracture-matrix flux ratio (qf/qf) output to file '*-equivalent-properties.text'
    void FullTensorAnalysis( const char* model_name );
  
    /**
        - Groups whatever (non-fracture) regions identified as matrix into region "MATRIX"
        - equivalent permeability analysis
    */
    void MatrixOnlyAnalysis( const char* model_name, bool with_vtk_output );

    /// by flow-based upscaling: works for multilayer matrix and any kind of boundary conditions
    double EquivalentPermeability( const char* region, double& flux_through_model );
  
    /// as above, but for csmp:Boundary objects with more flexibility
    double EquivalentPermeability( BOX_BOUNDARY inflow_boundary,
                                   BOX_BOUNDARY outflow_boundary,
                                   double& flux_through_model );

    /// ratio between fracture- and pore void space
    double FractureMatrixVoidRatio( bool verbose=false ) const;
  
    double FractureVolume( bool verbose=false ) const;
  
    /// initialize element variable that shall identify unique model regions
    void InitializeRegionIdentifiers( const char* region_identifying_variable );
  
  public:
    /// testing: move these methods later to high-level utilities
    bool TestBoundaryIntegrity();
    void OutputBoxBoundaryFlagsAsNumbers( const char* variable );
  
  private:
    void   PrintSpecifications() const;
    /// tolerance refers to how far out-of-plane nodes can be in terms of dim-max of boundary
    bool   DetectOutOfPlaneNodesAtBoxBoudaries( const char* model_name, double relative_tolerance=0.01 ) const;
    /// returns fracture area
    double IdentifyFractureAndMatrixRegions();
    double RegionVolume( const char* region ) const;
    void   ComputeTransmissivity();
    void   ComputeVelocityAndVolumeFlux( const char* region );
    void   OutputVelocityHistogramToMaple( size_t i ) const;
    void   PrintResultsToScreen( const char* model_name ) const;
    void   WriteResultsToFile( const char* file_name ) const;
    void   WriteResultsToVTK( const char* model_name, size_t i );
  
    /// integration of fluxes across box-model boundaries using the FEM method and including lower-dimensional elements
    double BoxBoundaryFluxFEM( BOX_BOUNDARY ) const;
  
  private:
   std::string              matrix_region_, fracture_region_;
   static constexpr double  fluid_viscosity_      = 1.6e-3,            ///< Pa.s  water at room temperature
                            farfield_pf_gradient_ = ACC_GRAVITY * 1e3; ///< 1000 kg/m3 * g;
   double                   xsect_area_;
   csmp::Model<dim>                        model_;
   RegionMonitor<dim>                      monitor_;
   Standard_IO_Handler                     stdio_;
   StatisticalAnalyzer<dim>                flux_histogram_;
   std::vector<std::pair<double,double> >  velo_bins_;
   std::vector<double>                     model_dimensions_,
                                           k_equivalent_,
                                           flux_ratio_,
                                           matrix_permeability_;
   double                                  model_volume_,
                                           model_pore_volume_,
                                           fracture_matrix_void_ratio_,
                                           fracture_surface_area_;
   bool verbose_;
};


} // end csmp


#endif
