// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef  CSMP_EXAMPLES_H
#define  CSMP_EXAMPLES_H

/// @file
/// @brief Example driver includes for the CSMP example suite.
///
/// @details Each include corresponds to a self-contained example demonstrating
/// a specific CSMP capability. Add new examples here in alphabetical order.
///
/// @note ReadingBinaries_Example.h is deprecated and excluded.

#include "Averaging_Example.h"                                  ///< averaging of distributed variables
#include "CuriouslyRecurringTemplate_Example.h"                 ///< C++ CRTP technique
#include "DESAdvectionDiffusion2D_Example.h"                    ///< DES advection-diffusion in 2D
#include "DESAdvectionDiffusion3D_Example.h"                    ///< DES advection-diffusion in 3D
#include "DES2PhaseFlowWithSplitBoundary_Example.h"             ///< two-phase flow with split boundary conditions
#include "DES2PhaseSlightlyCompressibleFlow_Example.h"          ///< two-phase slightly compressible flow
#include "DenseMatrix_Example.h"                                ///< illustration of dense matrix features
#include "EclipseMeshInterface_Example.h"                       ///< Eclipse mesh import and export
#include "EffectiveStressDilatation2D_Example.h"                ///< effective stress and dilatation in 2D
#include "ErrorMetric_Example.h"                                ///< 3D fluid pressure error metric
#include "Experimental_Example.h"                               ///< experimental / work-in-progress features
#include "Geothermal_Example.h"                                 ///< geothermal flow simulation
#include "Interfaces_Example.h"                                 ///< interface element handling
#include "Intrepid_Example.h"                                   ///< Intrepid integration
#include "LinearElasticity_Example.h"                           ///< mechanical computations with regions
#include "LinearSolver_Example.h"                               ///< linear solver demonstration
#include "ModelANSYS_Example.h"                                 ///< ANSYS mesh import
#include "ParallelPlateFracture_Example.h"                      ///< fracture permeability, linear or quadratic pressure
#include "PassiveAdvectionOfTracer_Example.h"                   ///< single-phase passive tracer advection in 3D
#include "PermeabilityTensor_Example.h"                         ///< trivial application of CSMP tensor variables
#include "PETSc_Example.h"                                      ///< PETSc linear solver integration
#include "PolicyBased_Example.h"                                ///< policy-based C++ design example
#include "PoreRadiusAnd_Pc_Example.h"                           ///< pore radius and capillary pressure
#include "PoroElasticity_Example.h"                             ///< coupling between stress and pore pressure in porous medium
#include "PressureDependentBlackOilProperties_Example.h"        ///< pressure-dependent black-oil fluid properties
#include "PressureDiffusion_Example.h"                          ///< 2D transient fluid pressure diffusion
#include "QuadraticPressure_parallelPlatePermeability_Example.h"///< quadratic pressure parallel-plate permeability
#include "QuadrilateratorToCSMPbinary_Example.h"                ///< Quadrilaterator 2D mesher to CSMP binary
#include "Region_Example.h"                                     ///< computations with regions
#include "RegionMonitor_Example.h"                              ///< monitoring of region-averaged quantities
#include "RegionProperties_Example.h"                           ///< region property assignment and retrieval
#include "RhinoMesh_Example.h"                                  ///< DFN Rhino mesh import
#include "SKUA_Example.h"                                       ///< SKUA geomodel mesh import
#include "SlopeMechanics_Example.h"                             ///< pore-pressure and gravitational loading (UoM SMD)
#include "StaggeredGridStokesSolver_Example.h"                  ///< staggered-grid Stokes solver
#include "StatisticalAnalyzer_Example.h"                        ///< creation and analysis of histograms
#include "StokesDiscrepancyMeasure_Example.h"                   ///< Darcy velocity vs. Stokes flow discrepancy
#include "StokesDiscrepancyMeasureQuadratic_Example.h"          ///< as above, with quadratic pressure elements
#include "StreamFunction_Example.h"                             ///< post-processing of streamlines
#include "TDDUnitTest_Example.h"                                ///< unit test and TDD demonstration
#include "TemplatizedIndex_Example.h"                           ///< templatised index usage
#include "TemperatureDensityPressure_Example.h"                 ///< 1D H2O temperature-density-pressure coupling
#include "ThermalConvectionETHZ_Example.h"                      ///< thermal convection (ETHZ benchmark)
#include "TimeSteppingApproaches_Example.h"                     ///< 1D comparison of time-stepping schemes
#include "TopographyDrivenFlow_Example.h"                       ///< topography-driven steady-state flow
#include "Tractions_Example.h"                                  ///< traction boundary condition application
#include "TransientPressure_Example.h"                          ///< transient 2D pressure calculation with well
#include "Triangulator_Example.h"                               ///< Triangulator mesher and 2D steady-state pressure
#include "Tutorial1_Example.h"                                  ///< introductory tutorial 1
#include "Tutorial2_Example.h"                                  ///< introductory tutorial 2
#include "Tutorial3_Example.h"                                  ///< introductory tutorial 3
#include "Tutorial4_Example.h"                                  ///< introductory tutorial 4
#include "UG4_ProMeshOutput_Example.h"                          ///< UG4 ProMesh output interface
#include "VariableManagement_Example.h"                         ///< reading and writing variables to CSMP binary files
#include "Variables_Example.h"                                  ///< general variable operations
#include "VariablesBasic_Example.h"                             ///< basic variable operations
#include "Visitor_Example.h"                                    ///< use and implementation of the visitor pattern
#include "CVFEM_fault_well_lithium_example.h"
#include "CVFEM_topo_magma_air_example.h"
#include "CVFEM_geothermal_cooling_magma_chamber_example.h"
#endif
