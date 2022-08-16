#ifndef  CSMP_EXAMPLES_H
#define  CSMP_EXAMPLES_H

// this is where you declare your examples
// =======================================

#include "SteadyStatePressure_Example.h"       // steady state calculation 2D
#include "TransientPressure_Example.h"         // transient calculation 2D with well
#include "PassiveAdvectionOfTracer_Example.h"  // single phase advection in 3D
#include "VariablesBasic_Example.h"         // variable operations EXTEND!
#include "Region_Example.h"                 // computations with regions
#include "PermeabilityTensor_Example.h"     // trivial application of csmp variables
#include "TransientPressure_Example.h"      // use property placement as template arg
#include "StatisticalAnalyzer_Example.h"    // creation of histograms
#include "Visitor_Example.h"                // use & implementation of visitor
#include "DenseMatrix_Example.h"            // illustration of matrix features
#include "ParallelPlateFracture_Example.h"  // fracture permeability calculation, linear or quadratic pressure
#include "PressureDiffusion_Example.h"      //  2D transient calculation of fluid pressure distribution
#include "LinearElasticity_Example.h"       // mechanical computations with regions
#include "TopographyDrivenFlow_Example.h"   // steady state
#include "TemperatureDensityPressure_Example.h" // 1D H2O
#include "RhinoMesh_Example.h"              // DFN rhino mesh import
#include "StreamFunction_Example.h"         // post processing of streamlines
#include "TimeSteppingApproaches_Example.h" // 1D comparison
#include "ErrorMetric_Example.h"            // 3D fluid pressure
// deprecated: #include "ReadingBinaries_Example.h"        // read from binaries created in SteadyStatePressureToVset_Example
#include "StokesDiscrepancyMeasure_Example.h" // Darcy velocity, stokes flow
#include "StokesDiscrepancyMeasureQuadratic_Example.h" // as above, but quadratic
#include "EffectiveStressDilatation2D_Example.h" 
#include "PolicyBased_Example.h"            // C++ example
#include "CuriouslyRecurringTemplate_Example.h" // C++ technique
#include "TDDUnitTest_Example.h"            // C++ unit test / TDD
#include "RegionProperties_Example.h"
#include "QuadraticPressure_parallelPlatePermeability_Example.h"
#include "TemplatizedIndex_Example.h"
#include "ModelANSYS_Example.h"
#include "Experimental_Example.h"
#include "EffectiveStressDilatation2D_Example.h"
#include "Intrepid_Example.h"
#include "SKUA_Example.h"
#include "Experimental_Example.h"
#include "Averaging_Example.h"
#include "Tutorial1_Example.h"
#include "Tutorial2_Example.h"
#include "Tutorial3_Example.h"
#include "Tutorial4_Example.h"
#include "Geothermal_Example.h"
#include "LinearSolver_Example.h"
#include "EclipseMeshInterface_Example.h"
#include "DESAdvectionDiffusion2D_Example.h"
#include "DESAdvectionDiffusion3D_Example.h"
#include "DES2PhaseSlightlyCompressibleFlow2D_Example.h"
#include "DES2PhaseSlightlyCompressibleFlow3D_Example.h"
#include "DES2PhaseFlowWithSplitBoundary_Example.h"
#include "Variables_Example.h"
#include "SlopeMechanics_Example.h" // developed for UoM SMD course, pore-pressure and gravitational loading
#include "QuadrilateratorToCSMPbinary_Example.h" // demonstrating the Quadrilaterator 2D mesher

#endif
