//
//  DirichletPressureBoxModel_VVCase.cpp
//  CSMP_GitHub_UnitTests-Intel
//
//  Created by Stephan Matthai on 1/02/2018.
//  Copyright © 2018 Stephan Matthai. All rights reserved.
//

#include "DirichletPressureBoxModel_VVCase.h"
#include "ANSYS_Interface.h"
#include "ANSYS_Model3D.h"
#include "MeshDiagnostics.h"
#include "ModelTopology.h"
//#include "CSMP_highLevelUtilities.h"
#include "meshManagementUtilities.h"

#ifdef CSMP_WITH_SAMG_SOLVER
#include "SAMG_Settings.h"
#include "SAMG_Solver.h"
#else
#include "LinearSolver.h"
#endif

#include "Region.h"
#include "PDE_Integrator.h"
#include "NumIntegral_dNT_dN_dV.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_NT_op_N_dV.h"
#include "LHS_Integral_dNT_dN_dV.h"
#include "Integral_NT_op_N_dV.h"
#include "LinearSolver.h"

#include "VTK_Interface.h"

using namespace std;

namespace csmp {

#define DIM 3


DirichletPressureBoxModel_VVCase::DirichletPressureBoxModel_VVCase( const char* model )
 : model_name_(model), verbose_(true)
 {
 }
 
 
void DirichletPressureBoxModel_VVCase::run()
  {
     // testing ANSYS derived isoparametric model
     TestModelFromANSYS();
     // testing model with analytic element integration
     TestModelFromANSYS_AnalyticallyIntegrated();
  }


/**
    Input model 3D from ANSYS using the default constructor for such models.
    Standard isoparametric set of elements for the computation.
*/
void DirichletPressureBoxModel_VVCase::TestModelFromANSYS()
 {
    // standard ANSYS model case
    const bool binary_file(true);

    ANSYS_Model3D  model( model_name_.c_str(), "DirichletPressureBoxModel_VVCase-variables.txt", binary_file );
    printModelDimensions( model, true );
  
    // configuration and steady-state pressure computation
    // ----------------------------------------------------------------------
  
    // testing
    // ------------------------------
    // 1. duplicate elements
    Region<DIM> model_domain = model.Region("Model");
    size_t duplicates = detectDuplicateCells<DIM,Element>( model_domain.CellsBegin(), model_domain.CellsEnd(), false );
    _test( duplicates == 0 );

    MeshDiagnostics<3U>  diagnostics;
   
    // 2. disconnected nodes etc.
    // _test( diagnostics.ScrutinizeMesh(model) );
    _test( !diagnostics.DetectPotentiallyMisnumberedElements(model) );

    // computed pressure field
    // ------------------------
      model.InputPropertyValue( "conductivity", makeScalar(PLAIN,1.0e-12) );
      model.InputPropertyValue( "fluid volume source", makeScalar(PLAIN,0.) );
      // !5Pa drop only
      model.InputPropertyValue( "fluid pressure", makeScalar(PLAIN,0.) );
      model.InputBoundaryValue( LEFT, "fluid pressure",  makeScalar(DIRICH,5.) );
      model.InputBoundaryValue( RIGHT, "fluid pressure",  makeScalar(DIRICH,1.) );

      printRangeOfVariable( model, "fluid pressure" );

     // 3. different values at the nodes of a single element
     _test( !diagnostics.DetectConflictingDirichletConditions( model, "fluid pressure" ) );
   
     // 4. inactive elements (all nodes have constraint status
     _test( !diagnostics.DetectOverConstrainedElements( model, "fluid pressure"  ) );
 
#ifdef CSMP_WITH_SAMG_SOLVER
      SAMG_Settings settings;
      SAMG_Solver   solver( &settings );
      PDE_Integrator<3U,Element> pde_integrator( solver );
      // solver configuration
      settings.SetSolverInstance(1);
      settings.Set_iout1( 0 );
      settings.Set_iout2( 0 );
      if ( !verbose_ ) settings.Set_idmp( -1 );
      /// more expensive relaxation parameter
      settings.Set_nxtyp(1); // GOOD
      /// one-time solver set-up: nothing is remembered for next try
//      settings.Set_iswit(5);
      /// SAMG solution criteria
      settings.Set_eps(0.); // absolute criterion
      /// Agressive first level coarsening nredlev(1) for decreased setup time and reduced no. of cycles
//    fluid_pressure.GetSolverSettings().Set_nred(1);
      /// Pre-adjust SAMG coarse matrix size relative to original size, based on solver output
//      settings.Set_a_cmplx(2);
      /// Pre-adjust SAMG mesh complexity, based on solver output
//      settings.Set_g_cmplx(1.5);
//      settings.Set_w_avrge(2);
#else
      CSMP_DEFAULT_LINEAR_SOLVER  solver;
      PDE_Integrator<3U,Element>  pde_integrator( solver );
#endif

#ifdef SAMG_OUTPUT_TO_FILE
      // trigger output to file
      settings.Set_idmp( 8 );
      // // define SAMG file output format for reduced file size, idmp > 1 is required
      settings.Set_ioform( "f" );
      // set filename for SAMG file output other than default "level", idmp > 1 is required
      settings.Set_filnam_dump( "DirichletPressureBoxModel_VVCase" );
#endif

      const PropertyDatabase<DIM>&   p_ref = model.Database();
      NumIntegral_dNT_op_dN_dV<DIM>  laplacian( p_ref, "conductivity",  "fluid pressure", "fluid pressure" );
      NumIntegral_NT_op_N_dV<DIM>    rhs( p_ref, "fluid volume source", "fluid pressure" );

      pde_integrator.Add( &laplacian );
      pde_integrator.Add( &rhs );           // End the numIntegration
      
      pde_integrator.IntegrateOver( model_domain );
      printRangeOfVariable( model, "fluid pressure" );
      VTK_Interface<DIM>  vtk_out;
      vtk_out.OutputDataToVTK( model, "fluid-pressure", "fluid pressure", 0 );
   
      // 5. smoothness of the computed pressure field
      _test( !diagnostics.DetectNonMonotonicity( model, "fluid pressure" ) );
 }



void DirichletPressureBoxModel_VVCase::TestModelFromANSYS_AnalyticallyIntegrated()
 {
   const bool isoparametric(false);
 
   VSet<DIM>        mesh_container;
   ModelTopology    mesh_topology(isoparametric);
   ANSYS_Interface  mesh_interface(isoparametric);

   const bool binary_file( true );
   mesh_interface.Read_ANSYS_Mesh( model_name_.c_str(), mesh_container, mesh_topology, binary_file, true );

   mesh_topology.ReduceToDomains( model_name_.c_str() );
   map<size_t,size_t>  old_and_new_elmtids;
   mesh_topology.CreateNewCellNumbers( old_and_new_elmtids );
   mesh_container.ReduceTo( old_and_new_elmtids );

   Model<DIM> model( mesh_topology, mesh_container, "DirichletPressureBoxModel_VVCase-variables.txt", true );

   model.InputPropertyValue( "conductivity", makeScalar(PLAIN,1.0e-12) );
   model.InputPropertyValue( "fluid volume source", makeScalar(PLAIN,0.) );
   // !5Pa drop only
   model.InputPropertyValue( "fluid pressure", makeScalar(PLAIN,0.) );
   model.InputBoundaryValue( LEFT, "fluid pressure",  makeScalar(DIRICH,5.) );
   model.InputBoundaryValue( RIGHT, "fluid pressure", makeScalar(DIRICH,1.) );

   printRangeOfVariable( model, "fluid pressure" );

   Region<DIM>& model_domain(model.Region("Model"));
  
#ifdef CSMP_WITH_SAMG_SOLVER
   SAMG_Settings  settings;
   settings.SetSolverInstance(1);
   // iout
   settings.Set_iout1( 0 );
   settings.Set_iout2( 0 );
   if ( !verbose_ ) settings.Set_idmp( -1 );
   // more expensive relaxation parameter
//         settings.Set_nxtyp(1); 
   settings.Set_eps(0.); // absolute criterion

   SAMG_Solver solver(&settings);
#else
   CSMP_DEFAULT_LINEAR_SOLVER solver;
#endif
   PDE_Integrator<DIM,Element>  pde_integrator(solver);
//      pde_integrator.ScaleEssentialConditions( 1.0e-15 ); // because the model is so small

   const PropertyDatabase<DIM>&  p_ref = model.Database();
   LHS_Integral_dNT_dN_dV<DIM>   laplacian( p_ref, "fluid pressure", "fluid pressure" );
//        Integral_dNT_op_dN_dV<DIM>  laplacian( p_ref, "conductivity",  "fluid pressure", "fluid pressure" );
   Integral_NT_op_N_dV<DIM>   rhs( p_ref, "fluid volume source", "fluid pressure" );

   pde_integrator.Add( &laplacian );
   pde_integrator.Add( &rhs );           // End the numIntegration
  
   pde_integrator.IntegrateOver( model_domain );
   printRangeOfVariable( model, "fluid pressure" );

   VTK_Interface<DIM>  vtk_out;
   vtk_out.OutputDataToVTK( model, "fluid-pressure", "fluid pressure", 0 );
   
} // analytically integrated



} // end csmp


