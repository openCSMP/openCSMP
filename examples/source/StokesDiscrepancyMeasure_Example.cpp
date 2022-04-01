#include "StokesDiscrepancyMeasure_Example.h"

#include "Region.h"
#include "PDE_Integrator.h"
#include "ANSYS_Model3D.h"
#include "PropertyHandle.h"
#include "NumIntegral_dNT_dN_dV.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_NT_op_N_dV.h"
#include "VelocityAndVolumeFlux.h"
#include "CSMP_highLevelUtilities.h"

// Interrelations
#include "ExtractVectorVariableLength.h"

// outputting
#include "VTK_Interface.h"

// implicit scheme
#include "InputDataManager.h"

using namespace std;

namespace csmp{

void StokesDiscrepancyMeasure_Example::Specifications()
{
  SetTitle( "Stokes discrepancy measure" );
  SetDifficulty( 2 );
  SetCategory( "Simulation of Physical Processes" );
  AddAuthor( "SKM" );
  AddDescription( "source in: StokesDiscrepancyMeasure_Example.cpp" );
  AddDescription( "calculates error metric" );
  AddRequirement( "one_sphere_0.45_tetra" );
} 

void StokesDiscrepancyMeasure_Example::Run()
{
  // *****************************************************************************************
  /*
       Stokes_discrepancy_measure:  calculates the difference between taking the laplacian
                                    of the Darcy velocity when taking a parabolic scalar
                                    field as an operator and the Stokes flow approximation.
  */
  // *****************************************************************************************
  
  // ESTABLISHING OUTPUTSTREAM FROM BASECLASS
  //ostream &cout = *GetStream();

  const size_t   dim(3U);
  const string   model_name("one_sphere_0.45_tetra");
  ANSYS_Model3D  model( model_name.c_str(), "example25.txt");

  printModelDimensions( model, true );

  InputDataManager<3U>  model_configuration;

  model_configuration.ConfigureFromFile( model, model_name.c_str(),
                                         false, true, true, true, false );

  VTK_Interface<3U>  vtk_output;


 // -----------------------------------------------------------------------
 // 1. Computing the steady-state parabolic function
 // -----------------------------------------------------------------------
    PropertyHandle<dim>  src( model, "source term", SCALAR, ELEMENT ); src = 1.;
    printRangeOfVariable( model, "source term" );
    PropertyHandle<dim>  pbf( model, "parabolic function", SCALAR, NODE );

    #ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings   settings;
    settings.Set_nxtyp(0);
    settings.Set_ncgtyp(5);
    settings.Set_ndefault(40);

    PDE_Integrator<dim,Region>  parabolic_profile(new SAMG_Solver(&settings));
    #else
    PDE_Integrator<dim,Region>  parabolic_profile(new CSMP_DEFAULT_LINEAR_SOLVER());
    #endif
    // the maximum computed 'parabolic function' value is equivalent to the pore radius of the corresponding pore space segment
    // laplacian matrix [L] on the left-hand side
    NumIntegral_dNT_dN_dV<dim,Element<dim> >   laplacian( model.Database(), "parabolic function", "parabolic function" );
    // source vector {q} on the right-hand side = 1
    NumIntegral_NT_op_N_dV<dim,Element<dim> >  rhs( model.Database(), "source term", "parabolic function" );
    parabolic_profile.Add( &laplacian );
    parabolic_profile.Add( &rhs );

    // boundary conditions: we set value of parabolic function at walls of pore space to zero
    // and fix this value by treating it as a Dirichlet boundary condition
    const string  channel_region("Model");
    Region<dim>&  gref(model.Region(channel_region.c_str()));
    gref.InputPropertyValue( "parabolic function", makeScalar(DIRICH,0.), PERIMETER );
    // we must now drop this Dirichlet boundary condition at inlet and outlet (i.e. FRONT and BACK of box-shaped model)
    // this is tricky because we do not want to drop them at the boundary of the boundary, so we can't use
    // the standard functions, but we must loop over the boundary nodes explicitly
    for ( auto nit=gref.PerimeterNodesBegin(); nit!=gref.NodesEnd(); nit++ )
      if ( (*nit)->AtBoundary() == FRONT  or  (*nit)->AtBoundary() == BACK )
        (*nit)->Status( pbf.Key(), PLAIN );
    // calculation
    parabolic_profile.IntegrateOver( gref );
    printRangeOfVariable( model, channel_region.c_str(), "parabolic function" );

    vtk_output.OutputDataToVTK( model, "parabolic-function", "parabolic function", 1 );


 // -----------------------------------------------------------------------
 // 2. Computing the gradient of the parabolic function
 // -----------------------------------------------------------------------
    // the 'parabolic function' divided by viscosity is now mapped to the variable "element parabolic function"
    PropertyHandle<dim>  lapbc( model, "element parabolic function", SCALAR, ELEMENT );
    ScalarVariable  sc;
    for ( auto eit=gref.CellsBegin(); eit!=gref.CellsEnd(); eit++ ) {
         (*eit)->PropertyValueAtBaryCenter( pbf.Key(), sc );
         (*eit)->Store( lapbc.Key(), (sc /= 1.) );
      }
    printRangeOfVariable( model, channel_region.c_str(), "element parabolic function" );

    PropertyHandle<dim>  grad( model, "grad parabolic function", VECTOR, ELEMENT );
    model.CopyGradientOfProperty_A_To_B( "parabolic function", "grad parabolic function" );
    printRangeOfVariable( model, channel_region.c_str(), "grad parabolic function" );

    PropertyHandle<dim>  magn( model, "magnitude grad parabolic function", SCALAR, ELEMENT );
    ExtractVectorVariableLength<dim>  grad_magnitude( model.Database(), "grad parabolic function", "magnitude grad parabolic function" );
    gref.Apply( grad_magnitude );
    printRangeOfVariable( model, channel_region.c_str(), "magnitude grad parabolic function" );


 // -----------------------------------------------------------------------
 // 3. Computing the pressure in the pore space
 // -----------------------------------------------------------------------
    #ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Solver  samg_solver;
    PDE_Integrator<dim,Region>  steady_state_pressure(samg_solver);
    #else
    CSMP_DEFAULT_LINEAR_SOLVER  linear_solver;
    PDE_Integrator<dim,Region>  steady_state_pressure(linear_solver);
    #endif

    NumIntegral_dNT_op_dN_dV<dim,Element<dim> >  flow_resistance( model.Database(), "element parabolic function", "fluid pressure", "fluid pressure" );
    // source vector {q} on the right-hand side = 0
    NumIntegral_NT_op_N_dV<dim,Element<dim> >    fsrc( model.Database(), "fluid volume source", "fluid pressure" );
    steady_state_pressure.Add( &flow_resistance );
    steady_state_pressure.Add( &fsrc );

    model.InputPropertyValue( "fluid volume source", makeScalar(PLAIN,0.) );
    steady_state_pressure.IntegrateOver( gref );
    printRangeOfVariable( model, channel_region.c_str(), "fluid pressure" );
    vtk_output.OutputDataToVTK( model, "fluid-pressure", "fluid pressure", 1 );


 // -----------------------------------------------------------------------
 // 4. Computing the hessian for the FEM solution of fluid pressure
 // -----------------------------------------------------------------------
  model.CopyGradientOfProperty_A_To_B( "fluid pressure", "fluid pressure gradient" );
  model.ExtrapolateCellToNodeProperty("fluid pressure gradient","nodal fluid pressure gradient");
  model.CopyGradientOfProperty_A_To_B( "nodal fluid pressure gradient", "fluid pressure gradient2" );
  model.ExtrapolateCellToNodeProperty("fluid pressure gradient2","nodal fluid pressure gradient2");

  vtk_output.OutputDataToVTK( model, "fluid-pressure-gradient", "fluid pressure gradient", 1 );
  // the Hessian matrix
  vtk_output.OutputDataToVTK( model, "hessian", "fluid pressure gradient2", 1 );


  // -----------------------------------------------------------------------
  // 5. Computing the discrepancy between method and Stokes solution
  //    (as in Lateef's thesis)
  // -----------------------------------------------------------------------
  /*
       E     = (2*(grad psi div nabla) * grad p + psi * grad Laplacian p)   normalized by gradient of p
                term1                           term2

     NB: The result is a vector placed on the element but when the absolute value is taken
         we have a scalar (functional) = to vector length.
  */
  // computing term1: by post-multiplying Hessian with gradient of psi function
   PropertyHandle<dim>  hessian( model, "fluid pressure gradient2", TENSOR, ELEMENT );
   PropertyHandle<dim>  term1( model, "term1", VECTOR, ELEMENT );
   TensorVariable<dim>  ts;
   VectorVariable<dim>  vc, vc2;
   for ( auto eit=gref.CellsBegin(); eit!=gref.CellsEnd(); eit++ ) {
        // read hessian
        (*eit)->Read( hessian.Key(), ts );
        (*eit)->Read( grad.Key(), vc );
        vc2 = ts * vc*2.0;
        //vc2 *= 2.
        (*eit)->Store( term1.Key(), vc2 );
     }

  // computing term2: Laplacian of fluid pressure is found as sum of the diagonal terms of the Hessian (=its trace)
    PropertyHandle<dim>  lap_p( model, "laplacian of pressure", SCALAR, ELEMENT );
    for ( auto nit=gref.CellsBegin(); nit!=gref.CellsEnd(); nit++ ) {
         // read Hessian from the nodes
         (*nit)->Read( hessian.Key(), ts );
         // add up the diagonal terms
         sc = ts(0,0) + ts(1,1) + ts(2,2);
         (*nit)->Store( lap_p.Key(), sc );
      }
    printRangeOfVariable( model, channel_region.c_str(), "laplacian of pressure" );
    vtk_output.OutputDataToVTK( model, "laplacian-pressure", "laplacian of pressure", 1 );

    PropertyHandle<dim>  nlap_pf( model, "nodal laplacian of pressure", SCALAR, NODE );
    model.ExtrapolateCellToNodeProperty("laplacian of pressure", "nodal laplacian of pressure");

    PropertyHandle<dim>  term2( model, "gradient of laplacian of pressure", VECTOR, ELEMENT );
    model.CopyGradientOfProperty_A_To_B( "nodal laplacian of pressure", "gradient of laplacian of pressure" );
    printRangeOfVariable( model, channel_region.c_str(), "gradient of laplacian of pressure");
    vtk_output.OutputDataToVTK( model, "grad-laplacian-pressure", "gradient of laplacian of pressure", 1 );

    for ( auto nit=gref.CellsBegin(); nit!=gref.CellsEnd(); nit++ ) {
         // read parabolic function psi
         (*nit)->Read( lapbc.Key(), sc );
         // read gradient of laplacian
         (*nit)->Read( term2.Key(), vc );
         // add up the diagonal terms
         vc*= sc;

         (*nit)->Store( term2.Key(), vc );
      }

    // forming the numerator of E
    term1 += term2;

    // obtain the final result that is collected into term1
    PropertyHandle<dim>  gradp( model, "fluid pressure gradient", VECTOR, ELEMENT );
    PropertyHandle<dim>  E_vec( model, "Stokes Discrepancy vector", VECTOR, ELEMENT );
    PropertyHandle<dim>  E_mag( model, "Stokes Discrepancy Measure", SCALAR, ELEMENT );

    ScalarVariable  em;
    for ( auto nit=gref.CellsBegin(); nit!=gref.CellsEnd(); nit++ ) {
         // take norm of grad p
         (*nit)->Read( gradp.Key(), vc );
         sc = vc.Length();
         // take norm of term 1
         (*nit)->Read( term1.Key(), vc );
         (*nit)->Store( E_vec.Key(), vc );
         em = vc.Length() / sc();
         (*nit)->Store( E_mag.Key(), em );
      }

    printRangeOfVariable( model, channel_region.c_str(), "Stokes Discrepancy Measure" );
    vtk_output.OutputDataToVTK( model, "Stokes-discrepancy-vector", "Stokes Discrepancy vector", 1 );
    vtk_output.OutputDataToVTK( model, "Stokes-discrepancy-measure", "Stokes Discrepancy Measure", 1 );

} // Run()

} // csmp
