#include "ErrorMetric_Example.h"

#include "Region.h"
#include "PDE_Integrator.h"
#include "ANSYS_Model3D.h"
#include "SteadyStateDiffusor.h"
#include "VelocityAndVolumeFlux.h"
#include "CSMP_highLevelUtilities.h"

// Interrelations
#include "ConstantFactor.h"

// outputting
#include "VTK_Interface.h"

// implicit scheme
#include "InputDataManager.h"

using namespace std;

namespace csmp{

void ErrorMetric_Example::Specifications()
{
  SetTitle( "Error Metric" );
  SetDifficulty( 3 );
  SetCategory( "Numerical Methods" );
  AddAuthor( "SKM" );
  AddDescription( "fluid pressure computation and evaluation of the discretization error" );
  AddDescription( "source in: ErrorMetric_Example.cpp" );
  AddRequirement( "file set: 'prism_test'");
}


/** *****************************************************************************************

   3D fluid pressure computation and evaluation of the discretization error incurred.

   Use models 'prism_test' or 'fracs4' (.dat, .asc, -regions.txt, -configuration.txt)
   as input file suites.

  **************************************************************************************** */
void ErrorMetric_Example::Run()
{
  // ESTABLISHING OUTPUTSTREAM FROM BASECLASS
  //ostream &cout = *GetStream();

    const string model_name("prism_test");
    const bool irregular_mesh(false);
    const bool binary_file(true);
    const bool use_regions_file(true);
    ANSYS_Model3D  model( model_name.c_str(), "example25.txt", irregular_mesh, binary_file, use_regions_file );

    printModelDimensions( model, true );

    InputDataManager<3U>  model_configuration;

    model_configuration.ConfigureFromFile( model, model_name.c_str(),
                                           false, true, true, true, false );


   // -----------------------------------------------------------------------
   // 1. hydraulic conductivity and other interrelations
   // -----------------------------------------------------------------------
    const double fluid_viscosity(1.0e-03);
    ConstantFactor<3U,divides>  conductivity( model.Database(),
                                             "conductivity", "permeability",
                                              fluid_viscosity );
    model.Apply( conductivity );
    printRangeOfVariable( model, "conductivity" );


   // -----------------------------------------------------------------------
   // 2. steady-state fluid pressure
   // -----------------------------------------------------------------------
    SteadyStateDiffusor<3U,Region> steady_state_pressure( model, "conductivity", "fluid pressure",
                                                                        "fluid volume source" );

    VelocityAndVolumeFlux<3U,Element<3U> >  postpro0( model, "conductivity", "porosity", "fluid pressure" );

    steady_state_pressure.AddPostProcess( &postpro0 );
    steady_state_pressure.ComputeSteadyState( model.Region("Model") );

    printRangeOfVariable( model, "fluid pressure" );
    printRangeOfVariable( model, "velocity" );
    printRangeOfVariable( model, "pore velocity" );

    VTK_Interface<3U>  vtk_output;
    vtk_output.OutputDataToVTK( model, "fluid-pressure", "fluid pressure", 1 );
    vtk_output.OutputDataToVTK( model, "velocity",       "velocity",       1 );
    vtk_output.OutputDataToVTK( model, "volume-flux",    "volume flux",    1 );


   // -----------------------------------------------------------------------
   // 3. Computing the error metric for the FEM solution of fluid pressure
   // -----------------------------------------------------------------------
    model.CopyGradientOfProperty_A_To_B( "fluid pressure", "fluid pressure gradient" );
    model.ExtrapolateCellToNodeProperty("fluid pressure gradient","nodal fluid pressure gradient");
    model.CopyGradientOfProperty_A_To_B( "nodal fluid pressure gradient", "fluid pressure gradient2" );
    model.ExtrapolateCellToNodeProperty("fluid pressure gradient2","nodal fluid pressure gradient2");

    vtk_output.OutputDataToVTK( model, "hessian", "fluid pressure gradient2", 1 );
    vtk_output.OutputDataToVTK( model, "FRAC_VOLUMES", "hessian", "fluid pressure gradient2", 1, true );
    vtk_output.OutputDataToVTK( model, "MATRIX",       "hessian", "fluid pressure gradient2", 1, true );

    assignLargestEigenValueOfTo( model, "nodal fluid pressure gradient2", "maximum curl" );

    printRangeOfVariable( model, "maximum curl" );
    vtk_output.OutputDataToVTK( model, "curl", "maximum curl", 1, true );
    vtk_output.OutputDataToVTK( model, "FRAC_VOLUMES", "curl", "maximum curl", 1, true );

    discretizationError3D( model, "fluid pressure gradient2", "discretization error", "discretization error magnitude" );

    printRangeOfVariable( model, "discretization error magnitude" );
    vtk_output.OutputDataToVTK( model, "discretization-error", "discretization error", 1, true );
    vtk_output.OutputDataToVTK( model, "FRAC_VOLUMES", "discretization-error", "discretization error", 1, true );
    vtk_output.OutputDataToVTK( model, "discretization-error-magnitude", "discretization error magnitude", 1, true );
    vtk_output.OutputDataToVTK( model, "FRAC_VOLUMES", "discretization-error-magnitude", "discretization error magnitude", 1, true );

} // Run()



/// computes eigenvectors and values of tensorvariable and puts largest eigenvalue in scalar output variable
void ErrorMetric_Example::assignLargestEigenValueOfTo( Model<3U>& sg, const char* of_var, const char* to_var )
 {
    csmp::Index  hes_key = sg.Database().StorageKey(of_var);
    assert( hes_key.type  == TENSOR );
    assert( hes_key.place == NODE );
    csmp::Index  eig_key = sg.Database().StorageKey(to_var);
    assert( eig_key.type  == SCALAR );
    assert( eig_key.place == NODE );
    TensorVariable<3U>   ts, evecs;
    VectorVariable<3U>   evals;
    Region<3>&  sgref(sg.Region("Model"));

    for ( vector<Node<3U>*>::const_iterator
          nit=sgref.NodesBegin(); nit!=sgref.NodesEnd(); nit++ ) {
         (*nit)->Read( hes_key, ts );

         // computing eigenvalues and vectors, no normalization is applied
         ts.EigenValues( evals );

         // finding the largest eigenvalue be it negative or positive
         ScalarVariable  emax(PLAIN, fabs(evals[0]));
         emax() = std::max( emax(), fabs(evals[1]) );
         emax() = std::max( emax(), fabs(evals[2]) );

         // storing the result
         (*nit)->Store( eig_key, emax );
      }

 } // end assignLargestEigenValueOfTo


void ErrorMetric_Example::discretizationError3D( Model<3U>& sg, const char* hessian_var, const char* error_var, const char* sc_err_var )
 {
    csmp::Index  hes_key = sg.Database().StorageKey(hessian_var);
    assert( hes_key.type  == TENSOR );
    assert( hes_key.place == ELEMENT );
    csmp::Index  err_key = sg.Database().StorageKey(error_var);
    assert( err_key.type  == TENSOR );
    assert( err_key.place == ELEMENT );
    csmp::Index  srr_key = sg.Database().StorageKey(sc_err_var);
    assert( srr_key.type  == SCALAR );
    assert( srr_key.place == ELEMENT );
    TensorVariable<3U>   ts, evecs;
    VectorVariable<3U>   vc, evals;
    ScalarVariable       emag;
    Region<3>&  sgref(sg.Region("Model"));

    for ( auto eit=sgref.CellsBegin(); eit!=sgref.CellsEnd(); eit++ )
      {
         (*eit)->Read( hes_key, ts );

         // computing eigenvalues and vectors, normalization is applied
         ts.Eigen( evals, evecs, true );

         // computing the element length in the directions of the eigenvectors
        vc(0) = evecs(0,0); vc(1) = evecs(0,1); vc(2) = evecs(0,2);
         const double  d1 = (*eit)->LengthInDirection( vc );
        vc(0) = evecs(1,0); vc(1) = evecs(1,1); vc(2) = evecs(1,2);
         const double  d2 = (*eit)->LengthInDirection( vc );
        vc(0) = evecs(2,0); vc(1) = evecs(2,1); vc(2) = evecs(2,2);
         const double  d3 = (*eit)->LengthInDirection( vc );

         // computing e = vT |H| v
         vc(0) = d1; vc(1) = d2; vc(2) = d3;
         //    |H|    v
         ts = evecs * vc;
         //   vT  (|H| v)
         vc = vc * ts;

         // error metric for visualisation
         for ( auto j=0U; j<3U; j++ )  {
              evecs(0,j) *= d1 * evals(0);
              evecs(1,j) *= d2 * evals(1);
              evecs(2,j) *= d3 * evals(2);
           }

         // getting the largest discretization error
         emag() = std::max( d1 * evals(0), std::max( d2 * evals(1), d3 * evals(2)) );

         // storing the results
         (*eit)->Store( err_key, evecs ); // discretization error visualized
         (*eit)->Store( srr_key, emag ); // error magnitude
      }

 } // end discretizationError3D

} // csmp
