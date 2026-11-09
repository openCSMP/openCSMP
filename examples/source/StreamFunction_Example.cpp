#include "StreamFunction_Example.h"

#include "Model.h"
#include "Boundary.h"
#include "Region.h"
#include "PDE_Integrator.h"
#include "CSMP_highLevelUtilities.h"
#include "variableOperations.h"

// Mesh Import & Property Assignment
#include "TRIANGLE_Interface.h"
#include "VSet.h"
#include "VSetConverter.h"
//#include "BinaryFileInterface.h"
#include "InputDataManager.h"

// PDE Operators
#include "NumIntegral_NT_rhsop_N_dV.h"
#include "NumIntegral_NT_lhsop_N_dV.h"
#include "NumIntegral_dNT_lhsop_dN_dV.h"
#include "VelocityAndVolumeFlux.h"
#include "NumIntegral_op_NT_dN_orthogonal_dV.h"
#include "LinearSolver.h"

// Interrelations
#include "ConstantFactor.h"

// Other
#include "PropertyConstraints.h"
#include "VTK_Interface.h"
#include "Standard_IO_Handler.h"
#include "PropertyHandle.h"
#include "ModelTime.h"

// Analysis
#include "StatisticalAnalyzer.h"
#include "RegionMonitor.h"


using namespace std;

namespace csmp{

void StreamFunction_Example::Specifications()
{
  SetTitle( "Streamfunction from velocity field. Function can be contoured to get streamlines" );
  SetDifficulty( 3 );
  SetCategory( "Simulation of Physical Processes" );
  AddAuthor( "SKM" );
  AddDescription( "post-processing of the velocity field using the streamfunction that can be contoured" );
  AddDescription( "steady-state analysis of flow and fluid pressure in heterogeneous medium" );
  AddDescription( "source in StreamFunction_Example.cpp" );
  AddRequirement( "file set: 'example21.1', example21.1-configuration.txt");
  AddRequirement( "variable file (example21.txt), 'example21_histogram.bin'");
}



/** *****************************************************************************************

  Analysis of steady-state fluid pressure and velocity distributions in fractured medium
  using quadratic FEM.
 
  Post-processing of the streamfunction that can be contoured to get streamlines.
 
  Example also allows to vary the permeability of model subregions (addressed by 
  their name) and calculates its effective permeability.

  Try 'example21.1'  input file set.

***************************************************************************************** */
void StreamFunction_Example::Run()
{
   /*
   TRIANGLE_Interface  mesh_interface;
   VSet<2U>            mesh_container;

   char   file_name[200];
   cout <<"\nmain: Enter name of 'Triangle' input file set: ";
   cin >> file_name;
   string model_name(file_name);
   model_name += "simulation.txt";

   Standard_IO_Handler  stdio( model_name.c_str() );

   mesh_interface.ReadTriangle2DMesh( file_name, mesh_container );

   VSetConverter<2U>  converter;
   converter.ConvertLinearToQuadraticTriangles( mesh_container );

 // ------------------------------------------------------------------------------------
 // 0. building Region object
 // ------------------------------------------------------------------------------------
  // to get isoparametric elements
   VSetConverter<2U>().ConvertElementTypesToOnesUsingLocalCoordinateSystem( mesh_container );
   Model<2U>  model( mesh_container, "example21.txt" );
   mesh_container.Erase();
   */

  // ------------------------------------------------------------
  // 0. Load CSMP native format model
  // ------------------------------------------------------------
  string model_name;
  cout<< "\nPlease enter the name of input model, or press ENTER to use the default model 'example21.1':"<<endl;
  cin.ignore();
  getline(cin, model_name);
  if (model_name.length() == 0) model_name = "example21.1";

  //find the name of current example source file
  string file_name = GetExampleFileName(__FILE__);
  string variable_file = "example21.txt";
  string config_file = model_name;
  //create of directory with current example name, go into this directory, and copy input files into it.
  CreateWorkingDirectoryAndCopyInputModelFiles(file_name, model_name, variable_file, config_file);
  //reads model from CSMP's native binary files, but creating (additional) storage based on supplied variable file
  Model<2U>  model(model_name, variable_file);
  model.EstablishBoxBoundariesFromOrientation();
  model.RegionsOut();
  model.BoundariesOut();
  printModelDimensions( model );

  //copy in 'example21_histogram.bins' file to be used by StatisticalAnalyzer
  string path = "../../example_inputs/variables_and_configuration_files/";
  string name = "example21_histogram.bins";
  file_name = path + name;
  if (filesystem::exists(file_name)) filesystem::copy(file_name, "./");
  else {
    string error_message = "\n\nError: file '";
    string input_directory = (filesystem::current_path().parent_path().parent_path()).string();
    input_directory += "/example_inputs/variables_and_configuration_files/";
    error_message += (name + "' does not exist in directory "  + input_directory);
    error_message += (", example cannot run, please copy this file into this directory\n");
    throw std::runtime_error(error_message);
  }

  Standard_IO_Handler  stdio( model_name.c_str() );


 // ------------------------------------------------------------------------------------
 // 1. Material properties, groups etc., initial & essential conditions &
 //    interrelations
 // ------------------------------------------------------------------------------------
   InputDataManager<2U>  model_configuration;
   model_configuration.ConfigureFromFile( model, config_file.c_str() );

   printModelDimensions( model );


 // ------------------------------------------------------------------------------------
 // 2. Calculating hydraulic conductivity from permeability using Interrelation subclass
 // ------------------------------------------------------------------------------------
   constexpr double fluid_viscosity(1.0e-03);
   PropertyHandle<2U> k( model, "permeability", SCALAR, ELEMENT );
   PropertyHandle<2U> K( model, "conductivity", SCALAR, ELEMENT );
   K  = 1. / fluid_viscosity;
   K *= k;
   // computing 'resistivity' as the inverse of conductivy
   PropertyHandle<2U> R( model, "resistivity", SCALAR, ELEMENT );
   R  = 1.;
   R /= K;
   printRangeOfVariable( model, "conductivity" );
   printRangeOfVariable( model, "resistivity" );


 // ------------------------------------------------------------------------------------
 // 3. Steady-state fluid pressure computation [K]{p} = {Q}
 // ------------------------------------------------------------------------------------
   EigenSolver                 verification_solver;
   PDE_Integrator<2U,Element>  fluid_pressure(verification_solver);

   NumIntegral_dNT_lhsop_dN_dV<2U>  conductance( model.Database(), "conductivity",   "fluid pressure", "fluid pressure" );
   NumIntegral_NT_rhsop_N_dV<2U>    source( model.Database(), "fluid volume source", "fluid pressure" );
   VelocityAndVolumeFlux<2U>        velo( model, "conductivity", "porosity", "fluid pressure" ); // element-based

   fluid_pressure.Add( &conductance );
   fluid_pressure.Add( &source );
   fluid_pressure.AddPostProcess( &velo );
   
   model.Apply( fluid_pressure );


 // ------------------------------------------------------------------------------------
 // 3.1 Analysis of results
 // ------------------------------------------------------------------------------------
   StatisticalAnalyzer<2U>                                 flux_histogram( model );
   vector<pair<double,double> >                            bins;
   map<string,pair<vector<pair<double,double> >,size_t> >  results;

   flux_histogram.DefineBins( "example21_histogram.bins", bins );
   flux_histogram.RegionPropertyHistograms( "volume flux", bins, results );
   flux_histogram.OutputRegionPropertyAbundancePolygonsMaple( "volume flux", bins, results, true ); // log10 of data

   printRangeOfVariable( model, stdio, "fluid pressure", true );
   printRangeOfVariable( model, stdio, "velocity", true ); // barycentre
   printRangeOfVariable( model, stdio, "volume flux", true );

   list<string> integral_properties;  integral_properties.push_back("volume flux");
   list<string> range_properties;     range_properties.push_back("fluid pressure");

   RegionMonitor<2U>  monitor( model, integral_properties, range_properties );
   bool normalize_integrated_values = false;
 //	    stdio.RecordLogicalChoice("Would you like to normalize region integral by region volume");
   double& model_time( ModelTime::Instance().modelTime );  model_time = 0.;
   monitor.DivideIntegralPropertiesByRegionVolumes( normalize_integrated_values );
   monitor.ScalarPropertyIntegrals( model, model_time );
   monitor.ScalarPropertyRanges( model, model_time );
   monitor.Out( "CSMP-monitored-regions" );


 // ------------------------------------------------------------------------------------
 // 4. effective permeability
 // ------------------------------------------------------------------------------------
   cout <<"\nmain: Enter direction of flow (parallel-to-X-axis=1, Y-axis=2): ";
   size_t  flow_direction;
   cin >> flow_direction;
   Point<2U>  xyz_min, xyz_max;
   model.Region("Model").MinMaxCoordinates( xyz_min, xyz_max );
   const double  xsec_X = xyz_max[0] - xyz_min[0];
   const double  xsec_Y = xyz_max[1] - xyz_min[1];
   double  pmin, pmax, farfield_pf_gradient,
             // non signalling initialization to NAN value
             k_effective(std::numeric_limits<double>::quiet_NaN());

   double  model_throughput = integrateDomainBoundaryFlux( model );

   model.MinMaxOf( "fluid pressure", pmin, pmax );

   if ( flow_direction == 1U ) {
        farfield_pf_gradient = (pmax - pmin) / xsec_X;
        // from Darcy's law
        //        effective permeability = (q * mu) / grad_P
        k_effective  = model_throughput / farfield_pf_gradient;
        k_effective /= xsec_Y;
     }
   else if ( flow_direction == 2U ) {
        farfield_pf_gradient = (pmax - pmin) / xsec_Y;
        // from Darcy's law
        //        effective permeability = (q * mu) / grad_P
        k_effective  = model_throughput / farfield_pf_gradient;
        k_effective /= xsec_X;
     }
   k_effective *= fluid_viscosity;
   cout <<"\nmain: The (flux-integrated) effective permeability of the model is: "<< k_effective << endl;
   cout <<"\nmain: The permeability range of the model was:"<< endl;
   printRangeOfVariable( model, stdio, "permeability", true );

   // output results for visualization using VTK
   VTK_Interface<2U>  vtk_output;
   vtk_output.OutputDataToVTK( model, "conductivity", "conductivity", 0 );
   vtk_output.OutputDataToVTK( model, "fluid-pressure", "fluid pressure", 0 );
   vtk_output.OutputDataToVTK( model, "velocity",       "velocity",       0 );
   vtk_output.OutputDataToVTK( model, "volume-flux",    "volume flux", 0 );


 // ------------------------------------------------------------------------------------
 // 4. Compute the streamfunction which corresponds to solution
 // ------------------------------------------------------------------------------------
   if ( flow_direction == 1 ) computeStreamFunction( model, BOTTOM, TOP, model_throughput, "stream function" );
   else                       computeStreamFunction( model, LEFT, RIGHT, model_throughput, "stream function" );
   double max_stream = printRangeOfVariable( model, stdio, "stream function", true );
   cout <<"\nmain: Total throughput versus max indicated by streamfunction solution: ";
   cout << model_throughput <<" vs. "<< max_stream <<" (should be the same)."<< endl;

   // dependent variable ranges in groups
   map<string,Region<2> >::const_iterator  grit(model.RegionsBegin());
   while ( grit!=model.RegionsEnd() ) {
        printRangeOfVariable( model, stdio, (*grit).first.c_str(), "fluid pressure" );
        printRangeOfVariable( model, stdio, (*grit).first.c_str(), "velocity" );
        printRangeOfVariable( model, stdio, (*grit).first.c_str(), "volume flux" );
        grit++;
     }

   // output results for visualization using VTK
   vtk_output.OutputDataToVTK( model, "stream-function", "stream function",  0 );

 //  BinaryFileInterface<double,2U>  bin_output;
 //  bin_output.WriteConnectivityFile( model, "model_connectivity" );
 //  bin_output.WriteDataTo( model, "fluid-pressure", "fluid pressure", 0 );


 // ------------------------------------------------------------------------------------
 // 4. Sensitivity test for a single region
 // ------------------------------------------------------------------------------------
   string region_name;

   if ( stdio.RecordLogicalChoice("Would you like to vary permeability in a single region") ) {
        cout <<"\nmain: Enter name of the desired region: ";
        cin >> region_name;

        printRangeOfVariable( model, stdio, region_name.c_str(), "fluid pressure" );
        printRangeOfVariable( model, stdio, region_name.c_str(), "velocity" );
        printRangeOfVariable( model, stdio, region_name.c_str(), "volume flux" );

        ConstantFactor<2U,divides> conductivity( model.Database(), "conductivity", "permeability", fluid_viscosity );

        AnalyseSensitivity( model, region_name.c_str(), stdio, conductivity, fluid_pressure );

        if ( stdio.RecordLogicalChoice("Output last result from sensitivity analysis ?") ) {
             vtk_output.OutputDataToVTK( model, "fluid-pressure", "fluid pressure",    0 );
             vtk_output.OutputDataToVTK( model, "velocity",       "velocity",          0 );
             vtk_output.OutputDataToVTK( model, "volume-flux",    "volume flux",       0 );
          }
     }

   cout <<"\nmain: That's it..."<< endl;

   filesystem::current_path("../../example_inputs/");


} // Run()






// auxiliary methods

// this method assumes that permeability is a scalar
void  StreamFunction_Example::AnalyseSensitivity( Model<2U>& model, const char* region,
                                                  Standard_IO_Handler& io,
                                                  Interrelation<2U>& itr,
                                                  PDE_Integrator<2U,Element>& algo )
 {
    assert( model.Database().Type("permeability") == SCALAR );
    for ( ;; ) {
          double average_k = model.Region(region).Average( "permeability" );
          cout <<"\nanalyze_sensitivity: Current average permeability, k = "<< average_k;
          cout <<"; enter new k value (-1. to break loop): ";
          ScalarVariable  perm;
          cin >> perm();
          if ( perm() < 0. ) break;
          model.Region(region).InputPropertyValue( "permeability", perm, COMPLETE );
          model.Apply( itr );
          model.Apply( algo );
          printRangeOfVariable( model, io, region, "fluid pressure" );
          printRangeOfVariable( model, io, region, "velocity" );
          printRangeOfVariable( model, io, region, "volume flux" );
      }
 }

/**

bool integrateDomainBoundaryFlux( Region<3>& sg, double& inflow, double& outflow )

@section description Description

For all boundary Faces of the current mesh,
from the fluid pressures at the nodes, velocities were found at the
integration points. These velocities are extrapolated to the nodes
points of each face and then integrated
by the faces over the faces.

Method returns total throughput through the model.

@return The results of the computation are returned into second and third
method arguments. They represent the total inflow and outflow of the
model, respectively.

@section application Application

This method is used to calculate the effective hydraulic conductivity
of the model domain.

tested: */
template<uint32_t dim>
double StreamFunction_Example::integrateDomainBoundaryFlux( Model<dim>& sg )
 {
  // ESTABLISHING OUTPUTSTREAM FROM BASECLASS
  //ostream &cout = *GetStream();

  csmp::Index ve_key=sg.Database().StorageKey("velocity");

  if ( ve_key.place != ELEMENT || ve_key.type != VECTOR ) {
       throw csmp::Exception( ERROR, "integrateDomainBoundaryFlux",
                                     "velocity variable must be a vector placed on the element");
       return false;
    }

  // integrating fluxes over the faces and adding ensuing contributions to flux balances
  // -----------------------------------------------------------------------------------
  double  velocity, flux(0.);
  double  sum_influx(0.), sum_outflux(0.);

  for ( auto bit=sg.BoundariesBegin(); bit!=sg.BoundariesEnd(); bit++ )
    for ( auto fit=(*bit).second.CellsBegin(); fit!=(*bit).second.CellsEnd(); fit++ )
      {
         // verify that the Face is indeed on the outside of the model
         assert( (*fit)->Parent(INSIDE)  != nullptr );
         assert( (*fit)->Parent(OUTSIDE) == nullptr );
              
         // computing and storing normal to face
         VectorVariable<2U> normal;
         (*fit)->UnitNormal( normal );

         // finding the neighbor element, getting the velocity from it, and projecting it onto the facet normal
         auto velo = (*fit)->Parent(INSIDE)->ReadVector( ve_key );
         // projecting flux onto facet normal
         velocity = dotProduct(normal,velo);
         flux    += (*fit)->Area() * velocity;

         // recording the integrated fluxes
         // -------------------------------------------------------------------------------------------
         // an outflux is is aligned with the normal and therefore positive
         // else this is an influx
         if      ( flux > 0. ) sum_outflux += flux;
         else if ( flux < 0. ) sum_influx  += fabs(flux);

      } // end for if faces

  // total throughput must be 1/2 of fabs() of what goes in and what goes out
  cout <<"\nintegrateDomainBoundaryFlux: Total cross-sectional throughput: "<< (sum_influx+sum_outflux)/2. << endl;

   return (sum_influx + sum_outflux) / 2.;

} // end integrateDomainBoundaryFlux





// @test tested: O.K.
void StreamFunction_Example::computeStreamFunction( Model<2U>& sg,
                                                    BOX_BOUNDARY boundary0, BOX_BOUNDARY boundary1,
                                                    double total_flux, const char* stream_func_var )
 {
  // ESTABLISHING OUTPUTSTREAM FROM BASECLASS
  //ostream &cout = *GetStream();

    // assigning boundary conditions at one boundary which should be parallel to flow direction
    // from this boundary away the flux will be integrated by this method
    sg.InputBoundaryValue( boundary0, stream_func_var, makeScalar(DIRICH, 0.) );
    sg.InputBoundaryValue( boundary1, stream_func_var, makeScalar(DIRICH, total_flux) );

    EigenSolver                 verification_solver;
    PDE_Integrator<2U,Element>  stream_function(verification_solver);

    NumIntegral_dNT_lhsop_dN_dV<2U>  conductance( sg.Database(), "resistivity", stream_func_var, stream_func_var );

// RENAME THIS OPERATOR into NumIntegral_NT_op_dN_orthogonal_dV
    NumIntegral_op_NT_dN_orthogonal_dV<2U>  rhs( sg.Database(), "fluid pressure", stream_func_var );

    stream_function.Add( &conductance );
    stream_function.Add( &rhs );

    cout <<"\n\ncomputeStreamFunction: Computing stream function..."<< endl;
    sg.Apply( stream_function );

 } // end computeStreamFunction

} // csmp
