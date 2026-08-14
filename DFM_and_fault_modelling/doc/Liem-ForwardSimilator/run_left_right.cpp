#include <iostream>

// the CSMP model
#include "Model.h"

#include <chrono>

#include <filesystem>

#include "FlowAndTransport.h"
#include "MeasureDataForEnKF.h"

// Transport scheme
#include "DESAdvectionDiffusion.h"

// monitoring individual regions
#include "RegionMonitor.h"


using namespace csmp;
using namespace std;


int main( int argc, char* argv[] )
{
  // Measure time when simulation starts
  auto time_simulation_start = std::chrono::high_resolution_clock::now();

  // Get input arguments
  size_t model_nr, input_iteration,input_realisation;
  std::string tmp_folder = "./";
  bool use_tmp_folder(false);
  if (argc==4) {
    model_nr=          strtol(argv[1],nullptr,10);
    input_iteration=   strtol(argv[2],nullptr,10);
    input_realisation= strtol(argv[3],nullptr,10);
  } else if (argc==5) {
    model_nr=          strtol(argv[1],nullptr,10);
    input_iteration=   strtol(argv[2],nullptr,10);
    input_realisation= strtol(argv[3],nullptr,10);
    tmp_folder= argv[4];
    tmp_folder= tmp_folder + "/";
    use_tmp_folder= true;
  } else {
    model_nr = 7;
    input_iteration= 0;
    input_realisation= 0;
  }

  // Define some variables
  std::string set_string;
  bool suppress_console_output(true);
  bool save_vtk_files(false);

  double hours = 3600.;
  double days = 24.*hours;

  std::string model_name = "odling720x720_separateSet_manual_2024";
  double time_end = 40.*days;
  double time_save = .5*days;
  double time_step = .5*days;


  cout << "Run run_left_right.cpp with:" << endl;
  cout << "\tModel " << to_string(model_nr) << " : " << model_name << endl;
  cout << "\ttime_end= " << to_string(time_end/hours) << " h" << endl;
  cout << "\ttime_save= " << to_string(time_save/hours) << " h" << endl;
  cout << "\ttime_step= " << to_string(time_step/hours) << " h" << endl;
  cout << "\tinput_iteration= " << to_string(input_iteration) << endl;
  cout << "\tinput_realisation= " << to_string(input_realisation) << endl;
  cout << "\ttmp_folder= " << tmp_folder << endl;


  // Load the csmp model
  std::string model_path_and_name = "model_csmp/" + model_name + "/" + model_name;
  cout << "\nLoad csmp model from " << model_path_and_name << endl;
  if (suppress_console_output ){ std::cout.setstate(std::ios_base::failbit); }
  Model<2U> model( model_path_and_name );
  const PropertyDatabase<2>& p_ref = model.Database();
  if (suppress_console_output ){ std::cout.clear(); }


  // Load fracture apertures
  std::string realisation_name = to_string(input_iteration) + "_" + to_string(input_realisation);
  std::string input_path_and_name = "input_data/" + model_name + "/" + realisation_name + ".txt";

  FlowAndTransport::LoadInputData( model, input_path_and_name );



  // Output folder name
  cout << "\nOutput folders" << endl;
  std::string write_measurement_folder, write_output_data_folder;
  std::string target_measurement_folder= "left_right_measurement/" + model_name + "/";
  std::string target_output_data_folder= "left_right_data/"        + model_name + "/";
  if (use_tmp_folder){
    write_measurement_folder= tmp_folder + realisation_name + "/";
    write_output_data_folder= tmp_folder + realisation_name + "/";
  } else {
    write_measurement_folder=  target_measurement_folder + realisation_name + "/";
    write_output_data_folder=  target_output_data_folder + realisation_name + "/";
  }
  cout << "\twrite_measurement_folder:  " << write_measurement_folder << endl;
  cout << "\twrite_output_data_folder:  " << write_output_data_folder << endl;
  cout << "\ttarget_measurement_folder: " << target_measurement_folder << endl;
  cout << "\ttarget_output_data_folder: " << target_output_data_folder << endl;


  // Create output folders
  if (filesystem::exists( write_measurement_folder ) ){
    cout << "\tDelete existing folder " << write_measurement_folder << " and its content" << endl;
    filesystem::remove_all( write_measurement_folder );
  }
  cout << "\tCreate folder " << write_measurement_folder << endl;
  filesystem::create_directories( write_measurement_folder );


  if (save_vtk_files){
    cout << "\nCreate output folders" << endl;
    if (filesystem::exists( write_output_data_folder ) ){
      cout << "\tDelete existing folder " << write_output_data_folder << " and its content" << endl;
      filesystem::remove_all( write_output_data_folder );
    }
    cout << "\tCreate folder " << write_output_data_folder << endl;
    filesystem::create_directories( write_output_data_folder );
  }


  // Change boundary condition of the model to left-right case
  model.Boundary("LEFT").InputPropertyValue("fluid pressure",makeScalar(DIRICH,2.0e7));
  model.Boundary("RIGHT").InputPropertyValue("fluid pressure",makeScalar(DIRICH,1.0e7));
  model.Boundary("BOTTOM").ChangePropertyStatus("fluid pressure",NEUMANN);
  model.Boundary("TOP").ChangePropertyStatus("fluid pressure",NEUMANN);

  model.Boundary("LEFT").InputPropertyValue("concentration",makeScalar(DIRICH,1.0));

  model.Region("WELL_FRACTURE").InputPropertyValue("fluid volume source",makeScalar(ANY,0.0));
  model.Region("WELL_FRACTURE").InputPropertyValue("concentration",makeScalar(ANY,0.0));

  printRangeOfVariable( model, "fluid volume source");
  printRangeOfVariable( model, "concentration");
  printRangeOfVariable( model, "fluid pressure");





  // Create three linear measurement regions in the interior of the domain
  MeasureDataForEnKF::CreateMeasRegion( model, 4 );



  // Solve for pressure field with SolveForPressure()
  FlowAndTransport::SolveForPressure(model, p_ref);


  // Measure steady-state pressure field and initial concentration field
  MeasureDataForEnKF::MeasurePressure(      model, "MeasRegion",   write_measurement_folder );
  MeasureDataForEnKF::MeasureConcentration( model, "MeasRegion",   write_measurement_folder ,0);



  // Write some output files (if desired)
  if (save_vtk_files){
    MeasureDataForEnKF::WriteSteadyStateVTK(model, write_output_data_folder);
    MeasureDataForEnKF::WriteConcentrationVTK( model, write_output_data_folder, 0);
  }




  // CFL multiplier for DES
  double cfl_multiplier_des = 0.4;

  DESAdvectionDiffusion<2U> DES_transport(model, "Model", cfl_multiplier_des, 1.0, false);

  // Solve for concentration
  for (size_t i_save = 0; double(i_save) < time_end/time_save; ++i_save) {
    for (size_t i_step = 0; double(i_step) < time_save/time_step; ++i_step) {
      cout << "\ni_save = " << i_save << " , i_step = " << i_step << endl;

      // Current (new) time
      double new_time = double(i_save)*time_save + double(i_step+1)*time_step;

      // Advect concentration to current (new) time
      if (suppress_console_output ){ std::cout.setstate(std::ios_base::failbit); }
      DES_transport.AdvectVariable_DES( new_time, 1);
      if (suppress_console_output ){ std::cout.clear(); }
    }

    // Measure current concentration field
    MeasureDataForEnKF::MeasureConcentration( model,"MeasRegion"   , write_measurement_folder,i_save+1 );

    // Write current concentration to file
    if (save_vtk_files) {
      MeasureDataForEnKF::WriteConcentrationVTK(model, write_output_data_folder, i_save + 1);
    }

    printRangeOfVariable(model, "concentration");
  }


  if (use_tmp_folder){
    cout << "\nCopy files from tmp folder to actual folder " << endl;

    // Create root folder (if not already there)
    filesystem::create_directories( target_measurement_folder );

    // Copy
    cout << "\tSource: " << write_measurement_folder << endl;
    cout << "\tTarget: " << target_measurement_folder + realisation_name << endl;
    const auto copyOptions = filesystem::copy_options::update_existing | filesystem::copy_options::recursive;
    filesystem::copy(write_measurement_folder, target_measurement_folder + realisation_name, copyOptions);

    if (save_vtk_files){
      // Create root folder (if not already there)
      filesystem::create_directories( target_output_data_folder );

      // Copy
      cout << "\tSource: " << write_output_data_folder << endl;
      cout << "\tTarget: " << target_output_data_folder + realisation_name << endl;
      filesystem::copy(write_output_data_folder, target_output_data_folder + realisation_name, copyOptions);
    }
  }


  // Measure simulation runtime
  auto time_simulation_end = std::chrono::high_resolution_clock::now();
  auto runtime = std::chrono::duration_cast<std::chrono::seconds>(time_simulation_end-time_simulation_start).count();
  cout << "\n\nModel run in " << runtime << "s = ~" << runtime/60 << "min = ~" << runtime/3600 << "h" << endl;
  return 0;

} // end main