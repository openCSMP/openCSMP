//
// Created by Michael Liem on 08.06.22.
//

#include "FlowAndTransport.h"

using namespace std;

namespace csmp{

  void FlowAndTransport::LoadInputData( Model<2U>& model, const std::string& input_path_and_name ) {
    cout << "\nReading fracture apertures from input file " << input_path_and_name << " ..." << endl;

    // Check if file is there
    if (!doesFileExist(input_path_and_name)){
      cerr << "\nError: Input file " << input_path_and_name << " not found." << endl;
      throw std::runtime_error("Input file not found.");
    }

    // Open file
    std::ifstream ifs;
    openFile( ifs, input_path_and_name);

    // Make sure the file is open
    if(!ifs.is_open()) throw std::runtime_error("Could not open file");

    // Create helper variables
    std::string line;
    size_t current_fracture_set = 0, current_fracture_number_in_set = 0;
    double current_aperture = 0.0;

    // Create vectors to store results
    std::vector<size_t> all_fracture_set, all_fracture_number_in_set;
    std::vector<double> all_apertures;

    // Read aperture and store values in struct
    while (std::getline(ifs,line)){
      // Create a string stream of the current line
      std::stringstream ss(line);

      // Extract first element of string stream: fracture set
      ss >> current_fracture_set;
      all_fracture_set.push_back(current_fracture_set);

      // If the next token is a comma, ignore it and move on. Else throw and error.
      if(ss.peek() == ',') {
        ss.ignore();
      } else {
        cerr << "Error in enkf_run_model.cpp::LoadInputData(): Something wrong with line " << ss.str() << endl;
        throw std::runtime_error("Could not read input file");
      }

      // Extract second element of string stream: fracture number in set
      ss >> current_fracture_number_in_set;
      all_fracture_number_in_set.push_back(current_fracture_number_in_set);

      // If the next token is a comma, ignore it and move on. Else throw and error.
      if(ss.peek() == ',') {
        ss.ignore();
      } else {
        cerr << "Error in enkf_run_model.cpp::LoadInputData(): Something wrong with line " << ss.str() << endl;
        throw std::runtime_error("Could not read input file");
      }

      // Extract third element of string stream: fracture aperture
      ss >> current_aperture;
      all_apertures.push_back(current_aperture);
    }
    ifs.close();

    // Storage key to access variables
    csmp::Index TFN_key =  model.Database().StorageKey("total fracture number");
    csmp::Index t_key =    model.Database().StorageKey("thickness");
    csmp::Index k_key =    model.Database().StorageKey("permeability");

    auto total_fracture_number = (size_t)model.Read( TFN_key );

    // Print warning if size of vector does not match total number of fractures
    if (all_apertures.size() > total_fracture_number){
      cout << "Warning: " << all_apertures.size()-total_fracture_number << " spare apertures defined" << endl;
      cout << "\tNon-existing regions will be skipped" << endl;
    }
    if (all_apertures.size() < total_fracture_number){
      cout << "Warning: Less apertures defined than fractures in model" << endl;
      cout << "\tUse default fracture aperture for remaining " << total_fracture_number-all_apertures.size() << " fractures" << endl;
    }

    // Loop through input data
    for (size_t idx = 0; idx < all_apertures.size(); ++idx){
      std::string region_str = "FRACTURE_" + to_string(all_fracture_set[idx]) + "_" + to_string(all_fracture_number_in_set[idx]);
      if (model.ContainsRegion(region_str)){
        // Thickness and permeability of fracture
        double t_value = all_apertures[idx];
        double k_value = t_value * t_value / 12;

        if (idx < 10){
          cout << "\t" << region_str << " :  \tt= " << t_value << " , \tk= " << k_value << endl;
        } else if (idx == 10) {
          cout << "\t" << "..." << endl;
        }


        // Loop through all elements of fracture region and store thickness and permeability
        for (auto & eit : model.Region(region_str).CellVector() ){
          (*eit).Store(t_key, makeScalar(PLAIN, t_value));
          (*eit).Store(k_key, makeScalar(PLAIN, k_value));
        }
      } else {
        cout << "\t" << region_str << " : region does not exist" << endl;
      }
    }

    // Manually set matrix permeability (make sure this is commented out if you don't need it!)
    //double k_value= 3e-13;
    //for (auto & eit : model.Region("MATRIX").CellVector() ){
    //  (*eit).Store(k_key, makeScalar(PLAIN, k_value));
    //}
    //for (auto & eit : model.Region("DELETED_ELEMENTS").CellVector() ){
    //  (*eit).Store(k_key, makeScalar(PLAIN, k_value));
    //}

    // Interrelation to compute hydraulic conductivity from permeability and fixed viscosity
    ConstantFactor<2U,divides>  conductivity( model.Database(), "conductivity", "permeability", 0.001 ); // viscosity 1 cp = 0.001 Pa s
    model.Apply( conductivity );

    // Print range of variable
    cout << "\nRange of variables on model" << endl;
    printRangeOfVariable(model, "porosity");
    printRangeOfVariable(model, "thickness");
    printRangeOfVariable(model, "permeability");
    printRangeOfVariable(model, "conductivity");
    printRangeOfVariable(model, "fluid viscosity");
    printRangeOfVariable(model, "fluid pressure");
    printRangeOfVariable(model, "fluid volume source");
  }



  void FlowAndTransport::SolveForPressure( Model<2U>& model, const PropertyDatabase<2>& p_ref ) {

    // create the CSMP FE Algorithm with SAMG solver
//#ifdef CSMP_WITH_SAMG_SOLVER
//    // Adjust SMAG solver settings (reduce console output)
//    SAMG_Settings  samg_settings;
//    samg_settings.Set_iout1(0);
//    samg_settings.Set_iout2(0);
//
//    SAMG_Solver  samg_solver(&samg_settings);
//    PDE_Integrator<2U,Region>  fluid_pressure(samg_solver);
//#else
//    CSMP_DEFAULT_LINEAR_SOLVER  linear_solver;
//  PDE_Integrator<2U,Region>  fluid_pressure(linear_solver);
//#endif

    EigenSolver eigen_solver;
    PDE_Integrator<2U,Element>  fluid_pressure(eigen_solver);

    // LHS stiffness matrix                                 operand         basis function    test function
    NumIntegral_dNT_op_dN_dV<2U>  stiffness_matrix( p_ref, "conductivity", "fluid pressure", "fluid pressure" );

    // RHS mass matrix for integrating source term
    NumIntegral_NT_op_N_dV<2U>    source_term( p_ref, "fluid volume source", "fluid pressure" );

    // use lumped formulation for all mass matrices (i.e., diagonalise matrices)
    source_term.LumpedFormulation(true);

    // define a post-processing step that computes the velocity in each finite element by solving Darcy's law
    VelocityAndVolumeFlux<2U>     velo( model, "conductivity", "porosity", "fluid pressure", true );

    // now add each FE operation (i.e., PDE Operator) to the FE algorithm
    fluid_pressure.Add( &stiffness_matrix );
    fluid_pressure.Add( &source_term );
    fluid_pressure.AddPostProcess( &velo );

    // compute fluid pressure field
    model.Apply( fluid_pressure );

    // reset algorithm and return memory because it is no longer needed (important for large models)
    fluid_pressure.Reset();

    // show results
    printRangeOfVariable( model, "fluid pressure");
    printRangeOfVariable( model, "velocity");
    printRangeOfVariable( model, "pore velocity");
    printRangeOfVariable( model, "fluid volume source");

  } // end SolveForPressure



  void FlowAndTransport::SolveTransientPressure( Model<2U>& model, PDE_Integrator<2U,Element>&  transient_pressure ) {
    // After TransientPressure_Example.cpp
    // 11. Transient loop: Compute fluid pressure during each time-step and output the results for each time step

    // transient pressure
    model.Apply( transient_pressure );

    // output variables screen
    printRangeOfVariable( model, "fluid pressure" );
    printRangeOfVariable( model, "velocity" );
    printRangeOfVariable( model, "pore velocity" );
    printRangeOfVariable( model, "volume flux" );
    printRangeOfVariable( model, "fluid volume source");

    // pressure in the WELL_FRACTURE
    double well_pressure = model.Region("WELL_FRACTURE").Average( "fluid pressure" );
    cout << "\nAverage pressure in WELL_FRACTURE: " << well_pressure << " Pa " << endl;

    //transient_pressure.Reset();

  } // end SolveTransientPressure

}
