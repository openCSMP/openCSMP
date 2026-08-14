//
// Created by Michael Liem on 08.06.22.
//

#include "MeasureDataForEnKF.h"

using namespace std;

namespace csmp{

  void MeasureDataForEnKF::MeasurePressure( Model<2U>& model, const std::string& region_str, const std::string& output_path ) {
    cout << "\nMeasure pressure in region/boundary: " << region_str << endl;

    // Path and name of output file
    std::string output_path_and_name = output_path + "pressure_" + region_str;

    // Storage key to access variables
    csmp::Index P_key =  model.Database().StorageKey("fluid pressure");
    csmp::Index Nn_key = model.Database().StorageKey("node number");

    // Open output file
    ofstream myFile;
    myFile.open((output_path_and_name + ".txt"));

    // Write header
    myFile << "Node-ID" << "\tx" << "\ty" << "\tfluid pressure" << endl;

    // Check if region/boundary exist
    if (model.ContainsRegion(region_str)) {
      // Iterate over all nodes in region_str
      for (auto & nit : model.Region(region_str).NodeVector() ){
        // Coordinate of node
        Point<2> xy= (*nit).Coordinate();

        // Write ID, x-coord, y-coord and fluid pressure of node to file
        myFile << (*nit).Read(Nn_key) << "\t" << xy[0] << "\t" << xy[1] << "\t" << (*nit).Read(P_key) << "\n";
      }

      myFile << std::flush;

      // Close output file
      myFile.close();
    }
    else if (model.ContainsBoundary(region_str)) {
      // Iterate over all nodes in region_str
      for (auto & nit : model.Boundary(region_str).NodeVector() ){
        // Coordinate of node
        Point<2> xy= (*nit).Coordinate();

        // Write ID, x-coord, y-coord and fluid pressure of node to file
        myFile << (*nit).Read(Nn_key) << "\t" << xy[0] << "\t" << xy[1] << "\t" << (*nit).Read(P_key) << "\n";
      }

      myFile << std::flush;

      // Close output file
      myFile.close();
    }
    else {
      // Close output file
      myFile.close();

      cerr << "\nError in MeasurePressure: Region or boundary \"" << region_str << "\" does not exist. Nothing is done" << endl;
      throw std::runtime_error("Unknown region");
    }
  } // end MeasurePressure



  void MeasureDataForEnKF::MeasureVolumeFlowRate( Model<2U>& model, const std::string& region_str, const std::string& output_path ) {
    cout << "\nMeasure volume flow rate in region/boundary: " << region_str << endl;

    // Path and name of output file
    std::string output_path_and_name = output_path + "volumeFlowRate_" + region_str;

    // Storage key to access variables
    csmp::Index VD_key =  model.Database().StorageKey("velocity");
    csmp::Index NVD_key = model.Database().StorageKey("nodal velocity");
    csmp::Index Nn_key =  model.Database().StorageKey("node number");

    // Open output file
    ofstream myFile;
    myFile.open((output_path_and_name + ".txt"));

    // Write header
    myFile << "Node-Nr" << "\tx" << "\ty" << "\tflux_FV" << "\tnodal_velocity_x" << "\tnodal_velocity_y" << endl;

    // Check if region/boundary exist
    if (model.ContainsRegion(region_str)) {
      // Iterate over all nodes in region_str
      for (auto & nit : model.Region(region_str).NodeVector() ){
        // Coordinate of node
        Point<2> xy= (*nit).Coordinate();

        // Flux through finite volume around node
        double flux_fv = fluxThroughFiniteVolume(*nit,VD_key);

        // Nodal velocity at node
        VectorVariable<2> v;
        (*nit).Read(NVD_key,v);

        // Write ID, x-coord, y-coord, flux through FV, x- and y-component of velocity of node to file
        myFile << (*nit).Read(Nn_key) << "\t" << xy[0] << "\t" << xy[1] << "\t" << flux_fv << "\t" << v[0] << "\t" << v[1] << "\n";
      }

      myFile << std::flush;

      // Close output file
      myFile.close();
    }
    else if (model.ContainsBoundary(region_str)) {
      // Iterate over all nodes in region_str
      for (auto & nit : model.Boundary(region_str).NodeVector() ){
        // Coordinate of node
        Point<2> xy= (*nit).Coordinate();

        // Flux through finite volume around node
        double flux_fv = fluxThroughFiniteVolume(*nit,VD_key);

        // Nodal velocity at node
        VectorVariable<2> v;
        (*nit).Read(NVD_key,v);

        // Write ID, x-coord, y-coord, flux through FV, x- and y-component of velocity of node to file
        myFile << (*nit).Read(Nn_key) << "\t" << xy[0] << "\t" << xy[1] << "\t" << flux_fv << "\t" << v[0] << "\t" << v[1] << "\n";
      }

      myFile << std::flush;

      // Close output file
      myFile.close();
    }
    else {
      // Close output file
      myFile.close();

      cout << "\nError in MeasureConcentration: Region or boundary \"" << region_str << "\" does not exist. Nothing is done" << endl;
      throw std::runtime_error("Unknown region");
    }

  } // end MeasureVolumeFlowRate



  void MeasureDataForEnKF::MeasureVolumeFlowRateAlongFracture( Model<2U>& model, const std::string& output_path ) {
    cout << "\nMeasure mean and max volume flow rate in all fractures" << endl;

    // Path and name of output file
    std::string output_path_and_name = output_path + "flow_along_fracture";

    // Storage key to access variables
    csmp::Index TFS_key =  model.Database().StorageKey("total fracture set");
    csmp::Index TFNS_key = model.Database().StorageKey("total fracture number in set");
    csmp::Index VD_key =   model.Database().StorageKey("velocity");
    csmp::Index t_key =    model.Database().StorageKey("thickness");

    // Open output file
    ofstream myFile;
    myFile.open((output_path_and_name + ".txt"));

    // Write header
    myFile << "Set" << "\tFracture in set" << "\tMean" << "\tMax" << endl;

    // Loop through all fracture sets
    auto total_fracture_set = (size_t)model.Read( TFS_key );
    for (size_t current_set = 1; current_set <= total_fracture_set; ++current_set){
      // Check if first fracture in set exists
      string region_str= "FRACTURE_" + to_string(current_set) + "_1";
      if (model.ContainsRegion(region_str)) {
        // Loop through all fractures in current set
        auto total_fracture_number_in_set = (size_t)model.Region(region_str).Read( TFNS_key );
        for(size_t current_number = 1; current_number <= total_fracture_number_in_set; ++current_number ) {
          // Check if region exist
          region_str= "FRACTURE_" + to_string(current_set) + "_" + to_string(current_number);
          if (model.ContainsRegion(region_str)) {
            // Reset some variables
            double max_volume_flow_along_fracture = 0.0;
            double sum_volume_flow_along_fracture = 0.0;
            size_t count_elements= 0;

            // Iterate over all elements in region_str
            for (auto & eit : model.Region(region_str).CellVector() ){
              // Get coordinate of nodes of current (line) element
              std::vector<Point<2>> xy;
              for (auto & nit : (*eit).NodeVector() ) {
                xy.push_back((*nit).Coordinate());
              }

              // Calculate unit vector in fracture direction
              Point<2> fracture_direction;
              fracture_direction[0] = (xy[0][0] - xy[1][0]) / (*eit).SegmentLength(0);
              fracture_direction[1] = (xy[0][1] - xy[1][1]) / (*eit).SegmentLength(0);

              // Velocity at barycenter of element
              VectorVariable<2> v;
              (*eit).Read(VD_key, v);

              // Velocity along fracture
              double volume_flow_along_fracture =
                  std::abs(v[0] * fracture_direction[0] + v[1] * fracture_direction[1]) * (*eit).Read(t_key);
              if (volume_flow_along_fracture > max_volume_flow_along_fracture) {
                max_volume_flow_along_fracture = volume_flow_along_fracture;
              }
              sum_volume_flow_along_fracture += volume_flow_along_fracture;
              count_elements++;
            }

            double mean_volume_flow_along_fracture = sum_volume_flow_along_fracture / double(count_elements);
            myFile << current_set << "\t" << current_number << "\t" << mean_volume_flow_along_fracture << "\t" << max_volume_flow_along_fracture << endl;

          } else {
            cerr << "Error in enkf_run_steadystate.cpp::MeasureVolumeFlowRateAlongFracture(): Region does not exist: " << region_str << endl;
            throw std::runtime_error("Unknown region");
          }
        }
      } else {
        cerr << "Error in enkf_run_steadystate.cpp::MeasureVolumeFlowRateAlongFracture(): Region does not exist: " << region_str << endl;
        throw std::runtime_error("Unknown region");
      }
    }

    // Close output file
    myFile.close();

  } // end MeasureVolumeFlowRateAlongFracture



  void MeasureDataForEnKF::MeasureConcentration( Model<2U>& model, const std::string& region_str, const std::string& output_path, const size_t timestep ) {
    cout << "\nMeasure concentration at time-step " << to_string(timestep) << " in region/boundary: " << region_str << endl;

    // Storage key to access variables
    csmp::Index C_key =  model.Database().StorageKey("concentration");
    csmp::Index Nn_key = model.Database().StorageKey("node number");

    // Path and name of output file
    std::string output_path_and_name = output_path + "concentration_" + region_str; // + "_" + to_string(timestep);

    // Open output file
    ofstream myFile;
    myFile.open((output_path_and_name +  ".txt"), std::ios_base::app);

    // Write header
    if (timestep==0){
      myFile << "t" << "\tNode-Nr" << "\tx" << "\ty" << "\tconcentration" << endl;
    }

    // Check if region/boundary exist
    if (model.ContainsRegion(region_str)) {
      // Iterate over all nodes in region_str
      for (auto & nit : model.Region(region_str).NodeVector() ){
        // Coordinate of node
        Point<2> xy= (*nit).Coordinate();

        // Write ID, x-coord, y-coord and fluid pressure of node to file
        myFile << timestep << "\t" << (*nit).Read(Nn_key) << "\t" << xy[0] << "\t" << xy[1] << "\t" << (*nit).Read(C_key) << "\n";
      }

      myFile << std::flush;

      // Close output file
      myFile.close();
    }
    else if (model.ContainsBoundary(region_str)) {
      // Iterate over all nodes in region_str
      for (auto & nit : model.Boundary(region_str).NodeVector() ){
        // Coordinate of node
        Point<2> xy= (*nit).Coordinate();

        // Write ID, x-coord, y-coord and fluid pressure of node to file
        myFile << timestep << "\t" << (*nit).Read(Nn_key) << "\t" << xy[0] << "\t" << xy[1] << "\t" << (*nit).Read(C_key) << "\n";
      }

      myFile << std::flush;

      // Close output file
      myFile.close();
    }
    else {
      // Close output file
      myFile.close();

      cerr << "\nError in MeasureConcentration: Region or boundary \"" << region_str << "\" does not exist. Nothing is done" << endl;
      throw std::runtime_error("Unknown region");
    }
  } // end MeasureConcentration



  void MeasureDataForEnKF::MeasurePressureAndConcentration( Model<2U>& model, const std::string& region_str, const std::string& output_path, const size_t timestep ) {
    cout << "\nMeasure pressure and concentration at time-step " << to_string(timestep) << " in region/boundary: " << region_str << endl;

    // Storage key to access variables
    csmp::Index p_key =  model.Database().StorageKey("fluid pressure");
    csmp::Index C_key =  model.Database().StorageKey("concentration");
    csmp::Index Nn_key = model.Database().StorageKey("node number");

    // Path and name of output file
    std::string output_path_and_name = output_path + "p_and_c_" + region_str; // + "_" + to_string(timestep);

    // Open output file
    ofstream myFile;
    myFile.open((output_path_and_name +  ".txt"), std::ios_base::app);

    // Write header
    if (timestep==0){
      myFile << "t" << "\tNode-Nr" << "\tx" << "\ty" << "\tpressure" << "\tconcentration" << endl;
    }

    // Check if region/boundary exist
    if (model.ContainsRegion(region_str)) {
      // Iterate over all nodes in region_str
      for (auto & nit : model.Region(region_str).NodeVector() ){
        // Coordinate of node
        Point<2> xy= (*nit).Coordinate();

        // Write ID, x-coord, y-coord and fluid pressure of node to file
        myFile << timestep << "\t" << (*nit).Read(Nn_key) << "\t" << xy[0] << "\t" << xy[1] << "\t" << (*nit).Read(p_key) << "\t" << (*nit).Read(C_key) << "\n";
      }

      myFile << std::flush;

      // Close output file
      myFile.close();
    }
    else if (model.ContainsBoundary(region_str)) {
      // Iterate over all nodes in region_str
      for (auto & nit : model.Boundary(region_str).NodeVector() ){
        // Coordinate of node
        Point<2> xy= (*nit).Coordinate();

        // Write ID, x-coord, y-coord and fluid pressure of node to file
        myFile << timestep << "\t" << (*nit).Read(Nn_key) << "\t" << xy[0] << "\t" << xy[1] << "\t" << (*nit).Read(p_key) << "\t" << (*nit).Read(C_key) << "\n";
      }

      myFile << std::flush;

      // Close output file
      myFile.close();
    }
    else {
      // Close output file
      myFile.close();

      cerr << "\nError in MeasurePressureAndConcentration: Region or boundary \"" << region_str << "\" does not exist. Nothing is done" << endl;
      throw std::runtime_error("Unknown region");
    }
  } // end MeasurePressureAndConcentration



  void MeasureDataForEnKF::WriteSteadyStateVTK( Model<2U>& model, const std::string& output_path) {
    // Create VTK interface
    VTK_Interface<2U> vtk_output;

    // Write some variables fields (at initial time t0) to file
    vtk_output.OutputDataToVTK( model, (output_path + "fluid_pressure"), "fluid pressure", 0 );
    vtk_output.OutputDataToVTK( model, (output_path + "velocity"),       "velocity",       0 );
    vtk_output.OutputDataToVTK( model, (output_path + "volume_flux"),    "volume flux",    0 );
    vtk_output.OutputDataToVTK( model, (output_path + "permeability"),   "permeability",   0 );

  } // end WriteSteadyStateVTK


  void MeasureDataForEnKF::WriteConcentrationVTK( Model<2U>& model, const std::string& output_path, int time_step) {
    // Create VTK interface
    VTK_Interface<2U> vtk_output;

    // Write concentration field at time_step to file
    vtk_output.OutputDataToVTK( model, (output_path + "concentration"),  "concentration",  time_step );

  } // end WriteConcentrationVTK



  void MeasureDataForEnKF::WriteFractureNumbers( Model<2U>& model, const std::string& output_path) {
    // Create text interface
    TextInterface txt_output;

    // Write some region properties
    txt_output.OutputDataAsTextColumns( model, (output_path + "_FS.txt" ).c_str(), "fracture set" );
    txt_output.OutputDataAsTextColumns( model, (output_path + "_FNS.txt" ).c_str(), "fracture number in set" );
    txt_output.OutputDataAsTextColumns( model, (output_path + "_TFNS.txt").c_str(), "total fracture number in set" );
    //txt_output.OutputDataAsTextColumns( model, (output_folder_name + "volume.txt").c_str(), "volume" );
    //txt_output.OutputDataAsTextColumns( model, (output_folder_name + "area.txt").c_str(), "surface area" );

  } // end WriteFractureNumbers


  void MeasureDataForEnKF::CheckMatlabConsistency(Model<2U> &model, const std::string &output_path) {
    //
    ofstream myFile;
    myFile.open((output_path + "/CheckMatlabConsistency.txt"));

    csmp::Index FNS_key =  model.Database().StorageKey("fracture number in set");

    // Loop through all unique fracture region
    for (auto it = model.UniqueRegionsBegin(); it != model.UniqueRegionsEnd(); ++it) {
      // Name of region
      std::string region_str = (*it).first;

      if ((region_str.find("FRACTURE") != string::npos)) {
        Region<2U>& tmp_reg= model.Region(region_str);

        myFile << size_t(tmp_reg.Read(FNS_key)) << "\t";

        // Loop through all nodes
        for (auto nit = tmp_reg.PerimeterNodesBegin(); nit != tmp_reg.NodesEnd(); ++nit ){
          if ( (*nit) != nullptr )
          {
            myFile << (*nit)->x() << "\t" << (*nit)->y() << "\t";
          }
        }
        myFile << endl;
      }
    }



    myFile.close();

  } // end CheckMatlabConsistency


  void MeasureDataForEnKF::PrintOutOfRangeConcentration(Model<2U>& model, const std::string & region_str, const std::string &output_path, size_t timestep, const double c_min, const double c_max) {
    // Storage key to access variables
    csmp::Index C_key =  model.Database().StorageKey("concentration");
    csmp::Index Nn_key = model.Database().StorageKey("node number");

    // Path and name of output file
    std::string output_path_and_name = output_path + "OutOfRangeConcentration_" + region_str;

    // Open output file
    ofstream myFile;
    myFile.open((output_path_and_name +  ".txt"), std::ios_base::app);

    // Write header
    if (timestep==0){
      myFile << "t" << "\tNode-Nr" << "\tx" << "\ty" << "\tconcentration" << endl;
    }

    // Check if region/boundary exist
    if (model.ContainsRegion(region_str)) {
      // Iterate over all nodes in region_str
      for (auto & nit : model.Region(region_str).NodeVector() ){
        // Check if concentration values are outside the range
        if ( (*nit).Read(C_key)<c_min || (*nit).Read(C_key)>c_max ){
          // Coordinate of node
          Point<2> xy= (*nit).Coordinate();

          // Write ID, x-coord, y-coord and fluid pressure of node to file
          myFile << timestep << "\t" << (*nit).Read(Nn_key) << "\t" << xy[0] << "\t" << xy[1] << "\t" << (*nit).Read(C_key) << endl;
        }
      }

      // Close output file
      myFile.close();
    }
    else if (model.ContainsBoundary(region_str)) {
      // Iterate over all nodes in region_str
      for (auto & nit : model.Boundary(region_str).NodeVector() ){
        // Check if concentration values are outside the range
        if ( (*nit).Read(C_key)<c_min || (*nit).Read(C_key)>c_max ) {
          // Coordinate of node
          Point<2> xy = (*nit).Coordinate();

          // Write ID, x-coord, y-coord and fluid pressure of node to file
          myFile << timestep << "\t" << (*nit).Read(Nn_key) << "\t" << xy[0] << "\t" << xy[1] << "\t" << (*nit).Read(C_key) << endl;
        }
      }

      // Close output file
      myFile.close();
    }
    else {
      // Close output file
      myFile.close();

      cerr << "\nError in PrintOutOfRangeConcentration: Region or boundary \"" << region_str << "\" does not exist. Nothing is done" << endl;
      return;
    }

  } // end PrintOutOfRangeConcentration

  void MeasureDataForEnKF::CreateRectangularHollowRegion(Model<2U> &model, const std::string & region_str, Point<2U> corner1, Point<2U> corner2, double_t thickness) {
    Point<2U> xy_min_outer, xy_max_outer, xy_min_inner, xy_max_inner;

    xy_min_outer[0]= corner1[0] - thickness;
    xy_min_outer[1]= corner1[1] - thickness;
    xy_min_inner[0]= corner1[0] + thickness;
    xy_min_inner[1]= corner1[1] + thickness;

    xy_max_inner[0]= corner2[0] - thickness;
    xy_max_inner[1]= corner2[1] - thickness;
    xy_max_outer[0]= corner2[0] + thickness;
    xy_max_outer[1]= corner2[1] + thickness;

    // Create two rectangular regions
    model.FormRectangularRegion("TmpOuter",xy_min_outer,xy_max_outer);
    model.FormRectangularRegion("TmpInner",xy_min_inner,xy_max_inner);

    model.RegionDifference("TmpOuter","TmpInner",region_str.c_str() );

    // Remove temp regions (without deleting their elements)
    model.RemoveRegion("TmpOuter",false);
    model.RemoveRegion("TmpInner",false);

  }

  void MeasureDataForEnKF::CreateMeasRegion(Model<2U> &model, const int meas_region_case ) {
    switch(meas_region_case) {
      case 0: {
        throw std::runtime_error("This does not work yet...");
        std::set<std::string> all_boundaries = {"LEFT", "BOTTOM", "RIGHT", "TOP"};
        model.MergeRegions(all_boundaries, "MeasRegion");
      }
      case 1: {
        // One inner rectangular hollow region
        Point<2U> x1, x2;
        double_t thickness(5.0);
        double_t distance(120.);

        x1[0] = distance;
        x1[1] = distance;
        x2[0] = 720. - distance;
        x2[1] = 720. - distance;

        CreateRectangularHollowRegion(model, "MeasRegion", x1, x2, 5.0);

        cout << "\nCreated 1 inner rectangular hollow measurement region." << endl;

        break;
      }
      case 3: {
        // Three inner rectangular hollow regions
        Point<2U> x1, x2;
        double_t thickness(2.0);
        double_t distance;

        distance = 60.;
        x1[0] = distance;
        x1[1] = distance;
        x2[0] = 720. - distance;
        x2[1] = 720. - distance;

        CreateRectangularHollowRegion(model, "InnerMeasRegion1", x1, x2, thickness);

        distance = 120.;
        x1[0] = distance;
        x1[1] = distance;
        x2[0] = 720. - distance;
        x2[1] = 720. - distance;

        CreateRectangularHollowRegion(model, "InnerMeasRegion2", x1, x2, thickness);

        distance = 180.;
        x1[0] = distance;
        x1[1] = distance;
        x2[0] = 720. - distance;
        x2[1] = 720. - distance;

        CreateRectangularHollowRegion(model, "InnerMeasRegion3", x1, x2, thickness);

        std::set<std::string> all_inner_meas_regions = {"InnerMeasRegion1", "InnerMeasRegion2", "InnerMeasRegion3"};

        model.MergeRegions(all_inner_meas_regions, "MeasRegion");

        model.RemoveRegion("InnerMeasRegion1",false);
        model.RemoveRegion("InnerMeasRegion2",false);
        model.RemoveRegion("InnerMeasRegion3",false);

        cout << "\nCreated 3 inner rectangular hollow measurement regions with thickness " << to_string(thickness) << "." << endl;

        break;
      }
      case 4: {
        Point<2U> x1, x2;
        double_t thickness(4.0);
        double_t location(180.);

        x1[0] = location - 0.5*thickness;
        x1[1] = 0.;
        x2[0] = location + 0.5*thickness;
        x2[1] = 720.;

        model.FormRectangularRegion("TmpReg1",x1,x2);

        location= 360.;
        x1[0] = location - 0.5*thickness;
        x2[0] = location + 0.5*thickness;

        model.FormRectangularRegion("TmpReg2",x1,x2);

        location= 540.;
        x1[0] = location - 0.5*thickness;
        x2[0] = location + 0.5*thickness;

        model.FormRectangularRegion("TmpReg3",x1,x2);


        std::set<std::string> all_tmp_meas_regions = {"TmpReg1", "TmpReg2", "TmpReg3"};

        model.MergeRegions(all_tmp_meas_regions, "MeasRegion");

        model.RemoveRegion("TmpReg1",false);
        model.RemoveRegion("TmpReg2",false);
        model.RemoveRegion("TmpReg3",false);

        cout << "\nCreated 3 line measurement regions with thickness " << to_string(thickness) << "." << endl;

        break;
      }
    }

  } // end CreateMeasRegion
}
