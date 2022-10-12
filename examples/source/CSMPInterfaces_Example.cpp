#include "CSMPInterfaces_Example.h"

// CSMP Interfaces
#include "ANSYS_Model2D.h"
#include "ANSYS_Model3D.h"
#include "TRIANGLE_Interface.h"
#include "Triangulator.h"
#include "Quadrilaterator.h"
#include "VSetConverter.h"
#include "ModelTopology.h"
#include "EclipseModel.h"
#include "IntrepidInterface.h"
#include "RhinoSurfaceReader.h"

#include "Standard_IO_Handler.h"
#include "VTU_Interface.h"
#include "TextInterface.h"
#include "convertColorToPermeability.h"


using namespace std;

namespace csmp{

void CSMPInterfaces_Example::Specifications()
{
  SetTitle( "Build CSMP native format models using various interfaces");
  SetDifficulty( 1 );
  SetCategory( "Software Interfaces" );
  AddAuthor( "Qi Shao" );
  AddDescription( "Build CSMP native format models (binary files) using various interfaces, including ANSYS, SKUA, ECLIPSE, RHINO, TRIANGLE, etc.");
  AddDescription( "These CSMP native format models are to be used in other examples");
  AddRequirement( "source files: 'CSMPInterfaces_Example.cpp' and '*.h'" );
  AddRequirement( "input mesh files that can be read by csmp.");
} 

/** 
    Build CSMP native format models (binary files) using various interfaces, including ANSYS, SKUA, ECLIPSE, RHINO, TRIANGLE, etc.

    User can select which model to build or build all the models that shall be used in other examples.
*/

void CSMPInterfaces_Example::Run()
{
  //set working directory to 'input_meshes/'
  if(!fs::is_directory("input_meshes")) {
    string message("\nCSMPInterfaces_Example::Run(): 'input_meshes/' directory does not exists in current working directory: ");
    message += fs::current_path();
    message += ", please ensure the working directory is set correctly (refer to open-csmp/examples/README.txt)";
    throw std::runtime_error(message);
  } else fs::current_path("input_meshes");


  Standard_IO_Handler  stdio;
  output_vtu_ = stdio.YesNo("\nDo you want to output vtu files for visualisation?");

  do {
    cerr << "\nPlease select from the following options (-1 to get back):" << endl;
    cerr << "  Build CSMP native format models from: " << endl;
    cerr << "    1.) ANSYS 2D model" << endl;
    cerr << "    2.) ANSYS 3D model" << endl;
    cerr << "    3.) TRIANGLE" << endl;
    cerr << "    4.) Triangulator" << endl;
    cerr << "    5.) Quadrilaterator" << endl;
    cerr << "    6.) gOcad_SKUA" << endl;
    cerr << "    7.) Eclipse" << endl;
    cerr << "    8.) GeoModeller" << endl;
    cerr << "    9.) Rhinoceros" << endl;
    cerr << "    10.) Build all models" << endl;

    uint32_t option(0U);
    cin >> option;

    //ANSYS 2D models
    if (option == 1) {
      do {
        cerr << "\nPlease choose from the following ANSYS 2D models (-1 to get back)" << endl;
        cerr << "      1.) box2d_fault" << endl;
        cerr << "      2.) pores" << endl;
        cerr << "      3.) 2000x1000_mesh" << endl;
        cerr << "      4.) Jura-slope1" << endl;
        cerr << "      5.) LeftRight" << endl;
        cerr << "      6.) build all above models" << endl;
        cerr << "      7.) specify another model name" << endl;
        uint32_t sub_option(0U);
        cin >> sub_option;
        if (sub_option == 1) BuildFromANSYS2DModel("box2d_fault");
        else if (sub_option == 2) BuildFromANSYS2DModel("pores");
        else if (sub_option == 3) BuildFromANSYS2DModel("2000x1000_mesh");
        else if (sub_option == 4) BuildFromANSYS2DModel("Jura-slope1");
        else if (sub_option == 5) BuildFromANSYS2DModel("LeftRight");
        else if (sub_option == 6) {
          BuildFromANSYS2DModel("box2d_fault");
          BuildFromANSYS2DModel("pores");
          BuildFromANSYS2DModel("2000x1000_mesh");
          BuildFromANSYS2DModel("Jura-slope1");
          BuildFromANSYS2DModel("LeftRight");
        } else if (sub_option == 7) {
          cerr << "\n        Please type in a model name" << endl;
          string model_name;
          cin >> model_name;
          BuildFromANSYS2DModel(model_name);
        } else if (sub_option == -1) break;
        else {
          cerr << "Wrong number, please type in a right number from the following options" << endl;
          continue;
        }
      } while(true);
      continue;
    }
    //ANSYS 3D models
    else if (option == 2) {
      do {
        cerr << "\nPlease choose from the following ANSYS 3D models (-1 to get back)" << endl;
        cerr << "      1.) hex2_3" << endl;
        cerr << "      2.) prism_test" << endl;
        cerr << "      3.) fracs4" << endl;
        cerr << "      4.) one_sphere_0.45_tetra" << endl;
        cerr << "      5.) FracBox" << endl;
        cerr << "      6.) 3D_box" << endl;
        cerr << "      7.) build all above models" << endl;
        cerr << "      8.) specify another model name" << endl;
        uint32_t sub_option(0U);
        cin >> sub_option;
        if (sub_option == 1) BuildFromANSYS3DModel("hex2_3");
        else if (sub_option == 2) BuildFromANSYS3DModel("prism_test");
        else if (sub_option == 3) BuildFromANSYS3DModel("fracs4");
        else if (sub_option == 4) BuildFromANSYS3DModel("one_sphere_0.45_tetra");
        else if (sub_option == 5) BuildFromANSYS3DModel("FracBox");
        else if (sub_option == 6) BuildFromANSYS3DModel("3D_box");
        else if (sub_option == 7) {
          BuildFromANSYS3DModel("hex2_3");
          BuildFromANSYS3DModel("prism_test");
          BuildFromANSYS3DModel("fracs4");
          BuildFromANSYS3DModel("one_sphere_0.45_tetra");
          BuildFromANSYS3DModel("FracBox");
          BuildFromANSYS3DModel("3D_box");
        } else if (sub_option == 8) {
          cerr << "\n        Please type in a model name" << endl;
          string model_name;
          cin >> model_name;
          BuildFromANSYS3DModel(model_name);
        } else if (sub_option == -1) break;
        else {
          cerr << "Wrong number, please type in a right number from the following options" << endl;
          continue;
        }
      } while(true);
      continue;
    }
    //TRIANGLE models
    else if (option == 3) {
      do {
        cerr << "\nPlease choose from the following TRIANGLE models (-1 to get back)" << endl;
        cerr << "      1.) well.1" << endl;
        cerr << "      2.) veins_20k.1" << endl;
        cerr << "      3.) frac30.1" << endl;
        cerr << "      4.) blunt30deg.1" << endl;
        cerr << "      5.) topo.1" << endl;
        cerr << "      6.) example21.1" << endl;
        cerr << "      7.) build all above models" << endl;
        cerr << "      8.) specify another model name" << endl;
        uint32_t sub_option(0U);
        cin >> sub_option;
        if (sub_option == 1) BuildFromTRIANGLEModel("well.1");
        else if (sub_option == 2) BuildFromTRIANGLEModel("veins_20k.1");
        else if (sub_option == 3) BuildFromTRIANGLEModel("frac30.1");
        else if (sub_option == 4) BuildFromTRIANGLEModel("blunt30deg.1");
        else if (sub_option == 5) BuildFromTRIANGLEModel("topo.1");
        else if (sub_option == 6) BuildFromTRIANGLEModel("example21.1");
        else if (sub_option == 7) {
          BuildFromTRIANGLEModel("well.1");
          BuildFromTRIANGLEModel("veins_20k.1");
          BuildFromTRIANGLEModel("frac30.1");
          BuildFromTRIANGLEModel("blunt30deg.1");
          BuildFromTRIANGLEModel("topo.1");
          BuildFromTRIANGLEModel("example21.1");
        } else if (sub_option == 8) {
          cerr << "\n        Please type in a model name" << endl;
          string model_name;
          cin >> model_name;
          BuildFromTRIANGLEModel(model_name);
        } else if (sub_option == -1) break;
        else {
          cerr << "Wrong number, please type in a right number from the following options" << endl;
          continue;
        }
      } while(true);
      continue;
    }
    //Triangulator models
    else if (option == 4) {
      do {
        cerr << "\nPlease choose from the following Triangulator models (-1 to get back)" << endl;
        cerr << "      1.) tutorial1_input" << endl;
        cerr << "      2.) specify another model name" << endl;
        uint32_t sub_option(0U);
        cin >> sub_option;
        if (sub_option == 1) BuildFromTriangulatorModel("tutorial1_input", 100., 80.);
        else if (sub_option == 2) {
          cerr << "\n        Please type in a model name" << endl;
          string model_name;
          cin >> model_name;
          double x, y;
          cout << "\n        Please type in the x- and y-dimensions of your model (in m): " << endl;
          cin >> x;
          cin >> y;
          BuildFromTriangulatorModel(model_name, x, y);
        } else if (sub_option == -1) break;
        else {
          cerr << "Wrong number, please type in a right number from the following options" << endl;
          continue;
        }
      } while(true);
      continue;
    }
    //Quadrilaterator models
    else if (option == 5) {
      do {
        cerr << "\nPlease choose from the following Quadrilaterator models (-1 to get back)" << endl;
        cerr << "      1.) tutorial2_input" << endl;
        cerr << "      2.) tutorial2_input_20x20" << endl;
        cerr << "      3.) build all above models" << endl;
        cerr << "      4.) specify another model name" << endl;
        uint32_t sub_option(0U);
        cin >> sub_option;
        if (sub_option == 1) BuildFromQuadrilateratorModel("tutorial2_input", 10., 10.);
        else if (sub_option == 2) BuildFromQuadrilateratorModel("tutorial2_input_20x20", 20., 20.);
        else if (sub_option == 3) {
          BuildFromQuadrilateratorModel("tutorial2_input", 10., 10.);
          BuildFromQuadrilateratorModel("tutorial2_input_20x20", 20., 20.);
        }
        else if (sub_option == 4) {
          cerr << "\n        Please type in a model name" << endl;
          string model_name;
          cin >> model_name;
          double x, y;
          cout << "\n        Please type in the x- and y-dimensions of your model (in m): " << endl;
          cin >> x;
          cin >> y;
          BuildFromQuadrilateratorModel(model_name, x, y);
        } else if (sub_option == -1) break;
        else {
          cerr << "Wrong number, please type in a right number from the following options" << endl;
          continue;
        }
      } while(true);
      continue;
    }
    //gOcad_SKUA models
    else if (option == 6) {
      do {
        cerr << "\nPlease choose from the following SKUA models (-1 to get back)" << endl;
        cerr << "      1.) SKUA_boundary_and_split_boundary" << endl;
        cerr << "      2.) SKUA_box_shaped_with_boundary" << endl;
        cerr << "      3.) SKUA_box_shaped_with_edges" << endl;
        cerr << "      4.) SKUA_box_shaped_with_split_boundary" << endl;
        cerr << "      5.) SKUA_cross_bedded_xsmall" << endl;
        cerr << "      6.) SKUA_intersecting_split_boundary_surfaces" << endl;
        cerr << "      7.) SKUA_model" << endl;
        cerr << "      8.) SKUA_split_boundary_layer" << endl;
        cerr << "      9.) SKUA_split_boundary_layer_with_well" << endl;
        cerr << "      10.) SKUA_split_boundary_surface" << endl;
        cerr << "      11.) SKUA_box_snowflake_split_boundaries" << endl;
        cerr << "      12.) build all above models" << endl;
        cerr << "      13.) specify another model name" << endl;
        uint32_t sub_option(0U);
        cin >> sub_option;
        if (sub_option == 1) BuildFromSKUAModel<3U>("SKUA_boundary_and_split_boundary");
        else if (sub_option == 2) BuildFromSKUAModel<3U>("SKUA_box_shaped_with_boundary");
        else if (sub_option == 3) BuildFromSKUAModel<3U>("SKUA_box_shaped_with_edges");
        else if (sub_option == 4) BuildFromSKUAModel<3U>("SKUA_box_shaped_with_split_boundary");
        else if (sub_option == 5) BuildFromSKUAModel<3U>("SKUA_cross_bedded_xsmall");
        else if (sub_option == 6) BuildFromSKUAModel<3U>("SKUA_intersecting_split_boundary_surfaces");
        else if (sub_option == 7) BuildFromSKUAModel<3U>("SKUA_model");
        else if (sub_option == 8) BuildFromSKUAModel<3U>("SKUA_split_boundary_layer");
        else if (sub_option == 9) BuildFromSKUAModel<3U>("SKUA_split_boundary_layer_with_well");
        else if (sub_option == 10) BuildFromSKUAModel<3U>("SKUA_split_boundary_surface");
        else if (sub_option == 11) BuildFromSKUAModel<3U>("SKUA_box_snowflake_split_boundaries");
        else if (sub_option == 12) {
          BuildFromSKUAModel<3U>("SKUA_boundary_and_split_boundary");
          BuildFromSKUAModel<3U>("SKUA_box_shaped_with_boundary");
          BuildFromSKUAModel<3U>("SKUA_box_shaped_with_edges");
          BuildFromSKUAModel<3U>("SKUA_box_shaped_with_split_boundary");
          BuildFromSKUAModel<3U>("SKUA_cross_bedded_xsmall");
          BuildFromSKUAModel<3U>("SKUA_intersecting_split_boundary_surfaces");
          BuildFromSKUAModel<3U>("SKUA_model");
          BuildFromSKUAModel<3U>("SKUA_split_boundary_layers");
          BuildFromSKUAModel<3U>("SKUA_split_boundary_layer_with_well");
          BuildFromSKUAModel<3U>("SKUA_split_boundary_surface");
          BuildFromSKUAModel<3U>("SKUA_box_snowflake_split_boundaries");
        } else if (sub_option == 13) {
          cerr << "\n        Please type in a model name" << endl;
          string model_name;
          cin >> model_name;
          BuildFromSKUAModel<3U>(model_name);
        } else if (sub_option == -1) break;
        else {
          cerr << "Wrong number, please type in a right number from the following options" << endl;
          continue;
        }
      } while(true);
      continue;
    }
      //Eclipse model
    else if (option == 7) {
      do {
        cerr << "\nPlease choose from the following options (-1 to get back)" << endl;
        cerr << "      1.) create CSMP native model from Eclipse model 'NPD5.grdecl'" << endl;
        cerr << "      2.) specify another Eclipse model name" << endl;
        uint32_t sub_option(0U);
        cin >> sub_option;
        if (sub_option == 1) BuildFromEclipseModel("NPD5");
        else if (sub_option == 2) {
          cerr << "\n        Please type in a model name" << endl;
          string model_name;
          cin >> model_name;
          BuildFromEclipseModel(model_name);
        } else if (sub_option == -1) break;
        else {
          cerr << "Wrong number, please type in a right number from the following options" << endl;
          continue;
        }
      } while(true);
      continue;
    }
      //GeoModeller
    else if (option == 8) {
      do {
        cerr << "\nPlease choose from the following options (-1 to get back)" << endl;
        cerr << "      1.) create CSMP native model from GeoModeller model 'Mansfield_H8_NoOrphans.mesh'" << endl;
        cerr << "      2.) specify another GeoModeller model name" << endl;
        uint32_t sub_option(0U);
        cin >> sub_option;
        if (sub_option == 1) BuildFromGeoModellerModel("Mansfield_H8_NoOrphans.mesh");
        else if (sub_option == 2) {
          cerr << "\n        Please type in a model name" << endl;
          string model_name;
          cin >> model_name;
          BuildFromGeoModellerModel(model_name);
        } else if (sub_option == -1) break;
        else {
          cerr << "Wrong number, please type in a right number from the following options" << endl;
          continue;
        }
      } while(true);
      continue;
    }
      //Rhino model
    else if (option == 9) {
      do {
        cerr << "\nPlease choose from the following options (-1 to get back)" << endl;
        cerr << "      1.) create CSMP native model from Rhino model 'example20.raw'" << endl;
        cerr << "      2.) specify another Rhino model name" << endl;
        uint32_t sub_option(0U);
        cin >> sub_option;
        if (sub_option == 1) BuildFromRhinoModel("example20.rawh");
        else if (sub_option == 2) {
          cerr << "\n        Please type in a model name" << endl;
          string model_name;
          cin >> model_name;
          BuildFromRhinoModel(model_name);
        } else if (sub_option == -1) break;
        else {
          cerr << "Wrong number, please type in a right number from the following options" << endl;
          continue;
        }
      } while(true);
      continue;
    }
      //Build all models
    else if (option == 10) {
      BuildFromANSYS2DModel("box2d_fault");
      BuildFromANSYS2DModel("pores");
      BuildFromANSYS2DModel("2000x1000_mesh");
      BuildFromANSYS2DModel("Jura-slope1");
      BuildFromANSYS2DModel("LeftRight");

      BuildFromANSYS3DModel("hex2_3");
      BuildFromANSYS3DModel("prism_test");
      BuildFromANSYS3DModel("fracs4");
      BuildFromANSYS3DModel("one_sphere_0.45_tetra");
      BuildFromANSYS3DModel("FracBox");
      BuildFromANSYS3DModel("3D_box");

      BuildFromTRIANGLEModel("well.1");
      BuildFromTRIANGLEModel("veins_20k.1");
      BuildFromTRIANGLEModel("frac30.1");
      BuildFromTRIANGLEModel("blunt30deg.1");
      BuildFromTRIANGLEModel("topo.1");
      BuildFromTRIANGLEModel("example21.1");

      BuildFromTriangulatorModel("tutorial1_input", 100., 80.);

      BuildFromQuadrilateratorModel("tutorial2_input", 10., 10.);
      BuildFromQuadrilateratorModel("tutorial2_input_20x20", 20., 20.);

      BuildFromSKUAModel<3U>("SKUA_boundary_and_split_boundary");
      BuildFromSKUAModel<3U>("SKUA_box_shaped_with_boundary");
      BuildFromSKUAModel<3U>("SKUA_box_shaped_with_edges");
      BuildFromSKUAModel<3U>("SKUA_box_shaped_with_split_boundary");
      BuildFromSKUAModel<3U>("SKUA_cross_bedded_xsmall");
      BuildFromSKUAModel<3U>("SKUA_intersecting_split_boundary_surfaces");
      BuildFromSKUAModel<3U>("SKUA_model");
      BuildFromSKUAModel<3U>("SKUA_split_boundary_layers");
      BuildFromSKUAModel<3U>("SKUA_split_boundary_layer_with_well");
      BuildFromSKUAModel<3U>("SKUA_split_boundary_surface");
      BuildFromSKUAModel<3U>("SKUA_box_snowflake_split_boundaries");

      BuildFromEclipseModel("NPD5");

      BuildFromGeoModellerModel("Mansfield_H8_NoOrphans.mesh");

      BuildFromRhinoModel("example20.rawh");

      break;
    }
    //quit
    else if (option == -1) {fs::current_path("../"); break;}
    //re-select
    else {
      cerr << "Wrong number, please type in a right number from the following options" << endl;
      continue;
    }
  } while(true);

} // Run()



void CSMPInterfaces_Example::BuildFromANSYS2DModel(const string& model_name) {
  cout<<"\nStart building ANSYS 2D model '"<<model_name<<"'..."<<endl;
  const string variables_file = "initial-variables.txt";
  ANSYS_Model2D model( model_name.c_str(), variables_file.c_str() );
  fs::create_directory("csmp_native_format_models");
  fs::current_path("csmp_native_format_models");
  model.OutputToBinaryFile( model_name.c_str() );
  if(output_vtu_) OutputRegionIDToVTU(model);
  cout<<"\nANSYS 2D model '"<<model_name<<"' saved to disk"<<endl;
  fs::current_path("../");
}


void CSMPInterfaces_Example::BuildFromANSYS3DModel(const string& model_name) {
  cout<<"\nStart building ANSYS 3D model '"<<model_name<<"'..."<<endl;
  const string variables_file = "initial-variables.txt";
  ANSYS_Model3D model( model_name.c_str(), variables_file.c_str() );
  fs::create_directory("csmp_native_format_models");
  fs::current_path("csmp_native_format_models");
  model.OutputToBinaryFile( model_name.c_str() );
  if(output_vtu_) OutputRegionIDToVTU(model);
  cout<<"\nANSYS 3D model '"<<model_name<<"' saved to disk"<<endl;
  fs::current_path("../");
}


void CSMPInterfaces_Example::BuildFromTRIANGLEModel(const string& model_name) {
  cout<<"\nStart building TRIANGLE model '"<<model_name<<"'..."<<endl;
  TRIANGLE_Interface  mesh_interface;
  VSet<2U>            mesh_container;
  VSetConverter<2U>   mesh_converter;
  mesh_interface.ReadTriangle2DMesh( model_name.c_str(), mesh_container );
  mesh_converter.ConvertLinearToQuadraticTriangles( mesh_container );
  const string variables_file = "initial-variables.txt";
  Model<2U>   model( mesh_container, variables_file.c_str() );
  fs::create_directory("csmp_native_format_models");
  fs::current_path("csmp_native_format_models");
  model.OutputToBinaryFile( model_name.c_str() );
  model.Name(model_name.c_str());
  if(output_vtu_) OutputRegionIDToVTU(model);
  cout<<"\nTRIANGLE model '"<<model_name<<"' saved to disk"<<endl;
  fs::current_path("../");
}


void CSMPInterfaces_Example::BuildFromTriangulatorModel(const string& model_name, double extent1, double extent2) {
  cout << "\nStart building Triangulator model '" << model_name << "'..." << endl;

  // 0.1 Reading a pixelated permeability image ASCII file into the new Matrix 'pixelcolors'
  size_t m, n;
  TextInterface().SizeofPixelTextImage256( model_name.c_str(), m, n );
  cout <<"\nreadTextPixelData: The size (in pixels) of the input image is: "<< n <<"h x "<< m;
  cout <<"v"<< endl;
  Matrix  pixelcolors(m,n);
  TextInterface().ReadPixelTextImage256( model_name.c_str(), pixelcolors );

  // 0.2 Converting the 256-color values into permeabilities
  convertColorToPermeability( 1., pixelcolors );

  // 0.3 Building a 2d mesh of triangular elements, incorporating the
  //    permeability data and the boundary conditions. The mesh is
  //    stored in the VSet object 'vset'
  VSet<2U>  vset;
  Triangulator().TrianglesFromRegularGrid( pixelcolors, vset );

  // 0.4 Scaling the geometrical input object that will become the Region
  //    The origin of the object is assumed to be zero.
  double zero(0.);
  vset.ScaleCoordinateToRange( 'x', zero, extent1 );
  vset.CoordinateRange( 'x', zero, extent1 );
  cout <<"\nAssigned X range: "<< zero <<" to "<< extent1 << " meter." << endl;
  vset.ScaleCoordinateToRange( 'y', zero, extent2 );
  vset.CoordinateRange( 'y', zero, extent2 );
  cout <<"\nAssigned Y range: "<< zero <<" to "<< extent2 << " meter." << endl;

  const string variables_file = "initial-variables.txt";
  Model<2U> model( vset, variables_file.c_str() );
  fs::create_directory("csmp_native_format_models");
  fs::current_path("csmp_native_format_models");
  model.OutputToBinaryFile( model_name.c_str() );
  model.Name(model_name.c_str());
  if(output_vtu_) OutputRegionIDToVTU(model);
  cout<<"\nTriangulator model '"<<model_name<<"' saved to disk"<<endl;
  fs::current_path("../");
}


void CSMPInterfaces_Example::BuildFromQuadrilateratorModel(const string& model_name, double x, double y) {
  cout << "\nStart building Quadrilaterator model '" << model_name << "'..." << endl;
  Quadrilaterator    quadrilaterator; // simple FE mesher
  VSet<2U>           mesh_container;  // container to store the input mesh

  // read in file and generate mesh
  quadrilaterator.QuadrilateralsFromRegularGrid( mesh_container, model_name.c_str(), x, y );
  const string variables_file = "initial-variables.txt";
  Model<2U> model( mesh_container, variables_file.c_str() ); // Quadrilaterator makes isoparametric FEs
  fs::create_directory("csmp_native_format_models");
  fs::current_path("csmp_native_format_models");
  model.OutputToBinaryFile( model_name.c_str() );
  model.Name(model_name.c_str());
  if(output_vtu_) OutputRegionIDToVTU(model);
  cout<<"\nQuadrilaterator model '"<<model_name<<"' saved to disk"<<endl;
  fs::current_path("../");
}


template<uint32_t dim>
void CSMPInterfaces_Example::BuildFromSKUAModel(const string& model_name) {
  cout<<"\nStart building SKUA model '"<<model_name<<"'..."<<endl;
  VSet<dim> vset;
  double model_time = 0.0;
  vset.InputFrom( model_name.c_str(), model_time );

  // Read ModelTopology file generated by SKUA
  ModelTopology model_topology(true);
  model_topology.InputFromTextFile( model_name.c_str() );

  // Build Model using variables file generated by SKUA
  //const string variables_file = "initial-variables.txt";
  //const string variables_file = model_name + "-variables.txt";
  const string variables_file = "SKUA_model-variables.txt";
  if (!fs::exists(variables_file)) {
    string error_message = "\n\nError: file '";
    string path = fs::current_path();
    error_message += (variables_file + "' does not exist in directory " + path);
    error_message += (", example cannot run, please check.\n");
    throw std::runtime_error(error_message);
  }
  Model<dim>  model( model_topology, vset, variables_file.c_str(), false );
  fs::create_directory("csmp_native_format_models");
  fs::current_path("csmp_native_format_models");
  model.OutputToBinaryFile( model_name.c_str() );
  if(output_vtu_) OutputRegionIDToVTU(model);
  cout<<"\nSKUA model '"<<model_name<<"' saved to disk"<<endl;
  fs::current_path("../");

}
template void CSMPInterfaces_Example::BuildFromSKUAModel<1U>(const string&);
template void CSMPInterfaces_Example::BuildFromSKUAModel<2U>(const string&);
template void CSMPInterfaces_Example::BuildFromSKUAModel<3U>(const string&);


void CSMPInterfaces_Example::BuildFromEclipseModel(const string& model_name) {
  cout<<"\nStart building Eclipse model '"<<model_name<<"'..."<<endl;
  // Setup mesh
  std::string regions_file(model_name);
  bool exclude_inactive_cells(true);
  bool tetra_mesh(false);
  const bool create_boundaries(false);
  EclipseModelSettings settings(model_name);
  settings.MeshSetup(regions_file,
                     exclude_inactive_cells,
                     tetra_mesh,
                     create_boundaries);
  // Setup properties
  // porosity
  std::string   porosity_name("porosity");
  VARIABLE_TYPE porosity_type(SCALAR);
  PLACEMENT     porosity_place(ELEMENT);
  settings.PoroPropertySetup(porosity_name,porosity_type,porosity_place);
  // permeability
  std::string   permeability_name("tensor permeability");
  VARIABLE_TYPE permeability_type(TENSOR);
  PLACEMENT     permeability_place(ELEMENT);
  std::string   permeability_unit("mD"); // other options m2,D
  settings.PermPropertySetup(permeability_name,permeability_type,
                             permeability_place,permeability_unit);
  // rock type
  std::string   rocktype_name("rock type");
  VARIABLE_TYPE rocktype_type(SCALAR);
  PLACEMENT     rocktype_place(ELEMENT);
  settings.RockNumPropertySetup(rocktype_name,rocktype_type,rocktype_place);

  const string variables_file = "initial-variables.txt";
  EclipseModel modelOut(settings, model_name, variables_file);

  // Get Fault Regions
  std::vector<std::string> faults;
  modelOut.GetFaults(faults);

  // Get Well Regions
  std::vector<std::string> wells;
  modelOut.GetWells(wells);

  // We are in the process of rewriting the Eclipse interface, and the
  // following part is not yet fully ported.  - AJB

  auto&  model_domain(modelOut.Region("Model"));
  modelOut.Region("Model").UpdateMemberIndexes();

  // 3. eliminating any potentially disfunctional elements / cells from the model
  // ----------------------------------------------------------------------------
  // elements that have a negative Jacobian determinant are assumed to be degenerate and flagged for deletion
  vector<uint32_t> degenerate_elements;
  int volume_e_removed(0U), surface_e_removed(0U), line_e_removed(0U);
  for (auto it = model_domain.CellsBegin(); it != model_domain.CellsEnd(); ++it) {
    // find broken elements
    // (an element is regarded as broken if the determinant of its Jacobian inverse is negative at least
    //  at one of the integration points
    bool broken_elmt(false);
    for (size_t ipoint = 0U; ipoint<(*it)->IntegrationPoints(); ++ipoint)
      if ((*it)->det_JINV_AtIntegrationPoint(ipoint) <= 0.) {
        broken_elmt = true;
        break;
      }
    if (broken_elmt) {
      if ((*it)->IsVolume()) volume_e_removed++;
      else if ((*it)->IsSurface()) surface_e_removed++;
      else if ((*it)->IsLine()) line_e_removed++;
      degenerate_elements.push_back((*it)->Idx());
      auto eclipseCoord = modelOut.EclipseCoordinates(*it);
      std::cerr << "Broken element at " << eclipseCoord.i << ' ' << eclipseCoord.j << ' ' << eclipseCoord.k << ' '
                << parseFiniteElementType((*it)->FE()->ElementType())
                << '\n';
#if 0
      for (auto nit = (*it)->NodesBegin(); nit != (*it)->NodesEnd(); ++nit) {
					std::cerr << (*nit)->Coordinate() << '\n';
				}
				std::cerr << "Done\n";
#endif
    }
  }
  if (volume_e_removed > 0 || surface_e_removed > 0 || line_e_removed > 0) {
    cout << "\nread_and_configure_ECLIPSE_model: removing degenerate elements:\n";
    cout << "\n\tvolume elements removed:  " << volume_e_removed;
    cout << "\n\tsurface elements removed: " << surface_e_removed;
    cout << "\n\tline elements removed:    " << line_e_removed;
  }

  fs::create_directory("csmp_native_format_models");
  fs::current_path("csmp_native_format_models");
  modelOut.OutputToBinaryFile( model_name.c_str() );
  if(output_vtu_) OutputRegionIDToVTU(modelOut);
  cout<<"\nEclipse model '"<<model_name<<"' saved to disk"<<endl;
  fs::current_path("../");
}


void CSMPInterfaces_Example::BuildFromGeoModellerModel(const string& model_name) {
  cout<<"\nStart building GeoModeller model '"<<model_name<<"'..."<<endl;
  IntrepidInterface  geomodel;
  VSet<3U>           vset;
  ModelTopology      model_topology(true);
  geomodel.Read( model_name.c_str(), vset, model_topology );

  const string variables_file = "initial-variables.txt";
  Model<3U>  model(model_topology, vset, variables_file.c_str(), false );

  fs::create_directory("csmp_native_format_models");
  fs::current_path("csmp_native_format_models");
  model.OutputToBinaryFile( model_name.c_str() );
  if(output_vtu_) OutputRegionIDToVTU(model);
  cout<<"\nGeoModeller model '"<<model_name<<"' saved to disk"<<endl;
  fs::current_path("../");

}

void CSMPInterfaces_Example::BuildFromRhinoModel(const string& model_name) {
  cout<<"\nStart building Rhino model '"<<model_name<<"'..."<<endl;
  //  reading in '.raw' meshes that were written to file as labeledentities
  SKM_RhinoSurfaceReader  rhino_surface( model_name.c_str() );

  VSet<3U>  mesh_container;
  rhino_surface.OutputObjectTo( "fault1", mesh_container );
  // setting to a numerically integrated element as this is later needed for the transport calculation
  mesh_container.SingleElementType(ISOPARAMETRIC_LINEAR_TRIANGLE);

  PropertyData  permdata( ELEMENT, SCALAR, 3U );
  permdata.Reserve( mesh_container.Elements() );
  for ( auto i{0}; i<mesh_container.Elements(); ++i )
    pushBack( permdata, makeScalar( PLAIN, 1.0e-12 ) );
  mesh_container.AddData( "permeability", permdata );

  const string variables_file = "initial-variables.txt";
  Model<3U>  model3D( mesh_container, variables_file.c_str() );
  fs::create_directory("csmp_native_format_models");
  fs::current_path("csmp_native_format_models");
  model3D.OutputToBinaryFile( model_name.c_str() );
  if(output_vtu_) OutputRegionIDToVTU(model3D);
  cout<<"\nRhino model '"<<model_name<<"' saved to disk"<<endl;
  fs::current_path("../");
}


template<uint32_t dim>
void CSMPInterfaces_Example::OutputRegionIDToVTU (Model<dim>& model)
{
  if(!model.Database().IsDefined("region id"))
    model.CreateProperty( "region id", "rid", "none", SCALAR, ELEMENT, 1, 0, 1000);

  model.InputPropertyValue("region id", makeScalar( PLAIN,0.0));

  uint32_t region_id(1U);
  for(auto rit=model.UniqueRegionsBegin(); rit != model.UniqueRegionsEnd(); rit++) {
    (*rit).second.InputPropertyValue("region id", makeScalar( PLAIN, region_id));
    region_id++;
  }
  VTU_Interface<dim>  vtu(model);
  string output_file = model.Name();
  output_file += "-region_id";
  vtu.OutputDataToVTU( output_file, "region id", "Model", 0 );
  model.DeleteProperty("region id");
}
template void CSMPInterfaces_Example::OutputRegionIDToVTU (Model<1U>& model);
template void CSMPInterfaces_Example::OutputRegionIDToVTU (Model<2U>& model);
template void CSMPInterfaces_Example::OutputRegionIDToVTU (Model<3U>& model);

} // csmp
