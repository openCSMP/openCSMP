//
//  Experimental_Example.cpp
//  CSMP_API_library2014
//
//  Created by Stephan Matthai on 1/27/14.
//  Copyright (c) 2014 Stephan Matthai. All rights reserved.
//

#include "Experimental_Example.h"
#include "compareFloats.h"

#ifdef CSMP_WITH_SAMG_SOLVER
#include "SAMG_Solver.h"
#include "SAMG_Settings.h"
#endif
#ifdef CSMP_WITH_MESCHACH
#include "Gauss_Solver.h"
#endif

#include "Boundary.h"
#include "DES2PhaseSlightlyCompressibleTransport.h"
#include "ErrorHandler.h"
#include "FlowFunctionsModule.h"
#include "InputDataManager.h"
#include "ModelTopology.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_dNT_op_dV.h"
#include "NumIntegral_NT_op_N_dV.h"
#include "PDE_Integrator.h"
#include "LinearSolver.h"
#include "Region.h"
#include "RegionMonitor.h"
#include "SplitBoundary.h"
#include "Standard_IO_Handler.h"
#include "VTK_Interface.h"
#include "VTU_Interface.h"
#include "DynamicArray2D.h"

#include "ANSYS_Model3D.h"
#include "UG4_UGX_FileExport.h"
#include "vsetMakers.h"


using namespace std;

namespace csmp {

void Experimental_Example::Specifications()
  {
     SetTitle( "Experimental_Example" );
     SetDifficulty( 1 );
     SetCategory( "Software Functionality" );
     AddAuthor( "You!" );
     AddDescription( "source in: Experimental_Example.cpp" );
     AddDescription( "Empty example for the user to experiment with" );
     AddRequirement( "none" );
     AddRequirement( "no predefined model or variables file" );
  }




void Experimental_Example::Run()
 {
    // a little demonstration of the pointer concept used extensively in C and C++
    int a{3};
    // pointer to the integer a
    int* ptr = &a;
    int* ptr2{ nullptr }; // safe initalisation of pointer
    delete ptr2;  // de-allocate the memory pointed to by a (in this case none;
    int ptr3[] = { 0, 1, 2, 3, 4, 5 }; // C-style fixed array based on the pointer concept
    // accessing the array's fourth element (counting from 0..arraysize-1
    *(ptr3 + 3) = 24;
    // assignment of values
    ptr3[3] = 24;
    
    // dynamic array (which lives in heap memory)
    int* array_ptr = new int[5];
    // size of array pointer (not the array)
    cout <<"\nmain: "<< sizeof(array_ptr) <<" (size of pointer = machine word) vs. "<< sizeof(int) <<" (integer)."<< endl;
    // what is the size of the array? - some work to calculate for a dynamic array
    cout <<"\nmain: array size: "<< 5 * sizeof(int) + sizeof(array_ptr) << endl;
    // now the array has to be de-allocated, otherwise we have a memory leak
    delete[] array_ptr;
    // IN CONTEMPORARY PRACTICE DO NOT USE RAW POINTERS, use autoptr instead
    std::shared_ptr<int> aptr{&a};
    // aptr will get deleted automatically when the program finishes, but its management incurs some overhead
    
    cout <<"main: value of a: ";
    cout << a <<" or "<< *ptr;
    // check the validity of the pointer before using (derefencing) it
    if ( ptr2 != nullptr ) cout << *ptr2;
    // shorthand
    if ( ptr2 ) cout << *ptr2;
 
 
    // testing numberToString and matrix output
    // ----------------------------------------
    double dvalue{ 1.0e-15 };
    float  fvalue{ 1.0e-01 };
    cout <<"\nexample: floating point numbers: "<< dvalue <<" vs. "<< fvalue << endl;

    cout <<"\n"<<"experimental example: sum: "<< numberToString( fvalue + dvalue );

    cout <<"\n"<<"experimental example: float with size "<< sizeof(dvalue) <<": "<< numberToString( dvalue );
    cout <<"\n"<<"experimental example: float with size "<< sizeof(fvalue) <<": "<< numberToString( fvalue );

 
    // comparing 1D array with vector of vectors
    DynamicArray2D<double> mat(3,5);
    // unsupported 2D initialisation: DynamicArray2D<double> mat2{ {1.,0}, {0.,2.} };
    // should work: DynamicArray2D<double> mat2{ 1., 0., 0., 2. };
    // would not know m, n: DynamicArray2D<double> mat2( vector<double>{1., 0., 0., 2.} );
//    DynamicArray2D<double> mat3{ 1., 0., 0., 2. }; // square matrix
    DynamicArray2D<double> mat4{ {1., 0.}, {0., 2.} };
    // assignments
//    DynamicArray2D<double> mat5 = { 1., 0., 0., 2. };
  
    // 2D Test case without SplitBoundary objects
    // ------------------------------------------
    VSet<2U> mesh;
    ModelTopology topo = test_Create_MeshPatchWithLineElements_VSet( mesh );
    const bool treat_domains_as_regions{true}; // model does not contain any Face objects!
    Model<2U> model2D( topo, mesh, "UG4_UGX_FileExport-variables.txt", treat_domains_as_regions );

    // intialising the variables 'element variable' and 'element vector'
    //model.InputPropertyValue( "element variable", makeScalar(ANY,1.0) );
    model2D.InputPropertyValue( "element vector", makeVector(ANY,ANY,1.0,2.0) );
    
    // forcing the creation of another unique region called quadrilaterals because Promesh does not show properties
    const bool unique_region{ true };
    model2D.FormRegionFrom( "quads", "element variable", 3.5, 5.0, unique_region );

    VTU_Interface<2U>  vtu(model2D);
    list<string> output_props{ "node number", "element number", "element variable" };
    vtu.OutputDataToVTU( "test_Create_MeshPatchWithLineElements_VSet", output_props, "Model", static_cast<int>(0) );
    
    // OUTPUTS MODEL TO UG (name will be the model name)
    UG4_UGX_FileExport<2U> ug4_exporter2D( model2D );
    ug4_exporter2D.Write_UGX_FileASCII( model2D, model2D.Name() );
    
    
    // 3D 'prism_mesh' testcase
    // ------------------------
    ANSYS_Model3D  model3D( "prism_test", "UG4_UGX_FileExport-variables.txt",
                             true   /* binary_file */
                          );
    printRangeOfVariable( model3D, "element number" );
    printRangeOfVariable( model3D, "face number" );
    printRangeOfVariable( model3D, "node number" );
    
    // assigning the number of nodes to 'element variable' so that this array can be tested for correctness
    Region<3U>& model_domain = model3D.Region("Model");
    const csmp::Index evar_key = model3D.Database().StorageKey("element variable");
    for ( auto& it : model_domain.CellVector() )
      it->Store( evar_key, makeScalar(ANY,it->Nodes()) );
      
    model3D.InputPropertyValue( "element vector", makeVector(ANY,ANY,ANY,1.,2.,3.) );
  
    UG4_UGX_FileExport<3U> ug4_exporter3D( model3D );
    ug4_exporter3D.Write_UGX_FileASCII( model3D, model3D.Name() );

    cout << endl << endl << "Run() finished." << endl;
}





#if 0
/**
     ANLOR bubble migration example
*/
void Experimental_Example::Run() {

    string example_name = "capillary_seal";
    string model_name = "hemisphere_csmp_2_native";
    string variable_file = "hemisphere_csmp_2_native_variables.txt";
    string config_file = model_name;
    //create of directory with current example name, go into this directory, and copy input files into it.
    CreateWorkingDirectoryAndCopyInputModelFiles(example_name, model_name, variable_file, config_file);

//?
    class ScopedFunctor {
    public:
        ~ScopedFunctor() {
            fs::current_path("../../example_inputs/");
        };
    };
    ScopedFunctor scope;

    //if( !ImportModelAndRunChecks( model_name, variable_file ) ) {
    //    return;
    //}

// loads CSMP native model
    Model<3U>  model( model_name, variable_file );
    
    // OUTPUTS MODEL TO UG (name will be the model name)
    UG4_UGX_FileExport<3U> ug4_exporter( model );
    ug4_exporter.Write_UGX_FileASCII( model, "test_model" );
    
    RunSimulation( model );

    cout << endl << endl << "Run finished" << endl;
} // end Run



bool Experimental_Example::ImportModelAndRunChecks(
    const string& model_name, const string& variables_file
) const {
  
  // Build Model using variables file generated by SKUA
  Model<3U>  model( model_name, variables_file );

  if( !RunChecks( model ) ) {
      return false;
  }

  cout << endl << "Checks successful!" << endl;
  return true;
} // end ImportModelAndRunChecks



bool Experimental_Example::RunChecks( const Model<3U>& model ) const {

  // Checks...
  printModelDimensions( model, true );

  size_t nb_things = model.Regions();
  if( nb_things == 0 ) {
      cout << endl << "Model does not have any Regions" << endl;
  } else {
      cout << endl << "Model has " << nb_things << " Regions:" << endl;
      for(
          std::map<std::string,Region<3U> >::const_iterator
          it  = model.UniqueRegionsBegin();
          it != model.UniqueRegionsEnd();
          it++
      ) {
          cout <<"\n\t"<< (*it).first;
      }
  }

  nb_things = model.Boundaries();
  if( nb_things == 0 ) {
      cout << endl << "Model does not have any Boundaries" << endl;
  } else {
      cout << endl << "Model has " << nb_things << " Boundaries:" << endl;
      for(
          BoundaryInterface<3U,Boundary>::boundaryConstIterator
          it  = model.BoundariesBegin();
          it != model.BoundariesEnd();
          it++
      ) {
          (*it).second.Out();
      }
  }
  
  nb_things = model.SplitBoundaries();
  if( nb_things == 0 ) {
      cout << endl << "Model does not have any SplitBoundaries" << endl;
  } else {
      cout << endl << "Model has " << nb_things << " SplitBoundaries:" << endl;
      for(
          SplitBoundaryInterface<3U,SplitBoundary>::splitBoundaryConstIterator
          it  = model.SplitBoundariesBegin();
          it != model.SplitBoundariesEnd();
          it++
      ) {
          (*it).second.Out();
      }
  }

  // Testing NodeManifolds
  nb_things = model.Mesh().NodeManifolds();
  bool result = true;
  if( nb_things == 0 ) {
      cout << "Model does not have any NodeManifolds" << endl;
  } else {
      cout << "Model has " << nb_things << " NodeManifolds:" << endl;

      uint32_t matched(0), unmatched(0);
      for(auto mit = model.Mesh().NodeManifoldsBegin(); mit!=model.Mesh().NodeManifoldsEnd(); mit++) {
          const NodeManifold<3U>& md = *mit;
          const uint32_t branches = md.Branches();
          for(auto n{0U}; n<branches; n++) {
              const Node<3U>* node = md.N(n);
              const NodeManifold<3U>* retrieved_md = node->Manifold();
              if(retrieved_md == &md) {
                  matched++;
                  continue;
              }
              cout << "Node " << node->Idx() << " at coordinates " << node->Coordinate() << " is inconsistent" << endl;
              cout<<"  retrieved_md size = "<<retrieved_md->Branches();
              retrieved_md->Out();
              cout<<endl;
              cout<<"  original_md size = "<<md.Branches();
              md.Out();
              cout<<endl<<endl;
              unmatched++;
        }
      }
      if( unmatched != 0 ) {
          cout << endl << unmatched <<" nodes do not retrieve the manifolds "
              "that contain them" << endl;
          result = false;
      } else {
          cout << "All manifolds are consistent" << endl;
      }
  }

  cout << "RunChecks: That's it..." << endl;
  return result;
} // end RunChecks




/**
    Performs either discrete-event simulation (DES) or time-driven simulation (TDS) of 2-phase flow through a porous medium.
*/
<<<<<<< HEAD
void Experimental_Example::Run()

  {
      // illustrating Boundary creation for faults
      // -----------------------------------------
      const string  input_file("fault_boundary_test");
      const uint32_t dim(3);
=======
void Experimental_Example::RunSimulation( Model<3U>& model ) const {
>>>>>>> f1171f99af48ad8765b86abd2ccd22558a98877e

    string split_boundary_name = "outer_split_boundary";
    if( !model.ContainsSplitBoundary( split_boundary_name ) ) {
        cout << "Error: split boundary " << split_boundary_name << " not found." << endl;
        return;
    }
    split_boundary_name = "inner_split_boundary";
    if( !model.ContainsSplitBoundary( split_boundary_name ) ) {
        cout << "Error: split boundary " << split_boundary_name << " not found." << endl;
        return;
    }

    // 1. MODEL CONFIGURATION
    // ----------------------
    // simulation settings
    Standard_IO_Handler stdio;
    VTU_Interface<3U>  vtu(model);

    bool  DES = stdio.YesNo("Do you want to solve the transport equation with DES? (y=DES, n=TDS)"); 
    double Courant_multiplier, PEP_parameter;
    cout <<"\nEnter CFL multiplier (suggested value: 0.3 for TDS, 0.1 for DES) and PEP parameter (suggested value: 0.1)" << endl;
    cin >> Courant_multiplier >> PEP_parameter;  
     
    const bool  with_gravity_forces = stdio.YesNo("Do you want to include gravity effect (y/n)?");
    const bool  with_capillary_spreading = stdio.YesNo("Do you want to include capillary effect (y/n)?");

    // give the model dimensions
    printModelDimensions( model, true );

    // Configure the simulation from a file
    /*
    bool ConfigureFromFile( Model<dim>&, const char* fname,
                            bool block1,            ///< region name from parameter range
                            bool block2,            ///< default property values
                            bool block3,            ///< regional property values
                            bool block4,            ///< boundary conditions for box-shaped model
                            bool block5,            ///< essential conditions for regions
                            bool block6 = false );  ///< boundary conditions for arbitrary-shaped model
    */
    InputDataManager<3U>  model_configuration;
    string config_file = model.Name();
    model_configuration.ConfigureFromFile( model, config_file.c_str(),
                                           false,  ///< region name from parameter range
                                           true,   ///< default property values
                                           true,   ///< regional property values
                                           false,  ///< boundary conditions for box-shaped model
                                           true,   ///< essential conditions for regions
                                           true );  ///< boundary conditions for arbitrary-shaped model

    // check whether tensor k is used
    bool with_tensor_permeability(false);
    if(model.Database().IsDefined("tensor permeability")) {
      TensorVariable<3U> kk;
      static Index  key_kk(model.Database().StorageKey("tensor permeability"));
      model.Region("Model").E(0U)->Read(key_kk, kk);
      if(!isnan(kk(0U,0U))) with_tensor_permeability = true;
    }
    if(with_tensor_permeability) cout<<"\ntensor permeability is in use"<<endl;
    else  cout<<"\nscalar permeability is in use"<<endl;

 
    // 2. PROVISIONS FOR LOWER-DIMENSIONAL REPRESENTATION OF SAND HORIZONS
    // -------------------------------------------------------------------
    //  {
    //     computeGravityDipVectors( model );
    //     if ( model.ContainsRegion("SAND") ) {
    //          const string dim_minus1_region("SAND");
    //          // not needed in current approach: 27/8/22
    //          //computeFV_Diameter_Normal_VerticalExtent( model, dim_minus1_region );
    //          //computeSpillPointSaturation( model, dim_minus1_region );
    //          //lowerDimensionalLayerDiagnostics( model, dim_minus1_region );
    //       }
    //     else cout <<"\nRunSimulation: no provisions made for lower-dimensional sandbodies."<< endl;
    //  }
    //const double bcp{2.}, swr{0.15}, snr{0.};
    //SandPropertiesFor_VE_Model  sandProperties( model, "SAND", bcp, swr, snr );


    // 3. RELATIVE PERMEABILITY & CAPILLARY PRESSURE MODEL
    // ---------------------------------------------------
    // flow functions (Brooks Corey)
    FlowFunctionsModule1<3U> flowfunctions(model.Database(), model.Read( model.Database().StorageKey("acceleration gravity") ));

    // solve static pressure before split boundaries are created
    Compute2PhaseFlowProperties(model, flowfunctions, with_gravity_forces, with_tensor_permeability);
    
    //if constexpr (dim == 2U )
    //  if ( sandProperties.HasLineElementRepresentation() )
    //     sandProperties.Compute2PhaseFlowPropertiesForSandLayer(model);
    
    ComputeSteadyStatePressure(model, with_gravity_forces, with_tensor_permeability);
    printRangeOfVariable( model, stdio, "fluid pressure" );
    vtu.OutputDataToVTU( "steady_state_pressure", "fluid pressure", "Model", 0 );

    // initialise split boundary
    size_t n_split_boundaries = model.SplitBoundaries();
    cout<<"number of splitboundafies = "<<n_split_boundaries<<endl;
    model.SplitBoundariesOut();
    if( n_split_boundaries == 0 ) {
        return;
    }

    // using 'nodal variable' to visualise which nodes are manifolds
    Region<3U>  model_domain = model.Region("Model");
    if(!model.Database().IsDefined("number of collocated nodes")) {
        model.CreateProperty( "number of collocated nodes", "n_nd", "uint", SCALAR, NODE, 1, 0, 100);;
    }
    model.Region("Model").InputPropertyValue( "number of collocated nodes", makeScalar(PLAIN,1), COMPLETE);
    const csmp::Index nodes_key = model.Database().StorageKey("number of collocated nodes");
    for ( auto nit=model_domain.NodesBegin(); nit!=model_domain.NodesEnd(); ++nit ) {
        if ( (*nit)->IsManifold() ) {
            // there should be as many branches as materials come together
            (*nit)->Store( nodes_key, makeScalar(PLAIN,(*nit)->Manifold()->Branches()));
        }
    }

    vtu.OmitZeroInFileName( true );
    vtu.OutputDataToVTU( "Collocated_nodes", "number of collocated nodes", "Model", static_cast<int>(0) );

    // set entry pressure on split boundary
    //split_boundary "entry pressure" 194.133

    // end initialise split boundaries

    //sort manifold nodes by entry pressure
    const csmp::Index pd_key = model.Database().StorageKey("entry pressure");
    for(auto mit = model.Mesh().NodeManifoldsBegin();mit!=model.Mesh().NodeManifoldsEnd();mit++)
        (*mit).SortByVariableValue(pd_key);

    //check parent elements of duplicated nodes
    if(!model.Database().IsDefined("parent element id")) {
        model.CreateProperty( "parent element id", "eID", "uint", SCALAR, ELEMENT, 1, 0 ,1.00E+08);
        model.Region("Model").InputPropertyValue( "parent element id", makeScalar(PLAIN,0), COMPLETE);
    }
    const csmp::INDEX<SCALAR,ELEMENT> key_parent = csmp::INDEX<SCALAR,ELEMENT>( model.Database().StorageKey("parent element id") );
    size_t parent_id(1);
    for(auto mit = model.Mesh().NodeManifoldsBegin();mit!=model.Mesh().NodeManifoldsEnd();mit++) {
        parent_id = 1;
        auto md = (*mit);
        auto master_node = md.N(0);
        for( auto e{0U}; e < master_node->Parents(); e++)
          master_node->Parent(e)->Store(key_parent, makeScalar(PLAIN, parent_id));

        for(size_t n=1;n<md.Branches();n++) {
          auto slave_node = md.N(n);
          parent_id++;
          for( auto e{0U}; e < slave_node->Parents(); e++)
            slave_node->Parent(e)->Store(key_parent, makeScalar(PLAIN, parent_id));
        }
    }
    vtu.OutputDataToVTU( "parent_elment_id", "parent element id",  "Model", 0 );


    // construct DES transport
    auto DEStransport = new DES2PhaseSlightlyCompressibleTransport<3U,FlowFunctionsModule1> (model,
                                                                                             "Model",
                                                                                             with_gravity_forces,
                                                                                             with_capillary_spreading,
                                                                                             Courant_multiplier,
                                                                                             PEP_parameter,
                                                                                             1., //relaxation factor
                                                                                             with_tensor_permeability, //tensor k
                                                                                             false, //2nd order in space
                                                                                             flowfunctions);

    // solve static pressure after split boundaries are created
    /*ComputeSteadyStatePressure(model, with_gravity_forces, with_tensor_permeability);
    printRangeOfVariable(model, stdio, "fluid pressure");
    vtu.OutputDataToVTU("steady_state_pressure_with_splitboundaries", "fluid pressure", "Model", 0);
    //vtu.OutputDataToVTU( "gravity_term", "gravity term", "Model", 0 );
    */

    // defining input properties
    std::list<string> input_properties;
    input_properties.emplace_back( "porosity" );
    input_properties.emplace_back( "fluid pressure" );
    input_properties.emplace_back( "permeability" );
    input_properties.emplace_back( "vertical permeability" );
    input_properties.emplace_back( "rocktype" );
    input_properties.emplace_back( "fluid volume source");
    input_properties.emplace_back( "nodal fluid volume source" );
    input_properties.emplace_back( "saturation carbonic phase" );
    input_properties.emplace_back( "saturation aqueous phase" );
    input_properties.emplace_back( "thickness" );
    input_properties.emplace_back( "entry pressure" );
    input_properties.emplace_back("pressure continuity status");
    input_properties.emplace_back("breakthrough status");
      
      // for model with lower dimensional representation of the sand horizons inside the split boundaries
      //input_properties.emplace_back("finite volume diameter");
      //input_properties.emplace_back("finite volume normal");
      //input_properties.emplace_back("finite volume vertical extent");
      //input_properties.emplace_back("finite volume diagnostics");


    //defining output properties
    std::list<string> output_properties;
    output_properties.emplace_back( "rocktype" );
    output_properties.emplace_back( "entry pressure" );
    output_properties.emplace_back( "fluid pressure" );
    output_properties.emplace_back( "saturation aqueous phase" );
    output_properties.emplace_back( "saturation carbonic phase" );
    output_properties.emplace_back( "update count" );
    if(!model.Database().IsDefined("out range value count"))
      model.CreateProperty( "out range value count", "orvc", "uint", SCALAR, NODE, 1, 0.00E+00 ,1.00E+10);
    model.Region("Model").InputPropertyValue( "out range value count", makeScalar(PLAIN,0), COMPLETE);
    output_properties.emplace_back( "out range value count" );
    output_properties.emplace_back("pressure continuity status");
    output_properties.emplace_back("breakthrough status");

    //write input properties
    std::string model_name ( model.Name() );
    std::string file_name = model_name + "-input_data";
    vtu.OutputDataToVTU( file_name, input_properties, "Model", 0 );

    // global time for simulated runtime
    double	model_time(0.);

    // setting up region monitor
    RegionMonitor<3U>  monitor( model, "saturation carbonic phase", "saturation carbonic phase",  true );
    monitor.ScalarPropertyRanges( model, model_time );
    monitor.ScalarPropertyIntegrals( model, model_time );
    monitor.Out( (model_name + "-monitored_CO2_saturation").c_str() );

    RegionMonitor<3U>  monitor2( model, "fluid pressure", "fluid pressure",  true );
    monitor2.ScalarPropertyRanges( model, model_time );
    monitor2.ScalarPropertyIntegrals( model, model_time );
    monitor2.Out( (model_name + "-monitored_pressure").c_str() );

    // setting up time parameters (HARDWIRED PRESSURE STEPS!)
    const double day{86400.}; // year{ 86400. * 365. };
    //double max_time (5. * year);      // run for # years
    double max_time (60. * day);      // run for 60 days
    //double time_increment(3600.); // timestep
    double time_increment(0.5 * day); // timestep
    //double save_interval = 5. * day;  // save every # days
     double save_interval = 3. * day;  // save every # days

    size_t time;
    double end_time = model_time + max_time;
    double next_save_time = model_time + save_interval;

    // number of threads for use in transport
    size_t n_threads=1;

    // starting simulation
    if (DES) cerr <<"\n\nmain: Starting DES simulation: "<<endl;
    else cerr <<"\n\nmain: Starting TDS simulation: "<<endl;

    long output_count(1);
    while ( model_time < end_time )
    {
      Compute2PhaseFlowProperties( model, flowfunctions, with_gravity_forces, with_tensor_permeability );

      // setting multiphase flow properties for dim-1 elements representing sand
      //if constexpr ( dim == 2U )
      //  if ( sandProperties.HasLineElementRepresentation() )
      //    sandProperties.Compute2PhaseFlowPropertiesForSandLayer( model );

      // solve pressure
      ComputeSteadyStatePressure(model, with_gravity_forces, with_tensor_permeability);

      //transport steps
      cout << std::fixed << std::setprecision(3) <<"\nTransporting from "<<model_time<<" sec to "<<model_time+time_increment<<" sec, dt = "<<time_increment<<" sec\n"<<endl;
      if (DES) DEStransport->AdvectVariable_DES( time_increment, model_time+time_increment, n_threads );
      else DEStransport->AdvectVariable_TDS( time_increment, n_threads );

      //increment time
      model_time += time_increment;

      //output variables
      if ( model_time >= next_save_time ) {
          const csmp::Index p_key = model.Database().StorageKey("fluid pressure");
          const csmp::Index p_status_key = model.Database().StorageKey("pressure continuity status");
          for (auto mit = model.Mesh().NodeManifoldsBegin(); mit != model.Mesh().NodeManifoldsEnd(); mit++) {
            const auto n_branches{(*mit).Branches()};
            for (auto n{1U}; n < n_branches; n++) {
              auto node = (*mit).N(n);
              if (node->Status(p_key) == ROBIN) node->Store(p_status_key, makeScalar(PLAIN, 1));
              else node->Store(p_status_key, makeScalar(PLAIN, 0));
            }
        }

        time = static_cast<long>(model_time/day);
        if (DES) {
          std::string runtime_file_name(model_name + "-DES_runtime_output");
          vtu.OutputDataToVTU(runtime_file_name, output_properties, "Model", time);
          vtu.OutputDataToVTU(runtime_file_name, output_properties, "SAND", time);
        } else {
          std::string runtime_file_name(model_name + "-TDS_runtime_output");
          vtu.OutputDataToVTU(runtime_file_name, output_properties, "Model", output_count);
        }

        monitor.ScalarPropertyRanges(model, model_time);
        monitor.ScalarPropertyIntegrals(model, model_time);
        monitor.Out((model_name + "-monitored_CO2_saturation").c_str());
        monitor2.ScalarPropertyRanges(model, model_time);
        monitor2.ScalarPropertyIntegrals(model, model_time);
        monitor2.Out((model_name + "-monitored_pressure").c_str());

        next_save_time += save_interval;
        output_count++;

        model.Region("Model").InputPropertyValue("update count", makeScalar(PLAIN, 0), COMPLETE);

      }

      // runtime info
      cout <<"\n\nmain: RUNTIME (SEC): "<< model_time << endl << endl;

    }

    delete DEStransport;
    // terminate
    cerr << "\nmain: That's it..."<< endl;

  } // RunSimulation()

//  template void DES2PhaseFlowWithSplitBoundary_Example::RunSimulation(Model<2U>& model);
//  template void DES2PhaseFlowWithSplitBoundary_Example::RunSimulation(Model<3U>& model);





  template<class FLOW_FUNCTIONS>
  void Experimental_Example::Compute2PhaseFlowProperties(
      Model<3U>& mdl, FLOW_FUNCTIONS& flowfunctions, bool with_gravity, bool with_tensor_k
  ) const {
    const size_t v( 1u );

    auto mref = mdl.Region("Model");
    // keys to properties
    static Index  sw_key(mdl.Database().StorageKey("saturation aqueous phase"));
    static Index  snw_key(mdl.Database().StorageKey("saturation carbonic phase"));
    static Index  thi_key(mdl.Database().StorageKey("thickness"));
    static Index  rhow_key(mdl.Database().StorageKey("density aqueous phase"));
    static Index  rhon_key(mdl.Database().StorageKey("density carbonic phase"));
    static Index  gt_key(mdl.Database().StorageKey("gravity term"));
    static Index  dip_key(mdl.Database().StorageKey("dip vector"));
    static Index  cw_key(mdl.Database().StorageKey("compressibility aqueous phase"));
    static Index  cn_key(mdl.Database().StorageKey("compressibility carbonic phase"));
    static Index  cR_key(mdl.Database().StorageKey("compressibility rock"));
    static Index  ct_key(mdl.Database().StorageKey("total system compressibility"));
    static Index  phi_key(mdl.Database().StorageKey("porosity"));
    static Index  swr_key(mdl.Database().StorageKey("residual saturation aqueous phase"));
    static Index  snr_key(mdl.Database().StorageKey("residual saturation carbonic phase"));
    static Index  muw_key(mdl.Database().StorageKey("viscosity aqueous phase"));
    static Index  mun_key(mdl.Database().StorageKey("viscosity carbonic phase"));

    const double gravity_acceleration = mdl.Read( mdl.Database().StorageKey("acceleration gravity") );

    // 1. Computing the saturation of water = 1 - So
    for ( auto nit = mref.NodesBegin(); nit != mref.NodesEnd(); nit++ )
    {
      if( (*nit)->Status(sw_key) != DIRICH ) {
        double sw = 1. - (*nit)->Read(snw_key);
        (*nit)->Store( sw_key, makeScalar((*nit)->Status(sw_key),sw));
      }
    }

    // 2. Computing the multiphase flow properties
    if(!with_tensor_k) { //scalar k
      static Index  k_key(mdl.Database().StorageKey("permeability"));
      static Index  mobt_key(mdl.Database().StorageKey("total mobility permeability product"));
      //vector<Element<dim>* >::const_iterator eit;
      for ( auto eit = mref.CellsBegin(); eit!= mref.CellsEnd(); eit++ ) {
        Element<3U>* eptr = *eit;
        const size_t nodes(eptr->Nodes());
        vector<double> IPOL;
        eptr->N_AtBaryCenter( IPOL );
        double ipol_sum(0.), e_sw(0.), e_sn(0.), e_muw (0.), e_mun(0.), e_rhow(0.), e_rhon(0.), e_cw(0.), e_cn(0.);
        for ( auto i{0U}; i<nodes; ++i ) {
          e_sw += IPOL[i] * eptr->N(i)->Read( sw_key );
          e_muw += IPOL[i] * eptr->N(i)->Read( muw_key );
          e_mun += IPOL[i] * eptr->N(i)->Read( mun_key );
          e_rhow += IPOL[i] * eptr->N(i)->Read( rhow_key );
          e_rhon += IPOL[i] * eptr->N(i)->Read( rhon_key );
          e_cw += IPOL[i] * eptr->N(i)->Read( cw_key );
          e_cn += IPOL[i] * eptr->N(i)->Read( cn_key );
          ipol_sum += IPOL[i];
        }
        if(ipol_sum > 0.) {
          e_sw *= 1. / ipol_sum;
          e_muw *= 1. / ipol_sum;
          e_mun *= 1. / ipol_sum;
          e_rhow *= 1. / ipol_sum;
          e_rhon *= 1. / ipol_sum;
          e_cw *= 1. / ipol_sum;
          e_cn *= 1. / ipol_sum;
        }
        e_sn = 1.0 - e_sw;

        const double e_swr = eptr->Read(swr_key);
        const double e_snr = eptr->Read(snr_key);
        double thickness = eptr->Read(thi_key);
        if(isnan(thickness)) thickness = 1.;
        double e_k = eptr->Read(k_key);

        //compute element mobilities
        double e_lw(0.), e_ln(0.), e_lt(0.);
        if(e_sw > e_swr) e_lw = e_k * thickness * flowfunctions.krw_at(eptr, e_sw) / e_muw; //krw/muw
        if(e_sn > e_snr) e_ln = e_k * thickness * flowfunctions.krn_at(eptr, e_sw) / e_mun; //krn/mun
        e_lt = e_lw + e_ln;
        eptr->Store( mobt_key, makeScalar(PLAIN, e_lt) );

        //gravity term
        if ( with_gravity ) {
          VectorVariable<3U> gravity;
          eptr->Read( dip_key, gravity );
          if(isnan(gravity(1))) { //dip vector has not been initialised
            gravity(0u) = 0.; gravity(1u) = -1.;
            if(3U==3u) gravity(2u) = 0.;
            eptr->Store( dip_key, gravity );
          }
          double gravity_w(0.), gravity_n(0.);
          if(e_sw > e_swr) gravity_w = gravity_acceleration * e_k * thickness * flowfunctions.krw_at(eptr, e_sw) / e_muw * (e_rhow);
          if(e_sn > e_snr) gravity_n = gravity_acceleration * e_k * thickness * flowfunctions.krn_at(eptr, e_sw) / e_mun * (e_rhon);
          gravity_w *= fabs(gravity[v]);
          gravity_n *= fabs(gravity[v]);
          gravity *= (gravity_w + gravity_n);
          (*eit)->Store( gt_key, gravity );
        }

        
        // total system compressibility (weighted average approach)
        //const double e_cR = eptr->Read( cR_key );
        //double e_phi = eptr->Read(phi_key); //porosity
        //double e_ct = (1. - e_phi) * e_cR + e_phi * (e_sw * e_cw + (1.-e_sw) * e_cn);
        //e_ct *= thickness;
        //eptr->Store( ct_key, makeScalar(PLAIN, e_ct)  );

      }

    } else { //tensor k
      static Index  kk_key(mdl.Database().StorageKey("tensor permeability"));
      static Index  LT_key(mdl.Database().StorageKey("tensor total mobility permeability product"));
      TensorVariable<3U> kk;
      VectorVariable<3U> kV;
      //vector<Element<dim>* >::const_iterator eit;
      for ( auto eit = mref.CellsBegin(); eit!= mref.CellsEnd(); eit++ ) {
        Element<3U>* eptr = *eit;
        const size_t nodes(eptr->Nodes());
        vector<double> IPOL;
        eptr->N_AtBaryCenter( IPOL );
        double ipol_sum(0.), e_sw(0.), e_sn(0.), e_muw (0.), e_mun(0.), e_rhow(0.), e_rhon(0.), e_cw(0.), e_cn(0.);
        for ( auto i{0U}; i<nodes; ++i ) {
          e_sw += IPOL[i] * eptr->N(i)->Read( sw_key );
          e_muw += IPOL[i] * eptr->N(i)->Read( muw_key );
          e_mun += IPOL[i] * eptr->N(i)->Read( mun_key );
          e_rhow += IPOL[i] * eptr->N(i)->Read( rhow_key );
          e_rhon += IPOL[i] * eptr->N(i)->Read( rhon_key );
          e_cw += IPOL[i] * eptr->N(i)->Read( cw_key );
          e_cn += IPOL[i] * eptr->N(i)->Read( cn_key );
          ipol_sum += IPOL[i];
        }
        if(ipol_sum > 0.) {
          e_sw *= 1. / ipol_sum;
          e_muw *= 1. / ipol_sum;
          e_mun *= 1. / ipol_sum;
          e_rhow *= 1. / ipol_sum;
          e_rhon *= 1. / ipol_sum;
          e_cw *= 1. / ipol_sum;
          e_cn *= 1. / ipol_sum;
        }
        e_sn = 1.0 - e_sw;

        const double e_swr = eptr->Read(swr_key);
        const double e_snr = eptr->Read(snr_key);
        double thickness = eptr->Read(thi_key); //thickness
        if(isnan(thickness)) thickness = 1.;

        //compute element mobilities
        double e_lw(0.), e_ln(0.), e_lt(0.);
        if(e_sw > e_swr) e_lw = thickness * flowfunctions.krw_at(eptr, e_sw) / e_muw; //krw/muw
        if(e_sn > e_snr) e_ln = thickness * flowfunctions.krn_at(eptr, e_sw) / e_mun; //krn/mun
        e_lt = e_lw + e_ln;
        eptr->Read( kk_key, kk );
        kk *= e_lt;
        eptr->Store( LT_key, kk );

        //gravity term
        if ( with_gravity ) {
          VectorVariable<3U> gravity;
          eptr->Read( dip_key, gravity );
          if(isnan(gravity(1))) { //dip vector has not been initialised
            gravity(0u) = 0.; gravity(1u) = -1.;
            if(3U==3u) gravity(2u) = 0.;
            eptr->Store( dip_key, gravity );
          }
          // projection of tensor on the dip vector
          eptr->Read( kk_key, kk );
          kV = kk * gravity;
          double kV_magnitude = kV.Length();

          double gravity_w(0.), gravity_n(0.);
          if(e_sw > e_swr) gravity_w = gravity_acceleration * kV_magnitude * thickness * flowfunctions.krw_at(eptr, e_sw) / e_muw * (e_rhow);
          if(e_sn > e_snr) gravity_n = gravity_acceleration * kV_magnitude * thickness * flowfunctions.krn_at(eptr, e_sw) / e_mun * (e_rhon);
          gravity_w *= fabs(gravity[v]);
          gravity_n *= fabs(gravity[v]);
          gravity *= (gravity_w + gravity_n);
          (*eit)->Store( gt_key, gravity );
        }


        // total system compressibility (weighted average approach)
        //const double e_cR = eptr->Read( cR_key );
        //double e_phi = eptr->Read(phi_key); //porosity
        //double e_ct = (1. - e_phi) * e_cR + e_phi * (e_sw * e_cw + (1.-e_sw) * e_cn);
        //e_ct *= thickness;
        //eptr->Store( ct_key, makeScalar(PLAIN, e_ct)  );

      }
    }
  }

//template void DES2PhaseFlowWithSplitBoundary_Example::Compute2PhaseFlowProperties( Model<2U>&, FlowFunctionsModule1<2U>&, bool, bool );
template void Experimental_Example::Compute2PhaseFlowProperties( Model<3U>&, FlowFunctionsModule1<3U>&, bool, bool ) const;






  void Experimental_Example::ComputeSteadyStatePressure(
      Model<3U>& mdl, bool with_gravity, bool with_tensor_k
  ) const {
    bool verbose(false);

    if(verbose) {
      cout << "\nExperimental_Example::ComputeSteadyStatePressure: Input parameters: " << endl;
      printRangeOfVariable( mdl, "fluid pressure" );
      printRangeOfVariable( mdl, "total mobility permeability product" );
      printRangeOfVariable( mdl, "permeability" );
      printRangeOfVariable( mdl, "fluid volume source");
      printRangeOfVariable( mdl, "nodal fluid volume source" );
      if(with_gravity) printRangeOfVariable( mdl, "gravity term" );
    }

    std::string	conductance_operator;
    if(!with_tensor_k) conductance_operator = "total mobility permeability product";
    else conductance_operator = "tensor total mobility permeability product";

    CSMP_DEFAULT_LINEAR_SOLVER  solver;
    PDE_Integrator<3U,Element>  steady_pressure( solver );

    NumIntegral_dNT_op_dN_dV<3U>  conductance( mdl.Database(), conductance_operator.c_str(), "fluid pressure", "fluid pressure" );
    NumIntegral_NT_op_N_dV<3U>    elmt_volume_source( mdl.Database(), "fluid volume source", "fluid pressure" );
    //PointSource_rhsop<dim>          nodal_volume_source( mdl.Database(), "nodal fluid volume source", "fluid pressure" );
    steady_pressure.Add( &conductance );
    steady_pressure.Add( &elmt_volume_source );
    //steady_pressure.AddBoundaryIntegral( &influx ); //added
    steady_pressure.Verbose(verbose);

    NumIntegral_dNT_op_dV<3U>* gravity(nullptr);
    if ( with_gravity ) {
      gravity = new NumIntegral_dNT_op_dV<3U>( mdl.Database(), "gravity term", "fluid pressure" );
      steady_pressure.Add( gravity );
    }


    // iout
    //settings.ExplicitSecondary(true);
    //settings.Set_iout1( -1 );
    //settings.Set_iout2( -1 );
    //settings.Set_idmp( -1 );
    //settings.Set_mode_mess( -2 );

    if(verbose) {
      cout << "\n\n\nDES2PhaseFlowWithSplitBoundary_Example::ComputeSteadyStatePressure: ";
      cout << " Computing '" << "steady state fluid pressure" << "'" << endl;
    }

    mdl.Apply( steady_pressure );

    delete gravity;

  } // end SteadyStatePressure



/**
      loops over the finite volumes of the supplied lower-dimensional supplied region, computing the diameters of the finite volumes and their normals
      @todo the "dip vector" variables and thickness attributes need to be computed elsewhere
*/
/*template<uint32_t dim>
void computeFV_Diameter_Normal_VerticalExtent( Model<dim>& model, const std::string& region_name )
 {
    csmp::ErrorHandler& csmp_error( ErrorHandler::Instance() );
    if ( !model.ContainsRegion( region_name ) )
      csmp_error.Note( ERROR, "computeFV_Diameter_Normal_VerticalExtent",
                      "target region does not exist");
                      
    Region<dim>& subdomain = model.Region( region_name );
    // checking that the region is indeed lower dimensional
    if ( subdomain.SpatialDimensions().second != dim - 1 )
      csmp_error.Note( ERROR, "computeFV_Diameter_Normal_VerticalExtent",
                      "target region is not lower dimensional (dim region != dim-1)");
                      
    // computing FV diameter and FV-averaged normals to lower-dimensional finite volumes
    const csmp::Index diam_key = model.Database().StorageKey("finite volume diameter");
    const csmp::Index nrml_key = model.Database().StorageKey("finite volume normal");
    const csmp::Index vext_key = model.Database().StorageKey("finite volume vertical extent");
    
    for ( auto nit=subdomain.NodesBegin(); nit!=subdomain.NodesEnd(); ++nit ) {
         // FV diameter/vertical extent
         const auto FV_props = diameterAndVerticalExtentOfLowerDimensional_FV( (*nit) );
         assert( FV_props.first > 0. ); // diameter must be greater than zero
         (*nit)->Store( diam_key, makeScalar(PLAIN,FV_props.first) );
         (*nit)->Store( vext_key, makeScalar(PLAIN,FV_props.second) );
         // FV normal
         Point<dim> nrml;
         bool was_able_to_compute_normal = (*nit)->UnitNormal( nrml );
         assert( was_able_to_compute_normal );
         // if this normal is not upward pointing, it is flipped
         if constexpr ( dim == 3U ) if ( dotProduct( nrml, Point<3U>(0.,1.,0.) ) < 0. ) nrml *= -1.;
         if constexpr ( dim == 2U ) if ( dotProduct( nrml, Point<2U>(0.,1.) ) < 0. )    nrml *= -1.;
         (*nit)->Store( nrml_key, std::move( VectorVariable<dim>(nrml) ) );
      }
 
 } // end computeFV_Diameter_Normal_VerticalExtent


template void computeFV_Diameter_Normal_VerticalExtent( Model<3U>&, const string& );
template void computeFV_Diameter_Normal_VerticalExtent( Model<2U>&, const string& );

*/


/**
    compute gravity dip vectors borrowed from ACGS
*/
/*template<uint32_t dim>
void computeGravityDipVectors( Model<dim>& model )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  // 1. volumetric model region
  // --------------------------
  csmp::Region<dim>&   sgref = model.Region( "Model" );
  VectorVariable<dim>  gvt;
  Point<dim>           gravity_direction;
  if constexpr ( dim == 3U ) gravity_direction = { 0., -1., 0. };
  else                       gravity_direction = { 0., -1. };

  const csmp::Index key_dip  = model.Database().StorageKey("dip vector");
  const csmp::Index key_dipf = model.Database().StorageKey("face dip vector");

  for ( auto it = sgref.CellsBegin(); it != sgref.CellsEnd(); it++ )
  {
    // computing gravity direction vectors for lower dimensional elements
    // line elements
    if ( (*it)->IsLine() ) {
      Point<dim> line_vector( (*it)->N( 0 )->Coordinate() - (*it)->N( 1 )->Coordinate() );
      line_vector.NormalizeLengthTo( 1. );
      if ( line_vector[1] > 0. ) line_vector *= -1.;
      gvt( 0 ) = line_vector[0];
      gvt( 1 ) = line_vector[1];
      if constexpr (dim == 3U ) gvt( 2 ) = line_vector[2];
    }
    // surface elements
    else if ( (*it)->IsSurface() ) {
      Point<dim> line_vector1( (*it)->N( 0 )->Coordinate() - (*it)->N( 1 )->Coordinate() );
      Point<dim> line_vector2( (*it)->N( 1 )->Coordinate() - (*it)->N( 2 )->Coordinate() );

      // finding the normal to the surface element
      Point<dim> surface_normal( crossProduct( line_vector1, line_vector2 ) );

      // finding the gravity direction from the cross-product of the horizontal surface vector and the gravity direction
      Point<dim> gravity_projection_on_surface( crossProduct( surface_normal, crossProduct( surface_normal, gravity_direction ) ) );
      gravity_projection_on_surface.NormalizeLengthTo( 1. );

      if ( gravity_projection_on_surface[1] > 0. ) gravity_projection_on_surface *= -1.;

      gvt( 0 ) = gravity_projection_on_surface[0];
      gvt( 1 ) = gravity_projection_on_surface[1];
      if constexpr (dim == 3U ) gvt( 2 ) = gravity_projection_on_surface[2];
    }
    // volume elements
    else {
      gvt( 0 ) = 0.;
      gvt( 1 ) = -1.;
      if constexpr (dim == 3U ) gvt( 2 ) = 0.;
    }
    if(isnan(gvt( 0 )) || isnan(gvt( 1 )) || isnan(gvt( 2 ))) {
        gvt( 0 ) = 0.;
        gvt( 1 ) = -1.;
        if constexpr (dim == 3U ) gvt( 2 ) = 0.;
        cerr<<"gvt reset to [0, -1, 0]"<<endl;
    }
    (*it)->Store( key_dip, gvt );
  }

  // 2. for all model boundaries
  // ---------------------------
  for ( auto git = model.BoundariesBegin(); git != model.BoundariesEnd(); git++ )
  {
    //          cout <<"\n\tBoundary: "<< (*git).first.first << (*git).first.second;
    //          cout.flush();
    assert( (*git).second.Cells() > 0 );
    for ( auto it = (*git).second.CellsBegin(); it != (*git).second.CellsEnd(); it++ )
    {
      if ( (*it) != nullptr and (*it)->IsSurface() ) {
        // finding the normal to the surface element
        Point<dim> line_vector1( (*it)->N( 0 )->Coordinate() - (*it)->N( 1 )->Coordinate() );
        Point<dim> line_vector2( (*it)->N( 1 )->Coordinate() - (*it)->N( 2 )->Coordinate() );
        Point<dim> surface_normal( crossProduct( line_vector1, line_vector2 ) );

        // finding the gravity direction from the cross-product of the hor izontal surface vector and the gravity direction
        Point<dim> gravity_projection_on_surface( crossProduct( surface_normal, crossProduct( surface_normal, gravity_direction ) ) );
        gravity_projection_on_surface.NormalizeLengthTo( 1. );

        if ( gravity_projection_on_surface[1] > 0. ) gravity_projection_on_surface *= -1.;

        gvt( 0 ) = gravity_projection_on_surface[0];
        gvt( 1 ) = gravity_projection_on_surface[1];
        if constexpr (dim == 3U ) gvt( 2 ) = gravity_projection_on_surface[2];

        (*it)->Store( key_dipf, gvt );
      }
      // line elements
      else if ( (*it) != nullptr and (*it)->IsLine() ) {
        Point<dim> gravity_projection_on_line( (*it)->N( 0 )->Coordinate() - (*it)->N( 1 )->Coordinate() );
        // flip gravity vector if it is upward pointing
        if ( (*it)->N( 0 )->y() > (*it)->N( 1 )->y() ) gravity_projection_on_line *= -1.;
        gravity_projection_on_line.NormalizeLengthTo( 1. );

        gvt( 0 ) = gravity_projection_on_line[0];
        gvt( 1 ) = gravity_projection_on_line[1];
        if constexpr (dim == 3U ) gvt( 2 ) = gravity_projection_on_line[2];

        (*it)->Store( key_dipf, gvt );
      }
      else {
        if ( (*it) != nullptr ) (*it)->Out();
        csmp_error.Note( WARNING, "CO2_GeoSequestrationSimulator::ComputeGravityDipVectors:",
                        "could not compute boundary normal for (?line?) element." );
      }
    }
  }
  
 } // end gravityDipVectors

template void computeGravityDipVectors( Model<3U>& );
template void computeGravityDipVectors( Model<2U>& );
*/


/**
    Computes the CO2 saturation value above which there a wedge of CO2  reaches from the lowest point of the FV to the highest point.
    Method uses 'vertical extent' and 'finite volume diameter' in this calculation.
    
    The scalar node variable 'spill-point saturation' is computed from the ratio of the (approximate) volume of CO2 required to create CO2 pool
    that spans its diameter in the dip direction and its pre-computed (exact total) 'finite volume'.
*/
/*template<uint32_t dim>
void computeSpillPointSaturation( Model<dim>& model, const std::string& region_name )
 {
    csmp::ErrorHandler& csmp_error( ErrorHandler::Instance() );
    if ( !model.ContainsRegion( region_name ) )
      csmp_error.Note( ERROR, "computeSpillPointSaturation",
                      "target region does not exist");
                      
    Region<dim>& subdomain = model.Region( region_name );
    // checking that the region is indeed lower dimensional
    if ( subdomain.SpatialDimensions().second != dim - 1 )
      csmp_error.Note( ERROR, "computeSpillPointSaturation",
                      "target region is not lower dimensional (dim region != dim-1)");
                      
    // computing FV diameter and FV-averaged normals to lower-dimensional finite volumes
    // input variables
    const csmp::Index fvol_key = model.Database().StorageKey("finite volume");
    const csmp::Index thi_key  = model.Database().StorageKey("thickness");
    const csmp::Index diam_key = model.Database().StorageKey("finite volume diameter");
    const csmp::Index vext_key = model.Database().StorageKey("finite volume vertical extent");
    const csmp::Index dip_key  = model.Database().StorageKey("dip vector");
    // output variables
    const csmp::Index spill_key = model.Database().StorageKey("spill-point saturation");
    
    // computation: finding the ratio between the approximate spill-point volume and the actual finite volume
    // (porosity is ignored in this geometric analysis)
    // (the approximate diameter of the finite volume is used 2R), see notes in MS Word.
    // ---------------------------------------------------------------------------------
    // V_wedge = 1/2 Pi r^2 (h1 + h2) -> h1=0, h2 = 'vertical extent'
    VectorVariable<dim> dip_vec, flat_vec;
    
    for ( auto nit=subdomain.NodesBegin(); nit!=subdomain.NodesEnd(); ++nit ) {
         const double h2 = (*nit)->Read( vext_key );
         // if the FV lies in the horizontal plane
         if ( fabs(h2) <= numeric_limits<double>::epsilon() * 100. ) {
              (*nit)->Store( spill_key, makeScalar(ANY,0.) );
              continue;
           }
         // if the FV is tilted
         const double r          = (*nit)->Read( diam_key ) / 2.;
         const double V_wedge    = 0.5 * PI * (r*r) * h2;
         const double sCO2_spill = V_wedge / (*nit)->Read( fvol_key );
         if ( sCO2_spill > 1. ) {
              cerr <<"\nNode "<< (*nit)->Idx() <<": sCO2_spill: "<< sCO2_spill <<", location: "<< (*nit)->Coordinate();
              cerr <<", dip of FV: ";
              // computing average dip of the lower-dim elmts making up the FV sectors
              double   avg_dip{ 0. };
              uint32_t n_surf_elmts{ 0u };
              for ( auto i{0U}; i<(*nit)->Parents(); i++ ) {
                   if constexpr ( dim == 3U ) {
                        if ( (*nit)->Parent(i)->IsSurface() ) {
                             (*nit)->Parent(i)->Read( dip_key, dip_vec );
                             // calculating dip
                             flat_vec    = dip_vec;
                             flat_vec(1) = 0.;
                             double dip = dip_vec.AngleTo( flat_vec );
                             avg_dip += dip;
                             n_surf_elmts++;
                          }
                     }
                   else if constexpr ( dim == 2U ) {
                        if ( (*nit)->Parent(i)->IsLine() ) {
                             (*nit)->Parent(i)->Read( dip_key, dip_vec );
                             // calculating dip
                             flat_vec    = dip_vec;
                             flat_vec(1) = 0.;
                             double dip = dip_vec.AngleTo( flat_vec );
                             avg_dip += dip;
                             n_surf_elmts++;
                          }
                     }
                   static_assert( dim != 1U, "computeSpillPointSaturation: method does not work in 1D." );
                   avg_dip /= static_cast<double>(n_surf_elmts);
                }
              cerr << avg_dip;
              csmp_error.Note( WARNING, "computeSpillPointSaturation:", "computed value greater than 1.");
           }
         else (*nit)->Store( spill_key, makeScalar(ANY,sCO2_spill) );
      }
 
  } // end computeSpillPointSaturation
*/


#endif



} // csmp
