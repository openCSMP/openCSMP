//
//  VariableManagement_Example.cpp
//  Open CSMP++
//
//  Created by Stephan Matthai on 27/9/2022.
//

#include "VariableManagement_Example.h"
#include "vsetMakers.h"
#include "ModelTopology.h"
#include "Model.h"

using namespace std;

namespace csmp {

void VariableManagement_Example::Specifications()
  {
    SetTitle( "Management of variables in CSMP native binary files" );
    SetDifficulty( 1 );
    SetCategory( "Software Functionality" );
    AddAuthor( "SKM" );
    AddDescription( "How to selectively read and store variables in existing CSMP binary" );
    AddDescription( "source in: VariableManagement_Example.cpp" );
    AddRequirement( "file set: 'LeftRight'");
    AddRequirement( "CSMP-1phase-variables.txt");
  }

// TODO: also use to demonstrate how to create header file


void VariableManagement_Example::Run()
 {
    // 1. building a complete model with boundaries and the variables 'node number' and 'element number'
    //   (NOTE: only those variables for which there is a previous definition in the PropertyDatabase
    //    can be read because the specification of variables in the VSet is incomplete)
    {
      VSet<2U>      vset;
      ModelTopology topo = test_Create_MeshPatchWithLineElements_VSet( vset );
      const bool convert_bdry_elmts_to_faces{true};
      Model<2U>     original_model( topo, vset, "Minimum-variables.txt", convert_bdry_elmts_to_faces );
      // saving model to CSMP binary
      original_model.OutputToBinaryFile("VariableManagement_Example_model");
    }
   
   // 2. building the model from binary file with only the single variable 'node number' in it
   //   (also see console output from variables database when model is built)
   {
      set<string> variable_subset{"node number"};
      Model<2U> model( string("VariableManagement_Example_model"), variable_subset );
      model.Database().FlushToScreen();
   }
   
   // 3. building the model from binary file as is, then adding two variables at runtime and saving it
   {
      Model<2U> model( string("VariableManagement_Example_model") );
      model.Database().FlushToScreen();
      // creating new variables
      model.CreateProperty("my element variable", "mev", "none" );
      // full specifications
      /*
      Model<dim>::CreateProperty( const char* new_prop, const char* unit,
                               VARIABLE_TYPE type = SCALAR, PLACEMENT place = NODE,
                               uint32_t vsize = 1, double vmin = -1.0e+30, double vmax = 1.0e+30,
                               std::string usage = "???" );
      */
      model.CreateProperty("integration-point variable", "ipv", "none", SCALAR, ELEMENT_INTEGRATION_POINT );
      model.Database().FlushToScreen();
      // saving model including the new variables to CSMP binary
      model.OutputToBinaryFile("VariableManagement_Example_model1");
   }
   
   // 4. building the model from binary file, but with a new variable set. Any non-null existing variables get included
   {
      Model<2U> model( string("VariableManagement_Example_model"), "CSMP-1phase-variables.txt" );
      model.Database().FlushToScreen();
      model.OutputToBinaryFile("VariableManagement_Example_model2");
   }
   
   // 5. Writing a header file with templatized indices for compile-time inclusion of variables of interest
   {
      Model<2U> model( string("VariableManagement_Example_model2"), "CSMP-1phase-variables.txt" );
      model.Database().FlushToScreen();
      const string header_file{"VariableManagement_Example_variable_set"}, variable_set_name{"ExampleVariableSet"};
      model.Database().WriteVariableSetToHeaderFile( header_file.c_str(), variable_set_name.c_str() );
   }

  // 6. TODO: show how to write a subset of variables back to CSMP binary file

 } // end run




} // end csmp
