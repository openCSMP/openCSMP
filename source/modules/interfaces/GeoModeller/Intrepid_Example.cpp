//
//  Experimental_Example.cpp
//  CSMP_API_library2014
//
//  Created by Stephan Matthai on 1/27/14.
//  Copyright (c) 2014 Stephan Matthai. All rights reserved.
//

#include "CSMP_definitions.h"
#include "Experimental_Example.h"
#include "CompareFloats.h"
#include "ANSYS_Model3D.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_NT_op_N_dV.h"
#include "PDE_Integrator.h"
#include "Face.h"
#include "NumIntegral_NT_op_N_dS.h"
#ifdef CSMP_WITH_SAMG_SOLVER
#include "SAMG_Solver.h"
#endif
#ifdef CSMP_WITH_MESCHACH
#include "Gauss_Solver.h"
#endif

#include "InputDataManager.h"
#include "VTK_Interface.h"
#include "VTU_Interface.h"
#include "IntrepidInterface.h"

// include any header files that you need here...

using namespace std;

namespace csmp {

void Experimental_Example::Specifications()
{
   SetTitle( "Experimental_Example" );
   SetDifficulty( 1 );
   SetCategory( "Software Functionality" );
   AddAuthor( "Stephan Matthai" );
   AddDescription( "source in: Experimental_Example.cpp" );
   AddDescription( "Empty example for the user to experiment with" );
   AddRequirement( "none" );
   AddRequirement( "no predefined model or variables file" );
}


/** 
    This is the documentation format that can be read by DOxygen automatically.
*/
void Experimental_Example::Run()
{
  std::string model_name("Mansfield_H8_NoOrphans.mesh");
  cout<< "\nPlease input the name of the model: ";
//  cin >> model_name;
 
  IntrepidInterface  geomodel;
  VSet<3U>           vset;
  ModelTopology      model_topology(true);
  geomodel.Read( model_name.c_str(), vset, model_topology );
  
  Model<3U>  model(model_topology, vset, "Geomodeller_test-variables.txt", false );

  printModelDimensions( model, true );
  
  cout<<"\ncurrent model boundaries:\n";
  for ( BoundaryInterface<3U,Boundary>::boundaryConstIterator it=model.BoundariesBegin(); it!=model.BoundariesEnd(); it++ )
    cout <<"\n\t"<< (*it).first;
  
  cout<<"\ncurrent model regions:\n";
  for ( map<string,Region<3U> >::const_iterator it=model.UniqueRegionsBegin(); it!=model.UniqueRegionsEnd(); it++ )
    cout <<"\n\t"<< (*it).first;
  
  InputDataManager<3U>().ConfigureFromFile( model, model_name.c_str(),
                                            false, true, true, true, true, true );

  printRangeOfVariable(model, "thermal conductivity");
  VTK_Interface<3U>  vtk_output;
  vtk_output.OutputRegionByRegionToVTK( model, model_name, "element volume", 0., true );
  
#ifdef CSMP_WITH_SAMG_SOLVER
  PDE_Integrator<3U,Region>  heat_conductor(new SAMG_Solver());
#else
  PDE_Integrator<3U,Region>  heat_conductor(new CSMP_DEFAULT_LINEAR_SOLVER());
#endif

  NumIntegral_dNT_op_dN_dV<3U> conductance (model.Database(), "thermal conductivity", "temperature", "temperature");
  
  NumIntegral_NT_op_N_dV<3U>  heat_source( model.Database(), "heat source", "temperature" );
  
//  NumIntegral_NT_op_N_dS<3U,Face >  basal_hfu( model.Database(), "basal heat flow", "temperature" );
//  basal_hfu.LumpedFormulation();
  
  heat_conductor.Add( &conductance );
  heat_conductor.Add( &heat_source );
  // heat_conductor.AddBoundaryIntegrals( &basal_hfu );
  
  Region<3U>& region(model.Region("Model"));
  
  heat_conductor.IntegrateOver( region, false );
  
  printRangeOfVariable(model, "temperature");

  list<string>  output_variables;
  output_variables.push_back("temperature");
  output_variables.push_back("thermal conductivity");
//  output_variables.push_back("basal heat flow");
  
  VTU_Interface<3U> vtu_output( model );
  vtu_output.OutputDataToVTU( model_name.c_str(), output_variables, "Model", static_cast<int>(0) );

  cout <<"\nExperimental_Example: That's it!\n";
  
} // end Run

} // csmp
