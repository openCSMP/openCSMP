// CSP Files
#include "SuperGroup.h"
#include "CSP_TRIANGLE_Interface.h"
#include "CSP_ICEM_Interface.h"
#include "CSP_VSetConverter.h"

using namespace std;

namespace csp
 {
    csp::ErrorHandler skm_err;
    csp_float         global_time;
    void              programInfo();
 }

using namespace csp;

// A simple interface to read in ICEM (2D or 3D) and Trianlge (2D) FE meshes 
// and convert them into a VSet file for partitioning and parallel simulations
// Verion 1.1: April 24, 2008 (SG)

int main()
{
  // Some info on program
  programInfo();

  // Set global variables
  global_time = 0.0;

  // Yes/No input
  Standard_IO_Handler  choice;

  // Strings for file input
  std::string physical_variables, input_mesh, output_vset;
  
  // Get the name of the input finite element mesh
  cout <<"\nPlease enter the full name of the finite element mesh"<< endl;
  cin  >> input_mesh;

  // Get the name for the output file
  cout <<"\nPlease enter the full name for the output VSet file";
  cout <<"\nto which the mesh will be saved (excluding *.vset extension)"<< endl;
  cin  >> output_vset;

  // Get the name of the CSMP physical variables text file
  cout <<"\nPlease enter the full name of the CSMP physical variables file (excluding *.txt extension)"<< endl;
  cin  >> physical_variables;
  physical_variables += ".txt";

  // Query if a 2D or 3D mesh should be read in 
  bool three_d_mesh = choice.RecordLogicalChoice("Is the input finite element mesh in 3D (n = 2D)?");

  // two-dimensional FE meshes
  if(!three_d_mesh) {
      // the SuperGroup
      SuperGroup<csp_float,2>* model;
     
      // Queary if an ICEM or Triangle mesh is used
      bool icem_mesh = choice.RecordLogicalChoice("Was the input finite element mesh generated with ICEM (n = Triangle mesher)?");

      if ( icem_mesh ) {
          // containers for reading in and storing the mesh
          VSet<csp_float,2>   mesh_container;
	      ICEM_Interface      mesh_interface_ICEM(true);
	      ModelTopology       mesh_topology(true);

	      // read mesh
	      mesh_interface_ICEM.ReadTetraSurfaceMeshBinary( input_mesh.c_str(), mesh_container, mesh_topology );

	      // selecting a subset of the total domains of the model
	      mesh_topology.ReduceToRegions( input_mesh.c_str() );
	      map<stl_index,stl_index>  old_and_new_elmtids;
	      mesh_topology.RenumberElementsConsecutively( old_and_new_elmtids );
	      mesh_container.ReduceTo( old_and_new_elmtids );

	      // create SuperGroup
	      model = new SuperGroup<csp_float,2>( mesh_topology, mesh_container, physical_variables.c_str() );

          // erase containers
	      mesh_container.Erase();
	      mesh_topology.Erase();
        }
      else {
          // containers for reading in and storing the mesh
          VSet<csp_float,2>          mesh_container;
	      VSetConverter<csp_float,2> converter;
	      TRIANGLE_Interface         mesh_interface;
 
		  // Query if isoparametric elements are used
		  bool isoparametric = choice.RecordLogicalChoice("Is the input finite element mesh isoparametric (n = analytical)?");

 	      // read mesh
	      mesh_interface.ReadTriangle2DMesh( input_mesh.c_str(), mesh_container, isoparametric );
 
 	      // create SuperGroup
	      model = new SuperGroup<csp_float,2>( mesh_container, physical_variables.c_str(), isoparametric );
         
          // erase containers
          mesh_container.Erase();

        }
     
      // save the model to a VSet file
      VSet<csp_float,2>  saved_model;
      model->OutputTo(saved_model);
      saved_model.OutputTo(output_vset.c_str(),global_time);
    }
    
  // three-dimensional FE meshes  
  else {
      // containers for reading in and storing the mesh
      VSet<csp_float,3>   mesh_container;
      ICEM_Interface      mesh_interface_ICEM(true);
      ModelTopology       mesh_topology(true);

      // read mesh
      mesh_interface_ICEM.ReadTetraMeshBinary( input_mesh.c_str(), mesh_container, mesh_topology );

      // selecting a subset of the total domains of the model
      mesh_topology.ReduceToRegions( input_mesh.c_str() );
      map<stl_index,stl_index>  old_and_new_elmtids;
      mesh_topology.RenumberElementsConsecutively( old_and_new_elmtids );
      mesh_container.ReduceTo( old_and_new_elmtids );

      // create SuperGroup
      SuperGroup<csp_float,3>* model;
      model = new SuperGroup<csp_float,3>( mesh_topology, mesh_container, physical_variables.c_str() );

      // erase containers
      mesh_container.Erase();
      mesh_topology.Erase();
      
      // save the model to a VSet file
      VSet<csp_float,3>   saved_model;
      model->OutputTo(saved_model);
      saved_model.OutputTo(output_vset.c_str(),global_time);
   }
  
  cout << "\nmain: Input file '" << input_mesh << "' successfully written to VSet '" << output_vset << "'";
  cout << "\nusing physical variables file '" <<  physical_variables << "'" << endl;
  cout << "\nmain: That's it..."<< endl;

  return 0;

} // end main

namespace csp {

void programInfo() {
  
  cout << "\n*******************************************************";
  cout << "\n*                                                     *";
  cout << "\n*  mesh_to_vset: A small program that translates 2D   *";
  cout << "\n*  and 3D finite element meshes into a CSMP VSet      *";
  cout << "\n*  file that can be used for partitioning the mesh    *";
  cout << "\n*  and parallel simulations.                          *";
  cout << "\n*                                                     *";
  cout << "\n*  Currently, binary ICEM meshes (2D and 3D) and      *";
  cout << "\n*  Triangle (2D) generated finite element meshes can  *";
  cout << "\n*  be converted. ICEM meshes will always be isopara-  *";
  cout << "\n*  matrix, Triangle meshes can be isoparametric or    *";
  cout << "\n*  with analytical finite elements.                   *";
  cout << "\n*                                                     *";
  cout << "\n*  The user must supply the file names for            *";
  cout << "\n*  (1) the input finite element mesh                  *";
  cout << "\n*  (2) the output VSet file                           *";
  cout << "\n*  (3) the CSMP physical variables text file          *";
  cout << "\n*                                                     *";
  cout << "\n*  NB: All file names must be supplied without the    *";
  cout << "\n*  respective extensions (e.g. *.txt, *.asc, etc.)    *";
  cout << "\n*                                                     *";
  cout << "\n*  Version 1.1 (April 24, 2008 SG)                    *";
  cout << "\n*                                                     *";
  cout << "\n*******************************************************";
  cout << endl << endl;
  cout.flush();
}

}

