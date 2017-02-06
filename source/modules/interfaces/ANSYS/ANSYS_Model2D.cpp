#include "ANSYS_Model2D.h"
#include "ANSYS_Interface.h"
#include "Region.h"
#include "Model.h"
#include "Exception.h"
#include "ModelTopology.h"

using namespace std;

namespace csmp {

/** Builds Model after it was constructed with the default constructor.
 */
void ANSYS_Model2D::Initialize( const char* mesh_file_set,
                                const char* regions_file_prefix,
                                bool irregular_mesh,
                                bool binary_input_file,
                                bool use_regions_file,
                                bool create_boundaries )
{
  // -------------------------------------------------  
  // initializing the empty Model from the ANSYS
  // data imported into a vset.  
  // -------------------------------------------------
  try {

       VSet<2U>         vset;
       bool isoparametric_elements( true );

       ModelTopology    mesh_topology (isoparametric_elements);
       ANSYS_Interface  mesh_interface(isoparametric_elements);

       // 0. reading the mesh from ANSYS CSMP-input files
       //    and eliminating the unwanted line/surface element regions
       mesh_interface.Read_ANSYS_Mesh( std::string(mesh_file_set), vset, mesh_topology, binary_input_file, irregular_mesh );

       // 1. construct model based on obtained model topology and vset
       //    ansys neighbor info will be overwritten later since it includes neighbor information
       //    of elements of different dimensionality (i.e. e volumetric element has a surface element neighbor)
       if ( use_regions_file )
         Model<2U>::Initialize( regions_file_prefix,
                                mesh_topology,
                                vset,
                                create_boundaries,
                                irregular_mesh );
       else
         Model<2U>::Initialize( mesh_topology,
                                vset,
                                create_boundaries,
                                irregular_mesh );
    }
    
  // --------------------------------------_-----------  
  // catching all possible standard and CSMP exceptions
  // ---------------------------------------_----------
  catch( bad_alloc& ba ) {
       cout <<"\nbad_alloc: Memory allocation error caused by: "<< ba.what() << endl;
    }
  catch( bad_cast& ba ) {
       cout <<"\nbad_cast: Type casting error caused by: "<< ba.what() << endl;
    }
  catch( bad_exception& ba ) {
       cout <<"\nbad_exception: Exception error caused by: "<< ba.what() << endl;
    }
  catch( bad_typeid& ba ) {
       cout <<"\nbad_typeid: Type ID error caused by: "<< ba.what() << endl;
    }
  catch( ios_base::failure& ba ) {
       cout <<"\nios_base::failure: Probable I/O error caused by: "<< ba.what() << endl;
    }
  // standard logic errors
  catch( domain_error& ba ) {
       cout <<"\ndomain_error: Logic error caused by: "<< ba.what() << endl;
    }
  catch( invalid_argument& ba ) {
       cout <<"\ninvalid_argument: Logic error caused by: "<< ba.what() << endl;
    }
  catch( length_error& ba ) {
       cout <<"\nlength_error: Logic error caused by: "<< ba.what() << endl;
     }
  catch( out_of_range& ba ) {
       cout <<"\nout_of_range: Logic error caused by: "<< ba.what() << endl;
    }
  // runtime errors
  catch( overflow_error& ba ) {
       cout <<"\noverflow_error: Runtime error caused by: "<< ba.what() << endl;
    }
  catch( range_error& ba ) {
       cout <<"\nrange_error: Runtime error caused by: "<< ba.what() << endl;
    }
  catch( underflow_error& ba ) {
       cout <<"\nunderflow_error: Runtime error caused by: "<< ba.what() << endl;
    }
  catch( Exception& ba ) {
       cout <<"\nException: Exception raised by: "<< ba.What() << endl;
       cout <<"\nDiagnostics:"<< endl;
       ba.Out(cout);
    }
         
 } // end Initialize




/**

Default constructor of Model is called. Then the model is build
from the 'icem_file_set' '*.asc' and '*.dat' files using the method
Initialize.
*/
ANSYS_Model2D::ANSYS_Model2D( const char* icem_file_set,
                              const char* regions_file_prefix,
                              const char* variable_file,
                              bool irregular_mesh,
                              bool binary_file,
                              bool use_regions_file,
                              bool create_boundaries)
 : Model<2U>( variable_file, false )
 {
    this->Name( icem_file_set );
    Initialize( icem_file_set,
                regions_file_prefix,
                irregular_mesh,
                binary_file,
                use_regions_file,
                create_boundaries);
 }


/**

Default constructor of Model is called. Then the model is build
from the 'icem_file_set' '*.asc' and '*.dat' files using the method
Initialize.   
*/
ANSYS_Model2D::ANSYS_Model2D( const char* icem_file_set, 
                              const char* variable_file,
                              bool irregular_mesh,
                              bool binary_file,
                              bool use_regions_file,
                              bool create_boundaries)
 : Model<2U>( variable_file, false )
 {
    this->Name( icem_file_set );
    Initialize( icem_file_set,
                icem_file_set,
                irregular_mesh,
                binary_file,
                use_regions_file,
                create_boundaries);
 }

ANSYS_Model2D::ANSYS_Model2D( const char* icem_file_set,
                              bool irregular_mesh,
                              bool binary_file,
                              bool use_regions_file,
                              bool create_boundaries)
 {
    this->Name(icem_file_set);
    Initialize( icem_file_set,
                icem_file_set,
                irregular_mesh,
                binary_file,
                use_regions_file,
                create_boundaries);
 }

 
// nothing to do here 
ANSYS_Model2D::~ANSYS_Model2D()
 {
 }


} // end csmp
