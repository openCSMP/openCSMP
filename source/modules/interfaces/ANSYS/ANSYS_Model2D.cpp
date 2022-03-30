#include "ANSYS_Model2D.h"
#include "ANSYS_Interface.h"
#include "Region.h"
#include "Model.h"
#include "Exception.h"
#include "ModelTopology.h"
#include "NodeManifold.h"

using namespace std;

namespace csmp {

void ANSYS_Model2D::InitializeANSYS( bool isoparametric,
                                      const char* mesh_file_set,
                                      const char* regions_file_prefix,
                                      bool irregular_mesh,
                                      bool binary_input_file,
                                      bool use_regions_file )
{
  // -------------------------------------------------  
  // initializing the empty Model from the ANSYS
  // data imported into a vset.  
  // -------------------------------------------------
  try {

    VSet<2U> vset;
    bool     isoparametric_elements( isoparametric );

    ModelTopology    mesh_topology( isoparametric_elements );
    ANSYS_Interface  mesh_interface( isoparametric_elements );

    // 0. reading mesh from ANSYS CSMP-input .asc and .dat files,
    //    eliminating unwanted line/surface element regions
    mesh_interface.Read_ANSYS_Mesh( std::string( mesh_file_set ), vset, mesh_topology, binary_input_file, true );
 
    // 1. recreating 'pfverts' information because ANSYS ICEM CFD does not get the neighbor connectivity right
    vset.RemovePfverts();
    vset.EstablishElementConnectivity2D();

    // 2. writing original element and node numbers to property data and storing them in VSet
    if ( Database().IsDefined( "element number" ) ) {
      // element numbers
      PropertyData elmt_nums( ELEMENT, SCALAR, 2U );
      elmt_nums.Reserve( vset.Elements() );
      for ( size_t i{0U}; i<vset.Elements(); ++i ) pushBack( elmt_nums, makeScalar( ANY, i ) );
      vset.AddData( "element number", elmt_nums );
    }
    if ( Database().IsDefined( "node number" ) ) {
      // node numbers
      PropertyData node_nums( NODE, SCALAR, 2U );
      node_nums.Reserve( vset.Vertices() );
      for ( size_t i{0U}; i<vset.Vertices(); ++i ) pushBack( node_nums, makeScalar( ANY, i ) );
      vset.AddData( "node number", node_nums );
    }

    // 3. constructing model from the polygonal data in the VSet and the region information in model topology
    // 3.1 using only the selected regions from the -regions.txt file
    if ( use_regions_file )
      Model<2U>::Initialize( regions_file_prefix,
                             mesh_topology,
                             vset );
    else
      // 3.2 using all regions from the ANSYS model
      Model<2U>::Initialize( mesh_topology,
                             vset );
  }

  // --------------------------------------_-----------  
  // catching all possible standard and CSMP exceptions
  // ---------------------------------------_----------
  catch ( bad_alloc& ba ) {
    cout << "\nbad_alloc: Memory allocation error caused by: " << ba.what() << endl;
  }
  catch ( bad_cast& ba ) {
    cout << "\nbad_cast: Type casting error caused by: " << ba.what() << endl;
  }
  catch ( bad_exception& ba ) {
    cout << "\nbad_exception: Exception error caused by: " << ba.what() << endl;
  }
  catch ( bad_typeid& ba ) {
    cout << "\nbad_typeid: Type ID error caused by: " << ba.what() << endl;
  }
  catch ( ios_base::failure& ba ) {
    cout << "\nios_base::failure: Probable I/O error caused by: " << ba.what() << endl;
  }
  // standard logic errors
  catch ( domain_error& ba ) {
    cout << "\ndomain_error: Logic error caused by: " << ba.what() << endl;
  }
  catch ( invalid_argument& ba ) {
    cout << "\ninvalid_argument: Logic error caused by: " << ba.what() << endl;
  }
  catch ( length_error& ba ) {
    cout << "\nlength_error: Logic error caused by: " << ba.what() << endl;
  }
  catch ( out_of_range& ba ) {
    cout << "\nout_of_range: Logic error caused by: " << ba.what() << endl;
  }
  // runtime errors
  catch ( overflow_error& ba ) {
    cout << "\noverflow_error: Runtime error caused by: " << ba.what() << endl;
  }
  catch ( range_error& ba ) {
    cout << "\nrange_error: Runtime error caused by: " << ba.what() << endl;
  }
  catch ( underflow_error& ba ) {
    cout << "\nunderflow_error: Runtime error caused by: " << ba.what() << endl;
  }
  catch ( Exception& ba ) {
    cout << "\nException: Exception raised by: " << ba.What() << endl;
    cout << "\nDiagnostics:" << endl;
    ba.Out();
  }

} // end Initialize


/** Builds Model after it was constructed with the default constructor.
*/
void ANSYS_Model2D::InitializeANSYS( const char* mesh_file_set,
                                     const char* regions_file_prefix,
                                     bool irregular_mesh,
                                     bool binary_input_file,
                                     bool use_regions_file )
{
  // -------------------------------------------------  
  // initializing the empty Model from the ANSYS
  // data imported into a vset.  
  // -------------------------------------------------
  try {

    VSet<2U> vset;
    bool isoparametric_elements( true );

    ModelTopology    mesh_topology( isoparametric_elements );
    ANSYS_Interface  mesh_interface( isoparametric_elements );

    // 0. reading the mesh from ANSYS CSMP-input files, eliminating the unwanted line/surface element regions
    //    node numbers of triangles and quadrilaterals are reversed if they are in clockwise order (establishing right-hand coordinate compliance)
    const bool recreate_node_boundary_flags{true}; // does this using the lower-dimensional boundary regions 
    mesh_interface.Read_ANSYS_Mesh( std::string( mesh_file_set ), vset, mesh_topology, binary_input_file, recreate_node_boundary_flags );
    
    // create 'pfverts' information because the one ANSYS does not get the line element orientations right
    vset.RemovePfverts();
    vset.EstablishElementConnectivity2D();
    
    // 1. writing element and node numbers to property data and storing them in the VSet
    if ( Database().IsDefined( "element number" ) ) {
      // element numbers
      const uint32_t dim{2};
      PropertyData elmt_nums( ELEMENT, SCALAR, dim );
      elmt_nums.Reserve( vset.Elements() );
      for ( size_t i{0U}; i<vset.Elements(); ++i ) pushBack( elmt_nums, makeScalar( ANY, i ) );
      vset.AddData( "element number", elmt_nums );
    }
    if ( Database().IsDefined( "node number" ) ) {
      // node numbers
      const uint32_t dim{2};
      PropertyData node_nums( NODE, SCALAR, dim );
      node_nums.Reserve( vset.Vertices() );
      for ( size_t i{0U}; i<vset.Vertices(); ++i ) pushBack( node_nums, makeScalar( ANY, i ) );
      vset.AddData( "node number", node_nums );
    }

    // 2. construct model based on obtained model topology and vset
    //    ansys neighbor info will be overwritten later since it includes neighbor information
    //    of elements of different dimensionality (i.e. e volumetric element has a surface element neighbor)
    if ( use_regions_file )
      Model<2U>::Initialize( regions_file_prefix,
                             mesh_topology,
                             vset );
    else
      Model<2U>::Initialize( mesh_topology,
                             vset );
  } // end try


  // --------------------------------------_-----------  
  // catching all possible standard and CSMP exceptions
  // ---------------------------------------_----------
  catch ( bad_alloc& ba ) {
    cout << "\nbad_alloc: Memory allocation error caused by: " << ba.what() << endl;
  }
  catch ( bad_cast& ba ) {
    cout << "\nbad_cast: Type casting error caused by: " << ba.what() << endl;
  }
  catch ( bad_exception& ba ) {
    cout << "\nbad_exception: Exception error caused by: " << ba.what() << endl;
  }
  catch ( bad_typeid& ba ) {
    cout << "\nbad_typeid: Type ID error caused by: " << ba.what() << endl;
  }
  catch ( ios_base::failure& ba ) {
    cout << "\nios_base::failure: Probable I/O error caused by: " << ba.what() << endl;
  }
  // standard logic errors
  catch ( domain_error& ba ) {
    cout << "\ndomain_error: Logic error caused by: " << ba.what() << endl;
  }
  catch ( invalid_argument& ba ) {
    cout << "\ninvalid_argument: Logic error caused by: " << ba.what() << endl;
  }
  catch ( length_error& ba ) {
    cout << "\nlength_error: Logic error caused by: " << ba.what() << endl;
  }
  catch ( out_of_range& ba ) {
    cout << "\nout_of_range: Logic error caused by: " << ba.what() << endl;
  }
  // runtime errors
  catch ( overflow_error& ba ) {
    cout << "\noverflow_error: Runtime error caused by: " << ba.what() << endl;
  }
  catch ( range_error& ba ) {
    cout << "\nrange_error: Runtime error caused by: " << ba.what() << endl;
  }
  catch ( underflow_error& ba ) {
    cout << "\nunderflow_error: Runtime error caused by: " << ba.what() << endl;
  }
  catch ( Exception& ba ) {
    cout << "\nException: Exception raised by: " << ba.What() << endl;
    cout << "\nDiagnostics:" << endl;
    ba.Out();
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
                              bool use_regions_file )
  : Model<2U>( variable_file )
{
  this->Name( icem_file_set );
  InitializeANSYS( icem_file_set,
                   regions_file_prefix,
                   irregular_mesh,
                   binary_file,
                   use_regions_file );
}


/**

Default constructor of Model is called. Then the model is build
from the 'icem_file_set' '*.asc' and '*.dat' files using the method
Initialize.
*/
ANSYS_Model2D::ANSYS_Model2D( bool isoparametric,
                              const char* icem_file_set,
                              const char* variable_file,
                              bool irregular_mesh,
                              bool binary_file,
                              bool use_regions_file )
  : Model<2U>( variable_file )
{
  this->Name( icem_file_set );
  InitializeANSYS( isoparametric,
                   icem_file_set,
                   icem_file_set,
                   irregular_mesh,
                   binary_file,
                   use_regions_file );
}

ANSYS_Model2D::ANSYS_Model2D( const char* icem_file_set,
                              const char* variable_file,
                              bool irregular_mesh,
                              bool binary_file,
                              bool use_regions_file )
  : Model<2U>( variable_file )
{
  this->Name( icem_file_set );
  InitializeANSYS( icem_file_set,
                   icem_file_set,
                   irregular_mesh,
                   binary_file,
                   use_regions_file );
}

ANSYS_Model2D::ANSYS_Model2D( const char* icem_file_set,
                              bool irregular_mesh,
                              bool binary_file,
                              bool use_regions_file )
{
  this->Name( icem_file_set );
  InitializeANSYS( icem_file_set,
                   icem_file_set,
                   irregular_mesh,
                   binary_file,
                   use_regions_file );
}


// nothing to do here 
ANSYS_Model2D::~ANSYS_Model2D()
{
}


} // end csmp
