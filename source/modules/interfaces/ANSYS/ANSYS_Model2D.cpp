#include "ANSYS_Model2D.h"
#include "ANSYS_Interface.h"
#include "Region.h"
#include "Model.h"
#include "Exception.h"
#include "ModelTopology.h"

using namespace std;

namespace csmp {

void ANSYS_Model2D::Initialize( bool isoparametric,
                                const char* mesh_file_set,
                                const char* regions_file_prefix,
                                bool irregular_mesh,
                                bool binary_input_file,
                                bool use_regions_file,
                                bool create_boundaries,
                                bool create_splitboundaries )
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

    // 0. reading the mesh from ANSYS CSMP-input files
    //    and eliminating the unwanted line/surface element regions
    mesh_interface.Read_ANSYS_Mesh( std::string( mesh_file_set ), vset, mesh_topology, binary_input_file, irregular_mesh );
    // create 'pfverts' information because the one ANSYS does not get the line element orientations right
    vset.EstablishElementConnectivity2D();
    vset.CreateConsistentLineElementOrientations2D();

    // 1. writing element and node numbers to property data and storing them in the VSet
    if ( Database().IsDefined( "element number" ) ) {
      // element numbers
      PropertyData elmt_nums( ELEMENT, SCALAR, 2U );
      elmt_nums.Reserve( vset.Elements() );
      for ( size_t i = 0U; i<vset.Elements(); ++i ) pushBack( elmt_nums, makeScalar( ANY, i ) );
      vset.AddData( "element number", elmt_nums );
    }
    if ( Database().IsDefined( "node number" ) ) {
      // node numbers
      PropertyData node_nums( NODE, SCALAR, 2U );
      node_nums.Reserve( vset.Vertices() );
      for ( size_t i = 0U; i<vset.Vertices(); ++i ) pushBack( node_nums, makeScalar( ANY, i ) );
      vset.AddData( "node number", node_nums );
    }

    // 2. construct model based on obtained model topology and vset
    //    ansys neighbor info will be overwritten later since it includes neighbor information
    //    of elements of different dimensionality (i.e. e volumetric element has a surface element neighbor)
    if ( use_regions_file )
      Model<2U>::Initialize( regions_file_prefix,
                             mesh_topology,
                             vset,
                             create_boundaries,
                             create_splitboundaries,
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
void ANSYS_Model2D::Initialize( const char* mesh_file_set,
                                const char* regions_file_prefix,
                                bool irregular_mesh,
                                bool binary_input_file,
                                bool use_regions_file,
                                bool create_boundaries,
                                bool create_splitboundaries )
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

    // 0. reading the mesh from ANSYS CSMP-input files
    //    and eliminating the unwanted line/surface element regions
    mesh_interface.Read_ANSYS_Mesh( std::string( mesh_file_set ), vset, mesh_topology, binary_input_file, irregular_mesh );
    // create 'pfverts' information because the one ANSYS does not get the line element orientations right
    vset.EstablishElementConnectivity2D();
    vset.VData::Out();
    vset.CreateConsistentLineElementOrientations2D();
    
    // 1. writing element and node numbers to property data and storing them in the VSet
    if ( Database().IsDefined( "element number" ) ) {
      // element numbers
      PropertyData elmt_nums( ELEMENT, SCALAR, 2U );
      elmt_nums.Reserve( vset.Elements() );
      for ( size_t i = 0U; i<vset.Elements(); ++i ) pushBack( elmt_nums, makeScalar( ANY, i ) );
      vset.AddData( "element number", elmt_nums );
    }
    if ( Database().IsDefined( "node number" ) ) {
      // node numbers
      PropertyData node_nums( NODE, SCALAR, 2U );
      node_nums.Reserve( vset.Vertices() );
      for ( size_t i = 0U; i<vset.Vertices(); ++i ) pushBack( node_nums, makeScalar( ANY, i ) );
      vset.AddData( "node number", node_nums );
    }

    // 2. construct model based on obtained model topology and vset
    //    ansys neighbor info will be overwritten later since it includes neighbor information
    //    of elements of different dimensionality (i.e. e volumetric element has a surface element neighbor)
    if ( use_regions_file )
      Model<2U>::Initialize( regions_file_prefix,
                             mesh_topology,
                             vset,
                             create_boundaries,
                             create_splitboundaries,
                             irregular_mesh );
    else
      Model<2U>::Initialize( mesh_topology,
                             vset,
                             create_boundaries,
                             create_splitboundaries,
                             irregular_mesh );
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
                              bool create_boundaries,
                              bool create_splitboundaries )
  : Model<2U>( variable_file )
{
  this->Name( icem_file_set );
  Initialize( icem_file_set,
              regions_file_prefix,
              irregular_mesh,
              binary_file,
              use_regions_file,
              create_boundaries,
              create_splitboundaries );
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
                              bool use_regions_file,
                              bool create_boundaries,
                              bool create_splitboundaries )
  : Model<2U>( variable_file )
{
  this->Name( icem_file_set );
  Initialize( isoparametric,
              icem_file_set,
              icem_file_set,
              irregular_mesh,
              binary_file,
              use_regions_file,
              create_boundaries,
              create_splitboundaries );
}

ANSYS_Model2D::ANSYS_Model2D( const char* icem_file_set,
                              const char* variable_file,
                              bool irregular_mesh,
                              bool binary_file,
                              bool use_regions_file,
                              bool create_boundaries,
                              bool create_splitboundaries )
  : Model<2U>( variable_file )
{
  this->Name( icem_file_set );
  Initialize( icem_file_set,
              icem_file_set,
              irregular_mesh,
              binary_file,
              use_regions_file,
              create_boundaries,
              create_splitboundaries );
}

ANSYS_Model2D::ANSYS_Model2D( const char* icem_file_set,
                              bool irregular_mesh,
                              bool binary_file,
                              bool use_regions_file,
                              bool create_boundaries,
                              bool create_splitboundaries )
{
  this->Name( icem_file_set );
  Initialize( icem_file_set,
              icem_file_set,
              irregular_mesh,
              binary_file,
              use_regions_file,
              create_boundaries,
              create_splitboundaries );
}


// nothing to do here 
ANSYS_Model2D::~ANSYS_Model2D()
{
}


} // end csmp
