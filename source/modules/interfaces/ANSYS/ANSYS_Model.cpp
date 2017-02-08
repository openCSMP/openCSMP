#include "ANSYS_Model.h"
#include "ModelTime.h"

using namespace std;

namespace csmp {


/**
Default constructor of Model is called.
The model is built from the 'icem_file_set' '*.asc' and '*.dat', variables file,
and the regions file prefix is used to read the regions file.
Or the model is build from *-1D-mesh.txt
in case if ansys interface is set to false
Uses the method Initialize.

*/

  template<size_t dim>
  ANSYS_Model<dim>::ANSYS_Model( const std::string& mesh_file_set,
                                 const std::string& regions_file_prefix,
                                 const std::string& variable_file,
                                 bool irregular_mesh,
                                 bool binary_file,
                                 bool use_regions_file,
                                 bool create_boundaries,
                                 bool ansys_interface
                               )
   : Model<dim>( variable_file.c_str(), false )
   {
      this->Name( mesh_file_set.c_str() );
      Initialize( mesh_file_set,
                  regions_file_prefix,
                  irregular_mesh,
                  binary_file,
                  use_regions_file,
                  create_boundaries,
                  ansys_interface );
   }


  template<size_t dim>
  ANSYS_Model<dim>::ANSYS_Model( const std::string& mesh_file_set,
                                 const std::string& variable_file,
                                 bool irregular_mesh,
                                 bool binary_file,
                                 bool use_regions_file,
                                 bool create_boundaries,
                                 bool ansys_interface
                               )
   : Model<dim>( variable_file.c_str(), false )
   {
      this->Name( mesh_file_set.c_str() );
      std::string regions_file_prefix( mesh_file_set );
      Initialize( mesh_file_set,
                  regions_file_prefix,
                  irregular_mesh,
                  binary_file,
                  use_regions_file,
                  create_boundaries,
                  ansys_interface );
   }


  template<size_t dim>
  ANSYS_Model<dim>::ANSYS_Model( const std::string& mesh_file_set,
                                 bool irregular_mesh,
                                 bool binary_file,
                                 bool use_regions_file,
                                 bool create_boundaries,
                                 bool ansys_interface
                               )
   {
      this->Name(mesh_file_set.c_str());
      std::string regions_file_prefix( mesh_file_set );
      Initialize( mesh_file_set,
                  regions_file_prefix,
                  irregular_mesh,
                  binary_file,
                  use_regions_file,
                  create_boundaries,
                  ansys_interface );
   }

  template<size_t dim>
  ANSYS_Model<dim>::~ANSYS_Model()
    {
    }

  /// Builds Model after it was constructed with the default constructor.
  template<size_t dim>
  void ANSYS_Model<dim>::Initialize( const std::string& mesh_file_set,
                                     const std::string& regions_file_prefix,
                                     bool irregular_mesh,
                                     bool binary_input_file,
                                     bool use_regions_file,
                                     bool create_boundaries,
                                     bool ansys_interface )
  {
    double64& model_time( ModelTime::Instance().modelTime );
    model_time = 0.;

    // -------------------------------------------------
    // initializing the empty Model from the ANSYS
    // data imported into a vset.
    // -------------------------------------------------
    try {

         VSet<dim>       vset;
         bool isoparametric_elements( true );

         if( ansys_interface && ( dim != 1U ) )
         {
             ModelTopology   mesh_topology ( isoparametric_elements );
             ANSYS_Interface mesh_interface( isoparametric_elements );

             // 0. reading the mesh from ANSYS CSMP-input files
             mesh_interface.Read_ANSYS_Mesh( mesh_file_set,
                                             vset,
                                             mesh_topology,
                                             binary_input_file,
                                             irregular_mesh );

             // 1. construct model based on obtained model topology and vset
             if ( use_regions_file )
               Model<dim>::Initialize( regions_file_prefix.c_str(),
                                       mesh_topology,
                                       vset,
                                       create_boundaries,
                                       irregular_mesh );
             else
               Model<dim>::Initialize( mesh_topology,
                                       vset,
                                       create_boundaries,
                                       irregular_mesh );
         }else{

             throw csmp::Exception( CSMP_ERROR,
                                    "ANSYS_Model<dim>::Initialize()",
                                    "Reading of ANSYS mesh file is not yet implemented in 1D! " );
         }

      }

    // -------------------------------------------------
    // catching all possible standard and csmp::Exceptions
    // -------------------------------------------------
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


  template class ANSYS_Model<1>;
  template class ANSYS_Model<2>;
  template class ANSYS_Model<3>;

} // end csmp
