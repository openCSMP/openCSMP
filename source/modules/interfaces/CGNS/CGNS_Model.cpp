#include "CGNS_Model.h"

using namespace std;

namespace csmp {


/**
Default constructor of Model is called.
The model is built from the 'mesh_file_set', variables file,
and the regions file prefix is used to read the regions file.
Uses the method Initialize.

*/

  template<size_t dim>
  CGNS_Model<dim>::CGNS_Model( const std::string& mesh_file_set,
                               const std::string& regions_file_prefix,
                               const std::string& variable_file,
                               bool use_regions_file,
                               bool create_boundaries
                             )
   : Model<dim>( variable_file.c_str(), false )
   {
      this->Name( mesh_file_set.c_str() );
      Initialize( mesh_file_set,
                  regions_file_prefix,
                  use_regions_file,
                  create_boundaries );
   }


  template<size_t dim>
  CGNS_Model<dim>::CGNS_Model( const std::string& mesh_file_set,
                               const std::string& variable_file,
                               bool use_regions_file,
                               bool create_boundaries
                             )
   : Model<dim>( variable_file.c_str(), false )
   {
      this->Name( mesh_file_set.c_str() );
      std::string regions_file_prefix( mesh_file_set );
      Initialize( mesh_file_set,
                  regions_file_prefix,
                  use_regions_file,
                  create_boundaries );
   }


  template<size_t dim>
  CGNS_Model<dim>::CGNS_Model( const std::string& mesh_file_set,
                               bool use_regions_file,
                               bool create_boundaries
                             )
   : Model<dim>()
   {
      this->Name( mesh_file_set.c_str() );
      std::string regions_file_prefix( mesh_file_set );
      Initialize( mesh_file_set,
                  regions_file_prefix,
                  use_regions_file,
                  create_boundaries );
   }

  template<size_t dim>
  CGNS_Model<dim>::CGNS_Model( const std::string& vset_dat_files,
                               const std::string& variable_file )
   : Model<dim>( vset_dat_files.c_str(), variable_file.c_str() )
   {
      this->Name(vset_dat_files.c_str());
   }

  template<size_t dim>
  CGNS_Model<dim>::CGNS_Model( const std::string& vset_dat_files )
    : Model<dim>( vset_dat_files.c_str() )
    {
      this->Name(vset_dat_files.c_str());
    }

  template<size_t dim>
  CGNS_Model<dim>::~CGNS_Model()
    {
    }

  template<size_t dim>
  void CGNS_Model<dim>::Write_CGNS_Mesh( const std::string& mesh_file_set )
  {
      bool isoparametric_elements( true );
      CGNS_Interface mesh_interface( isoparametric_elements );
      mesh_interface.Write_CGNS_Mesh( mesh_file_set, *this );
  }

  /// Builds Model after it was constructed with the default constructor.
  template<size_t dim>
  void CGNS_Model<dim>::Initialize( const std::string& mesh_file_set,
                                    const std::string& regions_file_prefix,
                                    bool use_regions_file,
                                    bool create_boundaries )
  {
    double64& model_time( ModelTime::Instance().modelTime );
    model_time = 0.;

    // -------------------------------------------------
    // initializing the empty Model from the CGNS
    // data imported into a vset.
    // -------------------------------------------------
    try {
         VSet<dim>       vset;

         bool isoparametric_elements( true );
         ModelTopology  mesh_topology ( isoparametric_elements );
         CGNS_Interface mesh_interface( isoparametric_elements );

         // 0. reading the mesh from CGNS CSMP-input files
         mesh_interface.Read_CGNS_Mesh( mesh_file_set,
                                        vset,
                                        mesh_topology );

         // 1. construct model based on obtained model topology and vset
         //    we won't use ansys neighbor info since it includes neighbor information
         //    of elements of different dimensionality (i.e. e volumetric element has a surface element neighbor)
         bool irregular_mesh( true );

         if ( use_regions_file )
           Model<dim>::Initialize( regions_file_prefix.c_str(),
                                   mesh_topology, vset,
                                   create_boundaries,
                                   irregular_mesh );
         else
           Model<dim>::Initialize( mesh_topology, vset,
                                   create_boundaries,
                                   irregular_mesh );
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
         ba.Out();
      }

   } // end Initialize

template class CGNS_Model<1U>;
template class CGNS_Model<2U>;
template class CGNS_Model<3U>;

} // end csmp
