  #include "EclipseModel.h"
#include "ModelTime.h"

namespace csmp {

/// Model constructor with provided "variables_file.txt" file is used
template<size_t dim>
EclipseModel<dim>::EclipseModel( const std::string& model_name,
                            const std::string& variables_file )
: csmp::Model<dim>( variables_file.c_str(), false ),
  eclipse_model_settings_( model_name )
{
   this->Name( model_name.c_str() );
}


/// Default Model constructor with empty property data base is called
template<size_t dim>
EclipseModel<dim>::EclipseModel( const std::string& model_name )
: csmp::Model<dim>(),
  eclipse_model_settings_( model_name )
{
   this->Name( model_name.c_str() );
}

template<size_t dim>
EclipseModel<dim>::~EclipseModel()
{
}

// MODEL NAME
template<size_t dim>
const char* EclipseModel<dim>::Name() const
{
    return model_name_.c_str();
}

template<size_t dim>
void EclipseModel<dim>::Name( const char* new_name )
{
    model_name_ = new_name;
}

// ACCESS TO MESH INTERFACE
template<size_t dim>
EclipseModelSettings& EclipseModel<dim>::EclipseModelSetup( )
{
    return eclipse_model_settings_;
}

/**  
    Master function that BUILDS CSMP MODEL FROM ECLIPSE DATA
*/
template<size_t dim>
void EclipseModel<dim>::BuildModel()
{
    csmp::ErrorHandler& error_handler(csmp::ErrorHandler::Instance());

    double64& model_time( csmp::ModelTime::Instance().modelTime );
    model_time = 0.;

    try {
         csmp::VSet<dim>  vset;
         bool isoparametric_elements( true );

         csmp::ModelTopology   mesh_topology( isoparametric_elements );
         EclipseInterface<dim> mesh_interface;

         // =====================================================================
         // 0. reads grid from ECLIPSE input files and converts into CSMP mesh
         // =====================================================================
         mesh_interface.SetProperties( eclipse_model_settings_.properties_ );
      
         // KEY METHOD here
         mesh_interface.ReadFile( vset,
                                  mesh_topology,
                                  eclipse_model_settings_.mesh_file_prefix_,
                                  eclipse_model_settings_.exclude_inactive_cells_,
                                  eclipse_model_settings_.tetra_mesh_ );

         // if there is no REGIONS section, i.e. no FIPNUM, SATNUM, EQLNUM or PVTNUM cell specifiers in the Eclipse input deck
         // NB: a regions file must be present for the region to be preserved
         if ( eclipse_model_settings_.regions_.empty() ) {
               // the sets of strings will be empty if no regions, faults or wells were detected by the mesh interface
               mesh_interface.GetRegions( regions_ );
               mesh_interface.GetFaults( faults_ );
               mesh_interface.GetWells( wells_ );
           }
         else {
               std::set<std::string>& desired_regions( eclipse_model_settings_.regions_ );
               std::set<std::string>  regions;

               /// assign regions
               regions.clear();
               mesh_interface.GetRegions( regions );
               std::set_intersection( desired_regions.begin(), desired_regions.end(),
                                      regions.begin(), regions.end(),
                                      std::inserter( regions_, regions_.begin() ) );
               /// assign faults
               regions.clear();
               mesh_interface.GetFaults( regions );
               std::set_intersection( desired_regions.begin(), desired_regions.end(),
                                      regions.begin(), regions.end(),
                                      std::inserter( faults_, faults_.begin() ) );

               /// assign wells
               regions.clear();
               mesh_interface.GetWells( regions );
               std::set_intersection( desired_regions.begin(), desired_regions.end(),
                                      regions.begin(), regions.end(),
                                      std::inserter( wells_, wells_.begin() ) );
           }

         // =====================================================================
         // 1. selectively read properties of interest, adding them to VSET
         // =====================================================================
         // porosity, permeability, saturations
         std::set<std::string> vset_props;
         std::set<std::string>::iterator prop_it;
         for( std::map<int,csmp::Parameter>::const_iterator
              epit = eclipse_model_settings_.properties_.begin(); epit != eclipse_model_settings_.properties_.end(); ++epit )
           {
               const csmp::Parameter& prop = (*epit).second;
               prop_it = vset_props.find( prop.name );
               if( prop_it != vset_props.end() )
               {
                   if( !this->Database().IsDefined( prop.name.c_str() ) )
                       // SKM FIX 
                       this->Database().AddProperty( prop.name.c_str(), prop.unit.c_str(),
                                                     prop.key.type, prop.key.place, this->Database().VariableCount( prop.key.place, prop.key.type), prop.min, prop.max,
                                                     prop.usage.c_str() );
                   vset_props.erase( prop_it );
               }
           }

         /// checking whether any undefined properties are left
         std::set<std::string> undefined_props;
         for( prop_it = vset_props.begin(); prop_it != vset_props.end(); ++prop_it )
             if( !this->Database().IsDefined( (*prop_it).c_str() ) )
                 undefined_props.insert( *prop_it );
         if( !undefined_props.empty() )
         {
             if( error_handler.Verbose() )
             {
                 std::string message = "VSet contains data for undefined properties: ";
                 for( prop_it = undefined_props.begin(); prop_it != undefined_props.end(); ++prop_it )
                 {
                     if( prop_it != undefined_props.begin() )
                         message += ", ";
                     message += *prop_it;
                 }
                 message += " !!!";
                 error_handler.notice( csmp::INFO, "EclipseModel<dim>::BuildModel", message.c_str() );
             }
         }

         // =====================================================================
         // 2. construct CSMP model from obtained topology and mesh in vset
         // =====================================================================
         //    we won't use eclipse neighbor info since it includes neighbor information
         //    of elements of different dimensionality (i.e. e volumetric element has a surface element neighbors )
         const bool non_box_shaped_model( true );

         // with regions
         if ( !eclipse_model_settings_.regions_.empty() ) {
             // performing a ckeck whether element numbers in the VSet and the model topology match; else something went wrong
             // and user is given the possibility to call the subsequent method or not.
             const bool require_unique_names_for_vol_surf_lines(true);
             const bool correct_orientation_of_surface_elements(false);
             const bool non_box_boundary(true);
             if ( !mesh_topology.CheckTopology( vset, require_unique_names_for_vol_surf_lines,
                                                correct_orientation_of_surface_elements, non_box_boundary ) );
             else
               error_handler.notice( csmp::INFO, "EclipseModel<dim>::BuildModel", "ModelTopology=subdivision into regions is broken.");
           
             csmp::Model<dim>::Initialize( eclipse_model_settings_.regions_file_prefix_.c_str(),
                                           mesh_topology,
                                           vset,
                                            ( ( dim != 1U ) ? eclipse_model_settings_.create_boundaries_ : false ),
                                           non_box_shaped_model );
           }
         // no regions
         else {
             // all cells are lumped into the region "Eclipse Model" that is stored in the model topology
             const bool isoparametric(true);
             csmp::Model<dim>::Initialize( isoparametric, vset,
                                           ( ( dim != 1U ) ? eclipse_model_settings_.create_boundaries_ : false ),
                                           non_box_shaped_model );
           }
      }

    // -------------------------------------------------
    // catching all possible standard and csmp::Exceptions
    // -------------------------------------------------
    catch( std::bad_alloc& ba ) {
         std::cout <<"\nbad_alloc: Memory allocation error caused by: "<< ba.what() << std::endl;
      }
    catch( std::bad_cast& ba ) {
         std::cout <<"\nbad_cast: Type casting error caused by: "<< ba.what() << std::endl;
      }
    catch( std::bad_exception& ba ) {
         std::cout <<"\nbad_exception: Exception error caused by: "<< ba.what() << std::endl;
      }
    catch( std::bad_typeid& ba ) {
         std::cout <<"\nbad_typeid: Type ID error caused by: "<< ba.what() << std::endl;
      }
    catch( std::ios_base::failure& ba ) {
         std::cout <<"\nios_base::failure: Probable I/O error caused by: "<< ba.what() << std::endl;
      }
    // standard logic errors
    catch( std::domain_error& ba ) {
         std::cout <<"\ndomain_error: Logic error caused by: "<< ba.what() << std::endl;
      }
    catch( std::invalid_argument& ba ) {
         std::cout <<"\ninvalid_argument: Logic error caused by: "<< ba.what() << std::endl;
      }
    catch( std::length_error& ba ) {
         std::cout <<"\nlength_error: Logic error caused by: "<< ba.what() << std::endl;
       }
    catch( std::out_of_range& ba ) {
         std::cout <<"\nout_of_range: Logic error caused by: "<< ba.what() << std::endl;
      }
    // runtime errors
    catch( std::overflow_error& ba ) {
         std::cout <<"\noverflow_error: Runtime error caused by: "<< ba.what() << std::endl;
      }
    catch( std::range_error& ba ) {
         std::cout <<"\nrange_error: Runtime error caused by: "<< ba.what() << std::endl;
      }
    catch( std::underflow_error& ba ) {
         std::cout <<"\nunderflow_error: Runtime error caused by: "<< ba.what() << std::endl;
      }
    catch( csmp::Exception& ba ) {
         std::cout <<"\nException: Exception raised by: "<< ba.What() << std::endl;
         std::cout <<"\nDiagnostics:"<< std::endl;
         ba.Out();
      }

} // end BuildModel

/// existing special regions


template<size_t dim>
template<class Container>
void EclipseModel<dim>
::GetRegions( Container& data )
{
    data.clear();
    typename Container::iterator dit = data.begin();
    std::copy( regions_.begin(), regions_.end(), std::inserter( data, dit ) );
}

template void EclipseModel<1U>::GetRegions( std::vector<std::string>& );
template void EclipseModel<1U>::GetRegions( std::list<std::string>& );
template void EclipseModel<1U>::GetRegions( std::set<std::string>& );

template void EclipseModel<2U>::GetRegions( std::vector<std::string>& );
template void EclipseModel<2U>::GetRegions( std::list<std::string>& );
template void EclipseModel<2U>::GetRegions( std::set<std::string>& );

template void EclipseModel<3U>::GetRegions( std::vector<std::string>& );
template void EclipseModel<3U>::GetRegions( std::list<std::string>& );
template void EclipseModel<3U>::GetRegions( std::set<std::string>& );

template<size_t dim>
template<class Container>
void EclipseModel<dim>
::GetFaults( Container& data )
{
    data.clear();
    typename Container::iterator dit = data.begin();
    std::copy( faults_.begin(), faults_.end(), std::inserter( data, dit ) );
}

template void EclipseModel<1U>::GetFaults( std::vector<std::string>& );
template void EclipseModel<1U>::GetFaults( std::list<std::string>& );
template void EclipseModel<1U>::GetFaults( std::set<std::string>& );

template void EclipseModel<2U>::GetFaults( std::vector<std::string>& );
template void EclipseModel<2U>::GetFaults( std::list<std::string>& );
template void EclipseModel<2U>::GetFaults( std::set<std::string>& );

template void EclipseModel<3U>::GetFaults( std::vector<std::string>& );
template void EclipseModel<3U>::GetFaults( std::list<std::string>& );
template void EclipseModel<3U>::GetFaults( std::set<std::string>& );


template<size_t dim>
template<class Container>
void EclipseModel<dim>
::GetWells( Container& data )
{
    data.clear();
    typename Container::iterator dit = data.begin();
    std::copy( wells_.begin(), wells_.end(), std::inserter( data, dit ) );
}

template void EclipseModel<1U>::GetWells( std::vector<std::string>& );
template void EclipseModel<1U>::GetWells( std::list<std::string>& );
template void EclipseModel<1U>::GetWells( std::set<std::string>& );

template void EclipseModel<2U>::GetWells( std::vector<std::string>& );
template void EclipseModel<2U>::GetWells( std::list<std::string>& );
template void EclipseModel<2U>::GetWells( std::set<std::string>& );

template void EclipseModel<3U>::GetWells( std::vector<std::string>& );
template void EclipseModel<3U>::GetWells( std::list<std::string>& );
template void EclipseModel<3U>::GetWells( std::set<std::string>& );

/// processing special regions
template<size_t dim>
void EclipseModel<dim>
::CreateBoundariesAroundFaults( bool keep_fault_regions )
{
    this->MergeRegions(faults_,"FAULTS");
    faults_.insert("FAULTS");
    this->InsertBoundary("FAULTS",csmp::IRREGULAR,keep_fault_regions);

    // create boundaries
    //for( std::set<std::string>::const_iterator
    //     rit = faults_.begin(); rit != faults_.end(); ++rit )
    //    this->InsertBoundary( (*rit).c_str(), csmp::IRREGULAR, keep_fault_regions );
}

template<size_t dim>
void EclipseModel<dim>
::CreateSplitBoundariesAroundFaults( bool delete_fault_regions )
{
    this->MergeRegions(faults_,"FAULTS");
    faults_.insert("FAULTS");
    this->InsertSplitBoundary("FAULTS",delete_fault_regions);

    // create splitboundaries
    //for( std::set<std::string>::const_iterator
    //     rit = faults_.begin(); rit != faults_.end(); ++rit )
    //    this->InsertSplitBoundary( (*rit), delete_fault_regions );
    //faults_.insert("FAULTS");

}


template class EclipseModel<1U>;
template class EclipseModel<2U>;
template class EclipseModel<3U>;

} // end csmp
