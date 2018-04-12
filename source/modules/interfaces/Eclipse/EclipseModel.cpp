#include "EclipseModel.h"
#include "ModelTime.h"
#include "variableOperations.h"

using namespace std;

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
         const bool non_box_shaped_model( !mesh_topology.BoxShapedModel() );

         // with regions
         if ( !eclipse_model_settings_.regions_.empty() ) {
             // performing a ckeck whether element numbers in the VSet and the model topology match; else something went wrong
             // and user is given the possibility to call the subsequent method or not.
             const bool require_unique_names_for_vol_surf_lines(true);
             const bool correct_orientation_of_surface_elements(false);
             if ( !mesh_topology.CheckTopology( vset, require_unique_names_for_vol_surf_lines,
                                                correct_orientation_of_surface_elements, non_box_shaped_model ) )
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
      
      
         // SKM FIX - retain critical information in the EclipseModel object
         // TODO: deal with other critical data as well
         // ----------------------------------------------------------------
         const CornerPointGrid<dim>& cnr_grid_ref = mesh_interface.GetCornerPointGrid();
      
         // store grid dimensions int the EclipseModel
         grid_dim_I_ = cnr_grid_ref.DimensionI();
         grid_dim_J_ = cnr_grid_ref.DimensionJ();
         grid_dim_K_ = cnr_grid_ref.DimensionK();
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



// trial versions

BOX_BOUNDARY  whichEdge( BOX_BOUNDARY side1, BOX_BOUNDARY side2 )
 {
    // dealing with the case where the 2 flags are the same so that this is not an EDGE
    if ( side1 == side2 ) return side1;
 
    if ( side1 == LEFT and side2 == BOTTOM ) return EDGE1;
    if ( side1 == LEFT and side2 == RIGHT )  return EDGE2;
    if ( side1 == LEFT and side2 == TOP )    return EDGE3;
    if ( side1 == LEFT and side2 == FRONT )  return EDGE4;
   
    if ( side1 == BOTTOM and side2 == FRONT ) return EDGE5;
    if ( side1 == BOTTOM and side2 == LEFT ) return EDGE1;
    if ( side1 == BOTTOM and side2 == BACK ) return EDGE6;
    if ( side1 == BOTTOM and side2 == TOP ) return EDGE7;
   
    if ( side1 == TOP and side2 == FRONT ) return EDGE8;
    if ( side1 == TOP and side2 == LEFT ) return EDGE3;
    if ( side1 == TOP and side2 == BACK ) return EDGE7;
    if ( side1 == TOP and side2 == RIGHT ) return EDGE12;

    if ( side1 == FRONT and side2 == BOTTOM ) return EDGE5;
    if ( side1 == FRONT and side2 == LEFT ) return EDGE4;
    if ( side1 == FRONT and side2 == TOP ) return EDGE8;
    if ( side1 == FRONT and side2 == RIGHT ) return EDGE12;

    if ( side1 == RIGHT and side2 == BOTTOM ) return EDGE9;
    if ( side1 == RIGHT and side2 == BACK ) return EDGE10;
    if ( side1 == RIGHT and side2 == TOP ) return EDGE11;
    if ( side1 == RIGHT and side2 == FRONT ) return EDGE12;
    return NOT;
 }

BOX_BOUNDARY  whichCorner( BOX_BOUNDARY side1, BOX_BOUNDARY side2, BOX_BOUNDARY side3 )
 {
    set<BOX_BOUNDARY> sides({side1,side2,side3});
    auto s3 = (*sides.begin());
    auto s2 = (*(next(sides.begin(),1)));
    auto s1 = (*sides.rbegin());
   
    if ( sides.empty() ) return NOT;
   
    if ( sides.size() == 1U ) return (*sides.begin());
   
    // dealing with duplicates (there will only be 1 or 2 entries so this can only be a corner if the entries are edges)
    if ( sides.size() == 2U ) {
          // 2 edges EDGE12, EDGE11... EDGE1
          if ( side1 == EDGE12 and side2 == EDGE11 ) return CNR8;
          if ( side1 == EDGE12 and side2 == EDGE9 ) return CNR5;
 
          if ( side1 == EDGE11 and side2 == EDGE10 ) return CNR7;
 
          if ( side1 == EDGE10 and side2 == EDGE9 ) return CNR6;
          if ( side1 == EDGE10 and side2 == EDGE6 ) return CNR6;
 
          if ( side1 == EDGE9 and side2 == EDGE6 ) return CNR6;
          if ( side1 == EDGE9 and side2 == EDGE5 ) return CNR5;

          if ( side1 == EDGE8 and side2 == EDGE4 ) return CNR4;
          if ( side1 == EDGE8 and side2 == EDGE3 ) return CNR4;
     
          if ( side1 == EDGE7 and side2 == EDGE3 ) return CNR3;
          if ( side1 == EDGE7 and side2 == EDGE2 ) return CNR3;

          if ( side1 == EDGE6 and side2 == EDGE2 ) return CNR2;
          if ( side1 == EDGE6 and side2 == EDGE1 ) return CNR2;

          if ( side1 == EDGE5 and side2 == EDGE4 ) return CNR1;
          if ( side1 == EDGE5 and side2 == EDGE1 ) return CNR1;

          if ( side1 == EDGE4 and side2 == EDGE3 ) return CNR4;
          if ( side1 == EDGE4 and side2 == EDGE1 ) return CNR1;

          if ( side1 == EDGE3 and side2 == EDGE2 ) return CNR3;
      
          if ( side1 == EDGE2 and side2 == EDGE1 ) return CNR2;

          // not a corner
          if ( side1 == BACK and side2 == TOP ) return EDGE7;
          if ( side1 == BACK and side2 == BOTTOM ) return EDGE6;
          if ( side1 == BACK and side2 == RIGHT ) return EDGE10;
          if ( side1 == BACK and side2 == LEFT ) return EDGE2;

          if ( side1 == FRONT and side2 == TOP ) return EDGE8;
          if ( side1 == FRONT and side2 == BOTTOM ) return EDGE5;
          if ( side1 == FRONT and side2 == RIGHT ) return EDGE12;
          if ( side1 == FRONT and side2 == LEFT ) return EDGE4;

          if ( side1 == TOP and side2 == RIGHT ) return EDGE11;
          if ( side1 == TOP and side2 == LEFT ) return EDGE3;

          if ( side1 == BOTTOM and side2 == LEFT ) return EDGE1;
          if ( side1 == BOTTOM and side2 == RIGHT ) return EDGE9;
      }
   
    // obeying increasing value constraint: BACK, FRONT, TOP, BOTTOM, RIGHT, LEFT
    if ( s1 == FRONT and s2 == BOTTOM and s3 == LEFT ) return CNR1;
    if ( s1 == BACK  and s2 == BOTTOM and s3 == LEFT ) return CNR2;
    if ( s1 == BACK  and s2 == TOP    and s3 == BOTTOM ) return CNR3;
    if ( s1 == FRONT and s2 == TOP    and s3 == LEFT ) return CNR4;

    if ( s1 == FRONT and s2 == BOTTOM and s3 == RIGHT ) return CNR5;
    if ( s1 == BACK  and s2 == BOTTOM and s3 == RIGHT ) return CNR6;
    if ( s1 == BACK  and s2 == TOP    and s3 == RIGHT ) return CNR7;
    if ( s1 == FRONT and s2 == TOP    and s3 == RIGHT ) return CNR5;
   
    return NOT;
 }


/**
    BOX flag nodes and elements of volumetric target region
*/
template<size_t dim>
void EclipseModel<dim>::AssignBoxBoundaryFlagsWherePossible( const char* target_region )
 {
    csmp::ErrorHandler& error_handler(csmp::ErrorHandler::Instance());

    Region<dim>&          domain( this->Region(target_region));
    vector<size_t>        fnids;
    multimap<size_t,pair<BOX_BOUNDARY,Node<dim>*> >  boundary_nodes;
    vector<double64>      nrml, nrml_right, nrml_left, nrml_top, nrml_bottom, nrml_front, nrml_back;
    Box                   box;
    double64              minLength(0.71); // dot-product of 2 unit vectors at an angle >=45 degrees
    BOX_BOUNDARY          bflag(NOT);
   
    box.UnitNormalTo( BOTTOM, dim, nrml_bottom );
    box.UnitNormalTo( TOP,    dim, nrml_top );
    box.UnitNormalTo( LEFT,   dim, nrml_left );
    box.UnitNormalTo( RIGHT,  dim, nrml_right );
    box.UnitNormalTo( FRONT,  dim, nrml_front );
    box.UnitNormalTo( BACK,   dim, nrml_back );
   
    cout <<"\nAssignBoxBoundaryFlagsWherePossible: scanning hexahedral elements for boundary adffiliation...";
    for ( auto it=domain.ElementsBegin(); it!=domain.ElementsEnd(); ++it )
      {
         // ignore elements that are not hexahedra
         if ( (*it)->FE_Type() != ISOPARAMETRIC_LINEAR_HEXAHEDRON ) {
              cout <<"\n\tignored: "<< parseFiniteElementType( (*it)->FE_Type() );
              continue;
           }
         // idea: loop over the faces of the cell and where there is no neighbor
         // check where the face is facing, assign boundary flags accordingly
         // if the element has more than one face idenfify edges and corners
         for ( size_t i=0U; i<(*it)->Faces(); ++i )
           if ( (*it)->Neighbor(i) == nullptr ) {
                // determining in which direction the face normal points
                (*it)->UnitNormalToFace( i, nrml );
                // projecting: perfect alignment would give dot-product equal 1, inclinations up to 37 degrees cos(37)~0.8 are tolerated
                if      ( dotProduct<dim>(nrml,nrml_bottom) >= minLength ) bflag = BOTTOM;
                else if ( dotProduct<dim>(nrml,nrml_top)    >= minLength ) bflag = TOP;
                else if ( dotProduct<dim>(nrml,nrml_left)   >= minLength ) bflag = LEFT;
                else if ( dotProduct<dim>(nrml,nrml_right)  >= minLength ) bflag = RIGHT;
                else if ( dotProduct<dim>(nrml,nrml_front)  >= minLength ) bflag = FRONT;
                else if ( dotProduct<dim>(nrml,nrml_back)   >= minLength ) bflag = BACK;
                // getting the nodes for flagging the faces
                (*it)->FE()->NodesOfFace( i, fnids );
                for ( auto j=0U; j<fnids.size(); ++j ) {
                    (*it)->N(fnids[j])->AtBoundary( bflag );
                    // storing the nodes to determine which ones lie on EDGES (duplicates) or even corners (triplicates)
                    //                           local node #              flag   local node #
                    boundary_nodes.insert( make_pair( fnids[j], make_pair( bflag, (*it)->N(fnids[j]) ) ) );
                 }
             }
         // flagging duplicate and triplicate nodes accordingly
         for ( size_t n=0U; n<(*it)->Nodes(); ++n ) {
              // dealing with any cases where there are multiple boundary flags
              // duplicates = edges
              if ( boundary_nodes.count(n) == 2U ) {
                   // figuring out which edge we are on
                   auto range = boundary_nodes.equal_range(n);
                   bflag = whichEdge( (*range.first).second.first, (*range.second).second.first );
                   (*it)->N(n)->AtBoundary( bflag );
                }
              // triplicates = corners
              else if ( boundary_nodes.count(n) == 3U ) {
                   // figuring out which corner we have found
                   auto range = boundary_nodes.equal_range(n);
                   auto it_2nd(range.first); it_2nd++;
                   bflag = whichCorner( (*range.first).second.first, (*it_2nd).second.first, (*range.second).second.first );
                   (*it)->N(n)->AtBoundary( bflag );
               }
              else {
                  (*it)->Out();
                  cerr <<"\n\tdetected "<< boundary_nodes.count(n) <<" boundary flags for element "<< (*it)->Idx();
                  error_handler.notice( WARNING, "EclipseModel<dim>::AssignBoxBoundaryFlagsWherePossible:",
                                        "this may be a completely disconnected element.");
               }
           }
 
 // testing
 cerr <<"\n("<< (*it)->Idx() <<"): ";
 for ( size_t n=0U; n<(*it)->Nodes(); ++n )
   cerr << parseBoundary( (*it)->N(n)->AtBoundary() ) <<" ";
        
         // resetting
         boundary_nodes.clear();
         bflag = NOT;
      }
 
 } // end AssignBoxBoundaryFlagsWherePossible
 




template class EclipseModel<1U>;
template class EclipseModel<2U>;
template class EclipseModel<3U>;

} // end csmp
