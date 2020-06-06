//
//  BoundaryInterface_Test.cpp
//
//  Created by Stephan Matthai on 23/03/2016.
//  Copyright © 2016 Stephan Matthai. All rights reserved.
//

#include <tuple>
#include "BoundaryInterface_Test.h"
#include "Boundary.h"
#include "Region.h"
#include "ErrorHandler.h"
#include "ANSYS_Model3D.h"
#include "VTK_Interface.h"
#include "VTU_Interface.h"
#include "FaceConstructionData.h"


using namespace std;

namespace csmp {

size_t labelRegionPatches( Model<3U>&, const char* dim_1_region, const char* diagnostic_elmt_variable,
                          const char* patch_variable, const std::vector<std::string>& region_names );


/**
     Testing 2016-17 refactored functionality of BoundaryInterface and Boundary classes:
     
     CreateInternalBoundaryFrom()
*/
void BoundaryInterface_Test::run()
  {
      TestBoxShapedModel();
  
      string  input_file("fault_boundary_test");	  
      //input_file = "hex2_s"; // hexahedral element test model
      //input_file = "prism_test"; // hexahedral element test model	  
      const size_t dim(3);

      // ------------------------------
      // 1. Building the CSMP Model
      // ------------------------------
      const bool irregular_mesh(true);
      const bool binary_file(true);
      const bool use_regions_file(true);
      const bool create_boundaries(true);

      ANSYS_Model3D model( input_file.c_str(), "CSMP-variables.txt", irregular_mesh, binary_file, use_regions_file, create_boundaries );
      printModelDimensions(model, true);

      /// assuming a dim-1 region, label and count material juxtaposition relationships
      const string    region_tag("region identifier");
      vector<string>  region_names;
      size_t regions = model.CountAndLabelRegions( region_tag.c_str(), region_names );
      //out( region_names );

      VTU_Interface<dim>  vtu(model);
      if ( verbose_ ) vtu.OutputDataToVTU( "BoundaryInterface_Test_", region_tag, string("Model"), 0 );

      const string patch_tag("patch identifier");
      // creating visual output that illustrates what the boundary should look like for testing
      size_t subregions = labelRegionPatches( model, "NORMAL_FAULT", region_tag.c_str(), patch_tag.c_str(), region_names );

      if ( verbose_ ) cout <<"\nBoundaryInterface_Test::run: identified "<< subregions <<" region patches in region NORMAL_FAULT touching "<< regions <<" model regions.\n";
      if ( verbose_ ) vtu.OutputDataToVTU( "test", patch_tag, string("NORMAL_FAULT"), 0 );
    
	  VTU_Interface<3> vtu_boundary(model);
	  {
		  std::string boundary_name("LEFT");
		  std::string variableName("face variable");
		  Boundary<dim>& boundary(model.Boundary(boundary_name));
		  Index areaKey(model.Database().StorageKey(variableName.c_str()));
		  assert(areaKey.place == FACE); // this needs to be done for Face::Read bc Read is inherited from Element and can
										 // either be Face, Element or InterFace
		  ScalarVariable area(PLAIN, 0.);
		  size_t surfaceElementCount(0);
		  const typename std::vector<Face<dim>*>::iterator domainElementsEnd(boundary.ElementsEnd());
		  for (typename std::vector<Face<dim>*>::iterator it = boundary.ElementsBegin(); it != domainElementsEnd; ++it)
		  {
			  area = (*it)->Area();
			  (*it)->Store(areaKey, area);
			  ++surfaceElementCount;
		  }
		  if(model.ContainsBoundary(boundary_name))
			vtu.OutputDataToVTU("test", variableName, boundary, 0);
	  }	  
	  {
		  std::string boundary_name("RIGHT");
		  std::string variableName("face variable");
		  Boundary<dim>& boundary(model.Boundary(boundary_name));
		  Index areaKey(model.Database().StorageKey(variableName.c_str()));
		  assert(areaKey.place == FACE); // this needs to be done for Face::Read bc Read is inherited from Element and can
										 // either be Face, Element or InterFace
		  ScalarVariable area(PLAIN, 0.);
		  size_t surfaceElementCount(0);
		  const typename std::vector<Face<dim>*>::iterator domainElementsEnd(boundary.ElementsEnd());
		  for (typename std::vector<Face<dim>*>::iterator it = boundary.ElementsBegin(); it != domainElementsEnd; ++it)
		  {
			  area = (*it)->Area();
			  (*it)->Store(areaKey, area);
			  ++surfaceElementCount;
		  }
		  if (model.ContainsBoundary(boundary_name))
			vtu.OutputDataToVTU("test", variableName, boundary, 0);
	  }
	  {
		  std::string boundary_name("FRONT");
		  std::string variableName("face variable");
		  Boundary<dim>& boundary(model.Boundary(boundary_name));
		  Index areaKey(model.Database().StorageKey(variableName.c_str()));
		  assert(areaKey.place == FACE); // this needs to be done for Face::Read bc Read is inherited from Element and can
										 // either be Face, Element or InterFace
		  ScalarVariable area(PLAIN, 0.);
		  size_t surfaceElementCount(0);
		  const typename std::vector<Face<dim>*>::iterator domainElementsEnd(boundary.ElementsEnd());
		  for (typename std::vector<Face<dim>*>::iterator it = boundary.ElementsBegin(); it != domainElementsEnd; ++it)
		  {
			  area = (*it)->Area();
			  (*it)->Store(areaKey, area);
			  ++surfaceElementCount;
		  }
		  if (model.ContainsBoundary(boundary_name))
			vtu.OutputDataToVTU("test", variableName, boundary, 0);
	  }
	  {
		  std::string boundary_name("BACK");
		  std::string variableName("face variable");
		  Boundary<dim>& boundary(model.Boundary(boundary_name));
		  Index areaKey(model.Database().StorageKey(variableName.c_str()));
		  assert(areaKey.place == FACE); // this needs to be done for Face::Read bc Read is inherited from Element and can
										 // either be Face, Element or InterFace
		  ScalarVariable area(PLAIN, 0.);
		  size_t surfaceElementCount(0);
		  const typename std::vector<Face<dim>*>::iterator domainElementsEnd(boundary.ElementsEnd());
		  for (typename std::vector<Face<dim>*>::iterator it = boundary.ElementsBegin(); it != domainElementsEnd; ++it)
		  {
			  area = (*it)->Area();
			  (*it)->Store(areaKey, area);
			  ++surfaceElementCount;
		  }
		  if (model.ContainsBoundary(boundary_name))
			vtu.OutputDataToVTU("test", variableName, boundary, 0);
	  }
	  {
		  std::string boundary_name("TOP");
		  std::string variableName("face variable");
		  Boundary<dim>& boundary(model.Boundary(boundary_name));
		  Index areaKey(model.Database().StorageKey(variableName.c_str()));
		  assert(areaKey.place == FACE); // this needs to be done for Face::Read bc Read is inherited from Element and can
										 // either be Face, Element or InterFace
		  ScalarVariable area(PLAIN, 0.);
		  size_t surfaceElementCount(0);
		  const typename std::vector<Face<dim>*>::iterator domainElementsEnd(boundary.ElementsEnd());
		  for (typename std::vector<Face<dim>*>::iterator it = boundary.ElementsBegin(); it != domainElementsEnd; ++it)
		  {
			  area = (*it)->Area();
			  (*it)->Store(areaKey, area);
			  ++surfaceElementCount;
		  }
		  if (model.ContainsBoundary(boundary_name))
			vtu.OutputDataToVTU("test", variableName, boundary, 0);
	  }
	  {
		  std::string boundary_name("BOTTOM");		  
		  std::string variableName("face variable");
		  Boundary<dim>& boundary(model.Boundary(boundary_name));
		  Index areaKey(model.Database().StorageKey(variableName.c_str()));
		  assert(areaKey.place == FACE); // this needs to be done for Face::Read bc Read is inherited from Element and can
										 // either be Face, Element or InterFace
		  ScalarVariable area(PLAIN, 0.);
		  size_t surfaceElementCount(0);
		  const typename std::vector<Face<dim>*>::iterator domainElementsEnd(boundary.ElementsEnd());
		  for (typename std::vector<Face<dim>*>::iterator it = boundary.ElementsBegin(); it != domainElementsEnd; ++it)
		  {
			  area = (*it)->Area();
			  (*it)->Store(areaKey, area);
			  ++surfaceElementCount;
		  }
		  if (model.ContainsBoundary(boundary_name))
			vtu.OutputDataToVTU("test", variableName, boundary, 0);
	  }

	  

      /// discerning patches by values for the region in terms of the diagnostic element variable
      // -------------------------------------------------------------
      // 2. Creating a Boundary from an internal region "NORMAL_FAULT"
      // -------------------------------------------------------------
      const string test_region("NORMAL_FAULT");
      Region<3U>&  test_subdomain(model.Region(test_region.c_str()));
      const size_t elmts_original_region(test_subdomain.Elements());
      // will remove the original region
      const size_t model_faces_before(model.Mesh().Faces());
      std::pair<std::set<std::string>,bool> boundaries = model.CreateInternalBoundaryFrom( "NORMAL_FAULT" );
      _test( boundaries.first.size() == subregions );
      _test( model.ContainsRegion("NORMAL_FAULT") == false );
      const size_t model_faces_after(model.Mesh().Faces());
      _test( (model_faces_after - model_faces_before) == elmts_original_region );

      // testing for existance of the new boundary patches
      if ( verbose_ ) {
           cout << "\nBoundaryInterface_Test::run: Printing the name of the boundaries in the model:";
           for ( auto it = model.BoundariesBegin(); it != model.BoundariesEnd(); ++it )
             cout << "\n\tBoundary: " << (*it).first <<" ("<< (*it).second.Elements() <<" faces)";
        }
      _test( model.Boundaries() == 12 );

      // creating property values on the boundaries and outputting these to VTU
      model.InputPropertyValue("nodal variable", makeScalar(PLAIN, 0.) );
      double64  bvalue(1.3e5);
      for ( auto it = model.BoundariesBegin(); it != model.BoundariesEnd(); ++it ) {
           if ( verbose_ ) cout << "\n Boundary: " << (*it).first;
           (*it).second.InputPropertyValue("nodal variable", makeScalar(PLAIN, bvalue) );
           // testing VTU output
           if ( verbose_ ) vtu.OutputDataToVTU( "patch", string("nodal variable"), (*it).second, 0 );
           // creating different pressure values for each boundary patch
           bvalue += 1.0e5;
        }
      if ( verbose_ ) cout << endl;

// TODO: JC: check this!!!
      // testing whether boundary segments can be found by combined search criteria
      //const set<string> intersected_regions{ "BOUNDARY", "BOTTOM", "TOP" };
      //string patch_name = findBoundary( model, intersected_regions );	  
      //_test( patch_name == "NORMAL_FAULT_BOUNDARY3_LAYER_BOTTOM_LAYER_TOP" );
      const set<string> search_strings{ "BOUNDARY", "NORMAL", "FAULT" };
      set<string> region_patches_found;
      const size_t patches_found = model.FindBoundaryNames( search_strings, region_patches_found );
      // the boundary was decomposed into 6 patches
      _test( patches_found == 6 );
    
      // -----------------------------------
      // 3. testing supporting functionality
      // -----------------------------------
      _test( TestRegionContactDetection(model) );
        
  } // end run
  
  
  
/**
   Tests:
    
     1. whether Box-shaped model creation be accomplished,
     
     2. whether the regions and their perimeters are correctly constructured
     
     3. how fast model creation is
*/
void BoundaryInterface_Test::TestBoxShapedModel()
 {
      // ----------------------------------------------------
      // 1. Building a CSMP Model from an ANSYS input dataset
      // ----------------------------------------------------
      string  input_file("cube_flag");
   
      const bool irregular_mesh(true);
      const bool binary_file(true);
      const bool use_regions_file(true);
      const bool create_boundaries(true);

      ANSYS_Model3D model( input_file.c_str(), "CSMP-variables.txt", irregular_mesh, binary_file, use_regions_file, create_boundaries );

 } // end TestBoxShapedModel





/**
   checks whether the contact surface of 2 contacting regions is recovered correctly
   
   For the model 'fault_boundary_test'  this test looks at whether the horizons juxtaposed 
   by the fault are correctly detected.
   
   TODO: develop this into a test for CreateBetween()
*/
bool BoundaryInterface_Test::TestRegionContactDetection( const Model<3U>& model )
 {
    assert( string(model.Name()) == "fault_boundary_test" );
    bool all_tests_passed(true);
   
    vector<tuple<Element<3U>*,Element<3U>*,size_t,size_t> > shared;
   
    if ( model.SharedPerimeterFaces( "LAYER_TOP", "LAYER_RESERVOIR", shared )  == 0 )
      all_tests_passed = false;

    if ( model.SharedPerimeterFaces( "LAYER_TOP", "LAYER_BOTTOM", shared )  == 0 )
      all_tests_passed = false;

    if ( model.SharedPerimeterFaces( "LAYER_BOTTOM", "LAYER_RESERVOIR", shared )  == 0 )
      all_tests_passed = false;
   
    return all_tests_passed;
   
 } // end TestRegionContactDetection









 
  
  
// TRANSFERRED
/**
    numbers the unique regions of the models, labeling their elements with the region number
    as "region identifier".
    
    @parameter vector of region names to retrieve them from the integer keys
*/
template<size_t dim>
size_t countAndLabelRegions( Model<dim>& model, const char* region_identifier, std::vector<std::string>& region_names )
 {
    // setting element variable up to identify all unique regions - uniquely
    const string rvariable(region_identifier);
    if ( !model.Database().IsDefined(rvariable.c_str()) )
       model.CreateProperty( rvariable.c_str(), "X", SCALAR, ELEMENT );
   
    // counting the regions and initialising them with the unique identifiers
    region_names.resize(model.UniqueRegions());
    size_t regions(0);
    for ( auto it=model.UniqueRegionsBegin(); it!=model.UniqueRegionsEnd(); it++ ) {
         (*it).second.InputPropertyValue( rvariable.c_str(), makeScalar(PLAIN,static_cast<double64>(regions)) );
         region_names[regions] = (*it).first;
         regions++;
      }
   
   assert( regions == region_names.size() );
   return regions;
   
 } // end countAndLabelRegions

template size_t countAndLabelRegions( Model<1U>&, const char*, std::vector<std::string>& );
template size_t countAndLabelRegions( Model<2U>&, const char*, std::vector<std::string>& );
template size_t countAndLabelRegions( Model<3U>&, const char*, std::vector<std::string>& );



//TRANSFERRED
/**
    loops over the surface cells of the model subdomain, 
    checking whether any of the projections of the normals of the neighbor 
    cells onto the normal of the current element are negative.
    
    @return if any of the projections is negative, method returns false.
    A negative projection will mean that the normals of the neighboring cells 
    are at > to 90^o to the current cell; this should not be the case
    unless the surface has cusps with in it.
    
    @return if the region does not consist entirely of surface elements
    an error is reported.
    
    @attention even if the normals have a consistent orientation, this method
    may return false if the surface contains a cusp (>90^o kink).
 
*/
bool checkNeighborNormalsForConsistentOrientation( const Region<3U>&  subdomain )
 {
    vector<double64> normal(3U), nbor_normal(3U);
   
    size_t non_surface_elements(0U);
    for ( auto it=subdomain.ElementsBegin(); it!=subdomain.ElementsEnd(); ++it )
      // this method only  considers surface elements
      if ( (*it)->IsSurfaceElement() ) {
           (*it)->UnitNormal( normal );
           const size_t neighbors((*it)->Neighbors());
           for ( size_t i=0U; i<neighbors; ++i )
             // only valid neighbor elements are considered
             if ( (*it)->Neighbor(i) != nullptr ) {
                  (*it)->Neighbor(i)->UnitNormal( nbor_normal );
                  // projection
                  double64 result(0.);
                  for ( size_t j=0U; j<3U; ++ j )
                    result += normal[j] * nbor_normal[j];
                  if ( result < 0. )
                    return false;
               }
         }
       else non_surface_elements++;
   
    if ( non_surface_elements > 0U )
      ErrorHandler::Instance().notice( ERROR, "checkNeighborNormalsForConsistentOrientation:",
                                       subdomain.Name(), "region contained not only surface elements." );
    return true;
   
 } // end checkNeighborNormalsForConsistentOrientation




// TRANSFERRED
/**
     creates underscore-separated unique names for the region patches based on the juxtapositions relationships
     across the lower-dimensional regions, the names are composed of:
 
     1. the name of the master region
     
     2. "BOUNDARY"
     
     3. the patch identifier number attached to boundary
     
     3. the name of the inner region, i.e. the region that the lower-dimensional element normals point away from
     
     4. the name of the outer region, i.e. that into which the normals point
     
     @attention  where the boundary just intersects a layer (same material on either side), the layer name appears
     only once. The second instance is replaced by INTERSECTION.
*/
std::string internalBoundaryNameFrom( const FaceConstructionData& fdata, const std::vector<std::string>& region_names )
 {
     assert( fdata.ElementMaterial() < region_names.size() );
     std::string boundary_name( region_names[ fdata.ElementMaterial() ] );
     boundary_name += "_BOUNDARY";
     boundary_name += to_string(fdata.PatchNumber());
     boundary_name += '_';
     pair<long,long> materials(fdata.Materials());
     assert( materials.first  < region_names.size() );
     assert( materials.second < region_names.size() );
     boundary_name += region_names[ materials.first ];
     boundary_name += '_';
     if ( materials.first == materials.second ) boundary_name +="INTERSECTION";
     else boundary_name += region_names[ materials.second ];
   
     return boundary_name;
 }



// TRANSFERRED
/**
     Returns name of the first Boundary (patch) that contains the supplied
     region name strings or any other other strings in an arbitrary order.
     
     if the boundary does not exist, method returns \O string.
     
     @author SKM 
     @date March 2016
*/
std::string  findBoundary( const Model<3U>& model, const set<string>& intersected_regions )
 {
    // if the substring set is empty
    if ( intersected_regions.empty() ) {
         ErrorHandler::Instance().notice( WARNING, "findBoundary:", "supplied set of substrings is empty; returning '\0'." );
         return std::string("\0");
      }
    // if the model has no boundaries
    if ( model.Boundaries() == 0 ) {
         ErrorHandler::Instance().notice( WARNING, "findBoundary:", "model has no boundaries; returning '\0'." );
         return std::string("\0");
      }
     // making a set of boundary names
     const size_t substrings_used_in_search(intersected_regions.size());
     for ( auto it=model.BoundariesBegin(); it!=model.BoundariesEnd(); ++it ) {
          size_t substrings_found(0U);
          for ( auto ir=intersected_regions.begin(); ir!=intersected_regions.end(); ++ir )
            // if the substring is found
            if ( (*it).first.find(*ir) !=std::string::npos ) substrings_found++;
            // when all substrings are contained in the boundary name, it is returned
          if ( substrings_found == substrings_used_in_search )
            return (*it).first;
       }
    return std::string("\0");
 }



// TRANSFERRED
/**
    Searches the model for boundaries the name of which contains the search strings
    provided via the first set. The results are returnd into the second set.
    
    @note Use this method, for example, to retrieve multiple boundary patches that were generated
    from a single lower-dimensional regon, like a fault surface.

     @author SKM 
     @date March 2016
*/
size_t  findBoundaries( const Model<3U>& model, const set<string>& intersected_regions,
                        set<string>& region_patches_found )
 {
    // if the substring set is empty
    if ( intersected_regions.empty() ) {
         ErrorHandler::Instance().notice( WARNING, "findBoundaries:", "supplied set of substrings is empty; returning '\0'." );
         return 0U;
      }
    // if the model has no boundaries
    if ( model.Boundaries() == 0 ) {
         ErrorHandler::Instance().notice( WARNING, "findBoundaries:", "model has no boundaries; returning '\0'." );
         return 0U;
      }
    region_patches_found.clear();
   
     // making a set of boundary names
     const size_t substrings_used_in_search(intersected_regions.size());
     for ( auto it=model.BoundariesBegin(); it!=model.BoundariesEnd(); ++it ) {
          size_t substrings_found(0U);
          for ( auto ir=intersected_regions.begin(); ir!=intersected_regions.end(); ++ir )
            // if the substring is found
            if ( (*it).first.find(*ir) !=std::string::npos ) substrings_found++;
          // when all substrings are contained in the boundary name, it is returned
          if ( substrings_found == substrings_used_in_search )
            region_patches_found.insert( (*it).first );
       }
    // return how many region patches contain the search string(s)
    return region_patches_found.size();
 }







// TRANSFERRED
/**
     higherDimensionalNeighbors()
 
     - finds the IDs of the higher dimensional neibors of the current element
 
     - identifies which of the neighbors is on the inside and which on the outside
       as indicated by the normal direction of the lower dimensional element
 
     - identifies the materials on either side
     
     @return two inside-outside pairs of element idx numbers and corresponding materials on either side
     all the data are stored in the returnd FaceConstructionData object.
     
     assumptions
     - assumes that the nodes and elements in the entire model domain are numbered continuously
     
     application
     - use this function for finding neighbors of a surface element that sits on the inside of another region
     
     TODO: output the numbers of matching faces of the discovered the elements so that they can later be connected
 
*/
FaceConstructionData  higherDimensionalNeighbors( const Element<3U>& e, const csmp::Index& mtrl_key )
 {
     assert( e.IsSurfaceElement() );

     // 1. looping over the parent elements of the nodes searching for the faces which are shared with the lower dimensional element
     // -----------------------------------------------------------------------------------------------------------------------------
     // making a set of element nodes to later identify faces by comparison
     set<size_t>   node_set, test_set;
     const size_t  nodes(e.Nodes());
     for ( size_t i=0U; i<nodes; ++i ) node_set.insert(e.N(i)->Idx());
     map<const Element<3U>*,size_t>  nbor_elmts;
     vector<size_t> fnids;
     for ( size_t i=0U; i<nodes; i++ ) {
          const size_t parents(e.N(i)->Parents());
          for ( size_t j=0U; j<parents; ++j ) {
               const Element<3U>* const eptr(e.N(i)->Parent(j));
               const size_t faces(eptr->Faces());
               for ( size_t k=0U; k<faces; ++k ) {
                     eptr->FE()->NodesOfFace( k, fnids );
                     size_t fnodes(fnids.size());
                     for ( size_t l=0U; l<fnodes; ++l )
                       test_set.insert( eptr->N( fnids[l])->Idx() );
                     // if the face is shared the element and its face are recorded
                     if ( node_set == test_set ) {
                          // storing a pointer to this element and its local face number
                          // making sure that no duplicate is received
                          nbor_elmts.insert( make_pair(eptr,k) );
                       }
                     test_set.clear();
                 }
            }
       }
    assert( nbor_elmts.size() == 2U );
   
    // 2. finding which of the neighbors is the inside one by projecting face normals onto lower dim element normal
    // -------------------------------------------------------------------------------------------------------------
    vector<double64>  enrml, fnrml;
    e.UnitNormal( enrml );
    map<const Element<3U>*,size_t>::const_iterator  nbit(nbor_elmts.begin());
    bool inside_elmt_found(false);
    bool outside_elmt_found(false);
    pair<size_t,size_t> nbors;
    pair<size_t,size_t> faces;
    pair<long,long>     materials;

    // first element
    // -------------
    assert( (*nbit).first != nullptr );
    faces.first = (*nbit).second;
    (*nbit).first->UnitNormalToFace( faces.first, fnrml );
    double64 dotproduct = enrml[0] * fnrml[0] + enrml[1] * fnrml[1] + enrml[2] * fnrml[2];
   
    // if the projection is negative, the first element lies on the outside
    if ( dotproduct < 0. ) {
         nbors.second       = (*nbit).first->Idx();
         materials.second   = static_cast<long>((*nbit).first->Read( mtrl_key ));
         outside_elmt_found = true;
      }
    else {
         nbors.first        = (*nbit).first->Idx();
         materials.first    = static_cast<long>((*nbit).first->Read( mtrl_key ));
         inside_elmt_found  = true;
      }
    nbit++;


// TODO: SKM: this could be done without 2 normal projections for speedup

    // second element
    // --------------
    assert( (*nbit).first != nullptr );
    faces.second = (*nbit).second;
    (*nbit).first->UnitNormalToFace( faces.second, fnrml );
    dotproduct = enrml[0] * fnrml[0] + enrml[1] * fnrml[1] + enrml[2] * fnrml[2];

   // if the projection is negative, the second element lies on the outside
    if ( dotproduct < 0. ) {
         // checking that we have no duplication here
         assert( outside_elmt_found == false );
         nbors.second     = (*nbit).first->Idx();
         materials.second = static_cast<long>((*nbit).first->Read( mtrl_key ));
      }
    else {
         assert( inside_elmt_found == false );
         nbors.first      = (*nbit).first->Idx();
         materials.first  = static_cast<long>((*nbit).first->Read( mtrl_key ));
      }
   
    // initialise with nbors, their faces, adjacent materials, and patch numbers
    return FaceConstructionData( e.Idx(), nbors,  faces, materials,
                                 static_cast<long>(e.Read( mtrl_key)) );
   
 } // end higherDimensionalNeighbors






/**
    Assuming a dim-1 = 2 dimensional region, this function,
    loops over lower dimensional region and records and labels the different intersection configurations
*/
size_t  labelRegionPatches( Model<3U>& model, const char* dim_1_region, const char* diagnostic_elmt_variable,
                            const char* patch_variable, const std::vector<std::string>& region_names )
 {
    // 1. verification of input to function
    // ------------------------------------
    if ( model.ContainsRegion(dim_1_region) == false ) {
          ErrorHandler::Instance().notice( ERROR, "labelRegionPatches:", dim_1_region, "does not exist; nothing was done." );
          return 0;
      }
    if ( model.UniqueRegions() <= 1 ) {
          ErrorHandler::Instance().notice( WARNING, "labelRegionPatches:", "model contains only a single unique region; so there is only one patch." );
          return 1;
      }
    // verifying that we are indeed dealing with a region of surface elements only
    Region<3U>&  subdomain(model.Region(dim_1_region));
    pair<int32,int32>  dimensionality = subdomain.ElementSpatialDimensions();
    //   number of dims in region      dimension of contained elements
    if ( dimensionality.first != 1 and dimensionality.second != 2 ) {
         ErrorHandler::Instance().notice( ERROR, "labelRegionPatches:", dim_1_region, "region does not consist of surface elements only; nothing was done." );
         return 0;
      }
    // verifying that the region lies inside of the model
    size_t boundary_elements(0);
    for ( auto eit=subdomain.ElementsBegin(); eit!=subdomain.ElementsEnd(); ++eit )
      {
         // TODO: improve these diagnostics
         if ( (*eit)->AtBoundary() != NOT and
              (*eit)->AtBoundary() != INTERNAL and
              (*eit)->AtBoundary() != IRREGULAR ) {
              cerr << parseBoundary( (*eit)->AtBoundary() ) <<" ";
              boundary_elements++;
           }
      }
    if ( boundary_elements > 0 ) {
         ErrorHandler::Instance().notice( ERROR, "labelRegionPatches:", dim_1_region, "region appears to lie at the model boundary; nothing was done." );
         return 0;
      }
    if ( !model.Database().IsDefined(diagnostic_elmt_variable) ) {
         ErrorHandler::Instance().notice( ERROR, "labelRegionPatches:", diagnostic_elmt_variable, "variable to discern regions is not defined; nothing was done." );
         return 0;
      }
      {  // check whether there are multiple region identifiers
         double64 rmin, rmax;
         model.MinMaxOf( diagnostic_elmt_variable, rmin, rmax );
         if ( fabs(rmax - rmin) <= numeric_limits<double64>::epsilon() ) {
              ErrorHandler::Instance().notice( WARNING, "labelRegionPatches:", diagnostic_elmt_variable, "is single valued; so there is only one patch." );
              return 1;
           }
      }

    // 2. setting up an element variable to show the partitions of the lower dimensional region
    // ----------------------------------------------------------------------------------------
    csmp::Index mtrl_key = model.Database().StorageKey(diagnostic_elmt_variable);
   
    if ( !model.Database().IsDefined(patch_variable) )
      model.CreateProperty( patch_variable, "X", SCALAR, ELEMENT );
    // patch-discerning variable
    csmp::Index pvar_key = model.Database().StorageKey(patch_variable);
   
    // renumbering nodes and elements
    model.Region("Model").UpdateMemberIndexes();
   
    // looping over the region, identifying and recording the juxtaposition relationships
    map<pair<long,long>,size_t>   patches;
    vector<FaceConstructionData>  face_construction_data;
    map<size_t,string>            patch_names;
    string                        patch_name;
    size_t                        n_juxtapositions(0);

    // 3. looping over lower dimensional region identifying juxtaposition relationships
    // ----------------------------------------------------------------------------------------
    face_construction_data.reserve(subdomain.Elements());
    for ( auto eit=subdomain.ElementsBegin(); eit!=subdomain.ElementsEnd(); ++eit )
      {
          // 3.1 identifying neighbors, facing relations, and juxtaposed materials for current element
          FaceConstructionData  fdata(higherDimensionalNeighbors( *(*eit), mtrl_key ));
        
          // 3.2 recording which category of juxtaposition element fall into, naming it and assigning a patch number
          pair<map<pair<long,long>,size_t>::iterator,bool>  it=patches.insert( make_pair(fdata.Materials(),n_juxtapositions) );
          // incrementing number of juxtapositions and corresponding patch names
          if ( it.second == true ) {
               fdata.PatchNumber( (*it.first).second );
               patch_name = internalBoundaryNameFrom( fdata, region_names );
               cout <<"\n\nlabelRegionPatches: found new patch: "<< patch_name <<"\n";
               patch_names.insert( make_pair(n_juxtapositions,patch_name) );
               n_juxtapositions++;
            }
          fdata.PatchNumber( (*it.first).second );

          // 3.3 recording the patch number as variable value
          (*eit)->Store( pvar_key, makeScalar(PLAIN,static_cast<double64>((*it.first).second)) );
        
          // 3.4 recording the data for the element that will later be used to construct the face from
          face_construction_data.push_back( fdata );
      }
   
    // 4. creating labeled boundary patches from the face-defining data
    // ----------------------------------------------------------------
    // 4.1 making a map 'patch_numbers' from 'patch_names' to search for patch identifiers
    map<string,size_t>  patch_numbers;
    for ( auto it=patch_names.begin(); it!=patch_names.end(); ++it )
      patch_numbers.insert( make_pair( (*it).second, (*it).first ) );
   
    // 4.2 building new map where the patch faces are organised by patch names
    map<string,vector<FaceConstructionData> > patch_simplexes;
    vector<FaceConstructionData>              empty_vec;
    for ( auto it=patch_names.begin(); it!=patch_names.end(); ++it )
      patch_simplexes.insert( make_pair( (*it).second, empty_vec ) );
   
    // 4.3 inserting the patch identifiers into it
    for ( auto pit=patch_simplexes.begin(); pit!=patch_simplexes.end(); ++pit )
      {
         assert( patch_numbers.find((*pit).first) != patch_numbers.end() );
         const size_t patch_number((*patch_numbers.find((*pit).first)).second);
         // reserving storage
         (*pit).second.reserve(face_construction_data.size());
         // looping over all face data assigning the ones that are suitable
         for ( auto it=face_construction_data.begin(); it!=face_construction_data.end(); ++it )
           if ( (*it).PatchNumber() == patch_number )
             (*pit).second.push_back( (*it) );
      }
    // 4.4 trimming excess storage of the face-data vectors
    for ( auto pit=patch_simplexes.begin(); pit!=patch_simplexes.end(); ++pit )
      //std::vector<size_t>(elementIds).swap(elementIds);
      vector<FaceConstructionData>( (*pit).second ).swap( (*pit).second );

/* TESTING
cout <<"\n\n\n\nBoundary patches and their face data:\n";
for ( auto it=patch_simplexes.begin(); it!=patch_simplexes.end(); ++it) {
     cout <<"\n\n\n"<< (*it).first <<"\n";
     for ( auto fit=(*it).second.begin(); fit!=(*it).second.end(); ++fit )
       (*fit).Out();
  }
*/
    return patches.size();
   
 } // end labelRegionPatches





} // end csmp
