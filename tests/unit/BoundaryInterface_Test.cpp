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
#include "meshManagementUtilities.h"
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
      const uint32_t dim(3);

      // ------------------------------
      // 1. Building the CSMP Model
      // ------------------------------
      const bool binary_file(true);

      ANSYS_Model3D  model( input_file.c_str(), "CSMP-variables.txt", binary_file );
      printModelDimensions(model, true);
      if ( verbose_ ) model.RegionsOut();

      const Region<3>& model_domain = model.Region("Model");
      set<TOPOTYPE> B_flags_model_before = nodeTopologyFlags<3,Element>( model_domain.NodesBegin(), model_domain.NodesEnd() );
      _test( B_flags_model_before.size() == 7 ); // mesh vertex, exterior points, lines & surfaces and interior lines and surfaces
      const Region<3>& fault = model.Region("NORMAL_FAULT");
      set<TOPOTYPE> B_flags_before = nodeTopologyFlags<3,Element>( fault.NodesBegin(), fault.NodesEnd() );
      _test( B_flags_before.size() == 4 ); // INTERIOR_SURFACE and INTERIOR_LINE, PERIMETER_LINE

      // testing VData::InitialiseTopoTypes()
      if ( verbose_ ) {
          TopoTypeToVTU( model, "LAYER_RESERVOIR" );
          TopoTypeToVTU( model, "NORMAL_FAULT" );
        }

      initialise_BREP_TopologyFlags( model ); // function defined in Model
 //     initialise_BREP_TopologyFlags_vs2<3>( model ); // Claude refactor, still gives wrong results

      set<TOPOTYPE> B_flags_after = nodeTopologyFlags<3,Element>( fault.NodesBegin(), fault.NodesEnd() );
      set<TOPOTYPE> B_flags_model_after = nodeTopologyFlags<3,Element>( model_domain.NodesBegin(), model_domain.NodesEnd() );

     // testing initialise_BREP_TopologyFlags()
      if ( verbose_ ) {
          TopoTypeToVTU( model, "LAYER_RESERVOIR" );
          TopoTypeToVTU( model, "NORMAL_FAULT" );
        }

      /// assuming a dim-1 region, label and count material juxtaposition relationships
      const string    region_tag("region identifier");
      vector<string>  region_names;
      size_t regions = model.CountAndLabelUniqueRegions( region_tag.c_str(), region_names );
      //out( region_names );
        
      const string patch_tag("patch identifier");
      // creating visual output that illustrates what the boundary should look like for testing
      size_t subregions = labelRegionPatches( model, "NORMAL_FAULT", region_tag.c_str(), patch_tag.c_str(), region_names );

      if ( verbose_ ) {
          cout <<"\nBoundaryInterface_Test::run: identified "<< subregions <<" region patches in region NORMAL_FAULT touching "<< regions <<" model regions.\n";
          VTU_Interface<dim>  vtu(model);
          vtu.OutputDataToVTU( "BoundaryInterface_Test_", region_tag, string("Model"), 0 );
          vtu.OutputDataToVTU( "test", patch_tag, string("NORMAL_FAULT"), 0 );
        }
	    VTU_Interface<3> vtu(model);
      {
        std::string boundary_name("LEFT");
        std::string variableName("face variable");
        Boundary<dim>& boundary(model.Boundary(boundary_name));
        Index areaKey(model.Database().StorageKey(variableName.c_str()));
        assert(areaKey.place == FACE); // this needs to be done for Face::Read bc Read is inherited from Element and can
                       // either be Face, Element or InterFace
        ScalarVariable area(PLAIN, 0.);
        size_t surfaceElementCount(0);
        const typename std::vector<Face<dim>*>::const_iterator domainElementsEnd(boundary.CellsEnd());
        for (typename std::vector<Face<dim>*>::const_iterator it = boundary.CellsBegin(); it != domainElementsEnd; ++it)
        {
          area = (*it)->Area();
          (*it)->Store(areaKey, area);
          ++surfaceElementCount;
        }
        if (model.ContainsBoundary(boundary_name))
        if ( verbose_ ) vtu.OutputDataToVTU("test", variableName, boundary, 0);
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
        const typename std::vector<Face<dim>*>::const_iterator domainElementsEnd(boundary.CellsEnd());
        for (typename std::vector<Face<dim>*>::const_iterator it = boundary.CellsBegin(); it != domainElementsEnd; ++it)
        {
          area = (*it)->Area();
          (*it)->Store(areaKey, area);
          ++surfaceElementCount;
        }
        if (model.ContainsBoundary(boundary_name))
        if ( verbose_ ) vtu.OutputDataToVTU("test", variableName, boundary, 0);
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
        const typename std::vector<Face<dim>*>::const_iterator domainElementsEnd(boundary.CellsEnd());
        for (typename std::vector<Face<dim>*>::const_iterator it = boundary.CellsBegin(); it != domainElementsEnd; ++it)
        {
          area = (*it)->Area();
          (*it)->Store(areaKey, area);
          ++surfaceElementCount;
        }
        if (model.ContainsBoundary(boundary_name))
        if ( verbose_ ) vtu.OutputDataToVTU("test", variableName, boundary, 0);
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
        const typename std::vector<Face<dim>*>::const_iterator domainElementsEnd(boundary.CellsEnd());
        for (typename std::vector<Face<dim>*>::const_iterator it = boundary.CellsBegin(); it != domainElementsEnd; ++it)
        {
          area = (*it)->Area();
          (*it)->Store(areaKey, area);
          ++surfaceElementCount;
        }
        if (model.ContainsBoundary(boundary_name))
        if ( verbose_ ) vtu.OutputDataToVTU("test", variableName, boundary, 0);
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
        const typename std::vector<Face<dim>*>::const_iterator domainElementsEnd(boundary.CellsEnd());
        for (typename std::vector<Face<dim>*>::const_iterator it = boundary.CellsBegin(); it != domainElementsEnd; ++it)
        {
          area = (*it)->Area();
          (*it)->Store(areaKey, area);
          ++surfaceElementCount;
        }
        if (model.ContainsBoundary(boundary_name))
        if ( verbose_ ) vtu.OutputDataToVTU("test", variableName, boundary, 0);
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
        const typename std::vector<Face<dim>*>::const_iterator domainElementsEnd(boundary.CellsEnd());
        for (typename std::vector<Face<dim>*>::const_iterator it = boundary.CellsBegin(); it != domainElementsEnd; ++it)
        {
          area = (*it)->Area();
          (*it)->Store(areaKey, area);
          ++surfaceElementCount;
        }
        if (model.ContainsBoundary(boundary_name))
        if ( verbose_ ) vtu.OutputDataToVTU("test", variableName, boundary, 0);
        if ( verbose_ ) model.BoundariesOut();
      }

	   if ( verbose_ ) TestTopoTypeIdentifiers( model );

      /// discerning patches by values for the region in terms of the diagnostic element variable
      // -------------------------------------------------------------
      // 2. Creating a Boundary from an internal region "NORMAL_FAULT"
      // -------------------------------------------------------------
      const bool connectivity_is_intact = integrityCheck<3,Element>( model.Mesh().ElementsBegin(), model.Mesh().ElementsEnd() );
      _test( connectivity_is_intact == true );

      const string test_region("NORMAL_FAULT");
      Region<3U>&  test_subdomain(model.Region(test_region.c_str()));
      const size_t elmts_original_region(test_subdomain.Cells());
      // will remove the original region
      const size_t model_faces_before(model.Mesh().Faces());
      std::pair<std::set<std::string>,bool> boundaries = model.CreateInternalBoundaryFrom( "NORMAL_FAULT" );
      //                                                       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      _test( boundaries.first.size() == subregions );
      _test( model.ContainsRegion("NORMAL_FAULT") == false );
      const size_t model_faces_after(model.Mesh().Faces());
      _test( (model_faces_after - model_faces_before) == elmts_original_region );

      // testing for existance of the new boundary patches
      if ( verbose_ ) {
           cout << "\nBoundaryInterface_Test::run: Printing the name of the boundaries in the model:";
           for ( auto it = model.BoundariesBegin(); it != model.BoundariesEnd(); ++it )
             cout << "\n\tBoundary: " << (*it).first <<" ("<< (*it).second.Cells() <<" faces)";
        }
      _test( model.Boundaries() == 12 );
      B_flags_model_after = nodeTopologyFlags<3,Element>( model_domain.NodesBegin(), model_domain.NodesEnd() );

      if ( verbose_ ) TestBoundaryAndTopoTypeIdentifiers( model );

      // creating property values on the boundaries and outputting these to VTU
      model.InputPropertyValue("nodal variable", makeScalar(PLAIN, 0.) );
      double  bvalue(1.3e5);
      for ( auto it = model.BoundariesBegin(); it != model.BoundariesEnd(); ++it ) {
           if ( verbose_ ) cout << "\n Boundary: " << (*it).first;
           (*it).second.InputPropertyValue("nodal variable", makeScalar(PLAIN, bvalue) );
           // testing VTU output
           if ( verbose_ ) vtu.OutputDataToVTU( (*it).first, string("nodal variable"), (*it).second, 0 );
           // creating different pressure values for each boundary patch
           bvalue += 1.0e5;
        }
      if ( verbose_ ) cout << endl;

      // testing whether boundary segments can be found by combined search criteria
      const set<string> intersected_regions{ "BOUNDARY", "LAYER_BOTTOM", "LAYER_TOP" };
      string patch_name = findBoundary( model, intersected_regions );
      _test( patch_name == "NORMAL_FAULT_BOUNDARY3_LAYER_BOTTOM_LAYER_TOP" );
      const set<string> search_strings{ "BOUNDARY", "NORMAL", "FAULT" };
      set<string> region_patches_found;
      const size_t patches_found = model.FindBoundaryByNames( search_strings, region_patches_found );
      // the boundary was decomposed into 6 patches
      _test( patches_found == 6 );
    
      // -----------------------------------
      // 3. testing supporting functionality
      // -----------------------------------
      _test( TestRegionContactDetection(model) );
      
      if ( verbose_ ) model.BoundariesOut();

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
      // box with a single volumetric region 'PORES' and all box boundaries and edges
      string  input_file("cube_flag");
   
      const bool binary_file(true);

      ANSYS_Model3D model( input_file.c_str(), "CSMP-variables.txt", binary_file );
      if ( verbose_ ) model.RegionsOut();
      
      const Region<3>& model_domain = model.Region("Model");
      set<TOPOTYPE> B_flags_model = nodeTopologyFlags<3,Element>( model_domain.NodesBegin(), model_domain.NodesEnd() );
      _test( B_flags_model.size() == 4 ); // mesh vertex, exterior point, line and surface
      
      if ( verbose_ ) cout <<"\nBoundaryInterface_Test::TestBoxShapedModel: done."<< endl;

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




void BoundaryInterface_Test::TestTopoTypeIdentifiers( Model<3U>& model )
 {
    model.CreateProperty( "topo type", "TT", "none");
    model.CreateProperty( "node BOX_BOUNDARY flag", "BBF", "none");
    model.CreateProperty( "element BOX_BOUNDARY flag", "EBF", "none", SCALAR, ELEMENT );
    
    // 0. getting diagnostics
    boxFlagsToVariable( model, "node BOX_BOUNDARY flag", "element BOX_BOUNDARY flag" );
    topoTypeToNumber( model, "topo type" );
 
    if ( verbose_ ) {
         VTU_Interface<3>  vtu(model);
         list<string> outvars{ "topo type", "node BOX_BOUNDARY flag","element BOX_BOUNDARY flag" };
         // the model as a whole
         vtu.OutputDataToVTU( "TestTopoTypeIdentifiers_regions_", outvars, model.Region("Model"), 0 );
         // all unique boundaries
         for ( auto it=model.UniqueRegionsBegin(); it!=model.UniqueRegionsEnd(); ++it )
           if ( !isDiagnosticBoxBoundaryClassifier( (*it).first ) )
             vtu.OutputDataToVTU( "TestTopoTypeIdentifiers_regions_", outvars, (*it).first, 0 );
      }
    
    model.DeleteProperty( "topo type" );
    model.DeleteProperty( "node BOX_BOUNDARY flag" );
    model.DeleteProperty( "element BOX_BOUNDARY flag" );
 }



/**
  Prints VTU file where topotype identifiers have been converted to number.
 */
void BoundaryInterface_Test::TopoTypeToVTU( Model<3U>& model, const char* region )
 {
    model.CreateProperty( "topo type", "TT", "none");
    
    // 0. getting diagnostics
    topoTypeToNumber( model, "topo type" );
 
    if ( verbose_ ) {
         VTU_Interface<3>  vtu(model);
         list<string> outvars{ "topo type" };
         vtu.OutputDataToVTU( "TestTopoTypes_of_region", outvars, model.Region(region), 0 );
      }
    
    model.DeleteProperty( "topo type" );
 }



/**
    Checks whether TOPOTYPE  node information is correct after the new boundaries have been inserted
    
    Step 1 - visualise the TOPOTYPE values using VTU
 */
void BoundaryInterface_Test::TestBoundaryAndTopoTypeIdentifiers( Model<3U>& model )
 {
    model.CreateProperty( "topo type", "TT", "none");
    model.CreateProperty( "node BOX_BOUNDARY flag", "BBF", "none");
    model.CreateProperty( "element BOX_BOUNDARY flag", "EBF", "none", SCALAR, ELEMENT );
    
    // 0. getting diagnostics
    boxFlagsToVariable( model, "node BOX_BOUNDARY flag", "element BOX_BOUNDARY flag" );
    topoTypeToNumber( model, "topo type" );
 
    if ( verbose_ ) {
         VTU_Interface<3>  vtu(model);
         // form a region of all boundaries
         list<string> outvars{ "topo type", "node BOX_BOUNDARY flag","element BOX_BOUNDARY flag" };
         for ( auto it=model.BoundariesBegin(); it!=model.BoundariesEnd(); ++it )
           if ( !isDiagnosticBoxBoundaryClassifier( (*it).first ) )
             vtu.OutputDataToVTU( "TestTopoTypeIdentifiers_boundaries_", outvars, (*it).second, 0 );
      }
    
    model.DeleteProperty( "topo type" );
    model.DeleteProperty( "node BOX_BOUNDARY flag" );
    model.DeleteProperty( "element BOX_BOUNDARY flag" );
 }






/**
    numbers the unique regions of the models, labeling their elements with the region number
    as "region identifier".
    
    @parameter vector of region names to retrieve them from the integer keys
*/
template<uint32_t dim>
size_t countAndLabelRegions( Model<dim>& model, const char* region_identifier, std::vector<std::string>& region_names )
 {
    // setting element variable up to identify all unique regions - uniquely
    const string rvariable(region_identifier);
    if ( !model.Database().IsDefined(rvariable.c_str()) )
       model.CreateProperty( rvariable.c_str(), "rid", "none", SCALAR, ELEMENT );
   
    // counting the regions and initialising them with the unique identifiers
    region_names.resize(model.UniqueRegions());
    size_t regions(0);
    for ( auto it=model.UniqueRegionsBegin(); it!=model.UniqueRegionsEnd(); it++ ) {
         (*it).second.InputPropertyValue( rvariable.c_str(), makeScalar(PLAIN,static_cast<double>(regions)) );
         region_names[regions] = (*it).first;
         regions++;
      }
   
   assert( regions == region_names.size() );
   return regions;
   
 } // end countAndLabelRegions

template size_t countAndLabelRegions( Model<1U>&, const char*, std::vector<std::string>& );
template size_t countAndLabelRegions( Model<2U>&, const char*, std::vector<std::string>& );
template size_t countAndLabelRegions( Model<3U>&, const char*, std::vector<std::string>& );





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
static bool checkNeighborNormalsForConsistentOrientation( const Region<3U>&  subdomain )
 {
    vector<double> normal(3U), nbor_normal(3U);
   
    size_t non_surface_elements(0U);
    for ( auto it=subdomain.CellsBegin(); it!=subdomain.CellsEnd(); ++it )
      // this method only  considers surface elements
      if ( (*it)->IsSurface() ) {
           (*it)->UnitNormal( normal );
           const size_t neighbors((*it)->Neighbors());
           for ( uint32_t i{0}; i<neighbors; ++i )
             // only valid neighbor elements are considered
             if ( (*it)->Neighbor(i) != nullptr ) {
                  (*it)->Neighbor(i)->UnitNormal( nbor_normal );
                  // projection
                  double result(0.);
                  for ( size_t j=0U; j<3U; ++ j )
                    result += normal[j] * nbor_normal[j];
                  if ( result < 0. )
                    return false;
               }
         }
       else non_surface_elements++;
   
    if ( non_surface_elements > 0U )
      ErrorHandler::Instance().Note( ERROR, "checkNeighborNormalsForConsistentOrientation:",
                                       subdomain.Name(), "region contained not only surface elements." );
    return true;
   
 } // end checkNeighborNormalsForConsistentOrientation






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
template<uint32_t dim>
static std::string internalBoundaryNameFrom( const FaceConstructionData<dim>& fdata, const std::vector<std::string>& region_names )
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
         ErrorHandler::Instance().Note( WARNING, "findBoundary:", "supplied set of substrings is empty; returning '\0'." );
         return std::string("\0");
      }
    // if the model has no boundaries
    if ( model.Boundaries() == 0 ) {
         ErrorHandler::Instance().Note( WARNING, "findBoundary:", "model has no boundaries; returning '\0'." );
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
          ErrorHandler::Instance().Note( ERROR, "labelRegionPatches:", dim_1_region, "does not exist; nothing was done." );
          return 0;
      }
    if ( model.UniqueRegions() <= 1 ) {
          ErrorHandler::Instance().Note( WARNING, "labelRegionPatches:", "model contains only a single unique region; so there is only one patch." );
          return 1;
      }
    // verifying that we are indeed dealing with a region of surface elements only
    Region<3U>&  subdomain(model.Region(dim_1_region));
    pair<int32_t,int32_t>  dimensionality = subdomain.ElementSpatialDimensions();
    //   number of dims in region      dimension of contained elements
    if ( dimensionality.first != 1 and dimensionality.second != 2 ) {
         ErrorHandler::Instance().Note( ERROR, "labelRegionPatches:", dim_1_region, "region does not consist of surface elements only; nothing was done." );
         return 0;
      }
    // verifying that the region lies inside of the model
    size_t boundary_elements(0);
    for ( auto eit=subdomain.CellsBegin(); eit!=subdomain.CellsEnd(); ++eit )
      {
         for ( uint32_t i{0}; i<(*eit)->Nodes(); ++i )
           if ( (*eit)->N(i)->AtBoundary() != NOT && (*eit)->N(i)->AtBoundary() != INTERNAL ) {
                cerr <<"\nElement is a boundary element:";
                (*eit)->Out();
                boundary_elements++;
                break;
             }
      }
    if ( boundary_elements > 0 ) {
         ErrorHandler::Instance().Note( ERROR, "labelRegionPatches:", dim_1_region, "region appears to lie at the model boundary; nothing was done." );
         return 0;
      }
    if ( !model.Database().IsDefined(diagnostic_elmt_variable) ) {
         ErrorHandler::Instance().Note( ERROR, "labelRegionPatches:", diagnostic_elmt_variable, "variable to discern regions is not defined; nothing was done." );
         return 0;
      }
      {  // check whether there are multiple region identifiers
         double rmin, rmax;
         model.MinMaxOf( diagnostic_elmt_variable, rmin, rmax );
         if ( fabs(rmax - rmin) <= numeric_limits<double>::epsilon() ) {
              ErrorHandler::Instance().Note( WARNING, "labelRegionPatches:", diagnostic_elmt_variable, "is single valued; so there is only one patch." );
              return 1;
           }
      }

    // 2. setting up an element variable to show the partitions of the lower dimensional region
    // ----------------------------------------------------------------------------------------
    csmp::Index mtrl_key = model.Database().StorageKey(diagnostic_elmt_variable);
   
    if ( !model.Database().IsDefined(patch_variable) )
      model.CreateProperty( patch_variable, "pvar", "node", SCALAR, ELEMENT );
    // patch-discerning variable
    csmp::Index pvar_key = model.Database().StorageKey(patch_variable);
   
    // renumbering nodes and elements
    model.Region("Model").UpdateMemberIndexes();
   
    // looping over the region, identifying and recording the juxtaposition relationships
    map<pair<long,long>,uint32_t>     patches;
    vector<FaceConstructionData<3U>>  face_construction_data;
    map<uint32_t,string>              patch_names;
    string                            patch_name;
    uint32_t                          n_juxtapositions(0);

    // 3. looping over lower dimensional region identifying juxtaposition relationships
    // ----------------------------------------------------------------------------------------
    face_construction_data.reserve(subdomain.Cells());
    for ( auto eit=subdomain.CellsBegin(); eit!=subdomain.CellsEnd(); ++eit )
      {
          // 3.1 identifying neighbors, facing relations, and juxtaposed materials for current element
          FaceConstructionData  fdata(higherDimensionalNeighbors( *(*eit), mtrl_key ));
        
          // 3.2 recording which category of juxtaposition element fall into, naming it and assigning a patch number
          pair<map<pair<long,long>,uint32_t>::iterator,bool>  it=patches.insert( make_pair(fdata.Materials(),n_juxtapositions) );
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
          (*eit)->Store( pvar_key, makeScalar(PLAIN,static_cast<double>((*it.first).second)) );
        
          // 3.4 recording the data for the element that will later be used to construct the face from
          face_construction_data.push_back( fdata );
      }
   
    // 4. creating labeled boundary patches from the face-defining data
    // ----------------------------------------------------------------
    // 4.1 making a map 'patch_numbers' from 'patch_names' to search for patch identifiers
    map<string,uint32_t>  patch_numbers;
    for ( auto it=patch_names.begin(); it!=patch_names.end(); ++it )
      patch_numbers.insert( make_pair( (*it).second, (*it).first ) );
   
    // 4.2 building new map where the patch faces are organised by patch names
    map<string,vector<FaceConstructionData<3U>> > patch_simplexes;
    vector<FaceConstructionData<3U>>              empty_vec;
    for ( auto it=patch_names.begin(); it!=patch_names.end(); ++it )
      patch_simplexes.insert( make_pair( (*it).second, empty_vec ) );
   
    // 4.3 inserting the patch identifiers into it
    for ( auto pit=patch_simplexes.begin(); pit!=patch_simplexes.end(); ++pit )
      {
         assert( patch_numbers.find((*pit).first) != patch_numbers.end() );
         const uint32_t patch_number((*patch_numbers.find((*pit).first)).second);
         // reserving storage
         (*pit).second.reserve(face_construction_data.size());
         // looping over all face data assigning the ones that are suitable
         for ( auto it=face_construction_data.begin(); it!=face_construction_data.end(); ++it )
           if ( (*it).PatchNumber() == patch_number )
             (*pit).second.push_back( (*it) );
      }
    // 4.4 trimming excess storage of the face-data vectors
    for ( auto pit=patch_simplexes.begin(); pit!=patch_simplexes.end(); ++pit )
      //std::vector<uint32_t>(elementIds).swap(elementIds);
      vector<FaceConstructionData<3U>>( (*pit).second ).swap( (*pit).second );

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
