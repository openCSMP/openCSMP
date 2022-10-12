//
//  Experimental_Example.cpp
//  CSMP_API_library2014
//
//  Created by Stephan Matthai on 1/27/14.
//  Copyright (c) 2014 Stephan Matthai. All rights reserved.
//

#include "Experimental_Example.h"
#include "compareFloats.h"

#ifdef CSMP_WITH_SAMG_SOLVER
#include "SAMG_Solver.h"
#include "SAMG_Settings.h"
#endif
#ifdef CSMP_WITH_MESCHACH
#include "Gauss_Solver.h"
#endif

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

void Experimental_Example::Specifications()
  {
     SetTitle( "Experimental_Example" );
     SetDifficulty( 1 );
     SetCategory( "Software Functionality" );
     AddAuthor( "You!" );
     AddDescription( "source in: Experimental_Example.cpp" );
     AddDescription( "Empty example for the user to experiment with" );
     AddRequirement( "none" );
     AddRequirement( "no predefined model or variables file" );
  }





/**
     Put CSMP code that you would like to test here and run it as part of the 
     example suite.

*/
void Experimental_Example::Run()
  {
      // illustrating Boundary creation for faults
      // -----------------------------------------
      const string  input_file("fault_boundary_test");
      const uint32_t dim(3);

      // ------------------------------
      // 1. Building Ansys model
      // ------------------------------
      const bool irregular_mesh(true);
      const bool binary_file(true);

      ANSYS_Model3D model( input_file.c_str(), "CSMP-variables.txt", irregular_mesh, binary_file );
      printModelDimensions(model, true);

      /// assuming a dim-1 region, label and count material juxtaposition relationships
      const string    region_tag("region identifier");
      vector<string>  region_names;
      size_t regions = model.CountAndLabelUniqueRegions( region_tag.c_str(), region_names );

      VTU_Interface<dim>  vtu(model);
      vtu.OutputDataToVTU( "BoundaryInterface_Test_", region_tag, string("Model"), 0 );

      const string patch_tag("patch identifier");
      // creating visual output that illustrates what the boundary should look like for testing
      size_t subregions; // = labelRegionPatches( model, "NORMAL_FAULT", region_tag.c_str(), patch_tag.c_str(), region_names );

      cout <<"\nBoundaryInterface_Test::run: identified "<< subregions <<" region patches in region NORMAL_FAULT touching "<< regions <<" model regions.\n";
      vtu.OutputDataToVTU( "test", patch_tag, string("NORMAL_FAULT"), 0 );
    
	  VTU_Interface<3> vtu_boundary(model);
      {
        string boundary_name("LEFT");
        string variableName("face variable");
        Boundary<dim>& boundary(model.Boundary(boundary_name));
        Index areaKey(model.Database().StorageKey(variableName.c_str()));
        assert(areaKey.place == FACE); // this needs to be done for Face::Read bc Read is inherited from Element and can
                       // either be Face, Element or InterFace
        ScalarVariable area(PLAIN, 0.);
        size_t surfaceElementCount(0);
        const auto domainElementsEnd(boundary.CellsEnd());
        for ( auto it = boundary.CellsBegin(); it != domainElementsEnd; ++it)
        {
          area = (*it)->Area();
          (*it)->Store(areaKey, area);
          ++surfaceElementCount;
        }
        if (model.ContainsBoundary(boundary_name))
          vtu.OutputDataToVTU("test", variableName, boundary, 0);
      }
      {
        string boundary_name("RIGHT");
        string variableName("face variable");
        Boundary<dim>& boundary(model.Boundary(boundary_name));
        Index areaKey(model.Database().StorageKey(variableName.c_str()));
        assert(areaKey.place == FACE); // this needs to be done for Face::Read bc Read is inherited from Element and can
                       // either be Face, Element or InterFace
        ScalarVariable area(PLAIN, 0.);
        size_t surfaceElementCount(0);
        const typename vector<Face<dim>*>::const_iterator domainElementsEnd(boundary.CellsEnd());
        for (typename vector<Face<dim>*>::const_iterator it = boundary.CellsBegin(); it != domainElementsEnd; ++it)
        {
          area = (*it)->Area();
          (*it)->Store(areaKey, area);
          ++surfaceElementCount;
        }
        if (model.ContainsBoundary(boundary_name))
        vtu.OutputDataToVTU("test", variableName, boundary, 0);
      }
      {
        string boundary_name("FRONT");
        string variableName("face variable");
        Boundary<dim>& boundary(model.Boundary(boundary_name));
        Index areaKey(model.Database().StorageKey(variableName.c_str()));
        assert(areaKey.place == FACE); // this needs to be done for Face::Read bc Read is inherited from Element and can
                       // either be Face, Element or InterFace
        ScalarVariable area(PLAIN, 0.);
        size_t surfaceElementCount(0);
        const typename vector<Face<dim>*>::const_iterator domainElementsEnd(boundary.CellsEnd());
        for (typename vector<Face<dim>*>::const_iterator it = boundary.CellsBegin(); it != domainElementsEnd; ++it)
          {
            area = (*it)->Area();
            (*it)->Store(areaKey, area);
            ++surfaceElementCount;
          }
        if (model.ContainsBoundary(boundary_name))
        vtu.OutputDataToVTU("test", variableName, boundary, 0);
      }
      {
        string boundary_name("BACK");
        string variableName("face variable");
        Boundary<dim>& boundary(model.Boundary(boundary_name));
        Index areaKey(model.Database().StorageKey(variableName.c_str()));
        assert(areaKey.place == FACE); // this needs to be done for Face::Read bc Read is inherited from Element and can
                       // either be Face, Element or InterFace
        ScalarVariable area(PLAIN, 0.);
        size_t surfaceElementCount(0);
        const typename vector<Face<dim>*>::const_iterator domainElementsEnd(boundary.CellsEnd());
        for (typename vector<Face<dim>*>::const_iterator it = boundary.CellsBegin(); it != domainElementsEnd; ++it)
        {
          area = (*it)->Area();
          (*it)->Store(areaKey, area);
          ++surfaceElementCount;
        }
        if (model.ContainsBoundary(boundary_name))
        vtu.OutputDataToVTU("test", variableName, boundary, 0);
      }
      {
        string boundary_name("TOP");
        string variableName("face variable");
        Boundary<dim>& boundary(model.Boundary(boundary_name));
        Index areaKey(model.Database().StorageKey(variableName.c_str()));
        assert(areaKey.place == FACE); // this needs to be done for Face::Read bc Read is inherited from Element and can
                       // either be Face, Element or InterFace
        ScalarVariable area(PLAIN, 0.);
        size_t surfaceElementCount(0);
        const typename vector<Face<dim>*>::const_iterator domainElementsEnd(boundary.CellsEnd());
        for (typename vector<Face<dim>*>::const_iterator it = boundary.CellsBegin(); it != domainElementsEnd; ++it)
          {
            area = (*it)->Area();
            (*it)->Store(areaKey, area);
            ++surfaceElementCount;
          }
        if (model.ContainsBoundary(boundary_name))
        vtu.OutputDataToVTU("test", variableName, boundary, 0);
      }
      {
        string boundary_name("BOTTOM");
        string variableName("face variable");
        Boundary<dim>& boundary(model.Boundary(boundary_name));
        Index areaKey(model.Database().StorageKey(variableName.c_str()));
        assert(areaKey.place == FACE); // this needs to be done for Face::Read bc Read is inherited from Element and can
                       // either be Face, Element or InterFace
        ScalarVariable area(PLAIN, 0.);
        size_t surfaceElementCount(0);
        const typename vector<Face<dim>*>::const_iterator domainElementsEnd(boundary.CellsEnd());
        for (typename vector<Face<dim>*>::const_iterator it = boundary.CellsBegin(); it != domainElementsEnd; ++it)
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
      pair<set<string>,bool> boundaries = model.CreateInternalBoundaryFrom( "NORMAL_FAULT" );

      // testing for existance of the new boundary patches
      cout << "\nBoundaryInterface_Test::run: Printing the name of the boundaries in the model:";
      for ( auto it = model.BoundariesBegin(); it != model.BoundariesEnd(); ++it )
        cout << "\n\tBoundary: " << (*it).first <<" ("<< (*it).second.Cells() <<" faces)";

      // creating property values on the boundaries and outputting these to VTU
      model.InputPropertyValue("nodal variable", makeScalar(PLAIN, 0.) );
      double  bvalue(1.3e5);
      for ( auto it = model.BoundariesBegin(); it != model.BoundariesEnd(); ++it ) {
           cout << "\n Boundary: " << (*it).first;
           (*it).second.InputPropertyValue("nodal variable", makeScalar(PLAIN, bvalue) );
           // testing VTU output
           vtu.OutputDataToVTU( (*it).first, string("nodal variable"), (*it).second, 0 );
           // creating different pressure values for each boundary patch
           bvalue += 1.0e5;
        }
      cout << endl;

      // testing whether boundary segments can be found by combined search criteria
      const set<string> intersected_regions{ "BOUNDARY", "LAYER_BOTTOM", "LAYER_TOP" };
      string patch_name = findBoundary( model, intersected_regions );
      const set<string> search_strings{ "BOUNDARY", "NORMAL", "FAULT" };
      set<string> region_patches_found;

} // end






} // csmp
