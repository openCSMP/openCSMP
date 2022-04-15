/*
 *  PropertyConstraints_Test.h
 *  csmp_core
 *
 *  Created by SKM 2/2/2022.
 *
 */

#include "PropertyConstraints_Test.h"

#include "vsetMakers.h"
#include "Model.h"
#include "Region.h"
#include "Boundary.h"
#include "PropertyConstraints.h"
#include "VSet.h"
#include "VTK_Interface.h"
#include "ANSYS_Model3D.h"
#include "MeshManagementUtilities.h"
#include "compareFloats.h"

using namespace std;

namespace csmp {

PropertyConstraints_Test::PropertyConstraints_Test( bool verbose )
    : verbose_(verbose)
    {
    }
    

PropertyConstraints_Test::~PropertyConstraints_Test()
    {
    }
    

/**
    Reading of all basic variable types and placements.
*/
void PropertyConstraints_Test::run()
  {
    // testing the repair of non-unique regions built using property constraints
    TestBuildRegionsFromPropertyConstraints();

 
   } // end run





 
 
 
bool PropertyConstraints_Test::TestBuildRegionsFromPropertyConstraints()
 {
     const bool    using_isoparametric_elements{true};
     ModelTopology topology( "FracBox", using_isoparametric_elements );
     VSet<3U>      vset;
     test_Create_FracBox( topology, vset );

     const bool create_boundaries_from_faces{true}, box_shaped{false};
     Model<3U>  model( topology, vset, "CSMP-1phase-variables.txt", false );
     
     // assigning some dummy values to verify functionality
     Region<3U>& model_domain = model.Region("Model");
     Region<3U>& matrix_domain = model.Region("MATRIX");
     const csmp::Index phi_key = model.Database().StorageKey("porosity");
     model_domain.InputPropertyValue( "permeability", makeScalar(ANY,1.0e-12) );
     model_domain.InputPropertyValue( "porosity", makeScalar(ANY,1.0) ); // fracture
     matrix_domain.InputPropertyValue( "porosity", makeScalar(ANY,0.4) ); // matrix
     // elevating the porosity of some extra elements
     for ( size_t eidx{50}; eidx<150; ++eidx )
       matrix_domain.E(eidx)->Store( phi_key, makeScalar(ANY,0.8) );
       
     // building new region from PropertyConstraints
     PropertyConstraints porosity_constraints( "porosity", 0.75, 0.85 );
     porosity_constraints.AddConstraint( "permeability", 1.0e-12, 1.0e-12 );
     
     model.FormRegionFrom( "medium porosity", porosity_constraints );
     const Region<3U>& medium_porosity_domain = model.Region("medium porosity");
     _test( medium_porosity_domain.Cells() == 100 );

     return true;
 }




 /**
   Tests pointInVolumeElement:
   
   Tests prism_test model because it contains elements of all
   types.
   */
  bool PropertyConstraints_Test::PointInVolumeElementTest()
  {
    ANSYS_Model3D model( "prism_test", "CSMP-variables.txt", true, true, true );
    
    Point<3u> query(2434.0f, -1510.0f, 5400.0f);
    
    auto& gref = model.Region("Model");
    
    auto eend = gref.CellsEnd();
    for (auto eit = gref.CellsBegin(); eit != eend; ++eit) {
      if (!(*eit)->IsVolume()) {
        continue;
      }
      const Point<3u> bctr = (*eit)->BaryCenter();
      
      Element<3u>* e = pointInVolumeElement(gref, bctr);
      _test(e == *eit);
    }
    return true;
  }
  






  
} // end csmp



