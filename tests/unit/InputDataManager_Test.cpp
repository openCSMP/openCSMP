/*
 *  InputDataManager_Test.cpp
 *  csmp_core
 *
 *  Created by Ali Tabatabaei on 11/4/10.
 *  Copyright 2010 CSMP. All rights reserved.
 *
 */

#include "InputDataManager_Test.h"

#include "InputDataManager.h"
#include "Model.h"
#include "Region.h"
#include "Boundary.h"
#include "VSet.h"
#include "ANSYS_Model3D.h"

using namespace std;

namespace csmp {

InputDataManager_Test::InputDataManager_Test( bool verbose )
    : verbose_(verbose)
    {
    }
    

InputDataManager_Test::~InputDataManager_Test()
    {
    }
    

/**
    Reading of all basic variable types and placements.
*/
void InputDataManager_Test::run()
  {
      ANSYS_Model3D  model( "InputDataManager_Test", "InputDataManager_Test-variables.txt" );
    
      // TESTING
      InputDataManager<3U>().ConfigureFromFile( model,"InputDataManager_Test", false, true, true, true, true, true );
      // element props
      const csmp::Index phi_key=model.Database().StorageKey("porosity");
      const csmp::Index velo_key=model.Database().StorageKey("velocity");
      const csmp::Index k_key=model.Database().StorageKey("tensor permeability");
      const csmp::Index aq_key=model.Database().StorageKey("aqueous species");
      // node props
      const csmp::Index pf_key=model.Database().StorageKey("fluid pressure");
      const csmp::Index nv_key=model.Database().StorageKey("nodal velocity");
    
      // testing the first element and some interior nodes of the respective model regions
      const Region<3U> matrix_domain(model.Region("MATRIX"));
      const Region<3U> fracture_domain(model.Region("FRAC_VOLUMES"));
      // testing
      const double tolerance(1.0e-10);
      // scalars
      _equal( matrix_domain.E(0)->Read(phi_key), 0.25, tolerance );
      _equal( fracture_domain.E(0)->Read(phi_key), 1.0, tolerance );
      // tensor
      TensorVariable<3U>  ts; // only upper triagonal
      matrix_domain.E(0)->Read( k_key, ts );
      _equal( ts(0,0), 1.0e-12, tolerance );
      _equal( ts(1,1), 1.0e-13, tolerance );
      _equal( ts(2,2), 1.0e-14, tolerance );
      // array variables
      ArrayVariable  ary(5U);
      matrix_domain.E(0)->Read( aq_key, ary );
      _equal( ary(0), 0.1, tolerance );
      _equal( ary(1), 0.2, tolerance );
      _equal( ary(2), 0.3, tolerance );
      _equal( ary(3), 0.4, tolerance );
      _equal( ary(4), 0.5, tolerance );
    
      // model props
      const csmp::Index t_key=model.Database().StorageKey("time");
      const csmp::Index xyz_key=model.Database().StorageKey("coordinate system");
      const csmp::Index sig_key=model.Database().StorageKey("stress");
      _equal( model.Read(t_key), 3600., tolerance );
      VectorVariable<3U>  vc;
      model.Read( xyz_key, vc );
      _equal( vc[0], 3., tolerance );
      _equal( vc[1], 6., tolerance );
      _equal( vc[2], 9., tolerance );
    
      // region props
      const csmp::Index fac_key=model.Database().StorageKey("facies");
      _equal( fracture_domain.Read( fac_key ), 2., tolerance );
    
      // boundary props
      const csmp::Index hfu_key=model.Database().StorageKey("heat flux");
      const csmp::Index bhfu_key=model.Database().StorageKey("boundary heat flux");
      const csmp::Index trac_key=model.Database().StorageKey("traction");
      const Boundary<3U>  bottom(model.Boundary("BOTTOM")), 
                          back(model.Boundary("BACK"));
      // testing
      _equal( back.E(0)->Read( hfu_key ), 2., tolerance );
      _equal( back.Read( bhfu_key ), 2., tolerance );
      bottom.Read( trac_key, vc );
      _equal( vc[0], 1.0e8,  tolerance );
      _equal( vc[1], 1.0e9,  tolerance );
      _equal( vc[2], 1.0e10, tolerance );
  }
  
} // end csmp



