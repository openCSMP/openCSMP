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
#include "VSet.h"
#include "ANSYS_Model3D.h"

namespace csmp {

    using namespace std;
  

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
        ANSYS_Model3D model("InputDataManager_Test","CSMP_DataInputManager_Test-variables.txt");

        InputDataManager<3U> idm;
      
        idm.ConfigureFromFile(model,"InputDataManager_Test",false, true, true, true, true );


        { // Create new no-name Scope
            csmp::Index porosity_idx=model.Database().StorageKey("porosity");
            csmp::Index permeability_idx=model.Database().StorageKey("permeability");

            double porosity(0.0);
            double permeability(0.0);
            for ( auto eit=model.Region("FRAC_VOLUMES").ElementsBegin();
                  eit!=model.Region("FRAC_VOLUMES").ElementsEnd(); ++eit )
            {
                porosity=(*eit)->Read(porosity_idx);
                permeability=(*eit)->Read(permeability_idx);
                _test(porosity==1.000e+00);
                _test(permeability==1.0e-10);
            }
            for ( auto eit=model.Region("MATRIX").ElementsBegin(); eit != model.Region("MATRIX").ElementsEnd();++eit)
            {
                porosity=(*eit)->Read(porosity_idx);
                permeability=(*eit)->Read(permeability_idx);
                _test(porosity==0.25);
                _test(permeability==1.0e-12);
            }
        } // End no-name Scope

    }
  
} // end csmp
