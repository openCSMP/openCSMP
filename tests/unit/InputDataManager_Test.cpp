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

namespace csmp
{
    using namespace std;
    InputDataManager_Test::InputDataManager_Test()
    {
    }
    
    InputDataManager_Test::~InputDataManager_Test()
    {
        if( _model!=NULL )
            delete _model;
    }
    
    void InputDataManager_Test::run()
    {
        _model=new csmp::ANSYS_Model3D("prism_test","CSMP_Variables.txt");

        InputDataManager<3U> idm;
        idm.ConfigureFromFile(*_model,"prism_test",false, true, true, true, true );


        { // Create new no-name Scope
            csmp::Index porosity_idx=_model->Database().StorageKey("porosity");
            csmp::Index permeability_idx=_model->Database().StorageKey("permeability");

            double porosity(0.0);
            double permeability(0.0);
            std::vector<Element<3U>*>::const_iterator eit(_model->Region("FRAC_VOLUMES").ElementsBegin());
            for (;eit != _model->Region("FRAC_VOLUMES").ElementsEnd();++eit)
            {
                porosity=(*eit)->Read(porosity_idx);
                permeability=(*eit)->Read(permeability_idx);
                _test(porosity==1.000e+00);
                _test(permeability==1.0e-10);
            }
            eit=_model->Region("MATRIX").ElementsBegin();
            for (;eit != _model->Region("MATRIX").ElementsEnd();++eit)
            {
                porosity=(*eit)->Read(porosity_idx);
                permeability=(*eit)->Read(permeability_idx);
                _test(porosity==0.25);
                _test(permeability==1.0e-12);
            }
        } // End no-name Scope

    }    
}
