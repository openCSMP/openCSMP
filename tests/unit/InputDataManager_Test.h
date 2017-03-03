/*
 *  InputDataManager_Test.h
 *  csmp_core
 *
 *  Created by Ali Tabatabaei on 11/4/10.
 *  Copyright 2010 csmp-project All rights reserved.
 *
 */

#ifndef INPUT_DATA_MANAGER_TEST_H
#define INPUT_DATA_MANAGER_TEST_H

#include "Test.h"


namespace csmp
{
    class ANSYS_Model3D;
    class InputDataManager_Test: public Test
    {
    public:
        InputDataManager_Test();
        ~InputDataManager_Test();
        void run();
    private:
        ANSYS_Model3D* _model;
    };
}

#endif