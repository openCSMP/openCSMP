// unit test
// Created by : Ali
// Date of Modification: 25.10.2010
// Passed All Tests? --> ( )

#include "Parameter_Test.h"

//CSMP Header File
#include "Parameter.h"
namespace csmp
{
    Parameter_Test::Parameter_Test()
    {
        t_param=new csmp::Parameter;
        t_param->name="absolute fluid pressure";
        t_param->notation="pf";
        t_param->unit="Pa";
        t_param->min=100000;
        t_param->max=1000000000;
        t_param->usage=NODE;
    }

    Parameter_Test::~Parameter_Test()
    {
        if( t_param!=NULL )
            delete t_param;
    }

    void Parameter_Test::RangeTest()
    {
        double vmin,vmax;

        t_param->Range(vmin,vmax);
        _test(vmin==100000);
        _test(vmax==1000000000);
    }

    void Parameter_Test::IsWithinRangeTest()
    {
        //test Range
        _test(t_param->IsWithinRange(99999)==false);
        _test(t_param->IsWithinRange(100000)==true);
        _test(t_param->IsWithinRange(1000000)==true);
        _test(t_param->IsWithinRange(10000000)==true);
        _test(t_param->IsWithinRange(100000000)==true);
        _test(t_param->IsWithinRange(1000000000)==true);
#ifndef COMPILER_IS_32_BIT
         _test(t_param->IsWithinRange(10000000001)==false);
#endif
    }

    void Parameter_Test::run()
    {
        _test(t_param->name=="absolute fluid pressure");
        IsWithinRangeTest();
        RangeTest();
    }
}

