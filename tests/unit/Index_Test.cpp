#include "Index_Test.h"
#include "PropertyDatabase.h"
#include <iostream>

using namespace std;

namespace csmp
{
    Index_Test::Index_Test()
    {

    }

    Index_Test::~Index_Test()
    {

    }

    void Index_Test::IndexOperatorsTest(Index &idx)
    {
        csmp::Index new_index;
        _test(new_index.IsDefined()==0);
        new_index=idx;
        _test(new_index.IsDefined()==1); // test = operator

        _test((new_index==idx)==1); //test == operator

        csmp::Index idx_tmp;
        PropertyDatabase<3> pdb;
        idx_tmp=pdb.AddProperty("test variable 2","m s",VECTOR,ELEMENT);
        _test((new_index==idx_tmp)==1); //test == operator
        _test((new_index!=idx_tmp)==0); //test != operator

        //_test((idx_tmp < new_index)==1); //test < operator

        _test(new_index.type==VECTOR);
        _test(new_index.place==ELEMENT);
    }

    void Index_Test::run()
    {
        PropertyDatabase<3> pdb;
        csmp::Index idx;
        _test(idx.IsDefined()==0);
        idx=pdb.AddProperty("test variable","m s-1",VECTOR,ELEMENT);
        _test(idx.IsDefined()==1);
        this->IndexOperatorsTest(idx);
    }
}
