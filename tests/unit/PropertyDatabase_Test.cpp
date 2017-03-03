
#include "PropertyDatabase_Test.h"

//CSMP Header Files
#include "PropertyDatabase.h"
#include "PropertyHandle.h"

#include <iostream>
using namespace std;
namespace csmp
{
PropertyDatabase_Test::PropertyDatabase_Test()
    :pdb(new PropertyDatabase<3>())
{

}

PropertyDatabase_Test::~PropertyDatabase_Test()
{
  if(pdb!=NULL)
    delete pdb;
}


void PropertyDatabase_Test::NameTest( const csmp::Index& idx,const char* name)
{
    _test(strcmp(pdb->Name(idx),name)==0);
}


void PropertyDatabase_Test::UnitTest( const csmp::Index& idx, const char* unit)
{
    _test(strcmp(pdb->Unit(pdb->Name(idx)),unit)==0);
}


/**
    At the beginning, there are variables of all types and placements,
    but no Index objects that need to be tracked.
    
    Two new properties are created here testing their addition and 
    registration with the IndexTracker.
    
    At the end of this text function, these 2 tracked Index objects
    are destructed again (leaving the function scope), triggering
    their removal from the tracked Index map inside of the 
    IndexTracker object.
*/
void PropertyDatabase_Test::AddPropertyTest()
{        
    csmp::Index idx = pdb->AddProperty("test variable","m s-1",VECTOR,REGION);
    this->UnitTest(idx,"m s-1");
    this->NameTest(idx,"test variable");
    _test(strcmp(csmp::parsePlacement(idx.place).c_str(),"REGION")==0); //Tests if the placement is REGION
    _test(strcmp(csmp::parsePlacement(idx.place).c_str(),"ELEMENT") !=0);
    _test(strcmp(csmp::parsePlacement(idx.place).c_str(),"MODEL")  !=0);
    _test(strcmp(csmp::parsePlacement(idx.place).c_str(),"NODE")  !=0);

    _test(strcmp(csmp::parseType(idx.type).c_str(),"VECTOR")==0); //Tests if the type is VECTOR
    _test(strcmp(csmp::parseType(idx.type).c_str(),"TENSOR")!=0);
    _test(strcmp(csmp::parseType(idx.type).c_str(),"SCALAR")!=0);

    csmp::Index idx2;
    idx2=pdb->AddProperty("test variable 2","m s",SCALAR,ELEMENT);
    this->UnitTest(idx2,"m s");
    this->NameTest(idx2,"test variable 2");
    double mintest(0.),maxtest(0.);
    pdb->RangeOf("test variable 2", mintest, maxtest);
    _test(mintest==-1.0e+30);
    _test(maxtest==1.0e+30);

    _test(strcmp(csmp::parsePlacement(idx2.place).c_str(),"ELEMENT")==0); //Tests if the placement is REGION
    _test(strcmp(csmp::parsePlacement(idx2.place).c_str(),"REGION") !=0);
    _test(strcmp(csmp::parsePlacement(idx2.place).c_str(),"MODEL")  !=0);
    _test(strcmp(csmp::parsePlacement(idx2.place).c_str(),"NODE")  !=0);

    _test(strcmp(csmp::parseType(idx2.type).c_str(),"SCALAR")==0); //Tests if the type is SCALAR
    _test(strcmp(csmp::parseType(idx2.type).c_str(),"TENSOR")!=0);
    _test(strcmp(csmp::parseType(idx2.type).c_str(),"VECTOR")!=0);

    _test( pdb->IsDefined(idx) == true );
}


void PropertyDatabase_Test::IsDefinedTest()
{
    _test(pdb->IsDefined("test variable 2"));
    _test(pdb->IsDefined("test variable"));
    _test(pdb->IsDefined("not defined")==0);
}


void PropertyDatabase_Test::PlacementTest()
{
    _test(strcmp(csmp::parsePlacement(pdb->Placement("test variable")).c_str(),"REGION")==0);
    _test(strcmp(csmp::parsePlacement(pdb->Placement("test variable 2")).c_str(),"ELEMENT")==0);
    _test(strcmp(csmp::parsePlacement(pdb->Placement("test variable 2")).c_str(),"REGION")!=0);
}


void PropertyDatabase_Test::TypeTest()
{
    _test(strcmp(csmp::parseType(pdb->Type("test variable 2")).c_str(),"SCALAR")==0);
    _test(strcmp(csmp::parseType(pdb->Type("test variable 2")).c_str(),"TENSOR")!=0);
}


void PropertyDatabase_Test::VariablesTotalTest(size_t s)
{
    _test(pdb->VariableCount()==s);
}


void PropertyDatabase_Test::VariablesRegionTest(size_t r)
{
    _test(pdb->VariableCount(REGION)==r);
}


void PropertyDatabase_Test::DeletePropertyTest()
{
    pdb->DeleteProperty("test variable 2");
    _test(pdb->IsDefined("test variable 2")==0);
    this->VariablesTotalTest(1);
}


void PropertyDatabase_Test::VariablesElementTest(size_t e)
{
    _test(pdb->VariableCount(ELEMENT)==e);
    csmp::Index idx=pdb->AddProperty("test variable 3","m s-1",VECTOR,ELEMENT);
    _test(pdb->VariableCount(ELEMENT)==(e+1));
    const char* pName = pdb->Name(idx);
    pdb->DeleteProperty(pName);
    _test(pdb->VariableCount(ELEMENT)==e);
}


void PropertyDatabase_Test::VariablesBoundaryTest(size_t b)
{
    _test(pdb->VariableCount(BOUNDARY)==b);

    csmp::Index idx=pdb->AddProperty("test variable 3","m s-1",VECTOR,BOUNDARY);
    _test(pdb->VariableCount(BOUNDARY)==(b+1));
    pdb->DeleteProperty(pdb->Name(idx));
    _test(pdb->VariableCount(BOUNDARY)==b);
}


void PropertyDatabase_Test::VariablesNodeTest(size_t n)
{
    _test(pdb->VariableCount(NODE)==n);
    csmp::Index idx=pdb->AddProperty("test variable 3","m s-1",VECTOR,NODE);
    _test(pdb->VariableCount(NODE)==(n+1));
    pdb->DeleteProperty(pdb->Name(idx));
    _test(pdb->VariableCount(NODE)==n);
}


void PropertyDatabase_Test::VariablesConstraintPointTest(size_t cp)
{
    _test(pdb->VariableCount(ELEMENT_INTEGRATION_POINT)==cp);
}


void PropertyDatabase_Test::VariablesFaceTest(size_t f)
{
    _test(pdb->VariableCount(FACE)==f);

    csmp::Index idx=pdb->AddProperty("test variable 3","m s-1",VECTOR,FACE);
    _test(pdb->VariableCount(FACE)==(f+1));
    pdb->DeleteProperty(pdb->Name(idx));
    _test(pdb->VariableCount(FACE)==f);

}

void PropertyDatabase_Test::VariablesInterFaceTest(size_t i)
{
    _test(pdb->VariableCount(INTER_FACE)==i);
}


void PropertyDatabase_Test::VariablesModelTest(size_t v)
{
    _test(pdb->VariableCount(MODEL)==v);

    csmp::Index idx=pdb->AddProperty("test variable 3","m s-1",VECTOR,MODEL);
    _test(pdb->VariableCount(MODEL)==(v+1));
    pdb->DeleteProperty(pdb->Name(idx));
    _test(pdb->VariableCount(MODEL)==v);
}


void PropertyDatabase_Test::run()
{
    FromFile();
    this->VariablesTotalTest(0);
    this->AddPropertyTest();
    this->VariablesElementTest(1);
    this->VariablesElementTest(1);
    this->VariablesBoundaryTest(0);
    this->VariablesConstraintPointTest(0);
    this->VariablesNodeTest(0);
    this->VariablesFaceTest(0);
    this->VariablesInterFaceTest(0);
    this->VariablesModelTest(0);
    this->VariablesTotalTest(2);
    this->IsDefinedTest();
    this->VariablesTotalTest(2);
    this->PlacementTest();
    this->TypeTest();
    this->DeletePropertyTest();
    pdb->Out();
}


<<<<<<< HEAD
void PropertyDatabase_Test::FromFile()
  {
    PropertyDatabase<3> pdb("CSMP-variables-vsTestLocked.txt");
  }
=======
    void PropertyDatabase_Test::PlacementTest()
    {
        _test(strcmp(csmp::parsePlacement(pdb->Placement("test variable")).c_str(),"REGION")==0);
        _test(strcmp(csmp::parsePlacement(pdb->Placement("test variable 2")).c_str(),"ELEMENT")==0);
        _test(strcmp(csmp::parsePlacement(pdb->Placement("test variable 2")).c_str(),"REGION")!=0);
    }


    void PropertyDatabase_Test::TypeTest()
    {
        _test(strcmp(csmp::parseType(pdb->Type("test variable 2")).c_str(),"SCALAR")==0);
        _test(strcmp(csmp::parseType(pdb->Type("test variable 2")).c_str(),"TENSOR")!=0);
    }


    void PropertyDatabase_Test::VariablesTotalTest(size_t s)
    {
        _test(pdb->VariableCount()==s);
    }


    void PropertyDatabase_Test::VariablesRegionTest(size_t r)
    {
        _test(pdb->VariableCount(REGION)==r);
    }


    void PropertyDatabase_Test::DeletePropertyTest()
    {
        pdb->DeleteProperty("test variable 2");
        _test(pdb->IsDefined("test variable 2")==0);
        this->VariablesTotalTest(1);
    }


    void PropertyDatabase_Test::VariablesElementTest(size_t e)
    {
        _test(pdb->VariableCount(ELEMENT)==e);
        csmp::Index idx=pdb->AddProperty("test variable 3","m s-1",VECTOR,ELEMENT);
        _test(pdb->VariableCount(ELEMENT)==(e+1));
        const char* pName = pdb->Name(idx);
        pdb->DeleteProperty(pName);
        _test(pdb->VariableCount(ELEMENT)==e);
    }


    void PropertyDatabase_Test::VariablesBoundaryTest(size_t b)
    {
        _test(pdb->VariableCount(BOUNDARY)==b);

        csmp::Index idx=pdb->AddProperty("test variable 3","m s-1",VECTOR,BOUNDARY);
        _test(pdb->VariableCount(BOUNDARY)==(b+1));
        pdb->DeleteProperty(pdb->Name(idx));
        _test(pdb->VariableCount(BOUNDARY)==b);
    }


    void PropertyDatabase_Test::VariablesNodeTest(size_t n)
    {
        _test(pdb->VariableCount(NODE)==n);
        csmp::Index idx=pdb->AddProperty("test variable 3","m s-1",VECTOR,NODE);
        _test(pdb->VariableCount(NODE)==(n+1));
        pdb->DeleteProperty(pdb->Name(idx));
        _test(pdb->VariableCount(NODE)==n);
    }


    void PropertyDatabase_Test::VariablesConstraintPointTest(size_t cp)
    {
        _test(pdb->VariableCount(ELEMENT_INTEGRATION_POINT)==cp);
    }


    void PropertyDatabase_Test::VariablesFaceTest(size_t f)
    {
        _test(pdb->VariableCount(FACE)==f);

        csmp::Index idx=pdb->AddProperty("test variable 3","m s-1",VECTOR,FACE);
        _test(pdb->VariableCount(FACE)==(f+1));
        pdb->DeleteProperty(pdb->Name(idx));
        _test(pdb->VariableCount(FACE)==f);

    }

    void PropertyDatabase_Test::VariablesInterFaceTest(size_t i)
    {
        _test(pdb->VariableCount(INTER_FACE)==i);
    }


    void PropertyDatabase_Test::VariablesModelTest(size_t v)
    {
        _test(pdb->VariableCount(MODEL)==v);

        csmp::Index idx=pdb->AddProperty("test variable 3","m s-1",VECTOR,MODEL);
        _test(pdb->VariableCount(MODEL)==(v+1));
        pdb->DeleteProperty(pdb->Name(idx));
        _test(pdb->VariableCount(MODEL)==v);
    }


    void PropertyDatabase_Test::run()
    {
        FromFile();
        this->VariablesTotalTest(0);
        this->AddPropertyTest();
        this->VariablesElementTest(1);
        this->VariablesElementTest(1);
        this->VariablesBoundaryTest(0);
        this->VariablesConstraintPointTest(0);
        this->VariablesNodeTest(0);
        this->VariablesFaceTest(0);
        this->VariablesInterFaceTest(0);
        this->VariablesModelTest(0);
        this->VariablesTotalTest(2);
        this->IsDefinedTest();
        this->VariablesTotalTest(2);
        this->PlacementTest();
        this->TypeTest();
        this->DeletePropertyTest();
        pdb->Out(getInfoStream());
    }


    void PropertyDatabase_Test::FromFile()
      {
        PropertyDatabase<3> pdb("CSMP-variables-vsTestLocked.txt");
      }
>>>>>>> b71117cb444b1539e747fa5855ab341062d3c4b3

}
