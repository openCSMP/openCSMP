#ifndef PROPERTY_DATABASE_TEST_H
#define PROPERTY_DATABASE_TEST_H
#include "Test.h"

namespace csmp
{
    class Index;
    template<size_t> class PropertyDatabase;
    class PropertyDatabase_Test : public Test
    {
    public:
        PropertyDatabase_Test();
        ~PropertyDatabase_Test();
        void AddPropertyNotFromCommandLineTest();
        //Tests
        void AddPropertyTest();
        void IsDefinedTest();
        void PlacementTest();
        void TypeTest();
        void VariablesTotalTest(size_t);
        void DeletePropertyTest();
        void NameTest(csmp::Index, const char*);
        void UnitTest(csmp::Index, const char*);
        void VariablesElementTest(size_t);
        void VariablesRegionTest(size_t);
        void VariablesBoundaryTest(size_t);
        void VariablesNodeTest(size_t);
        void VariablesConstraintPointTest(size_t);
        void VariablesFaceTest(size_t);
        void VariablesInterFaceTest(size_t);
        void VariablesModelTest(size_t);
        void FromFile();

        void run();
    private:
        PropertyDatabase<3>* pdb;
    };
}

#endif // PROPERTYDATABASE_TEST_H
