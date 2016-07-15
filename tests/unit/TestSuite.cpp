// TestSuite.cpp

#include "TestSuite.h"
#include <iostream>
#include <cassert>

using namespace std;

namespace csmp {

void TestSuite::addTest(Test* t) throw(TestSuiteError)
{
    // Make sure test has a stream:
    if (t == 0)
        throw TestSuiteError( "Null test in TestSuite::addTest" );
    else if (m_osptr != 0 && t->getStream() == 0)
        t->setStream(m_osptr);

    m_tests.push_back(t);
    t->reset();
}

void TestSuite::addTest(const char* test_name, Test* t) throw(TestSuiteError)
{
    // Make sure test has a stream:
    if (t == 0)
        throw TestSuiteError( "Null test in TestSuite::addTest" );
    else if (m_osptr != 0 && t->getStream() == 0)
        t->setStream(m_osptr);

    m_tests.push_back(t);
    t->reset();
    t->setName(test_name);
}

void TestSuite::addTestSuite(const TestSuite& s) throw(TestSuiteError)
{
    for (size_t i = 0; i < s.m_tests.size(); ++i)
        addTest(s.m_tests[i]);
}

void TestSuite::FreeAllButSpecificTest(const char* test_name)
{
    // This is not a destructor because tests
    // don't have to be on the heap.
    for (size_t i = 0; i < m_tests.size(); ++i)
    {
        if (!(string(test_name)==m_tests[i]->getName()))
        {
            delete m_tests[i];
            m_tests[i] = 0;
        }
    }
}

void TestSuite::free()
{
    // This is not a destructor because tests
    // don't have to be on the heap.
    for (size_t i = 0; i < m_tests.size(); ++i)
    {
        delete m_tests[i];
        m_tests[i] = 0;
    }
}


/// @todo (1-D) Put exception catch around each test and report as failed if throw, instead of throwing all the way out of this member, halting test routine
void TestSuite::run()
{
    int   status;
    char* realname;
    reset();
    for (size_t i = 0; i < m_tests.size(); ++i)
    {
        assert(m_tests[i]);
        if (m_tests[i]->hasName())
            cout<<"Running Test: "<<m_tests[i]->getName()<<endl;
        else {
          #ifdef __GNUC__
            //This is a fix for gcc name demangling.
            const std::type_info  &ti = typeid(*m_tests[i]);
            realname = abi::__cxa_demangle(ti.name(), 0, 0, &status);
            cout<<"Running Test: "<<realname<<endl;
          #else
            cout<<"Running Test: "<<typeid(*m_tests[i]).name()<<endl;
          #endif
        }
        m_tests[i]->run();
    }
}

void TestSuite::RunSpecificTest(const char* test_name)
{

    reset();
    for (size_t i = 0; i < m_tests.size(); ++i)
    {
        assert(m_tests[i]);  //just a check
        if (string(test_name)==m_tests[i]->getName())
        {
            m_tests[i]->run();
            break; // only one test with said name should run. Names should not be repeated!
        }
    }
}

bool TestSuite::IsInSuite( const char* test_name )
{
    for (size_t i = 0; i < m_tests.size(); ++i)
     {
        if (m_tests[i])
          {
            if (string(test_name)==m_tests[i]->getName())
              {
                return true;
              }
            else {
                return false;
              }
         }
     }
    return false;
}


long TestSuite::report() const
{
    if (m_osptr)
    {
        long totFail = 0;
        *m_osptr << "\nTestSuite \"" << m_name << "\"\n=======";
        size_t i;
        for (i = 0; i < m_name.size(); ++i)
            *m_osptr << '=';
        *m_osptr << "=\n";

        for (i = 0; i < m_tests.size(); ++i)
        {
            //assert(m_tests[i]);
            if (m_tests[i])
                totFail += m_tests[i]->report();
        }

        *m_osptr << "=======";
        for (i = 0; i < m_name.size(); ++i)
            *m_osptr << '=';
        *m_osptr << "=\n";
        return totFail;
    }
    else
        return getNumFailed();
}



long TestSuite::getNumPassed() const
{
    long totPass = 0;
    for (size_t i = 0; i < m_tests.size(); ++i)
    {
        assert(m_tests[i]);
        totPass += m_tests[i]->getNumPassed();
    }
    return totPass;
}

long TestSuite::getNumFailed() const
{
    long totFail = 0;
    for (size_t i = 0; i < m_tests.size(); ++i)
    {
        assert(m_tests[i]);
        totFail += m_tests[i]->getNumFailed();
    }
    return totFail;
}

void TestSuite::reset()
{
    for (size_t i = 0; i < m_tests.size(); ++i)
    {
        assert(m_tests[i]);
        m_tests[i]->reset();
    }
}

} // end namespace csmp
