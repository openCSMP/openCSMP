#ifndef TESTSUITE_H
#define TESTSUITE_H

#include "Test.h"   // includes <string>, <iosfwd>
#include <vector>
#include <stdexcept>
using std::string;
using std::ostream;
using std::vector;

namespace csmp {

class TestSuiteError : public std::logic_error 
{
public:
    TestSuiteError(const string& s = "")
        : logic_error(s)
    {}
    
};

class TestSuite
{
public:
    TestSuite(const string& name, ostream* osptr = 0);

    string getName() const;
    long getNumPassed() const;
    long getNumFailed() const;
    const ostream* getStream() const;
    void setStream(ostream* osptr);
    
    void addTest(Test* t);
    void addTest(const char* test_name, Test* t);
    void addTestSuite(const TestSuite&);

    void run();     // Calls Test::run() repeatedly
    void RunSpecificTest(const char*);
    void FreeAllButSpecificTest(const char* test_name);
    bool IsInSuite(const char* test_name);
    long report() const;
    void free();    // deletes tests

private:
    string m_name;
    ostream* m_osptr;
    vector<Test*> m_tests;
    void reset();

    // Disallowed ops:
    TestSuite(const TestSuite&);
    TestSuite& operator=(const TestSuite&);
};

inline
TestSuite::TestSuite(const string& name, ostream* osptr)
     : m_name(name)
{
    m_osptr = osptr;
}

inline
string TestSuite::getName() const
{
    return m_name;
}

inline
const ostream* TestSuite::getStream() const
{
    return m_osptr;
}

inline
void TestSuite::setStream(ostream* osptr)
{
    m_osptr = osptr;
}

} // end namespace csmp

#endif

