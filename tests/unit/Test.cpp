// Test.cpp

#include "Test.h"
#include <iostream>
#include <cstdlib>
#include <typeinfo>     // Visual Studio requires /GR""




#ifdef _MSC_VER
//Allow return-less mains:
#pragma warning(disable: 4541)
#endif

using namespace std;

namespace csmp {

Test::Test(std::ostream* osptr)
  {
    m_osptr = osptr;
    m_nPass = m_nFail = 0;
  }

long Test::getNumPassed() const
  {
    return m_nPass;
  }

long Test::getNumFailed() const
  {
    return m_nFail;
  }

const std::ostream* Test::getStream() const
  {
    return m_osptr;
  }

void Test::setStream(std::ostream* osptr)
  {
    m_osptr = osptr;
  }

void Test::do_succeed()
  {
    ++m_nPass;
  }

void Test::reset()
  {
    m_nPass = m_nFail = 0;
  }



void Test::do_test( bool cond, const std::string& lbl,
                    const char* fname, long lineno )
{
    if (!cond)
        do_fail(lbl, fname, lineno);
    else
        _succeed();
}


// SKM 8/2/17: introduced appropriate name de-mangling
void Test::do_fail( const std::string& lbl,
                    const char* fname, long lineno )
{
    int   status;
    char* realname(nullptr);

    ++m_nFail;
    if (m_osptr)
      {
        if( !hasName() ) {
#ifdef __GNUC__
             //This is a fix for gcc name demangling.
             const std::type_info& ti = typeid(*this);
             realname = abi::__cxa_demangle(ti.name(), 0, 0, &status);
             *m_osptr << realname;
#else
             *m_osptr << typeid(*this).name();
#endif
          }
        else
          *m_osptr << (*this).getName();

        *m_osptr << " failure: (" << lbl << ") , " << fname
                 << " (line " << lineno << ")\n";
    }
}


// SKM 8/2/17: introduced appropriate name de-mangling
void Test::do_equal( double expr, double value, double tol,
                     const string& lbl, const char* fname, long lineno ) 
{
    bool  cond = ((expr >= value - tol) && (expr <= value + tol));
    int   status;
    char* realname(nullptr);

    if ( !cond ) {
        if( !hasName() ) {
#ifdef __GNUC__
             //This is a fix for gcc name demangling.
             const std::type_info& ti = typeid(*this);
             realname = abi::__cxa_demangle(ti.name(), 0, 0, &status);
             *m_osptr << realname;
#else
             *m_osptr << typeid(*this).name();
#endif
          }
        else
          *m_osptr << (*this).getName();

        *m_osptr << " failure: ( expr = " << expr <<" , value =  " << value<< ", diff = "<< expr-value <<", tolerance = "<< tol << ") , "
                 << fname << " (line " << lineno << ")\n";
      }

    do_test(cond, lbl, fname, lineno);                
}


// SKM 8/2/17: introduced appropriate name de-mangling
long Test::report() const
{
  int   status;
  char* realname(nullptr);

  if (m_osptr)
  {
    if( hasName() )
    {
      *m_osptr << "Test \""
          << getName() << "\":\n"
          << "\tPassed: " << m_nPass
          << "\tFailed: " << m_nFail
          << endl;
    }
    else
    {
      #ifdef __GNUC__
        //This is a fix for gcc name demangling.
        const std::type_info& ti = typeid(*this);
        realname = abi::__cxa_demangle(ti.name(), 0, 0, &status);
        *m_osptr << "Test \""
            << realname << "\":\n"
            << "\tPassed: " << m_nPass
            << "\tFailed: " << m_nFail
            << endl;
        free(realname);
      #else
        *m_osptr << "Test \""
            << typeid(*this).name() << "\":\n"
            << "\tPassed: " << m_nPass
            << "\tFailed: " << m_nFail
            << endl;
      #endif
    }
  }
  return m_nFail;
}


} // end namespace csmp
