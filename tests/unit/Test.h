#ifndef TEST_H
#define TEST_H

#include <string>
#include <iosfwd>
#ifdef __GNUC__
#include "cxxabi.h"
#endif

using std::string;
using std::ostream;

// The following have underscores because they are macros
// (and it's impolite to usurp other users' functions!).

#ifdef RUNNING_UNDER_CATCH

#include "catch.hpp"

#define _test(cond) REQUIRE( (cond) )
#define _fail(e) FAIL( e )
#define _equal(expr,value,tol) REQUIRE( (expr) == Approx( (value) ).epsilon( (tol) ) )
#define _succeed()
#define _warn(e) WARN( e )
#define _info(e) INFO( e )

#else

#define _test(cond) do_test(cond, #cond, __FILE__, __LINE__)
#define _fail(str) do_fail(str, __FILE__, __LINE__)
#define _equal(expr,value,tol) do_equal(expr, value, tol, #expr " == " #value, __FILE__, __LINE__)
#define _succeed() do_succeed()
#define _warn(str) std::cout << e << '\n'
#define _info(str) std::cout << e << '\n'

#endif


namespace csmp {

/** Test  base class for unit testing.

    Usage example:
      1. define default constructor and run() function
      2. in run() function, use the macros provided to test whether desired conditions
         apply or not. 
      
      Examples:
         - check whether the floating point value, a, is equal to b within given tolerance, tol:
         
           _equal(a,b,tol);
        
         - check whether a boolean conditions applies:
         
           _test(isoparametric==false);
          
         - add a message to a failure report:
         
           _fail("for this input this did not work");
           
         - report that a test passed: if (...) _succeed();
    
     3. Obtain a fail/succeed report after running the test.
     
        test.report();
*/
class Test
  {
  public:
    Test(ostream* osptr = 0);
    virtual ~Test(){}
    virtual void run() = 0;

#ifdef RUNNING_UNDER_CATCH
    
    long getNumPassed() const;
    long getNumFailed() const;
    const ostream* getStream() const;
    void setStream(ostream* osptr);
    void setName( string testName ) { testName_ = testName; }
    string getName() const { return testName_; }
    bool hasName() const { return !testName_.empty(); }
    
    void do_succeed();
    long report() const;
    virtual void reset();


  protected:
    void do_equal( double expr, double value, double tol, 
                   const string& lbl, const char* fname, long lineno );
                  
    void do_test( bool cond, const string& lbl,
                  const char* fname, long lineno );
                 
    void do_fail( const string& lbl,
                  const char* fname, long lineno );
    const char* prefix_;
  private:
    ostream* m_osptr;
    string testName_;
    long m_nPass;
    long m_nFail;
    // Disallowed:
    Test(const Test&);
    Test& operator=(const Test&);
#endif
};

#ifdef RUNNING_UNDER_CATCH
  inline
      Test::Test(ostream* osptr)
  {
    m_osptr = osptr;
    m_nPass = m_nFail = 0;
  }

  inline
      long Test::getNumPassed() const
  {
    return m_nPass;
  }

  inline
      long Test::getNumFailed() const
  {
    return m_nFail;
  }

  inline
     const ostream* Test::getStream() const
  {
    return m_osptr;
  }

  inline
      void Test::setStream(ostream* osptr)
  {
    m_osptr = osptr;
  }

  inline
      void Test::do_succeed()
  {
    ++m_nPass;
  }

  inline
      void Test::reset()
  {
    m_nPass = m_nFail = 0;
  }
#endif

} // end namespace csmp

#endif

