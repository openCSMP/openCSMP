#ifndef TEST_H
#define TEST_H

#include <string>
#include <iosfwd>
#ifdef __GNUC__
#include "cxxabi.h"
#endif

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
#define _warn(str) std::cout << str << '\n'
#define _info(str) std::cout << str << '\n'

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
    Test( std::ostream* osptr = 0 );
    virtual ~Test(){}
    virtual void run() = 0;

    long getNumPassed() const;
    long getNumFailed() const;
    const std::ostream* getStream() const;
    void setStream(std::ostream* osptr);
    void setName( std::string testName ) { testName_ = testName; }
    std::string getName() const { return testName_; }
    bool hasName() const { return !testName_.empty(); }
    
    void do_succeed();
    long report() const;
    virtual void reset();


  protected:
    void do_equal( double expr, double value, double tol, 
                   const std::string& lbl, const char* fname, long lineno );
                  
    void do_test( bool cond, const std::string& lbl,
                  const char* fname, long lineno );
                 
    void do_fail( const std::string& lbl,
                  const char* fname, long lineno );
    const char* prefix_;
  private:
    std::ostream* m_osptr;
    std::string testName_;
    long m_nPass;
    long m_nFail;
    // Disallowed:
    Test(const Test&);
    Test& operator=(const Test&);
};

} // end namespace csmp

#endif

