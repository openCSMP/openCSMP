#ifndef GENERIC_SINGLETON_HPP
#define GENERIC_SINGLETON_HPP

namespace csmp {

/**
 
@brief references based generic singleton base class template

@author P. Lang
@date 2011

Class template that allows to inherit from with singleton(GoF) functionality.
The inherited class should look as outlined:

@code
class TestSingleton : public GenericSingleton<TestSingleton>
{
  friend class GenericSingleton<TestSingleton>;
 
  private:
    TestSingleton() : i(3) {}

  public:
    int i;
};
@endcode

and is consequently (globally!!!) available. 
For usage see:

@code
TestSingleton& single( TestSingleton::Instance() );
single.i = 4;
@endcode

alternatively

@code
TestSingleton::Instance().i = 4;
@endcode

@attention Strictly speaking singletons(and statics for that matter) are considered to be a violation of OOP and should be treated as such.

*/
template<class CT>
class GenericSingleton
{
  // type of inherited class(child type)
  typedef CT SingletonType;

  public:
    static SingletonType& Instance()
    {
      static SingletonType instance;
      return instance;
    }

  protected:
    GenericSingleton() {}
 
  private:
    GenericSingleton( const GenericSingleton& );
    GenericSingleton& operator = ( const GenericSingleton& ) { return *this; }
};

} // csmp
 
 #endif
