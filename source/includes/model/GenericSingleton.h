// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef GENERIC_SINGLETON_HPP
#define GENERIC_SINGLETON_HPP

#include <type_traits>

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

That the inherited class has a private constructor and cannot be copied is checked at compile time
when Instance() is used (static_assert).

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
      // Compile-time check of what makes the inherited class a singleton: the instance returned here must be
      // the only one. Objects of the class can therefore neither be created by other code nor be copied or moved.
      static_assert( !std::is_default_constructible_v<SingletonType>,
                     "a class inheriting from GenericSingleton must have a private constructor (and declare GenericSingleton<itself> a friend)" );
      static_assert( !std::is_copy_constructible_v<SingletonType> && !std::is_move_constructible_v<SingletonType>,
                     "a class inheriting from GenericSingleton must not be copy- or move-constructible" );
      static_assert( !std::is_copy_assignable_v<SingletonType> && !std::is_move_assignable_v<SingletonType>,
                     "a class inheriting from GenericSingleton must not be copy- or move-assignable" );

      // initialised once, on first use; this initialisation is thread-safe (since C++11) and the object
      // is shared by all translation units of a program
      static SingletonType instance;
      return instance;
    }

    // singletons are neither copied nor assigned; this also removes these operations from the inherited class
    GenericSingleton( const GenericSingleton& ) = delete;
    GenericSingleton& operator = ( const GenericSingleton& ) = delete;

  protected:
    GenericSingleton() = default;
};

} // csmp
 
 #endif
