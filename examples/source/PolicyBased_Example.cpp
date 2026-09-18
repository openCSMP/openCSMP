// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "PolicyBased_Example.h"

using namespace std;

namespace csmp {


class PolicyOne
  {
  public:
    void doSomething1(){}
  };


class PolicyTwo
  {
  public:
    void doSomething1(){}
    void doSomething2(){}
  };


template<class Policy>
class Base : public Policy {};


class AsYouKnowIt
  {
  public:
    void Interface1(){}
    void Interface2(){}
  };


void PolicyBased_Example::Specifications()
  {
    SetTitle( "Policy-based class design" );
    SetDifficulty( 2 );
    SetCategory( "C++" );
    AddAuthor( "P. Lang" );
    AddDescription( "source in: PolicyBased_Example.cpp" );
    AddDescription( "brief overview of the concept" );
  }


void PolicyBased_Example::Run()
  {
    AsYouKnowIt asYouKnowIt;
    asYouKnowIt.Interface1();
    // ...

    Base<PolicyOne> singleInterface;
    singleInterface.doSomething1();
    
    // this doesn't work, compiler would nicely spot
    // singleInterface.doSomething2();

    Base<PolicyTwo> doubleInterface;
    doubleInterface.doSomething1();
    doubleInterface.doSomething2();

    //ostream& cout = *GetStream();
    cout << "\nCODE EXAMPLE\n";

  } // Run

} // csmp
