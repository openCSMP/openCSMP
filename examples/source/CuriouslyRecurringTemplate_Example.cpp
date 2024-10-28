#include "CuriouslyRecurringTemplate_Example.h"

using namespace std;

namespace csmp {

/// base class
template<typename T>
class Base {
  public:
    void CallChildClass() { static_cast<T>(*this).ThisIsTheCatch(); }
};


/// sub class (demonstrating compile time polymorphism)
class Sub : public Base<Sub> {
  public:
    void ThisIsTheCatch() {}
};


void CuriouslyRecurringTemplate_Example::Specifications()
{
  SetTitle( "Curiously recurring template pattern" );
  SetDifficulty( 2 );
  SetCategory( "C++" );
  AddAuthor( "P. Lang" );
  AddDescription( "source in: CuriouslyRecurringTemplate_Example.cpp" );
  AddDescription( "combined with safe-by-design static cast" );
}


/**
     Step through this example with the debugger.
*/
void CuriouslyRecurringTemplate_Example::Run()
 {
    Sub  sub;
    //ostream& cout = *GetStream();
    cout <<"\n"<<"CuriouslyRecurringTemplate_Example::Run: object type is: " << typeid(sub).name() << endl;
    cout << "\nCODE EXAMPLE\n";

 } // Run

} // csmp
