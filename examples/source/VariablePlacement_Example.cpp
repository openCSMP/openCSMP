#include "VariablePlacement_Example.h"

#include "ANSYS_Model3D.h"
#include "Region.h"
#include "InputDataManager.h"
#include "PropertyHandle.h"

#define DIM 3U

using namespace std;

namespace csmp {


  // POLICIES
  template<uint32_t dim>
  class FromElementToElement{
  public:
    void Store( Element<dim>* e, const Index& i, double v )
      { cache() = v; e->Store( i,cache ); }
  private:
    ScalarVariable cache;
  };

  template<uint32_t dim>
  class FromElementToNodes{
  public:
    void Store( Element<dim>* e, const Index& i, double v )
      {
        cache() = v;
        for( typename vector<Node<dim>*>::const_iterator it = e->NodesBegin(); it != e->NodesEnd(); ++it )
          (*it)->Store( i,cache );
      }
  private:
    ScalarVariable cache;
  };

  // 'BASE' CLASS
  template<template<uint32_t> class StoragePolicy,template<uint32_t> class SIMPLEX, uint32_t dim = 3U>
  class GenericStore : public StoragePolicy<dim>
  {
  public:
    void operator () ( SIMPLEX<dim>* s, const Index& i, double v )
      { this->Store( s,i,v );  }
  };


void VariablePlacement_Example::Specifications()
{
  SetTitle( "Variable Placement" );
  SetDifficulty( 2 );
  SetCategory( "Simulation of Physical Processes" );
  AddAuthor( "P. Lang" );
  SetCategory( "Software Functionality" );
  AddDescription( "source in: VariablePlacement_Example.cpp" );
  AddDescription( "CSMP: properties on the fly, property access " );
  AddDescription( "C++: Template template parameters, Policy Based Classes" );
  AddDescription( "C++: Templatized Function Objects" );
  AddRequirement( "'fracs4' .dat, .asc, -regions & -configuration .txt" );
  AddRequirement( "CSMP-1phase-variables.txt" );
} 


void VariablePlacement_Example::Run()
{
  // ESTABLISHING OUTPUTSTREAM FROM BASECLASS
    //ostream &cout = *GetStream();

  // MODEL SETUP & AUXILIARIES
    // build model based on ansys .dat & .asc file...
    string modelName = "fracs4";
    ANSYS_Model3D  model( modelName.c_str(), "CSMP-1phase-variables.txt" );
    // ...and configuring it based on -configuration.txt & -regions.txt
    InputDataManager<DIM>().ConfigureFromFile( model, modelName.c_str(),
                                               false, true, true, true, false );

  // TASK OF THIS EXAMPLE
    // we create variables on the fly
    PropertyHandle<DIM>* elementProperty = new PropertyHandle<DIM>( model, "property element", SCALAR, ELEMENT );
    PropertyHandle<DIM>* nodalProperty = new PropertyHandle<DIM>( model, "property node", SCALAR, NODE );
    //                   ^^^^^^^^^^^^^^^^^^^^^ handles live on the heap. Why? Read on...

    // we instantiate indizes that later on allow for access to the given properties
    Index elementKey( model.Database().StorageKey( "property element") );
    Index nodeKey( model.Database().StorageKey( "property node") );


    GenericStore<FromElementToElement,Element> storeToElement;
    GenericStore<FromElementToNodes,Element> storeToNode;

    // we loop over all elements to access the created properties
    Region<DIM>& rref = model.Region( "Model" );
    const vector<Element<DIM>*>::const_iterator elementsEnd = rref.CellsEnd();
    for( vector<Element<DIM>*>::const_iterator it = rref.CellsBegin(); it != elementsEnd; ++ it)
    {
      // now we may use the function object like an ordinary function
      storeToElement( *it, elementKey, 5. );
      storeToNode( *it, nodeKey, 10. );
      // which does the same as accessing the inherited interface
      storeToElement.Store( *it, elementKey, 5. );
      storeToNode.Store( *it, nodeKey, 10. );
    }

    // since calling the destructor of a property handle causes the removal
    // of the created property from the database, we are able(within the same
    // scope, to remove the two properties and free some space during runtime
    delete elementProperty;
    delete nodalProperty;

    // NOW WHAT HAVE WE ACHIEVED?
    // in this very example, we have complicated a trivial task.
    // storing, or reading for that matter, may be achieved
    // with a lot less effort and much cleaner by e.g trivial function
    // overloading.
    // the intention was to provide a beginner's introduction into
    // policy based classes, their implementation and features
    // within a csmp related application. Also, we had a first glance
    // at function objects .


} // Run()



} // csmp
