#ifndef EXAMPLESUITE_H
#define EXAMPLESUITE_H

#include <map>
#include <vector>
#include <ostream>

#include "CSMP_definitions.h"
#include "Example.h"

namespace csmp{

/// A suite to contain a number of examples with UI
class ExampleSuite
{
public:
  ExampleSuite( const std::string &title, std::ostream *ostream = &std::cout );

  // Property interface
  std::string GetTitle() const { return title_; }
  const std::ostream* GetStream() const { return ostream_; }
  void SetStream( std::ostream* ostream) { ostream_ = ostream; }

  // User interface
  void Run();                                   // runs listing and user choice
  void Free();                                  // removes examples
  size_t RegisterExample( Example *example );   // register new example

  // GUI interface
  typedef std::map<std::string,std::vector<Example*>* >::iterator category_iterator;
  category_iterator GetCategoriesBegin() { return categories_.begin(); }
  category_iterator GetCategoriesEnd() { return categories_.end(); }
  std::string StripExampleTitle( std::string title ) const;
  void   OstreamDoubleUnderlined( std::string ) const;


private:
  // IO
  void   UIMenu() const;
  void   UICategoryMenu() const;
  void   UIExampleDetailsMenu( size_t exampleIndex,
                               std::map<std::string,std::vector<Example*>* >::const_iterator  ) const;
  void   UIExamplesMenu( std::map<std::string,std::vector<Example*>* >::const_iterator ) const;
  void   Initialize();
  void   CategoriesList() const;
  void   SuiteHeader() const;
  int    ChoiceWithinRange( double min, double max ) const;
  void   ExampleDetails( Example* example ) const;
  bool   ChooseExampleAction() const;
  void   SortExampleVector( std::vector<Example*>& examples );
  void   ExamplesList( std::map<std::string,std::vector<Example*>* >::const_iterator ) const;
  void   OstreamUnderlined( std::string ) const;
  size_t ChooseExample( const std::vector<Example*>* ) const;
  void   WriteExamplesFile();
  std::map<std::string,std::vector<Example*>* >::const_iterator ChooseCategory() const;
  // DB
  std::string title_;                                         ///< Title of suite
  std::ostream* ostream_;                                     ///< output stream suite writes messages to
  std::map<std::string,std::vector<Example*>* > categories_;  ///< registerd examples of different categories

  // Disabled
  ExampleSuite( const ExampleSuite& );
  ExampleSuite &operator = ( const ExampleSuite& );
};

/**
@class ExampleSuite ExampleSuite "examples/ExampleSuite.h"

@author P. Lang
@date 2010

@section implementation Implementation

Contains a map which holds one entry for each example category (Example.h --> category_)
This map entry holds a pointer to a vector which in turn contains pointers to all
examples registered to the suite of the respective category.

@section examples Application Examples

@code
...
// instantiating ExampleSuite
ExampleSuite exampleSuite( "CSMP EXAMPLES", &cout );

// register new example
exampleSuite.RegisterExample( new Example1ForSuite() );
...
exampleSuite.RegisterExample ( new YourExampleClass() );

// execute examples in order of registration
exampleSuite.Execute();
...
@endcode
*/

} // csmp


#endif // EXAMPLESUITE_H
