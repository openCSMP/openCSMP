#ifndef CSMP_EXAMPLE_SUITE_H
#define CSMP_EXAMPLE_SUITE_H

#include "CSMP_definitions.h"
#include "Example.h"

namespace csmp {

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
class ExampleSuite
{
public:
  explicit ExampleSuite( const std::string &title, std::ostream *ostream = &std::cout );
  ExampleSuite( const ExampleSuite& ) = delete;
  ExampleSuite &operator = ( const ExampleSuite& ) = delete;

  // Property interface
  [[nodiscard]] std::string GetTitle() const { return title_; }
  [[nodiscard]] const std::ostream* GetStream() const { return ostream_; }
  void SetStream( std::ostream* ostream) { ostream_ = ostream; }

  // User interface
  void Run();                                   // runs listing and user choice
  void Free();                                  // removes examples
  size_t RegisterExample( Example *example );   // register new example

  // GUI interface
  typedef std::map<std::string,std::vector<Example*>* >::iterator category_iterator;
  category_iterator GetCategoriesBegin() { return categories_.begin(); }
  category_iterator GetCategoriesEnd() { return categories_.end(); }
  [[nodiscard]] std::string StripExampleTitle( std::string title ) const;
  void   OstreamDoubleUnderlined( std::string ) const;

private:
  // IO
  void   UIMenu() const;
  void   UICategoryMenu() const;
  void   UIExampleDetailsMenu( int exampleIndex,
                               std::map<std::string,std::vector<Example*>* >::const_iterator  ) const;
  void   UIExamplesMenu( std::map<std::string,std::vector<Example*>* >::const_iterator ) const;
  void   Initialize();
  void   CategoriesList() const;
  void   SuiteHeader() const;
  [[nodiscard]] int ChoiceWithinRange( long min, long max ) const;
  void   ExampleDetails( Example* example ) const;
  [[nodiscard]] bool   ChooseExampleAction() const;
  void   SortExampleVector( std::vector<Example*>& examples );
  void   ExamplesList( std::map<std::string,std::vector<Example*>* >::const_iterator ) const;
  void   OstreamUnderlined( std::string ) const;
  int    ChooseExample( const std::vector<Example*>* ) const;
  void   WriteExamplesFile();
  [[nodiscard]] std::map<std::string,std::vector<Example*>* >::const_iterator ChooseCategory() const;
  // DB
  std::string title_;                                         ///< Title of suite
  std::ostream* ostream_;                                     ///< output stream suite writes messages to
  std::map<std::string,std::vector<Example*>* > categories_;  ///< registerd examples of different categories
};

} // csmp

#endif // CSMP_EXAMPLE_SUITE_H
