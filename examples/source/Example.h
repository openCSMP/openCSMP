#ifndef EXAMPLE_H
#define EXAMPLE_H

#include <string>
#include <list>
#include "CSMP_definitions.h"

namespace csmp{

class ExampleSuite;

/// Subclassing allows to run examples in ExamplesSuite
class Example
{

friend class ExampleSuite;

public:
  Example( std::ostream *osptr = 0 );
  virtual ~Example(){}

  // Run Function
  virtual void Run() = 0;
  // Describe your example
  virtual void Specifications() = 0;

  /// Reporting interface
  void SetTitle( const std::string& text ) { title_ = text; }
  void SetCategory( const std::string& text ) { category_ = text; }
  void SetDifficulty( const size_t& difficulty );
  void AddDescription( const std::string& text );
  void AddRequirement( const std::string& text ); // use this to tell which source file is used
  void AddAuthor( const std::string& text );

  // Return functions
  std::string  GetTitle() const;
  size_t       GetDifficulty() const;
  std::string  GetCategory() const;

protected:
  void Initialize();
  // Stream interface
  std::ostream *GetStream() { return ostream_; }
  void SetStream( std::ostream *ostream ) { ostream_ = ostream; }
  void CreateWorkingDirectoryAndCopyInputModelFiles(std::string& example_name, std::string& model_name,
                                                    std::string& variable_file, std::string config_file = "");
public:
  // DB
  std::list<std::string>::const_iterator GetDescriptionsBegin();
  std::list<std::string>::const_iterator GetDescriptionsEnd() const;
  std::list<std::string>::const_iterator GetAuthorsBegin();
  std::list<std::string>::const_iterator GetAuthorsEnd() const;
  std::list<std::string>::const_iterator GetRequirementsBegin();
  std::list<std::string>::const_iterator GetRequirementsEnd() const;

private:
  // output stream
  std::ostream *ostream_;
  // db
  std::string title_;
  std::string category_;
  std::list<std::string> descriptions_;
  std::list<std::string> requirements_;
  std::list<std::string> authors_;
  size_t difficulty_;                       ///< Example Difficulty: 1-5
  // disabled
  Example( const Example& );
  Example& operator = ( const Example& );
};

/**
@class Example Example "examples/Example.h"

@author P. Lang
@date 2010

@section examples Application Examples

@code
// DECLARATION
#include "Example.h"

namespace csmp {

class  YourExampleClass : public Example{
public:
  virtual void Run();
  virtual void Initialize();
};

} // csmp


// DEFINITION
#include "YourExampleClass.h"

#include "CSMP_definitions.h�
...

using namespace std;

namespace csmp{

void YourExampleClass::Initialize()
{
  SetTitle ( "CSMP Example 1" );
  SetDifficulty( 1 );
  AddAuthor ( "SKM" );
  AddDescription( "computes steady state pressure distribution" );
  AddDescription( "2D reservoir profile" );
  AddRequirement( "image(granite_model1)" );
  AddRequirement( "configuration(example1.txt)" );
} // Initialize()

void Example1ForSuite::Execute()
{
  ...your application code, or a simple function call for an external resource
} // Execute()

} // csmp

@endcode

*/

} // csmp

#endif // EXAMPLE_H
