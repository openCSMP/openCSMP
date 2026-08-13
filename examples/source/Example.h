#ifndef EXAMPLE_H
#define EXAMPLE_H

#include "CSMP_definitions.h"

namespace csmp {

class ExampleSuite;

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

#include "CSMP_definitions.h"
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

*/class Example
{

friend class ExampleSuite;

public:
  explicit Example( std::ostream *osptr = nullptr );
  Example( const Example& ) = delete;
  Example& operator=( const Example& ) = delete;
  virtual ~Example() = default;

  // Run Function
  virtual void Run() = 0;
  // Describe your example
  virtual void Specifications() = 0;

  /// Reporting interface
  void SetTitle( const std::string& text ) { title_ = text; }
  void SetCategory( const std::string& text ) { category_ = text; }
  void SetDifficulty( const int& difficulty );
  void AddDescription( const std::string& text );
  void AddRequirement( const std::string& text ); // use this to tell which source file is used
  void AddAuthor( const std::string& text );

  // Return functions
  [[nodiscard]] std::string  GetTitle() const;
  [[nodiscard]] int          GetDifficulty() const;
  [[nodiscard]] std::string  GetCategory() const;

protected:
  void Initialize();
  // Stream interface
  [[nodiscard]] std::ostream *GetStream() const { return ostream_; }
  void SetStream( std::ostream *ostream ) { ostream_ = ostream; }

static std::string GetExampleFileName(const char* path);
  static void CreateWorkingDirectoryAndCopyInputModelFiles( const std::string& example_name, const std::string& model_name,
                                                            const std::string& variable_file, std::string config_file ="" );
public:
  // DB
  std::list<std::string>::const_iterator GetDescriptionsBegin();
  [[nodiscard]] std::list<std::string>::const_iterator GetDescriptionsEnd() const;
  std::list<std::string>::const_iterator GetAuthorsBegin();
  [[nodiscard]] std::list<std::string>::const_iterator GetAuthorsEnd() const;
  std::list<std::string>::const_iterator GetRequirementsBegin();
  [[nodiscard]] std::list<std::string>::const_iterator GetRequirementsEnd() const;

private:
  // output stream
  std::ostream *ostream_;
  // db
  std::string title_;
  std::string category_;
  std::list<std::string> descriptions_;
  std::list<std::string> requirements_;
  std::list<std::string> authors_;
  int difficulty_;                       ///< Example Difficulty: 1-5
 };


} // csmp

#endif // EXAMPLE_H
