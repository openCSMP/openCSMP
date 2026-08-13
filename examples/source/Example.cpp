#include "Example.h"

// for creating and changing directories etc.
#include <filesystem>

#ifdef __GNUC__
#include "cxxabi.h"
#endif

using namespace std;

namespace csmp {

/// base class constructor, sets ostream for example output
Example::Example( ostream *osptr )
  : ostream_(osptr), difficulty_(0) {
}

/// set function, takes care that 0 < difficulty < 4
void Example::SetDifficulty( const int& difficulty )
{
  difficulty_ = difficulty;
  difficulty_ = difficulty_ < 1 ? 1 : difficulty_;
  difficulty_ = difficulty_ > 5 ? 5 : difficulty_;
}

/// allows to add a descrtiption to an example
void Example::AddDescription( const std::string& text )
{
  descriptions_.push_back( text );
}

/// allows to add a requirement to an example (files etc..)
void Example::AddRequirement( const std::string& text )
{
  requirements_.push_back( text );
}

/// allows to add an author to an example
void Example::AddAuthor( const std::string& text )
{
  authors_.push_back( text );
}


/// example title return function
string Example::GetTitle() const
{
  string text(title_);
  text +=" (";
  
#ifdef __GNUC__
  //This is a fix for gcc name demangling.
  int   status;
  const auto &ti = typeid(*this);
  const char* realname = abi::__cxa_demangle(ti.name(), nullptr, nullptr, &status);
  text += realname;
#else
  text += typeid(*this).name();
#endif
  text +=")";

  return text;
}

/// return function
int Example::GetDifficulty() const
{
  return difficulty_;
}

/// return function
string Example::GetCategory() const
{
  return category_;
}

/// assures db to contain data
void Example::Initialize()
{
  Specifications();
  difficulty_ = difficulty_ < 1 || difficulty_ > 5 ? 1 : difficulty_;
  if( category_.empty() )
    category_ = "General";
  if( descriptions_.empty() )
    descriptions_.emplace_back("none" );
  if( requirements_.empty() )
    requirements_.emplace_back("none" );
  if( authors_.empty() )
    authors_.emplace_back("unknown" );
}

string Example::GetExampleFileName(const char* path)
{
  const filesystem::path p(path);
  return p.stem().string();
}

/**
   Creates new directory "example_outputs/" side by side with the working directory "example_inputs", moving any file writing by the current example into there.

   @attention This function expects that the directory that contains the executable is    open-csmp/examples/example_inputs/

   Context: to run the examples, you need to create a new directory somewhere outside of your repository and copy the example_inputs/  directory from the repository in there,
   including its subdirectories. Then set this example_inputs/  directory as your working directory (XCode:  Product -> Scheme -> Edit Scheme -> custom working directory in the Options tab)

   When you create a CSMP native binary, it gets  placed into the directory  example_inputs/csmp_native_format_models/
   Then, when you run any of the examples, CSMP will create a new directory  example_outputs/
   in your new directory side-by-side with example_inputs/ and this is where you will find the results of your computations.
   All this is done by this method.
*/
void Example::CreateWorkingDirectoryAndCopyInputModelFiles( const string& example_name, const string& model_name,
                                                            const string& variable_file, string config_file )
{
  //create of directory with current example name and go into this directory
  string current_path = filesystem::current_path().parent_path().string();
  filesystem::current_path(filesystem::path(current_path));
  filesystem::create_directory(filesystem::path("example_outputs"));
  filesystem::current_path(filesystem::path("example_outputs"));
  if(filesystem::is_directory(filesystem::path(example_name))) filesystem::remove_all(filesystem::path(example_name)); //if directory already exists, delete it
  filesystem::create_directory(example_name);
  filesystem::current_path(filesystem::path(example_name));

  //copy input files into working directory
  if(!model_name.empty()) {
    string input_directory = filesystem::current_path().parent_path().parent_path().string();
    input_directory += "/example_inputs/csmp_native_format_models/";
    //.vset
    string path = "../../example_inputs/csmp_native_format_models/";
    string name = model_name + ".vset";
    string file_name = path + name;
    if (filesystem::exists(filesystem::path(file_name))) filesystem::copy(filesystem::path(file_name), filesystem::path("./"));
    else {
      string error_message = "\n\nError: file '";
      error_message += name + "' does not exist in directory "  + input_directory;
      error_message += ", example cannot run, please check.\n";
      throw std::runtime_error(error_message);
    }
    //boundaries.dat
    name = model_name + "_boundaries.dat";
    file_name = path + name;
    if (filesystem::exists(filesystem::path(file_name))) filesystem::copy(filesystem::path(file_name), filesystem::path("./"));
    else {
      string error_message = "\n\nError: file '";
      error_message += name + "' does not exist in directory "  + input_directory;
      error_message += ", example cannot run, please check.\n";
      throw std::runtime_error(error_message);
    }
    //regions.dat
    name = model_name +  "_regions.dat";
    file_name = path + name;
    if (filesystem::exists(filesystem::path(file_name))) filesystem::copy(filesystem::path(file_name), filesystem::path("./"));
    else {
      string error_message = "\n\nError: file '";
      error_message += name + "' does not exist in directory "  + input_directory;
      error_message += ", example cannot run, please check.\n";
      throw std::runtime_error(error_message);
    }
    //variables.dat
    name = model_name  + "_variables.dat";
    file_name = path + name;
    if (filesystem::exists(filesystem::path(file_name))) filesystem::copy(filesystem::path(file_name), filesystem::path("./"));
    else {
      string error_message = "\n\nError: file '";
      error_message += name + "' does not exist in directory "  + input_directory;
      error_message += ", example cannot run, please check.\n";
      throw std::runtime_error(error_message);
    }
    //splitboundaries.dat - optional
    name = model_name  + "_splitboundaries.dat";
    file_name = path + name;
    if (filesystem::exists(filesystem::path(file_name))) filesystem::copy(filesystem::path(file_name), filesystem::path("./"));
  }

  //copy variable file
  if(!variable_file.empty()) {
    string input_directory = filesystem::current_path().parent_path().parent_path().string();
    input_directory += "/example_inputs/variables_and_configuration_files/";

    string path = "../../example_inputs/variables_and_configuration_files/";
    if (string file_name = path + variable_file; filesystem::exists(filesystem::path(file_name))) filesystem::copy(filesystem::path(file_name), filesystem::path("./"));
    else {
      string error_message = "\n\nError: file '";
      error_message += variable_file + "' does not exist in directory "  + input_directory;
      error_message += ", example cannot run, please check.\n";
      throw std::runtime_error(error_message);
    }
  }

  //copy configuration file if required
  if(!config_file.empty()) {
    string input_directory = filesystem::current_path().parent_path().parent_path().string();
    input_directory += "/example_inputs/variables_and_configuration_files/";

    string path = "../../example_inputs/variables_and_configuration_files/";
    config_file += "-configuration.txt";
    if(string file_name = path + config_file; filesystem::exists(filesystem::path(file_name))) filesystem::copy(filesystem::path(file_name), filesystem::path("./"));
    else {
      string error_message = "\n\nError: file '";
      error_message += config_file + "' does not exist in directory "  + input_directory;
      error_message +=", example cannot run, please check.\n";
      throw std::runtime_error(error_message);
    }
  }
}

/// return function
list<string>::const_iterator Example::GetDescriptionsBegin()
{
  return descriptions_.begin();
}

/// return function
list<string>::const_iterator Example::GetDescriptionsEnd() const
{
  return descriptions_.end();
}

/// return function
list<string>::const_iterator Example::GetAuthorsBegin()
{
  return authors_.begin();
}

/// return function
list<string>::const_iterator Example::GetAuthorsEnd() const
{
  return authors_.end();
}


/// return function
list<string>::const_iterator Example::GetRequirementsBegin()
{
  return requirements_.begin();
}


/// return function
list<string>::const_iterator Example::GetRequirementsEnd() const
{
  return requirements_.end();
}


} // csmp
