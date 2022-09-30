#include <iostream>
#include <cstdlib>
#include <typeinfo>

#ifdef __GNUC__
#include "cxxabi.h"
#endif

#include "Example.h"

#include <filesystem>
namespace fs = std::filesystem;

using namespace std;

namespace csmp{

/// base class constructor, sets ostream for example output
Example::Example( ostream *ostream )
  : ostream_( ostream )
{
}

/// assures db to contain data
void Example::Initialize()
{
  Specifications();
  difficulty_ = (difficulty_ < 1 || difficulty_ > 5) ? 1 : difficulty_;
  if( category_.empty() )
    category_ = "General";
  if( descriptions_.empty() )
    descriptions_.push_back( "none" );
  if( requirements_.empty() )
    requirements_.push_back( "none" );
  if( authors_.empty() )
    authors_.push_back( "unknown" );
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
   const std::type_info& ti = typeid(*this);
   char* realname = abi::__cxa_demangle(ti.name(), 0, 0, &status);
   text += realname;
#else
  text += typeid(*this).name();
#endif
  text +=")";

  return text;
}

/// return function
size_t Example::GetDifficulty() const
{
  return difficulty_;
}

/// return function
string Example::GetCategory() const
{
  return category_;
}

/// set function, takes care that 0 < difficulty < 4
void Example::SetDifficulty( const size_t& difficulty )
{
  difficulty_ = difficulty;
  difficulty_ = difficulty_ < 1 ? 1 : difficulty_;
  difficulty_ = difficulty_ > 5 ? 5 : difficulty_;
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
list<string>::const_iterator Example::GetRequirementsBegin()
{
  return requirements_.begin();
}

/// return function
list<string>::const_iterator Example::GetRequirementsEnd() const
{
  return requirements_.end();
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


void Example::CreateWorkingDirectoryAndCopyInputModelFiles(std::string& example_name, std::string& model_name,
                                                           std::string& variable_file, std::string config_file)
{
  //create of directory with current example name and go into this directory
  if(fs::is_directory(example_name)) fs::remove_all(example_name); //if directory already exists, delete it
  fs::create_directory(example_name);
  fs::current_path(example_name);

  //copy input files into working directory
  string file_name_pre = "../";
  file_name_pre += model_name;
  string file_name = file_name_pre + ".vset";
  if(fs::exists(file_name)) fs::copy(file_name, "./");
  else {
    file_name.erase(0,3);
    cerr<<"\n\nError: file '"<<file_name<<"' does not exist in directory "<<fs::current_path().parent_path();
    cerr<<", example cannot run, please check. Did you run the CSMPInterfaces_example first?"<<endl;
    return;
  }
  file_name = file_name_pre + "_boundaries.dat";
  if(fs::exists(file_name)) fs::copy(file_name, "./");
  else {
    file_name.erase(0,3);
    cerr<<"\n\nError: file '"<<file_name<<"' does not exist in directory "<<fs::current_path().parent_path();
    cerr<<", example cannot run, please check. Did you run the CSMPInterfaces_example first?"<<endl;
    return;
  }
  file_name = file_name_pre + "_regions.dat";
  if(fs::exists(file_name)) fs::copy(file_name, "./");
  else {
    file_name.erase(0,3);
    cerr<<"\n\nError: file '"<<file_name<<"' does not exist in directory "<<fs::current_path().parent_path();
    cerr<<", example cannot run, please check. Did you run the CSMPInterfaces_example first?"<<endl;
    return;
  }
  file_name = file_name_pre + "_variables.dat";
  if(fs::exists(file_name)) fs::copy(file_name, "./");
  else {
    file_name.erase(0,3);
    cerr<<"\n\nError: file '"<<file_name<<"' does not exist in directory "<<fs::current_path().parent_path();
    cerr<<", example cannot run, please check. Did you run the CSMPInterfaces_example first?"<<endl;
    return;
  }
  file_name = file_name_pre + "_splitboundaries.dat";
  if(fs::exists(file_name)) fs::copy(file_name, "./");

  //copy variable file
  file_name = "../variable_and_config_files/";
  if(!fs::is_directory(file_name)) fs::create_directory(file_name);
  file_name += variable_file;
  if(fs::exists(file_name)) fs::copy(file_name, "./");
  else {
    string path = fs::current_path().parent_path();
    path += "/variable_and_config_files/";
    cerr<<"\n\nError: variable file '"<<variable_file<<"' does not exist in directory "<<path;
    cerr<<", example cannot run, please copy this variable file into this directory"<<endl;
    return;
  }

  //copy configuration file if required
  if(!config_file.empty()) {
    file_name = "../variable_and_config_files/";
    if(!fs::is_directory(file_name)) fs::create_directory(file_name);
    file_name += config_file;
    if(fs::exists(file_name)) fs::copy(file_name, "./");
    else {
      string path = fs::current_path().parent_path();
      path += "/variable_and_config_files/";
      cerr<<"\n\nError: configuration file '"<<config_file<<"' does not exist in directory "<<path;
      cerr<<", example cannot run, please copy this configuration file into this directory"<<endl;
      return;
    }
  }
}


} // csmp
