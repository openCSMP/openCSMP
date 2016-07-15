#include <iostream>
#include <cstdlib>
#include <typeinfo>

#include "Example.h"

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
  text += typeid(*this).name();
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



} // csmp
