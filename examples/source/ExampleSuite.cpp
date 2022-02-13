#include "Exception.h"
#include "ExampleSuite.h"


using namespace std;

namespace csmp{

// INTERFACE

/// public constructor
ExampleSuite::ExampleSuite( const string &title, ostream *ostream )
  : title_( title ),
    ostream_( ostream )
{
}

/// start the stdIO interface
void ExampleSuite::Run()
{
  Initialize();
  WriteExamplesFile();
  SuiteHeader();
  UIMenu();
}

/// clear suite of all examples, calls their destructors
void ExampleSuite::Free()
{
  // loop over all categories(map entries) and call the destructors of thei examples(vector entries)
  for( map<string,vector<Example*>* >::iterator it = categories_.begin();
       it != categories_.end(); ++it )
  {
    for( vector<Example*>::iterator iit = (*it).second->begin();
         iit != (*it).second->end(); ++iit )
    {
      delete (*iit);
    } // loop vector
  } // loop map
  categories_.clear();
}

/// add example to suite (use heap)
size_t ExampleSuite::RegisterExample( Example *example )
{
  // initialize example - sets database entries( Author, Category, Desrciption...)
  example->Initialize();
  // get pointer to example category entry in map. If not contained yet, create one.
  // either way, an iterator to the respective map entry is established
  pair<map<string,vector<Example*>* >::iterator,bool> inserter;
  inserter.first = categories_.find( example->GetCategory() );
  if( inserter.first == categories_.end() )
  {
    vector<Example*>* newExampleContainer = new vector<Example*>();
    inserter = categories_.insert( make_pair( example->GetCategory(), newExampleContainer ) );
  }
  // take the established iterator to the category(map), then push back example
  // in category example container(vector)
  inserter.first->second->push_back( example );
  size_t check( inserter.first->second->size() );
  // IO
  /*
  *ostream_ << "\nExampleSuite: '" << example->GetTitle();
  *ostream_ << "' added to '" << example->GetCategory() << "' No.";
  *ostream_ << check;
  */
  return check;
}

// INTERNALS

/// initialize suite db
void ExampleSuite::Initialize()
{

}

/// ostream all categories
void ExampleSuite::CategoriesList() const
{
  size_t counter( 1 );
  for( map<string,vector<Example*>* >::const_iterator it = categories_.begin();
       it != categories_.end(); ++it )
  {
    *ostream_ << "\n  " << counter++ << ".) " << (*it).first;
    *ostream_ << " (" << (*it).second->size() << ")";
  }
}

/// ostream all examples of a category
void ExampleSuite::ExamplesList( map<string,vector<Example*>* >::const_iterator map ) const
{
  size_t counter( 1 );
  for( vector<Example*>::const_iterator it = map->second->begin();
       it != map->second->end(); ++it )
  {
    *ostream_ << "\n  " << counter++ << ".) " << StripExampleTitle( (*it)->GetTitle() );
    *ostream_ << " (" << (*it)->GetDifficulty() << ")";
  }
}

/// starts UI menu
void ExampleSuite::UIMenu() const
{
  UICategoryMenu();
}

/// UI menu for Level Example Details
void ExampleSuite::UIExampleDetailsMenu( size_t exampleIndex,
                                         map<string,vector<Example*>* >::const_iterator category ) const
{
  Example* example = category->second->at( exampleIndex );
  ExampleDetails( example );

  if( ChooseExampleAction() == true )
  {
    OstreamDoubleUnderlined( "Example '" + example->GetTitle() + "' started..." );
    example->Run();
    OstreamDoubleUnderlined( "Example '" + example->GetTitle() + "' done..." );
    UIExamplesMenu( category );
  }
  else
    UIExamplesMenu( category );
}

/// UI menu for Level Examples
void ExampleSuite::UIExamplesMenu( map<string,vector<Example*>* >::const_iterator category ) const
{
  OstreamUnderlined( "Category '" + category->first + "'" );
  ExamplesList( category );
  size_t exampleIndex( ChooseExample( (category->second ) ) );
  if( exampleIndex == -1 )
    UICategoryMenu();
  else
  {
    UIExampleDetailsMenu( exampleIndex-1, category );
  }
}

/// UI menu for Level Category
void ExampleSuite::UICategoryMenu() const
{
  OstreamUnderlined( "Suite Categories:" );
  map<string,vector<Example*>* >::const_iterator category;
  CategoriesList();
  category = ChooseCategory();
  if( category == categories_.end() )
      return;
  UIExamplesMenu( category );
}


/// sorts examples in vector<Example*> depending on difficulty
void ExampleSuite::SortExampleVector( vector<Example*>& examples )
{
  // finding all occuring diffiulties and establish an example vector for each
  vector<uint32_t> difficulties;
  vector<vector<Example*>* > sortedExamples;
  bool newDifficulty;
  for( vector<Example*>::const_iterator it = examples.begin(); it != examples.end(); ++it )
  {
    newDifficulty = true;
    for( vector<uint32_t>::iterator iit = difficulties.begin(); iit != difficulties.end(); ++iit )
      if( (*it)->GetDifficulty() == (*iit) )
        newDifficulty = false;
    if( newDifficulty )
    {
      difficulties.push_back( (*it)->GetDifficulty() );
      sortedExamples.push_back( new vector<Example*>() );
    }
  }
  // putting examples in their respectiv difficulty vector
  size_t index( 0 );
  for( vector<Example*>::const_iterator it = examples.begin(); it != examples.end(); ++it )
  {
    index = 0;
    // STOPPED HERE
  }

}

/// ostreams example in detail
void ExampleSuite::ExampleDetails( Example* example ) const
{
  OstreamUnderlined( "Example '" + StripExampleTitle( example->GetTitle() ) + "'" );
  *ostream_ << "\n  Authors: ";
  for( list<string>::const_iterator it = example->GetAuthorsBegin();
       it != example->GetAuthorsEnd(); ++it )
    *ostream_ << (*it) << "  ";
  *ostream_ << "\n  Difficulty: " << example->GetDifficulty();
  *ostream_ << "\n  Description:";
  for( list<string>::const_iterator it = example->GetDescriptionsBegin();
       it != example->GetDescriptionsEnd(); ++it )
    *ostream_ << "\n  - " << (*it);
  *ostream_ << "\n  Requirements:";
  for( list<string>::const_iterator it = example->GetRequirementsBegin();
       it != example->GetRequirementsEnd(); ++it )
    *ostream_ << "\n  - " << (*it);
}

void ExampleSuite::OstreamUnderlined( string text ) const
{
  *ostream_ << endl << text << endl;
  for( auto i = 1; i <= text.length(); ++i )
    *ostream_ << "-";
  *ostream_ << endl;
}

void ExampleSuite::OstreamDoubleUnderlined( string text ) const
{
  *ostream_ << endl;
  for( auto i = 1; i <= text.length(); ++i )
    *ostream_ << "=";
  *ostream_ << endl << text << endl;
  for( auto i = 1; i <= text.length(); ++i )
    *ostream_ << "=";
  *ostream_ << endl;
}

/// ostream suite header
void ExampleSuite::SuiteHeader() const
{
  string header = "\nCSMP Example Suite '" + title_ + "'\n";
  OstreamDoubleUnderlined( header );
}

/// UI to choose example
size_t ExampleSuite::ChooseExample( const std::vector<Example*>* examples ) const
{
  *ostream_ << "\n\nPlease enter example number(-1 to get back): ";
  size_t examplesCount( examples->size() );
  size_t exampleIndex( ChoiceWithinRange( -1, examplesCount ) );
  return exampleIndex;
}

/// UI to choose what to do with given example
bool ExampleSuite::ChooseExampleAction() const
{
  *ostream_ << "\n\nHow do you want to proceed(-1 to get back, 1 to run example): ";
  size_t choice( ChoiceWithinRange( -1, 1 ) );
  if( choice == -1 )
    return false;
  else
    return true;
}

/// UI to choose category
map<string,vector<Example*>* >::const_iterator ExampleSuite::ChooseCategory() const
{
  *ostream_ << "\n\nPlease enter category number or q to exit: ";
  int categoryIndex( ChoiceWithinRange( 1, categories_.size() ) );

  // quit
  if( categoryIndex == -1 )
      return categories_.end();

  size_t counter( 1 );
  for( map<string,vector<Example*>* >::const_iterator it = categories_.begin();
       it != categories_.end(); ++it )
  {
    if( categoryIndex == counter++ )
      return it;
  }
  // should never get to here: maybe implement ERROR
  return categories_.end();
}

/// assure that user choice is within range, otherwise reprompt
int ExampleSuite::ChoiceWithinRange( double min, double max ) const
{
  string input;
  int choice;
  bool valid( true );
  do{
      cin >> input;
      if( input == std::string("q") || input == std::string("Q") )
      {
          choice = -1;
          return choice;
      }
      else
      {
          choice = atoi(input.c_str());
          if ( (!choice) || (atoi(input.c_str()) - choice) || choice < min || choice > max)
          {
            *ostream_ << "\n\nInvalid. Please reenter: ";
            valid = false;
          }
          else
            valid = true;
      }
  }while( !valid );
  return choice;
}

/// remove typeid from example name
string ExampleSuite::StripExampleTitle( string title ) const
{
  size_t end = title.find_first_of( "(" );
  title.erase( end - 1 );
  return title;
}

/// writes all registered examples to file, sorted by section headings
void ExampleSuite::WriteExamplesFile()
{
  // establishing output file
  ofstream toFile;
  toFile.open( "ExampleSuiteList.txt", ios::out );
  // backing up and redirecting ostream
  ostream* ostreamCache( ostream_ );
  ostream_ = &toFile;

  // ostream SuiteHeader
  SuiteHeader();

  // loop over all categories & ostream
  for( map<string,vector<Example*>* >::const_iterator it = categories_.begin();
       it != categories_.end(); ++it )
  {
    // ostream category title
    OstreamDoubleUnderlined( (*it).first );
    // ostream all examples in category
    ExamplesList( it );
    *ostream_ << endl;
  } // categories

  // reset ostream and close file
  ostream_ = ostreamCache;
  toFile.close();
}



} // csmp
