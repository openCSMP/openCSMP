#include "XML_Document.h"

#include <fstream>

namespace csmp{


/// default constructor
XML_Document::XML_Document()
 : currentLevel_( 0 )
{
}


/// opens a new xml node with name and label/settings
void XML_Document::OpenNode( const char* nodeLabel )
{
  // setting current line indentation
  BringToLevel();

  // establishing content
  sCache_.clear(); sCache_ += "<"; sCache_ += nodeLabel; sCache_ += ">\n";

  // writing content
  WriteToData( sCache_.c_str() ); sCache_.clear();

  // updating indentation level
  Up();
}


/// closes given node
void XML_Document::CloseNode( const char* nodeName )
{
  // updating indentation level
  Down();

  // setting current line indentation
  BringToLevel();

  // establishing content
  sCache_.clear(); sCache_ += "</"; sCache_ += nodeName; sCache_ += ">\n";

  // writing content
  WriteToData( sCache_.c_str() ); sCache_.clear();
}


/// inserts a xml comment
void XML_Document::AddComment( const char* comment )
{
  BringToLevel();
  sCache_.clear(); sCache_ += "<!--\n"; sCache_ += comment; sCache_ += "\n-->\n\n";
  WriteToData( sCache_.c_str() ); sCache_.clear();
}


/// inserts a xml system information
void XML_Document::AddInfo( const char* info )
{
  BringToLevel();
  sCache_.clear(); sCache_ += "<\?"; sCache_ += info; sCache_ += "\?>\n\n";
  WriteToData( sCache_.c_str() ); sCache_.clear();
}


/// inserts data at current position(not responsible for line breaks/levels!!)
void XML_Document::InsertData( const char* data )
{
  WriteToData( data );
}


/// writes data to file(no exception handling)
bool XML_Document::WriteToFile( const char* fileName ) const
{
  // open file to output to
  std::fstream xmlFile;
  xmlFile.open ( fileName, std::fstream::out );

  // return false if failed
  if( !xmlFile.is_open() )
    return false;

  // writing data to file
  xmlFile << data_;

  // closing xml file
  xmlFile.close();

  return true;
}


} // csmp
