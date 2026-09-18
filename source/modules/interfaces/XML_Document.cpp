// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "XML_Document.h"

#include <fstream>
#include <cassert>

using namespace std;

namespace csmp {


/// default constructor
XML_Document::XML_Document()
 : currentLevel_( 0 )
{
}


/// opens a new xml node with name and label/settings
void XML_Document::OpenNode( const char* nodeLabel )
{
  // extract tag name (everything before first space or end)
  std::string label( nodeLabel );
  nodeStack_.push_back( label.substr( 0, label.find(' ') ) );

  // setting current line indentation
  BringToLevel();

  // establishing content
  string tag;
  tag.reserve( std::strlen(nodeLabel) + 3 );
  tag += "<";
  tag += nodeLabel;
  tag += ">\n";

  // writing content
  WriteToData( tag.c_str() );

  // updating indentation level
  Up();
}


/// closes given node
void XML_Document::CloseNode( const char* nodeName )
{
  assert( !nodeStack_.empty() && nodeStack_.back() == nodeName );
  nodeStack_.pop_back();

  // updating indentation level
  Down();

  // setting current line indentation
  BringToLevel();

  // establishing content
  string tag;
  tag.reserve( std::strlen(nodeName) + 3 );
  tag += "</";
  tag += nodeName;
  tag += ">\n";

  // writing content
  WriteToData( tag.c_str() );
}


/// inserts a xml comment
void XML_Document::AddComment( const char* comment )
{
  BringToLevel();
  string tag;
  tag.reserve( std::strlen(comment) + 3 );
  tag += "<!--\n";
  tag += comment;
  tag += "\n-->\n\n";
  WriteToData( tag.c_str() );
}


/// inserts a xml system information
void XML_Document::AddInfo( const char* info )
{
  BringToLevel();
  string tag;
  tag.reserve( std::strlen(info) + 3 );
  tag += "<?";
  tag += info;
  tag += "?>";
  WriteToData( tag.c_str() );
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
  std::ofstream xmlFile;
  xmlFile.open ( fileName, std::fstream::out );

  // return false if failed
  if( !xmlFile.is_open() )
    return false;

  // writing data to file
  xmlFile << data_;

  return xmlFile.good();
}


} // csmp
