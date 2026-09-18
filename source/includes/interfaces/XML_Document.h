// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef XML_DOCUMENT_H
#define XML_DOCUMENT_H

#include <cstring>
#include <string>
#include <vector>

namespace csmp {

/**
@brief trivial class to provide an interface for creating and parsing XML documents.

@author P. Lang
@date Jan 2011

@todo (3) Range check current level up/down
@todo (3) Integrate reading/parsing features
@todo (3) Implement vector to hold node names per level to allow for auto closing etc...

*/
class XML_Document
{
public:
  XML_Document();

  XML_Document& operator=( const XML_Document& src ){ data_ = src.data_; currentLevel_ = src.currentLevel_; return *this; }

  void OpenNode( const char* nodeLabel );
  void CloseNode( const char* nodeName );

  void AddComment( const char* comment );
  void AddInfo( const char* comment );

  void InsertData( const char* data );
  void LineBreak( size_t n = 1 ) { for( size_t i = 0; i < n; ++i ) WriteToData( "\n" ); }
  void Tab( size_t n = 1 ) { for( size_t i = 0; i < n; ++i ) WriteToData( "\t" ); }
  void BringToLevel() { for( size_t i = 0; i < Level(); ++i ) WriteToData( "\t"); }

  /// returns the number of characters in current line
  size_t CurrentLineCharacterCount() const { return data_.size() - data_.rfind( "\n" ); }

  bool WriteToFile( const char* fileName ) const;

protected:
  // disabled
  XML_Document ( const XML_Document& src );

  // indentation level
  size_t Up( size_t n = 1 ) { currentLevel_ += n; return currentLevel_; }
  size_t Down( size_t n = 1 ) { currentLevel_ -= n; return currentLevel_; }
  void   Level( size_t level ) { currentLevel_ = level; }
  size_t Level() const { return currentLevel_; }

  // append to data
  void WriteToData( const char* content ) { data_ += content; }

private:
  std::vector<std::string> nodeStack_; ///< to make sure that open nodes are closed again
  std::string data_;                   ///< stl string to hold document data
  size_t currentLevel_;                ///< current level of indentation/node
};

} // csmp

#endif // XML_DOCUMENT_H
