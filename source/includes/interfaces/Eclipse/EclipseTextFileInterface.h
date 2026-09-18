// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef ECLIPSE_TEXT_FILE_INTERFACE_H
#define ECLIPSE_TEXT_FILE_INTERFACE_H

#include "Exception.h"
#include "ErrorHandler.h"

#include <fstream>
#include <iostream>
#include <string>
#include <set>
#include <map>

#include "TextFileIO.h"

namespace csmp {

struct GlobalFunction;

/**

@brief parses the ECLIPSE text file(s) into streams

@author R. Manasipov
@date 2014

Text file parser that separates out comments and finds data blocks on the 
basis of keywords.

*/
template<class STREAM>
class EclipseTextFileInterface {
  public:

    EclipseTextFileInterface( size_t max_line_length = 256 );
    ~EclipseTextFileInterface();

    bool ReadFile( std::ifstream& ifs,
                   STREAM*,
                   bool (STREAM::*IsCommentLineFunction)(char*),
                   bool (STREAM::*ReadFunction)( std::ifstream&, char*, size_t,
                   std::set<std::string>&, std::string&, std::vector<std::string>& ),
                   bool fail_for_unknown_keyword = true );

    bool ReadFile( std::ifstream& ifs,
                   STREAM*,
                   bool (*IsCommentLineFunction)(char*),
                   bool (STREAM::*ReadFunction)( std::ifstream&, char*, size_t,
                   std::set<std::string>&, std::string&, std::vector<std::string>& ),
                   bool fail_for_unknown_keyword = true );

    bool ReadFile( std::ifstream& ifs,
                   STREAM*,
                   bool (STREAM::*IsCommentLineFunction)(char*),
                   bool (STREAM::*ReadFunction)( std::ifstream&, char*, size_t ) );

    bool ReadFile( std::ifstream& ifs,
                   STREAM*,
                   bool (*IsCommentLineFunction)(char*),
                   bool (STREAM::*ReadFunction)( std::ifstream&, char*, size_t ) );

    bool FindKeyword   ( const std::string& );
    void AddKeyword    ( const std::string& );
    void AddKeywords   ( const std::set<std::string>& );
    void ClearKeywords();

  protected:

    std::set<std::string>       keywords_;
    std::string                 keyword_;
    std::vector<std::string>    keyword_parameters_;

    size_t                      line_length_;
    char*                       text_line_;

};


template<>
class EclipseTextFileInterface<GlobalFunction> {

  public:

    EclipseTextFileInterface( size_t max_line_length = 256 );
    ~EclipseTextFileInterface();

    bool ReadFile( std::ifstream& ifs,
                   bool (*IsCommentLineFunction)(char*),
                   bool (*ReadFunction)( std::ifstream&, char*, size_t,
                   std::set<std::string>&, std::string&, std::vector<std::string>& ),
                   bool fail_for_unknown_keyword = true );

    bool ReadFile( std::ifstream& ifs,
                   bool (*IsCommentLineFunction)(char*),
                   bool (*ReadFunction)( std::ifstream&, char*, size_t ) );

    bool FindKeyword   ( const std::string& );
    void AddKeyword    ( const std::string& );
    void AddKeywords   ( const std::set<std::string>& );
    void ClearKeywords ( );

  protected:

    std::set<std::string>       keywords_;
    std::string                 keyword_;
    std::vector<std::string>    keyword_parameters_;

    size_t                      line_length_;
    char*                       text_line_;

};









// TEXT FILE INTERFACE IN CLASS


template<class STREAM>
EclipseTextFileInterface<STREAM>
::EclipseTextFileInterface( size_t line_length )
  : line_length_( line_length )
{
    text_line_ = new char[ line_length ];
}

template<class STREAM>
EclipseTextFileInterface<STREAM>
::~EclipseTextFileInterface()
{
    delete text_line_;
}

template<class STREAM>
bool EclipseTextFileInterface<STREAM>
::FindKeyword( const std::string& keyword )
{
    return ( keywords_.find( keyword ) != keywords_.end() );
}

template<class STREAM>
void EclipseTextFileInterface<STREAM>
::AddKeyword( const std::string& keyword )
{
    keywords_.insert( keyword );
}

template<class STREAM>
void EclipseTextFileInterface<STREAM>
::AddKeywords( const std::set<std::string>& keywords )
{
    std::set<std::string>::const_iterator kitEnd = keywords.end();
    for( std::set<std::string>::const_iterator
         kit = keywords.begin(); kit != kitEnd; kit++ )
        keywords_.insert( *kit );
}

template<class STREAM>
void EclipseTextFileInterface<STREAM>
::ClearKeywords( )
{
    keywords_.clear();
}


/** Reading block data marked by keywords

    @author Roman, 2014

*/
template<class STREAM>
bool EclipseTextFileInterface<STREAM>::ReadFile( std::ifstream& ifs,
                                          STREAM* stream,
                                          bool (STREAM::*IsCommentLineFunction)(char*),
                                          bool (STREAM::*ReadFunction)( std::ifstream&, char*, size_t, std::set<std::string>&, std::string&, std::vector<std::string>& ),
                                          bool fail_for_unknown_keyword )
{
    csmp::ErrorHandler& error_handler( csmp::ErrorHandler::Instance() );

    // Read block data marked by keywords
    do{
        // read non-empty text line
        readNonBlankLine( ifs, text_line_, line_length_ );
        if ( ifs.eof() )
        {
            if( error_handler.Verbose() )
                std::cout <<"\nEclipseTextFileInterfaceInClass<dim>::ReadFile: Reading completed!" << std::endl;
            ifs.close();
            return true;
        }
        // analyze text line
        if ( !(stream->*IsCommentLineFunction)(text_line_) )
        {
            // read keyword and it's parameters
            readKeyword( ifs, text_line_, line_length_,
                         keyword_, keyword_parameters_ );
            // read first line in block if not blank
            csmp::readNonBlankLine( ifs, text_line_, line_length_ );
            if ( ifs.eof() )
            {
                if( error_handler.Verbose() )
                    std::cout <<"\nEclipseTextFileInterfaceInClass<dim>::ReadFile: Reading completed!" << std::endl;
                ifs.close();
                return true;
            }
            // check the keyword
            if( fail_for_unknown_keyword && keywords_.find( keyword_) == keywords_.end() )
            {
                std::cerr<< "EclipseTextFileInterfaceInClass<dim>::ReadFile: Can not read data in block marked by keyword: "<< keyword_ << std::endl;
                std::cerr<< "Please correct the keyword that you've used in your file" << std::endl;
                std::cerr<< "The list of the avaliable keywords: "<< std::endl;
                size_t i(1);
                for( std::set<std::string>::const_iterator it = keywords_.begin(); it != keywords_.end(); it++, i++)
                    std::cerr<<"Keyword [ "<<i<< " ] = "<<(*it)<<std::endl;
                std::cerr <<"\nEclipseTextFileInterfaceInClass<dim>::ReadFile: Reading was done with errors!" << std::endl;
                ifs.close();
                error_handler.Note( csmp::FATAL_ERROR,
                                   "EclipseTextFileInterfaceInClass<dim>::ReadFile:",
                                   "Undefined keyword:",
                                   keyword_ );
                return false;
            }
            // read data
            else if( !(stream->*ReadFunction)( ifs, text_line_, line_length_,
                                               keywords_, keyword_, keyword_parameters_ ) )
            {
                std::cerr <<"\nEclipseTextFileInterfaceInClass<dim>::ReadFile: Reading was done with errors!" << std::endl;
                ifs.close();
                error_handler.Note( csmp::FATAL_ERROR,
                                   "EclipseTextFileInterfaceInClass<dim>::ReadFile:",
                                   "Can not read data in block marked by keyword: ",
                                   keyword_ );
                return false;
            }
        }
    }
    while ( !ifs.eof() );

    if( error_handler.Verbose() )
        std::cout <<"\nEclipseTextFileInterfaceInClass<dim>::ReadFile: Reading completed!" << std::endl;

    ifs.close();

    return true;

} // end ReadInputFile






template<class STREAM>
bool EclipseTextFileInterface<STREAM>::ReadFile( std::ifstream& ifs,
                                          STREAM* stream,
                                          bool (*IsCommentLineFunction)(char*),
                                          bool (STREAM::*ReadFunction)( std::ifstream&, char*, size_t, std::set<std::string>&, std::string&, std::vector<std::string>& ),
                                          bool fail_for_unknown_keyword )
{
    csmp::ErrorHandler& error_handler( csmp::ErrorHandler::Instance() );

    // Read block data marked by keywords
    do{
        // read non-empty text line
        csmp::readNonBlankLine( ifs, text_line_, line_length_ );
        if ( ifs.eof() )
        {
            if( error_handler.Verbose() )
                std::cout <<"\nEclipseTextFileInterfaceInClass<dim>::ReadFile: Reading completed!" << std::endl;
            ifs.close();
            return true;
        }
        // analyze text line
        if ( !(*IsCommentLineFunction)(text_line_) )
        {
            // read keyword and it's parameters
            csmp::readKeyword( ifs, text_line_, line_length_,
                               keyword_, keyword_parameters_ );
            // read first line in block if not blank
            csmp::readLineTellIfBlank( ifs, text_line_, line_length_ );
            if ( ifs.eof() )
            {
                if( error_handler.Verbose() )
                    std::cout <<"\nEclipseTextFileInterfaceInClass<dim>::ReadFile: Reading completed!" << std::endl;
                ifs.close();
                return true;
            }
            // check the keyword
            if( fail_for_unknown_keyword && keywords_.find( keyword_) == keywords_.end() )
            {
                std::cerr<< "EclipseTextFileInterfaceInClass<dim>::ReadFile: Can not read data in block marked by keyword: "<< keyword_  << std::endl;
                std::cerr<< "Please correct the keyword that you've used in your file" << std::endl;
                std::cerr<< "The list of the avaliable keywords: "<< std::endl;
                size_t i(1);
                for( std::set<std::string>::const_iterator it = keywords_.begin(); it != keywords_.end(); it++, i++)
                    std::cerr<<"Keyword [ "<<i<< " ] = "<<(*it)<<std::endl;
                std::cerr <<"\nEclipseTextFileInterfaceInClass<dim>::ReadFile: Reading was done with errors!" << std::endl;
                ifs.close();
                error_handler.Note( csmp::FATAL_ERROR,
                                   "EclipseTextFileInterfaceInClass<dim>::ReadFile:",
                                   "Undefined keyword:",
                                   keyword_.c_str() );
                return false;
            }
            // read data
            else if( !(stream->*ReadFunction)( ifs, text_line_, line_length_,
                                               keywords_, keyword_, keyword_parameters_ ) )
            {
                std::cerr <<"\nEclipseTextFileInterfaceInClass<dim>::ReadFile: Reading was done with errors!" << std::endl;
                ifs.close();
                error_handler.Note( csmp::FATAL_ERROR,
                                   "EclipseTextFileInterfaceInClass<dim>::ReadFile:",
                                   "Can not read data in block marked by keyword: ",
                                   keyword_.c_str() );
                return false;
            }
        }
    }
    while ( !ifs.eof() );

    if( error_handler.Verbose() )
        std::cout <<"\nEclipseTextFileInterfaceInClass<dim>::ReadFile: Reading completed!" << std::endl;

    ifs.close();

    return true;

} // end ReadInputFile


template<class STREAM>
bool EclipseTextFileInterface<STREAM>
::ReadFile( std::ifstream& ifs,
            STREAM* stream,
            bool (STREAM::*IsCommentLineFunction)(char*),
            bool (STREAM::*ReadFunction)( std::ifstream&, char*, size_t ) )
{
    csmp::ErrorHandler& error_handler( csmp::ErrorHandler::Instance() );

    // Read block data marked by keywords
    do{
        // read non-empty text line
        readNonBlankLine( ifs, text_line_, line_length_ );
        if ( ifs.eof() )
        {
            if( error_handler.Verbose() )
                std::cout <<"\nEclipseTextFileInterfaceInClass<dim>::ReadFile: Reading completed!" << std::endl;
            ifs.close();
            return true;
        }
        // analyze text line
        if ( !(stream->*IsCommentLineFunction)(text_line_) )
        {
            if( !(stream->*ReadFunction)( ifs, text_line_, line_length_ ) )
            {
                std::cerr <<"\nEclipseTextFileInterfaceInClass<dim>::ReadFile: Reading was done with errors!" << std::endl;
                ifs.close();
                error_handler.Note( csmp::FATAL_ERROR,
                                   "EclipseTextFileInterfaceInClass<dim>::ReadFile:",
                                   "Can not read data in block marked by keyword: ",
                                   keyword_ );
                return false;
            }
        }
    }
    while ( !ifs.eof() );

    if( error_handler.Verbose() )
        std::cout <<"\nEclipseTextFileInterfaceInClass<dim>::ReadFile: Reading completed!" << std::endl;

    ifs.close();

    return true;

} // end ReadInputFile




template<class STREAM>
bool EclipseTextFileInterface<STREAM>
::ReadFile( std::ifstream& ifs,
            STREAM* stream,
            bool (*IsCommentLineFunction)(char*),
            bool (STREAM::*ReadFunction)( std::ifstream&, char*, size_t ) )
{
    csmp::ErrorHandler& error_handler( csmp::ErrorHandler::Instance() );

    // Read block data marked by keywords
    do{
        // read non-empty text line
        csmp::readNonBlankLine( ifs, text_line_, line_length_ );
        if ( ifs.eof() )
        {
            if( error_handler.Verbose() )
                std::cout <<"\nEclipseTextFileInterfaceInClass<dim>::ReadFile: Reading completed!" << std::endl;
            ifs.close();
            return true;
        }
        // analyze text line
        if ( !(*IsCommentLineFunction)(text_line_) )
        {
            if( !(stream->*ReadFunction)( ifs, text_line_, line_length_ ) )
            {
                std::cerr <<"\nEclipseTextFileInterfaceInClass<dim>::ReadFile: Reading was done with errors!" << std::endl;
                ifs.close();
                error_handler.Note( csmp::FATAL_ERROR,
                                   "EclipseTextFileInterfaceInClass<dim>::ReadFile:",
                                   "Cannot read data in block marked by keyword: ",
                                   keyword_.c_str() );
                return false;
            }
        }
    }
    while ( !ifs.eof() );

    if( error_handler.Verbose() )
        std::cout <<"\nEclipseTextFileInterfaceInClass<dim>::ReadFile: Reading completed!" << std::endl;

    ifs.close();

    return true;

} // end ReadInputFile



} // end csmp


#endif
