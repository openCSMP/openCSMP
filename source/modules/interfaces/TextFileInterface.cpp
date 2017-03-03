#include "TextFileInterface.h"

using namespace std;

namespace csmp {

// TEXT FILE INTERFACE IN GLOBAL SPACE

TextFileInterface<GlobalFunction>
::TextFileInterface( size_t line_length )
  : line_length_( line_length )
{
    text_line_ = new char[ line_length ];
}


TextFileInterface<GlobalFunction>
::~TextFileInterface()
{
    delete text_line_;
}

bool TextFileInterface<GlobalFunction>
::FindKeyword( const std::string& keyword )
{
    return ( keywords_.find( keyword ) != keywords_.end() );
}


void TextFileInterface<GlobalFunction>
::AddKeyword( const std::string& keyword )
{
    keywords_.insert( keyword );
}

void TextFileInterface<GlobalFunction>
::AddKeywords( const std::set<std::string>& keywords )
{
    std::set<std::string>::const_iterator kitEnd = keywords.end();
    for( std::set<std::string>::const_iterator
         kit = keywords.begin(); kit != kitEnd; kit++ )
        keywords_.insert( *kit );
}

void TextFileInterface<GlobalFunction>
::ClearKeywords( )
{
    keywords_.clear();
}


/** Reading block data marked by keywords

    @author Roman, 2014

*/
bool TextFileInterface<GlobalFunction>
::ReadFile( std::ifstream& ifs,
            bool (*IsCommentLineFunction)(char*),
            bool (*ReadFunction)( std::ifstream&, char*, size_t, std::set<std::string>&, std::string&, std::vector<std::string>& ),
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
                std::cout <<"\nTextFileInterface<dim>::ReadFile: Reading completed!" << std::endl;
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
            csmp::readNonBlankLine( ifs, text_line_, line_length_ );
            if ( ifs.eof() )
            {
                if( error_handler.Verbose() )
                    std::cout <<"\nTextFileInterface<dim>::ReadFile: Reading completed!" << std::endl;
                ifs.close();
                return true;
            }
            // check the keyword
            if( fail_for_unknown_keyword && keywords_.find( keyword_) == keywords_.end() )
            {
                std::cerr<< "TextFileInterface<dim>::ReadFile: Can not read data in block marked by keyword: "<< keyword_ << std::endl;
                std::cerr<< "Please correct the keyword that you've used in your file" << std::endl;
                std::cerr<< "The list of the avaliable keywords: "<< std::endl;
                size_t i(1);
                for( std::set<std::string>::const_iterator it = keywords_.begin(); it != keywords_.end(); it++, i++)
                    std::cerr<<"Keyword [ "<<i<< " ] = "<<(*it)<<std::endl;
                std::cerr <<"\nTextFileInterface<dim>::ReadFile: Reading was done with errors!" << std::endl;
                ifs.close();
                error_handler.notice( CSMP_FATAL_ERROR,
                                   "TextFileInterface<dim>::ReadFile:",
                                   "Undefined keyword:",
                                   keyword_.c_str() );
                return false;
            }
            // read data
            else if( !ReadFunction( ifs, text_line_, line_length_,
                                    keywords_, keyword_, keyword_parameters_ ) )
            {
                std::cerr <<"\nTextFileInterface<dim>::ReadFile: Reading was done with errors!" << std::endl;
                ifs.close();
                error_handler.notice( CSMP_FATAL_ERROR,
                                   "TextFileInterface<dim>::ReadFile:",
                                   "Can not read data in block marked by keyword: ",
                                   keyword_.c_str() );
                return false;
            }
        }
    }
    while ( !ifs.eof() );

    if( error_handler.Verbose() )
        std::cout <<"\nTextFileInterface<dim>::ReadFile: Reading completed!" << std::endl;

    ifs.close();

    return true;

} // end ReadInputFile


bool TextFileInterface<GlobalFunction>
::ReadFile( std::ifstream& ifs,
            bool (*IsCommentLineFunction)(char*),
            bool (*ReadFunction)( std::ifstream&, char*, size_t ) )
{
    csmp::ErrorHandler& error_handler( csmp::ErrorHandler::Instance() );

    // Read block data marked by keywords
    do{
        // read non-empty text line
        csmp::readNonBlankLine( ifs, text_line_, line_length_ );
        if ( ifs.eof() )
        {
            if( error_handler.Verbose() )
                std::cout <<"\nTextFileInterface<dim>::ReadFile: Reading completed!" << std::endl;
            ifs.close();
            return true;
        }
        // analyze text line
        if ( !(*IsCommentLineFunction)(text_line_) )
        {
            if( !ReadFunction( ifs, text_line_, line_length_ ) )
            {
                std::cerr <<"\nTextFileInterface<dim>::ReadFile: Reading was done with errors!" << std::endl;
                ifs.close();
                error_handler.notice( CSMP_FATAL_ERROR,
                                   "TextFileInterface<dim>::ReadFile:",
                                   "Can not read data in block marked by keyword: ",
                                   keyword_.c_str() );
                return false;
            }
        }
    }
    while ( !ifs.eof() );

    if( error_handler.Verbose() )
        std::cout <<"\nTextFileInterface<dim>::ReadFile: Reading completed!" << std::endl;

    ifs.close();

    return true;

} // end ReadInputFile



} // end csmp
