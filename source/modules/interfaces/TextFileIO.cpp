#include "TextFileIO.h"
#include "Boundary.h"

#include <sys/types.h>
#include <sys/stat.h>
#ifndef _WIN32
#include <unistd.h>
#endif

using namespace std;

namespace csmp {

// GLOBAL READING FUNCTIONS

bool isFileExist (const std::string& name)
{
    ifstream f(name.c_str());
    if (f.good()) {
        f.close();
        return true;
    } else {
        f.close();
        return false;
    }
}

bool isFileExistStat (const std::string& name)
{
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

/// checks whether a file with the given name exists in the current directory
bool checkExistance( const char* filename )
{
   ifstream Infield(filename);
   if ( !Infield.good() ) return false;
   return true;
}

/**

Appends file extension ( by default: '.txt') to the supplied file name before it attempts to open this file.
If successful the method returns true and initializes the argument file stream;
else it returns false.

@section arguments Input Arguments

The name of the file that shall be opened.

@return The file stream argument. This will be a reference to the beginning of
the input file if it could be opened successfully.

@section messages Messages

An error is reported if the input file cannot be opened.
*/

template<class FileStream>
bool openFile( FileStream& fs, const std::string& fname )
{
    char    filename[200];

    assert( !fname.empty() );
    strcpy( filename, fname.c_str() );

    fs.open( filename );

    if ( !fs.is_open() )
        throw csmp::Exception( ERROR,
                               "openInputFile",
                               "File seems to be missing",
                               filename );

    return true;

} // end openFile

template bool openFile(std::ifstream&,const std::string&);
template bool openFile(std::ofstream&,const std::string&);

template<class FileStream>
bool openFile( FileStream& fs, const std::string& fname, const std::string& file_extension )
{
    char    filename[200];

    assert( !fname.empty() );
    strcpy( filename, fname.c_str() );
    strcat( filename, file_extension.c_str() );

    fs.open( filename );

    if ( !fs.is_open() )
        throw csmp::Exception( ERROR,
                               "openInputFile",
                               "File seems to be missing",
                               filename );

    return true;

} // end openFile

template bool openFile(std::ifstream&,const std::string&,const std::string&);
template bool openFile(std::ofstream&,const std::string&,const std::string&);

template<class FileStream>
bool openFile( FileStream& fs, const std::string& fname, const std::vector<std::string>& file_extensions )
{
    bool success( false );
    char    filename[200];
    assert( !fname.empty() );

    strcpy( filename, fname.c_str() );
    strcat( filename, file_extensions[0].c_str() );
    fs.open( filename );

    if ( !fs.is_open() )
    {
        const size_t extensions( file_extensions.size() );
        for( size_t i =1; i < extensions; i++ )
        {
            strcpy( filename, fname.c_str() );
            strcat( filename, file_extensions[i].c_str() );
            fs.open( filename );
            if( fs.is_open( ) ){
                success = true;
                break;
            }
        }
    }
    else
        success = true;

    if(!success){
        throw csmp::Exception( ERROR,
                               "openInputFile",
                               "File seems to be missing",
                               filename );
        return false;
    }

    return true;

} // end openFile

template bool openFile(std::ifstream&,const std::string&,const std::vector<std::string>&);
template bool openFile(std::ofstream&,const std::string&,const std::vector<std::string>&);


/** The method reads the ASCII file headline and echoes it to the
screen. The correctly read string is returned.

@section arguments Input Arguments

The method reads from an initialized file stream.

@param header The file title is returned into the second method argument. If there is
a title line,

@return bool the method returns true, else false.

*/
bool readFileHeader( ifstream& ifs, string& header, bool verbose )
{
    char  text_line[256];

    // reading file header
    ifs.getline( text_line, 256 ); // title line

    if( verbose )
    {
        cout <<"\nreadFileHeader: File header: "<< endl;
        cout <<"\n\t"<< text_line << endl;
        header = text_line;
    }

    if ( header.empty() )
        return false;

    return true;

} // readFileHeader

bool readLineTellIfBlank( std::ifstream& ifs, char* text_line, size_t line_length )
{
    ifs.getline( text_line, line_length );
    return isBlankLine(text_line);
}

bool readNonBlankLine( std::ifstream& ifs, char* text_line, size_t line_length )
{
    ifs.getline( text_line, line_length );
    while ( isBlankLine(text_line) && !ifs.eof() )
        ifs.getline( text_line, line_length );
    return true;
}

int readFirstLineInBlock( std::ifstream& ifs, char* text_line, size_t line_length )
{
    // read first line in block
    if ( !readNonBlankLine( ifs, text_line, line_length ) )
        return 0;

    if ( ifs.eof() )
    {
        ifs.close();
        std::cout <<"\nreadFirstLineInBlock: reading completed!" << std::endl;
        return 2;
    }

    return 1;
}

int readIncludeFileName( std::ifstream& ifs, char* text_line, size_t line_length,
                         std::string& fname )
{
    // read keyword and it's parameters
    const char* filename_delims       = ":,\t,\n,\r'`/";

    int firstLine( readFirstLineInBlock( ifs, text_line, line_length ) );
    if( firstLine == 0 || firstLine == 2 )
        return firstLine;

    // read file name
    fname = strtok( text_line, filename_delims);

    return 1;
}

bool readKeyword( std::ifstream& ifs, char* text_line, size_t line_length,
                  std::string& keyword, std::vector<std::string>&  keyword_parameters )
{
    // read keyword and it's parameters
    char*       token( 0 );
    const char* keyword_delims       = ":,\t,\n,\r";
    const char* keyword_param_delims = " :,\t,\n,\r";

    // read keyword and make it upper case
    // in order to avoid case sensitive mistakes
    keyword = strtok( text_line, keyword_delims);
    std::transform( keyword.begin(), keyword.end(), keyword.begin(), ::toupper );

    // read keyword parameters
    keyword_parameters.clear();
    do
    {
        token = strtok( NULL, keyword_param_delims );
        if( token != NULL )
            keyword_parameters.push_back( token );
    }
    while( token != NULL );

    return true;
}

bool readKeywordsAndParameters( std::ifstream& ifs, char* text_line, size_t line_length, bool verbose,
                                std::vector<std::pair<std::string, std::vector<std::string> > >& keywords_and_parameters )
{
    if( verbose )
        cout <<"\n\nreadKeywordsAndParameters: reading keywords and it's parameters...\n\n";

    do {
        if ( !isCommentLine( text_line ) )
        {
            std::string  keyword("undefined");
            std::vector<std::string>  parameters;
            csmp::readKeyword( ifs, text_line, line_length,
                               keyword, parameters );
            keywords_and_parameters.push_back( std::make_pair(keyword,parameters) );
        }
        ifs.getline( text_line, line_length );
    }
    while ( !isBlankLine(text_line) && !ifs.eof() );

    if( verbose )
        cout <<"\n\nreadKeywordsAndParameters: Finished reading keywords and it's parameters!\n\n";

    return true;
}

/**
@section arguments Input Arguments

The method acts on the input character string.

@return The method returns true if the line is blank, else false.
 */
bool isBlankLine( const char* str )
 {
    if ( strlen(str) <= 1 || str == NULL ) return true;

    std::string s(str);
    std::string::iterator it=s.begin();

    do {
          if (it == s.end()) return true;
       }
    while (*it >= 0 && *it <= 0x7f && std::isspace(*(it++)));

    return false;
 }



/**

If string commences with # sign true is returned, else
false is returned but if there is a # sign anywhere within the
string reading is terminated there.

@section arguments Input Arguments

The method acts on the input character string.

@return The method returns true if the line is a comment, else false.
*/
bool isCommentLine( char* str )
 {
    if ( str    == NULL )
        return false;
    if ( str[0] == '#' || str[0] == '%' )
        return true;
    if ( ( str[0] == '-' ) && ( str[1] == '-' ) )
        return true;

    const size_t strlength( strlen(str) );
    for ( size_t i = 0U; i<strlength; i++ )
      if ( str[i] == '#' || str[i] == '%' ){
           str[i] = '\0';
           break;
      }
      else if ( str[i] == '-' ) {
         if( i != strlength-1 )
             if ( str[i+1] == '-' ){
                 str[i] = '\0';
                 break;
             }
      }
    return false;
 }



/** Advances the file stream to behind the comment line.
  */
void advancePastCommentLine( std::ifstream& ifs )
 {
    char  text_line[512];
    do {
          ifs.getline( text_line, 512 );
       }
    while ( (!isCommentLine(text_line) && !ifs.eof()) );

    std::cout <<"\n\tSkipped comment line: "<< text_line << std::endl;

 } // end

bool isYES( const std::string& text_line )
{
    // positive answer: yes, exist, use it, etc.
    std::set<std::string> yes;
    yes.insert("Y");
    yes.insert("y");
    yes.insert("Yes");
    yes.insert("YES");
    yes.insert("YEs");
    yes.insert("Yes");
    yes.insert("Yep");
    yes.insert("YEP");
    yes.insert("YEp");
    yes.insert("Yep");
    yes.insert("OK");
    yes.insert("Ok");
    yes.insert("ok");

    // it is definite yes
    if( yes.find( text_line ) != yes.end() )
        return true;
    return false;
}


bool isNO( const std::string& text_line )
{
    // negative answer: no, no answer, not applicable, not avaliable, etc.
    std::set<std::string> no;
    no.insert("N");
    no.insert("n");
    no.insert("NO");
    no.insert("No");
    no.insert("no");
    no.insert("NOT");
    no.insert("NOt");
    no.insert("Not");
    no.insert("not");
    no.insert("NA");
    no.insert("Na");
    no.insert("na");

    // it is definite no
    if( no.find( text_line )  != no.end() )
        return false;

    throw csmp::Exception( FATAL_ERROR,
                           "yes_or_no:",
                           "The answer is uncertain!",
                           text_line );

    return false;
}



int yes_or_no( const std::string& text_line )
{
    // it is definite yes
    if( isYES( text_line ) )
        return 1;

    // it is definite no
    if( isNO( text_line ) )
        return 0;

    // The answer is uncertain!
    return 2;
}


// NUMBERS

bool isRealNumber( const std::string& s, int len )
{
    const size_t slen( s.length() );
    size_t nlen = ( len > slen ? slen : len );
    std::string str = s.substr(0,nlen);
    std::istringstream iss( str );
    double num = 0.0;
    if( ( iss >> num ).fail() )
        return false;
    return true;
}

bool isRealNumber( const char* s, int len )
{
    const size_t slen( strlen( s ) );
    size_t nlen = ( len > slen ? slen : len );
    std::string str( s, nlen );
    std::istringstream iss( str );
    double num = 0.0;
    if( ( iss >> num ).fail() )
        return false;
    return true;
}

bool isIntegerNumber( const std::string& s, int len )
{
    const size_t slen( s.length() );
    size_t nlen = ( len > slen ? slen : len );
    std::string str = s.substr(0,nlen);
    std::istringstream iss( str );
    int num = 0;
    if( ( iss >> num ).fail() )
        return false;
    return true;
}

bool isIntegerNumber( const char* s, int len )
{
    const size_t slen( strlen( s ) );
    size_t nlen = ( len > slen ? slen : len );
    std::string str( s, nlen );
    std::istringstream iss( str );
    int num = 0;
    if( ( iss >> num ).fail() )
        return false;
    return true;
}




// READ PROPERTY VALUES



/**

reads any of the CSP basic variables from the current string stream
tokenizing it with the delimiters ' ,:,\t,\n,\r'. Depending on the
variable type and the dimension of the model a certain number of
tokens is expected to be contained in the stream. Thus, for a
vector and tensor variables in 2D, 2 and 4 variable values are
expected, respectively. and dim and dim x dim in any other dimension.

@param sc The correctluy parsed data is returned into the variable that forms
the method argument.

@section messages Messages

The method will issue an error if the variable cannot be parsed correctly.
*/
void readPropertyValue( ScalarVariable& sc )
{
    char*        token(0);
    const char*  delims =" ,:,\t,\n,\r";

    token = strtok( NULL, delims );
    if ( token == NULL ) {
        sc = std::numeric_limits<double64>::quiet_NaN();
        cout <<"\n"<< token << endl;
        throw csmp::Exception(ERROR,
                              "readPropertyValue(scalar)",
                              "Property value could not be read");
    }
    else if ( !isdigit(token[0]) && token[0] != '-' )
    {
        sc = std::numeric_limits<double64>::quiet_NaN();
        cout <<"\n"<< token << endl;
        throw csmp::Exception(ERROR,
                              "readPropertyValue(scalar)",
                              "Property value is not a digit");
    }
    else
        sc = atof(token);

} // end


template<size_t dim>
void readPropertyValue( VectorVariable<dim>& vc )
{
    char*        token(0);
    const char*  delims =" ,:,\t,\n,\r";

    size_t i;
    vc = 0.0;
    for ( i = 0; i < dim; ++i )
    {
        token = strtok( NULL, delims );
        if ( token != NULL )
            vc(i) = atof(token);
        else if( i != 0 )
        {
            for ( size_t j = i; j < dim; ++j )
                vc(j) = vc(i-1);
            break;
        }
        else
            break;
    }

    if( i == 0 )
    {
        vc = std::numeric_limits<double64>::quiet_NaN();
        throw csmp::Exception( ERROR,
                               "readPropertyValue(vector)",
                               "Property value could not be read properly");
    }
}


template<size_t dim>
void readPropertyValue( TensorVariable<dim>& ts )
{
    char*  token(0);
    const char*  delims =" ,:,\t,\n,\r";
  
    double64 values[dim*dim];
    std::memset((void*)values, 0, sizeof(values));
  
  size_t value_count = 0;
  for (value_count = 0; value_count < dim*dim; ++value_count) {
    token = strtok( NULL, delims );
    if ( token == NULL ) {
      break;
    }
    values[value_count] = atof(token);
  }
  
  ts = 0.0;

  if (value_count == dim*dim) {
    // Interpret dim*dim numbers as the full tensor
    const double64* val = &values[0];
    for ( size_t i = 0; i < dim; ++i )
    {
      for ( size_t j = 0; j < dim; ++j )
      {
        ts(i,j) = *val++;
      }
    }
  }
  else if (value_count == dim) {
    // Interpret dim numbers as the tensor diagonal
    const double64* val = &values[0];
    for ( size_t i = 0; i < dim; ++i )
    {
      ts(i,i) = *val++;
    }
  }
  else if (value_count == 1) {
    // Interpret 1 number as the tensor diagonal
    for ( size_t i = 0; i < dim; ++i )
    {
      ts(i,i) = values[0];
    }
  }
  else {
    ts = std::numeric_limits<double64>::quiet_NaN();
    throw csmp::Exception( ERROR,
                          "readPropertyValue(tensor)",
                          "Property value could not be read properly");
  }

}

void readPropertyValue( ArrayVariable& av )
{
    char*        token(0);
    const char*  delims =" ,:,\t,\n,\r";

    size_t i;
    const size_t depth( av.Size() );
    for ( i=0; i < depth; ++i )
    {
        token = strtok( NULL, delims );
        if ( token != NULL )
            av(i) = atof(token);
        else if( i != 0 )
        {
            for ( size_t j = i; j < depth; j++ )
                av(j) = av(i-1);
            break;
        }
        else
            break;
    }

    if ( i == 0 )
    {
        av = std::numeric_limits<double64>::quiet_NaN();
        throw csmp::Exception( ERROR,
                               "readPropertyValue(array)",
                               "Property value could not be read properly");
    }
}

void readPropertyValue( FlaggedArrayVariable& fv )
{
    char*        token(0);
    const char*  delims =" ,:,\t,\n,\r";

    size_t i;
    const size_t depth( fv.Size() );
    for ( i=0; i < depth; i++ )
    {
        token = strtok( NULL, delims );
        if ( token != NULL )
            fv(i) = atof(token);
        else if( i != 0 )
        {
            for ( size_t j = i; j < depth; j++ )
                fv(j) = fv(i-1);
            break;
        }
        else
            break;
    }

    if ( i == 0 )
    {
        fv  = std::numeric_limits<double64>::quiet_NaN();
        throw csmp::Exception( ERROR,
                               "readPropertyValue(flagged array)",
                               "Property value could not be read properly");
    }
}

// EXPLICIT INSTANTIATIONS


template void readPropertyValue<1U>( VectorVariable<1U>& );
template void readPropertyValue<2U>( VectorVariable<2U>& );
template void readPropertyValue<3U>( VectorVariable<3U>& );

template void readPropertyValue<1U>( TensorVariable<1U>& );
template void readPropertyValue<2U>( TensorVariable<2U>& );
template void readPropertyValue<3U>( TensorVariable<3U>& );


// READ PROPERTY STATUS

void readPropertyStatusOfScalar( size_t depth, bool& digit, size_t& position, vector<VARIABLE_FLAG>& flags )
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    char*        token(0);
    const char*  delims =" ,:,\t,\n,\r";
    std::string  flag_name;

    token = strtok( NULL, delims );
    if ( token == NULL )
    {
        cout <<"\n"<< token << endl;
        csmp_error.notice( ERROR,
                           "readPropertyStatusOfScalar:",
                           "Flag of 'scalar' property could not be read");
        terminate();
    }
    else
    {
        if ( !isdigit(token[0]) )
        {
            digit = false;
            flags.resize(1U,parseCondition( flag_name = token ) );
        }
        else
        {
            digit = true;
            flags.resize(1U);

            position = 0U;
            token = strtok( NULL, delims );
            if ( token == NULL )
            {
                cout <<"\nflag"<< position <<": "<< token << endl;
                csmp_error.notice( ERROR,
                                   "readPropertyStatusOfScalar:",
                                   "Flag of 'scalar' property could not be read");
                terminate();
            }
            else
                flags[0U] = parseCondition( flag_name = token );
        }
    }
} // end


void readPropertyStatusOfVector( size_t depth, bool& digit, size_t& position, vector<VARIABLE_FLAG>& flags )
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    char*        token(0);
    const char*  delims =" ,:,\t,\n,\r";
    std::string  flag_name;

    token = strtok( NULL, delims );
    if ( token == NULL )
    {
        cout <<"\n"<< token << endl;
        csmp_error.notice( ERROR,
                           "readPropertyStatusOfVector:",
                           "Flag of 'vector' property could not be read");
        terminate();
    }
    else
    {
        if ( !isdigit(token[0]) )
        {
            digit = false;
            flags.resize(depth,parseCondition( flag_name = token ) );

            flags[0] = parseCondition( flag_name = token );
            for ( size_t i=1U; i<depth; i++ )
            {
                token=strtok( NULL, delims );
                if ( token == NULL )
                {
                    for ( size_t j=i; j<depth; j++ )
                        flags[j] = flags[i-1];
                    break;
                }
                else
                    flags[i] = parseCondition( flag_name = token );
            }
        }
        else
        {
            digit = true;
            flags.resize(1U);

            position = atoi(token);
            token=strtok( NULL, delims );
            if ( token == NULL )
            {
                cout <<"\nflag"<< position <<": "<< token << endl;
                csmp_error.notice( ERROR,
                                   "readPropertyStatusOfVector:",
                                   "Flag of 'vector' property could not be read");
                terminate();
            }
            else
                flags[0U] = parseCondition( flag_name = token );
        }
    }
}


void readPropertyStatusOfTensor( size_t depth, bool& digit, size_t& position, vector<VARIABLE_FLAG>& flags )
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    char*  token(0);
    const char*  delims =" ,:,\t,\n,\r";
    std::string  flag_name;

    token = strtok( NULL, delims );
    if ( token == NULL )
    {
        cout <<"\n"<< token << endl;
        csmp_error.notice( ERROR,
                           "readPropertyStatusOfTensor:",
                           "Flag of 'tensor' property could not be read");
        terminate();
    }
    else
    {
        if ( !isdigit(token[0]) )
        {
            digit = false;
            flags.resize(depth,parseCondition( flag_name = token ) );

            flags[0] = parseCondition( flag_name = token );
            const size_t depth2( depth*depth );
            for ( size_t i=0U; i<depth; i++ )
            {
                for ( size_t j=0U; j<depth; j++ )
                {
                    if( !(i==0) || !(j==0) )
                    {
                        token=strtok( NULL, delims );
                        if ( token == NULL )
                        {
                            for ( size_t n=i*depth+j ; n<depth2; n++ )
                                flags[ n ] = flags[i*depth+j-1];
                            break;
                        }
                        else
                            flags[i*depth+j] = parseCondition( flag_name = token );
                    }
                }
            }
        }
        else
        {
            digit = true;
            flags.resize(1U);

            position = 0U;
            token = strtok( NULL, delims );
            if ( token == NULL )
            {
                cout <<"\nflag"<< position <<": "<< token << endl;
                csmp_error.notice( ERROR,
                                   "readPropertyStatusOfTensor:",
                                   "Flag of 'tensor' property could not be read");
                terminate();
            }
            else
                flags[0U] = parseCondition( flag_name = token );
        }
    }
}

void readPropertyStatusOfArray( size_t depth, bool& digit, size_t& position, vector<VARIABLE_FLAG>& flags )
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    char*        token(0);
    const char*  delims =" ,:,\t,\n,\r";
    std::string  flag_name;

    token = strtok( NULL, delims );
    if ( token == NULL )
    {
        cout <<"\n"<< token << endl;
        csmp_error.notice( ERROR,
                           "readPropertyStatusOfArray:",
                           "Flag of 'array' property could not be read");
        terminate();
    }
    else
    {
        if ( !isdigit(token[0]) )
        {
            digit = false;
            flags.resize(1U,parseCondition( flag_name = token ) );
        }
        else
        {
            digit = true;
            flags.resize(1U);

            position = 0U;
            token = strtok( NULL, delims );
            if ( token == NULL )
            {
                cout <<"\nflag"<< position <<": "<< token << endl;
                csmp_error.notice( ERROR,
                                   "readPropertyStatusOfArray:",
                                   "Flag of 'array' property could not be read");
                terminate();
            }
            else
                flags[0U] = parseCondition( flag_name = token );
        }
    }
}

void readPropertyStatusOfFlaggedArray( size_t depth, bool& digit, size_t& position, vector<VARIABLE_FLAG>& flags )
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    char*        token(0);
    const char*  delims =" ,:,\t,\n,\r";
    std::string  flag_name;

    token = strtok( NULL, delims );
    if ( token == NULL )
    {
        cout <<"\n"<< token << endl;
        csmp_error.notice( ERROR,
                           "readPropertyStatusOfFlaggedArray:",
                           "Flag of 'flagged array' property could not be read");
        terminate();
    }
    else
    {
        if ( !isdigit(token[0]) )
        {
            digit = false;
            flags.resize(depth,parseCondition( flag_name = token ) );
            flags[0] = parseCondition( flag_name = token );

            for ( size_t i=1U; i<depth; i++ )
            {
                token=strtok( NULL, delims );
                if ( token == NULL )
                {
                    for ( size_t j=i; j<depth; j++ )
                        flags[j] = flags[i-1];
                    break;
                }
                else
                    flags[i] = parseCondition( flag_name = token );
            }
        }
        else
        {
            digit = true;
            flags.resize(1U);

            position = 0U;
            token = strtok( NULL, delims );
            if ( token == NULL )
            {
                cout <<"\nflag"<< position <<": "<< token << endl;
                csmp_error.notice( ERROR,
                                   "readPropertyStatusOfFlaggedArray:",
                                   "Flag of 'flagged array' property could not be read");
                terminate();
            }
            else
                flags[0U] = parseCondition( flag_name = token );
        }
    }
}


// PRINTING PROPERTY VALUES

/**

Writes any basic CSMP variable to screen in a more compact format than
is achieved by the Out() interface.

@section arguments Input Arguments

A CSMP ScalarVariable, VectorVariable, TensorVariable, ArrayVariable, FlaggedArrayVariable objects.
*/

void printPropertyValue( const ScalarVariable& sc )
{
    cout << sc << endl;
}

template<size_t dim>
void printPropertyValue( const VectorVariable<dim>& vc )
{
    for ( size_t i=0; i<dim; i++ )
        cout << vc[i] <<", ";
    cout << endl;
}


template<size_t dim>
void printPropertyValue( const TensorVariable<dim>& ts )
{
    for ( size_t i=0; i<dim; i++ ) {
        cout <<"\t";
        for ( size_t j=0; j<dim; j++ )
            cout << ts(i,j) <<", ";
        cout << endl;
    }
    cout << endl;
}

void printPropertyValue( const ArrayVariable& av )
{
    for ( size_t i=0; i<av.Size(); i++ )
        cout << av[i] <<", ";
}

void printPropertyValue( const FlaggedArrayVariable& fv )
{
    for ( size_t i=0; i<fv.Size(); i++ )
        cout << fv[i] <<", ";
}

// EXPLICIT INSTANTIATIONS

template void printPropertyValue<1U>( const VectorVariable<1U>& );
template void printPropertyValue<2U>( const VectorVariable<2U>& );
template void printPropertyValue<3U>( const VectorVariable<3U>& );

template void printPropertyValue<1U>( const TensorVariable<1U>& );
template void printPropertyValue<2U>( const TensorVariable<2U>& );
template void printPropertyValue<3U>( const TensorVariable<3U>& );





// REPORTING ABOUT READ DATA

/**

Echoes the supplied arguments to screen so that the user gets a feedback
whether all boundary condition values have been read correctly.

@section arguments Input Arguments

Strings that specify the name of the property its SI unit, and its
lower and upper bound, i.e.  physically meaningful range as specified
in the '*-variables.txt' text file.

The last two (2D) or 4 (3D) method arguments specify the boundary
condition values that were assigned for linear or tri-linear interpolation
assuming that the boundaries are planar. These can be either of the
types ScalarVariable, VectorVariable, TensorVariable, ArrayVariable, FlaggedArrayVariable.
*/
void reportAssignment( const string& prop,
                       const string& unit,
                       const string& bound,
                       ScalarVariable& sc1,
                       ScalarVariable& sc2 )
{
    string  property(" '"); property += prop;  property +="' ";
    string  si_unit(" [");  si_unit  += unit;  si_unit  +="] ";
    string  target(" '");   target   += bound; target   +="' ";

    if ( sc1() == sc2() ) {
        cout <<"\tassigned fixed value to"<< setw(35) << property << setw(15) << si_unit <<"of  "<< sc1();
        cout <<"\tto '"<< bound <<"' model boundary."<< endl;
    }
    else
    {
        cout <<"\tassigned linear '"<< prop <<"' gradient ["<< unit <<"]\t\t\t"<< sc1() <<" - "<< sc2();
        cout <<"\tto '"<< bound <<"' model boundary."<< endl;
    }
}


void reportAssignment( const string& prop,
                       const string& unit,
                       const string& bound,
                       ScalarVariable& sc1,
                       ScalarVariable& sc2,
                       ScalarVariable& sc3,
                       ScalarVariable& sc4 )
{
    if ( (sc1() == sc2() && sc2() == sc3()) || sc3() == sc4() ) {
        cout <<"\tassigned fixed '"<< prop <<"' value ["<< unit <<"]\t\t\t"<< sc1();
        cout <<"\tto '"<< bound <<"' model boundary."<< endl;
    }
    else
    {
        cout <<"\tassigned bi-linear '"<< prop <<"' gradient field ["<< unit <<"]\t\t\t";
        cout << sc1() <<" - "<< sc2() <<" - "<< sc3() <<" - "<< sc4();
        cout <<"\tto '"<< bound <<"' model boundary."<< endl;
    }
}


template<size_t dim>
void reportAssignment( const string& prop,
                       const string& unit,
                       const string& bound,
                       VectorVariable<dim>& vc1,
                       VectorVariable<dim>& vc2 )
{
    if ( vc1 == vc2 ) {
        cout <<"\tassigned fixed vector '"<< prop <<"' value ["<< unit <<"] ";
        cout <<" to '"<< bound <<"' model boundary:\t\t"<< endl;
        printPropertyValue( vc1 ); cout << endl;
    }
    else
    {
        cout <<"\tassigned linear vector '"<< prop <<"' gradient ["<< unit <<"] ";
        cout <<" to '"<< bound <<"' model boundary:"<< endl;
        cout <<"\t", printPropertyValue( vc1 ), cout <<" to ", printPropertyValue( vc2 );
        cout << endl;
    }
}


template<size_t dim>
void reportAssignment( const string& prop,
                       const string& unit,
                       const string& bound,
                       VectorVariable<dim>& vc1,
                       VectorVariable<dim>& vc2,
                       VectorVariable<dim>& vc3,
                       VectorVariable<dim>& vc4 )
{
    if ( (vc1 == vc2 && vc2 == vc3) || vc3 == vc4 ) {
        cout <<"\tassigned fixed vector '"<< prop <<"' value ["<< unit <<"] ";
        cout <<" to '"<< bound <<"' model boundary:\t\t";
        printPropertyValue( vc1 );
    }
    else
    {
        cout <<"\tassigned bi-linear '"<< prop <<"' gradient field ["<< unit <<"] ";
        cout <<" to '"<< bound <<"' model boundary. Corner values:"<< endl;
        cout <<"\t", printPropertyValue( vc1 ), cout <<" to ", printPropertyValue( vc2 ), cout << endl;
        cout <<"\t", printPropertyValue( vc3 ), cout <<" to ", printPropertyValue( vc4 ), cout << endl;
        cout << endl;
    }
}

template<size_t dim>
void reportAssignment( const string& prop,
                       const string& unit,
                       const string& bound,
                       TensorVariable<dim>& ts1,
                       TensorVariable<dim>& ts2 )
{
    if ( ts1 == ts2 ) {
        cout <<"\tassigned fixed vector '"<< prop <<"' value ["<< unit <<"] ";
        cout <<" to '"<< bound <<"' model boundary:\t\t"<< endl;
        printPropertyValue( ts1 ); cout << endl;
    }
    else
    {
        cout <<"\tassigned linear vector '"<< prop <<"' gradient ["<< unit <<"] ";
        cout <<" to '"<< bound <<"' model boundary:"<< endl;
        cout <<"\t", printPropertyValue( ts1 ), cout <<" to ", printPropertyValue( ts2 );
        cout << endl;
    }
}


template<size_t dim>
void reportAssignment( const string& prop,
                       const string& unit,
                       const string& bound,
                       TensorVariable<dim>& ts1,
                       TensorVariable<dim>& ts2,
                       TensorVariable<dim>& ts3,
                       TensorVariable<dim>& ts4 )
{
    if ( (ts1 == ts2 && ts2 == ts3) || ts3 == ts4 ) {
        cout <<"\tassigned fixed tensor '"<< prop <<"' value ["<< unit <<"] ";
        cout <<" to '"<< bound <<"' model boundary:\t\t";
        printPropertyValue( ts1 );
    }
    else
    {
        cout <<"\tassigned bi-linear '"<< prop <<"' gradient field ["<< unit <<"] ";
        cout <<" to '"<< bound <<"' model boundary. Corner values:"<< endl;
        cout <<"\t", printPropertyValue( ts1 ), cout <<" to ", printPropertyValue( ts2 ), cout << endl;
        cout <<"\t", printPropertyValue( ts3 ), cout <<" to ", printPropertyValue( ts4 ), cout << endl;
        cout << endl;
    }
}

void reportAssignment( const string& prop,
                       const string& unit,
                       const string& bound,
                       ArrayVariable& av1,
                       ArrayVariable& av2 )
{
    string  property(" '"); property += prop;  property +="' ";
    string  si_unit(" [");  si_unit  += unit;  si_unit  +="] ";
    string  target(" '");   target   += bound; target   +="' ";

    if ( av1 == av2 ) {
        cout <<"\tassigned fixed array '"<< prop <<"' value ["<< unit <<"] ";
        cout <<" to '"<< bound <<"' model boundary:\t\t";
        printPropertyValue( av1 );
    }
    else
    {
        cout <<"\tassigned bi-linear '"<< prop <<"' gradient field ["<< unit <<"]\t\t\t";
        cout <<" to '"<< bound <<"' model boundary. Corner values:"<< endl;
        cout <<"\t", printPropertyValue( av1 ), cout <<" to ", printPropertyValue( av2 ), cout << endl;
        cout << endl;
    }
}


void reportAssignment( const string& prop,
                       const string& unit,
                       const string& bound,
                       ArrayVariable& av1,
                       ArrayVariable& av2,
                       ArrayVariable& av3,
                       ArrayVariable& av4 )
{
    if ( ( av1 == av2 && av2 == av3 ) || av3 == av4 ) {
        cout <<"\tassigned fixed array '"<< prop <<"' value ["<< unit <<"] ";
        cout <<" to '"<< bound <<"' model boundary:\t\t";
        printPropertyValue( av1 );
    }
    else
    {
        cout <<"\tassigned bi-linear '"<< prop <<"' gradient field ["<< unit <<"]\t\t\t";
        cout <<" to '"<< bound <<"' model boundary. Corner values:"<< endl;
        cout <<"\t", printPropertyValue( av1 ), cout <<" to ", printPropertyValue( av2 ), cout << endl;
        cout <<"\t", printPropertyValue( av3 ), cout <<" to ", printPropertyValue( av4 ), cout << endl;
        cout << endl;
    }
}


void reportAssignment( const string& prop,
                       const string& unit,
                       const string& bound,
                       FlaggedArrayVariable& fv1,
                       FlaggedArrayVariable& fv2 )
{
    string  property(" '"); property += prop;  property +="' ";
    string  si_unit(" [");  si_unit  += unit;  si_unit  +="] ";
    string  target(" '");   target   += bound; target   +="' ";

    if ( fv1 == fv2 ) {
        cout <<"\tassigned fixed array '"<< prop <<"' value ["<< unit <<"] ";
        cout <<" to '"<< bound <<"' model boundary:\t\t";
        printPropertyValue( fv1 );
    }
    else
    {
        cout <<"\tassigned bi-linear '"<< prop <<"' gradient field ["<< unit <<"]\t\t\t";
        cout <<" to '"<< bound <<"' model boundary. Corner values:"<< endl;
        cout <<"\t", printPropertyValue( fv1 ), cout <<" to ", printPropertyValue( fv2 ), cout << endl;
        cout << endl;
    }
}


void reportAssignment( const string& prop,
                       const string& unit,
                       const string& bound,
                       FlaggedArrayVariable& fv1,
                       FlaggedArrayVariable& fv2,
                       FlaggedArrayVariable& fv3,
                       FlaggedArrayVariable& fv4 )
{
    if ( ( fv1 == fv2 && fv2 == fv3 ) || fv3 == fv4 ) {
        cout <<"\tassigned fixed array '"<< prop <<"' value ["<< unit <<"] ";
        cout <<" to '"<< bound <<"' model boundary:\t\t";
        printPropertyValue( fv1 );
    }
    else
    {
        cout <<"\tassigned bi-linear '"<< prop <<"' gradient field ["<< unit <<"]\t\t\t";
        cout <<" to '"<< bound <<"' model boundary. Corner values:"<< endl;
        cout <<"\t", printPropertyValue( fv1 ), cout <<" to ", printPropertyValue( fv2 ), cout << endl;
        cout <<"\t", printPropertyValue( fv3 ), cout <<" to ", printPropertyValue( fv4 ), cout << endl;
        cout << endl;
    }
}

// EXPLICIT INSTANTIATIONS


template void reportAssignment<1U>( const std::string& prop, const std::string& unit, const std::string& bound,
                                    VectorVariable<1U>& vc1, VectorVariable<1U>& vc2 );

template void reportAssignment<1U>( const std::string& prop, const std::string& unit, const std::string& bound,
                                    VectorVariable<1U>& vc1, VectorVariable<1U>& vc2,
                                    VectorVariable<1U>& vc3, VectorVariable<1U>& vc4 );

template void reportAssignment<2U>( const std::string& prop, const std::string& unit, const std::string& bound,
                                    VectorVariable<2U>& vc1, VectorVariable<2U>& vc2 );

template void reportAssignment<2U>( const std::string& prop, const std::string& unit, const std::string& bound,
                                    VectorVariable<2U>& vc1, VectorVariable<2U>& vc2,
                                    VectorVariable<2U>& vc3, VectorVariable<2U>& vc4 );

template void reportAssignment<3U>( const std::string& prop, const std::string& unit, const std::string& bound,
                                    VectorVariable<3U>& vc1, VectorVariable<3U>& vc2 );

template void reportAssignment<3U>( const std::string& prop, const std::string& unit, const std::string& bound,
                                    VectorVariable<3U>& vc1, VectorVariable<3U>& vc2,
                                    VectorVariable<3U>& vc3, VectorVariable<3U>& vc4 );

template void reportAssignment<1U>( const std::string& prop, const std::string& unit, const std::string& bound,
                                    TensorVariable<1U>& ts1, TensorVariable<1U>& ts2 );

template void reportAssignment<1U>( const std::string& prop, const std::string& unit, const std::string& bound,
                                    TensorVariable<1U>& ts1, TensorVariable<1U>& ts2,
                                    TensorVariable<1U>& ts3, TensorVariable<1U>& ts4 );

template void reportAssignment<2U>( const std::string& prop, const std::string& unit, const std::string& bound,
                                    TensorVariable<2U>& ts1, TensorVariable<2U>& ts2 );

template void reportAssignment<2U>( const std::string& prop, const std::string& unit, const std::string& bound,
                                    TensorVariable<2U>& ts1, TensorVariable<2U>& ts2,
                                    TensorVariable<2U>& ts3, TensorVariable<2U>& ts4 );

template void reportAssignment<3U>( const std::string& prop, const std::string& unit, const std::string& bound,
                                    TensorVariable<3U>& ts1, TensorVariable<3U>& ts2 );

template void reportAssignment<3U>( const std::string& prop, const std::string& unit, const std::string& bound,
                                    TensorVariable<3U>& ts1, TensorVariable<3U>& ts2,
                                    TensorVariable<3U>& ts3, TensorVariable<3U>& ts4 );






// CONFIG FILE
bool isConfigFileExist( const std::string& filename_prefix )
{
    std::string config_filename = filename_prefix;
    config_filename += "-configuration.txt";
    if( !csmp::isFileExist( config_filename ) )
        return false;
    return true;
}

// VARIABLES FILE
bool isVariablesFileExist( const std::string& filename_prefix )
{
    std::string config_filename = filename_prefix;
    config_filename += "-variables.txt";
    if( !csmp::isFileExist( config_filename ) )
        return false;
    return true;
}




// READ POINT DATA


/**

Parses the supplied ASCII text line to read the point name, its
coordinates and the a users-specified number of scalar data.

@section arguments Input Arguments

The second argument gives the number of scalar data that shall be
retrieved from the text string (third method argument).

@return The point data record will be stored in the argument map to be retrieved
by the point name.

The method returns true if the reading process was successful, i.e.
the correct number of scalar data could be read.

@section messages Messages

The method reports problems during the reading process.

*/
template<size_t dim>
bool readPointData( std::map<std::string,std::vector<double64> >& pdata,
                    int ndata,
                    std::ifstream& ifs, char* text_line, size_t line_length, bool verbose )
{
    char*  token(0);
    const char*  delims =":,\t,\n,\r";
    string  point_name;

    cout <<"\n\nreading well names, locations, associated flow rates and saturations..."<< endl << endl;

    do {
        if ( !isCommentLine(text_line) )
        {
            // 1. tokenizing string and interpreting point name
            token = strtok( text_line, delims );
            if ( token == NULL ) {
                cout <<"\n"<< text_line << endl;
                throw csmp::Exception(ERROR,"readPointData",
                                      "point name could not be read from input file.");
                return false;
            }
            else if ( strcmp( "no point data", token ) == 0 )
                return false;
            else
                point_name = token;

            // 2. configuring the output map
            pair<typename map<string,vector<double64> >::iterator,bool>  it =
                    pdata.insert( make_pair( point_name, vector<double64>(dim+ndata)) );
            if ( !it.second )
                throw csmp::Exception(ERROR,"readPointData",
                                      "Duplicate record for", token );

            // 3. point coordinates and data
            for ( size_t i=0; i<(dim + ndata); i++ ) {
                token = strtok( NULL, delims );
                if ( token == NULL ) {
                    cout <<"\n"<< text_line << endl;
                    throw csmp::Exception(ERROR,"readPointData",
                                          "Point data read incompletely for", point_name.c_str() );
                    return false;
                }
                else
                    (*it.first).second[i] = atof(token);
            }

            return true;
        }
        ifs.getline( text_line, line_length );
    }
    while ( !isBlankLine(text_line) && !ifs.eof() );

   return true;

} // end readPointData

template bool readPointData<1U>(std::map<std::string,std::vector<double64> >&,int,std::ifstream&,char*,size_t,bool);
template bool readPointData<2U>(std::map<std::string,std::vector<double64> >&,int,std::ifstream&,char*,size_t,bool);
template bool readPointData<3U>(std::map<std::string,std::vector<double64> >&,int,std::ifstream&,char*,size_t,bool);



/**

This method allows the user to define regions on the basis of input
data ranges and is used if the supplied finite element mesh has no
region names associated with it.

The method reads Block 1 of the file structure outlined in the class
documentation further below. The text string that it interprets for
this purpose must have the form:

@code
joint	tab		permeability	tab		1.0e-8 1.0e-8
@endcode

Where the first string is that regionname that will be used, then comes
the name of the variable whose range will be tested and then the range
of variable values that shall discern the finite elements that will be
accumulated into the intended region. Importantly, the variable that
is used to form the region must be an element property because any
region must have at least one element.

@section arguments Input Arguments

The method takes a reference to the Model object in which the
region shall be created and the textline which contains the data on
the basis of which the region shall be formed.

@return The method returns the name of the newly created region into its
third argument.

The method also returns the boolean variable true if the attempt to build
the region was successful.

@section implementation Implementation

The region is created directly so that it can be accessed immediately
after the method is called provided that the call was successful.

@section messages Messages

The method provides feedback on the intended operation.
 */
template<size_t dim>
bool buildRegionsBasedOnPropertyRange( Model<dim>& model,
                                       std::set<std::string>& groups,
                                       std::ifstream& ifs, char* text_line, size_t line_length, bool verbose )
{
    char*  token(0);
    const char*  delims =":,\t,\n,\r";
    string  group_name, prop_name;
    double  prop_min, prop_max;

    if( verbose )
        cout <<"\n\nIdentifying model subregions as Regions based on particular property value range..."<< endl;

    do {
        if ( !isCommentLine(text_line) ){

            // 1. tokenizing string and interpreting region name
            token = strtok( text_line, delims );
            if ( token == NULL )
            {
                cout <<"\n"<< text_line << endl;
                throw csmp::Exception(ERROR,
                                      "buildRegionsBasedOnPropertyRange",
                                      "region name could not be read from input file.");
                return false;
            }
            else if ( strcmp( "no regions", token ) == 0 )
            {
                return true;
            }
            else
                group_name = token;

            // 2. region formation property
            token = strtok( NULL, delims );
            if ( token == NULL )
            {
                std::cerr <<"\n"<< text_line << endl;
                throw csmp::Exception(ERROR,
                                      "buildRegionsBasedOnPropertyRange",
                                      "property name could not be read from input file.");
                return false;
            }
            else
                prop_name = token;

            // 3. property range (value 1)
            token = strtok( NULL, delims );
            if ( token == NULL )
            {
                std::cerr <<"\n"<< text_line << endl;
                throw csmp::Exception(ERROR,
                                      "buildRegionsBasedOnPropertyRange",
                                      "property range value 1 could not be read from input file.");
                return false;
            }
            else
                prop_min = atof(token);

            // 4. property range (value 2)
            token = strtok( NULL, delims );
            if ( token == NULL )
            {
                std::cerr <<"\n"<< text_line << endl;
                throw csmp::Exception(ERROR,
                                      "buildRegionsBasedOnPropertyRange",
                                      "property range value 2 could not be read from input file.");
                return false;
            }
            else
                prop_max = atof(token);

            // 5. building region using the aquired input data
            if ( verbose )
            {
                cout <<"\nSimulationConfigurator<"<< dim <<">::buildRegionFrom: Forming region: '";
                cout << group_name <<"' using '"<< prop_name <<"' range: ";
                cout << prop_min <<" - "<< prop_max << endl;
            }
            model.FormRegionFrom( group_name.c_str(), prop_name.c_str(), prop_min, prop_max );

            if ( !model.ContainsRegion( group_name.c_str() ) )
            {
                throw csmp::Exception( FATAL_ERROR,
                                       "buildRegionsBasedOnPropertyRange",
                                       "region could not be build", text_line );
                return false;
            }
            else
                groups.insert( group_name );
        }
        ifs.getline( text_line, line_length );
    }
    while ( !isBlankLine(text_line) && !ifs.eof() );

    return true;

} // end buildRegionsBasedOnPropertyRange


template bool buildRegionsBasedOnPropertyRange<1U>(Model<1U>&,std::set<std::string>&,std::ifstream&,char*,size_t,bool);
template bool buildRegionsBasedOnPropertyRange<2U>(Model<2U>&,std::set<std::string>&,std::ifstream&,char*,size_t,bool);
template bool buildRegionsBasedOnPropertyRange<3U>(Model<3U>&,std::set<std::string>&,std::ifstream&,char*,size_t,bool);



// READ PROPERTIES AND BOUNDARY CONDITIONS
// NOTE: doesn't support all the latest functionality
// RECOMMENDATION: use instead functionality to read properties and boundary conditions separately from each other


/**

This method reads data of the type Block 4, see class documentation below.
It assigns boundary conditions to a the supplied box-shaped model.
Any of the BOX_BOUNDARY enumeration switches may be used to specify the
boundary conditions. The input should have the following form:

@code
# a uniform (scalar) fluid pressure value is applied to a line 2D model
# BOX_BOUNDARY   condition		variable name	variable values
TOP		tab		Dirichlet	tab	fluid pressure	tab	1.0 1.0

# 2D model, (vector) variable velocity with x, y components specified for
# the 2 endpoints of the linear boundary. The x component is fixed, the
# y component can be modified by the computation. The PLAIN condition
# flag is ignored since velocity is a vector property
# the flags after the variable values specifiy the desired condition
# flags for the x and y velocity components
LEFT    PLAIN   velocity    5.  0.  5.  0.  Dirichlet   Plain
@endcode

Importantly, behind the values that shall be assigned to end points of
the linear boundary (2D) or the corner points of the planar boundary
in 3D only two or three flags can / must be specified. These overwrite
the boundary condition flag specified earlier in the line in the
case where the property is a vector or tensor variable the individual
components of which should be treatable separately.

@section arguments Input Arguments

The method takes a reference to the Model representing the box
shaped model to which the boundary condition shall be applied to.

@section implementation Implementation

Thus far assignments can only be made for scalar and vector type
CSP basic variables.

@section application Application

Only for box-shaped models or models with an irregular top surface. In
this latter case TOP must not be used but the boundary can be accessed
if was properly assigned an IRREGULAR flag instead.

@section messages Messages

The method reports potential errors that may occur if the input textline
contains errors or imcomplete variable specifications.
*/
template<size_t dim>
bool readBoxBoundaryPropertyValuesAndConditions( Model<dim>& model,
                                                 std::ifstream& ifs, char* text_line, size_t line_length, bool verbose )
{
    char*               token(0);
    const char*         delims =":,\t,\n,\r";
    string              bound_name, cond_type, prop_name, flag1, flag2, flag3, unit;
    VARIABLE_TYPE       prop_type;
    size_t              length;
    ScalarVariable      sc1, sc2, sc3, sc4;
    VectorVariable<dim> vc1, vc2, vc3, vc4;
    TensorVariable<dim> ts1, ts2, ts3, ts4;
    double64            val;

    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    if( verbose )
        cout <<"\n\nAssigning essential property values and flags to Box-Boundaries..."<< endl << endl;

    do {
        if ( !isCommentLine(text_line) )
        {
            // 1. name of boundary
            bound_name = strtok( text_line, delims );

            // 2. type of boundary condition
            token = strtok( NULL, delims );
            if ( token == NULL )
            {
                std::cerr <<"\n"<< text_line << endl;
                csmp_error.notice( ERROR,
                                   "readBoxBoundaryPropertyValuesAndConditions",
                                   "Condition type specifier missing (Dirichlet, Neumann...");
                return false;
            }
            else
                cond_type = token;

            // 3. property name
            token = strtok( NULL, delims );
            if ( token == NULL )
            {
                std::cerr <<"\n"<< text_line << endl;
                csmp_error.notice( ERROR,
                                   "readBoxBoundaryPropertyValuesAndConditions",
                                   "Property name could not be read");
                return false;
            }
            else
                prop_name = token;

            // 3. property value & assignment to Model<dim>
            prop_type = model.Database().Type( prop_name.c_str() );
            length    = model.Database().Components( prop_name.c_str() );
            unit      = model.Database().Unit( prop_name.c_str() );

            if ( prop_type == SCALAR ) {
                vector<ScalarVariable >  bvalues;
                readPropertyValue( sc1 );
                model.Database().CheckRange( prop_name.c_str(), sc1() );
                readPropertyValue( sc2 );
                model.Database().CheckRange( prop_name.c_str(), sc2() );
                if ( dim == 1U || dim == 2U ) {
                    if ( sc1 == sc2 )
                        model.InputBoundaryValue( parseBoundary(bound_name), prop_name.c_str(),
                                               makeScalar( parseCondition(cond_type), sc1() ) );
                    else {
                        bvalues.reserve(2U);
                        bvalues.push_back( sc1 ); assert( parseCondition(cond_type) == sc1.Flag() );
                        bvalues.push_back( sc2 ); assert( parseCondition(cond_type) == sc2.Flag() );
                        model.InterpolateBoundaryValues( parseBoundary(bound_name), prop_name.c_str(),
                                                      bvalues );
                    }
                    reportAssignment( prop_name, unit, bound_name, sc1, sc2 );
                }
                else {
                    readPropertyValue( sc3 );
                    model.Database().CheckRange( prop_name.c_str(), sc3() );
                    readPropertyValue( sc4 );
                    model.Database().CheckRange( prop_name.c_str(), sc4() );
                    if ( sc1 == sc2 and sc2 == sc3 and sc3 == sc4 )
                        model.InputBoundaryValue( parseBoundary(bound_name), prop_name.c_str(),
                                               makeScalar( parseCondition(cond_type), sc1() ) );
                    else {
                        bvalues.reserve(4U);
                        bvalues.push_back( sc1 ); assert( parseCondition(cond_type) == sc1.Flag() );
                        bvalues.push_back( sc2 ); assert( parseCondition(cond_type) == sc2.Flag() );
                        bvalues.push_back( sc3 ); assert( parseCondition(cond_type) == sc3.Flag() );
                        bvalues.push_back( sc4 ); assert( parseCondition(cond_type) == sc4.Flag() );
                        model.InterpolateBoundaryValues( parseBoundary(bound_name), prop_name.c_str(),
                                                      bvalues );
                    }
                    reportAssignment( prop_name, unit, bound_name, sc1, sc2, sc3, sc4 );
                }
            }
            else if ( prop_type == VECTOR ) {
                vector<VectorVariable<dim> >  bvalues;
                readPropertyValue( vc1 );
                val = vc1.Length(); model.Database().CheckRange( prop_name.c_str(), val );
                readPropertyValue( vc2 );
                val = vc2.Length(); model.Database().CheckRange( prop_name.c_str(), val );
                if ( dim == 1U || dim == 2U ) {
                    // reading setting condition flags for assignment
                    if ( (token=strtok( NULL, delims )) == NULL ) { // condition flag 1
                        if ( token != nullptr ) std::cerr <<"\n"<< token << endl;
                        csmp_error.notice( ERROR,
                                           "readBoxBoundaryPropertyValuesAndConditions",
                                           "Vector variable entry 1 VARIABLE_FLAG could not be read");
                        return false;
                    }
                    else flag1 = token;
                    if ( (token=strtok( NULL, delims )) == NULL ) { // condition flag 2
                        if ( token != nullptr ) std::cerr <<"\n"<< token << endl;
                        csmp_error.notice( ERROR,
                                           "readBoxBoundaryPropertyValuesAndConditions",
                                           "Vector variable entry 2 VARIABLE_FLAG could not be read");
                        return false;
                    }
                    else flag2 = token;
                    vc1.Flag(0) = vc2.Flag(0) = parseCondition(flag1);
                    vc1.Flag(1) = vc2.Flag(1) = parseCondition(flag2);
                    if ( vc1 == vc2 )
                        model.InputBoundaryValue( parseBoundary(bound_name), prop_name.c_str(), vc1 );
                    else {
                        bvalues.reserve(2U);
                        bvalues.push_back( vc1 );
                        bvalues.push_back( vc2 );
                        model.InterpolateBoundaryValues( parseBoundary(bound_name), prop_name.c_str(), bvalues );
                    }
                    reportAssignment( prop_name, unit, bound_name, vc1, vc2 );
                }
                else {
                    readPropertyValue( vc3 );
                    val = vc3.Length(); model.Database().CheckRange( prop_name.c_str(), val );
                    readPropertyValue( vc4 );
                    val = vc4.Length(); model.Database().CheckRange( prop_name.c_str(), val );
                    // reading and setting condition flags for assignment
                    if ( (token=strtok( NULL, delims )) == NULL ) { // condition flag 1
                        std::cerr <<"\n"<< token << endl;
                        csmp_error.notice( ERROR,
                                           "readBoxBoundaryPropertyValuesAndConditions",
                                           "Vector variable entry 1 VARIABLE_FLAG could not be read");
                        return false;
                    }
                    else flag1 = token;
                    if ( (token=strtok( NULL, delims )) == NULL ) { // condition flag 2
                        std::cerr <<"\n"<< token << endl;
                        csmp_error.notice( ERROR,
                                           "readBoxBoundaryPropertyValuesAndConditions",
                                           "Vector variable entry 2 VARIABLE_FLAG could not be read");
                        return false;
                    }
                    else flag2 = token;
                    if ( (token=strtok( NULL, delims )) == NULL ) { // condition flag 3
                        std::cerr <<"\n"<< token << endl;
                        csmp_error.notice( ERROR,
                                           "readBoxBoundaryPropertyValuesAndConditions",
                                           "Vector variable entry 3 VARIABLE_FLAG could not be read");
                        return false;
                    }
                    else flag3 = token;
                    vc1.Flag(0) = vc2.Flag(0) = vc3.Flag(0) = vc4.Flag(0) = parseCondition(flag1);
                    vc1.Flag(1) = vc2.Flag(1) = vc3.Flag(1) = vc4.Flag(1) = parseCondition(flag2);
                    vc1.Flag(2) = vc2.Flag(2) = vc3.Flag(2) = vc4.Flag(2) = parseCondition(flag3);
                    if ( vc1 == vc2 and vc2 == vc3 and vc3 == vc4 )
                        model.InputBoundaryValue( parseBoundary(bound_name), prop_name.c_str(), vc1 );
                    else {
                        bvalues.reserve(4U);
                        bvalues.push_back( vc1 );
                        bvalues.push_back( vc2 );
                        bvalues.push_back( vc3 );
                        bvalues.push_back( vc4 );
                        model.InterpolateBoundaryValues( parseBoundary(bound_name), prop_name.c_str(), bvalues );
                    }
                    reportAssignment( prop_name, unit, bound_name, vc1, vc2, vc3, vc4 );
                }
            }
            else if ( prop_type == TENSOR ) {
                parseCondition(cond_type);
                readPropertyValue( ts1 );
                val = ts1.MinElement();
                model.Database().CheckRange( prop_name.c_str(), val );
                val = ts1.MaxElement();
                model.Database().CheckRange( prop_name.c_str(), val );
                readPropertyValue( ts2 );
                val = ts2.MinElement();
                model.Database().CheckRange( prop_name.c_str(), val );
                val = ts2.MaxElement();
                model.Database().CheckRange( prop_name.c_str(), val );
                if ( dim == 1U || dim == 2U )
                    ; // model.AssignBoundaryValues( parseBoundary(bound_name), prop_name.c_str(), ts1, ts2 );
                else {
                    readPropertyValue( ts3 );
                    val = ts3.MinElement();
                    model.Database().CheckRange( prop_name.c_str(), val );
                    val = ts3.MaxElement();
                    model.Database().CheckRange( prop_name.c_str(), val );
                    readPropertyValue( ts4 );
                    val = ts4.MinElement();
                    model.Database().CheckRange( prop_name.c_str(), val );
                    val = ts4.MaxElement();
                    model.Database().CheckRange( prop_name.c_str(), val );
                }
            }
            else if( prop_type == ARRAY )
            {
                throw csmp::Exception( FATAL_ERROR,"readBoxBoundaryPropertyValuesAndConditions",
                                       "Array box-boundary condition assignment not implemented yet");
                return false;
            }
            else if( prop_type == FLAGGEDARRAY )
            {
                throw csmp::Exception( FATAL_ERROR,"readBoxBoundaryPropertyValuesAndConditions",
                                       "Flagged Array box-boundary condition assignment not implemented yet");
                return false;
            }
        }
        ifs.getline( text_line, line_length );
    }
    while ( !isBlankLine(text_line) && !ifs.eof() );

    return true;

} // end readBoxBoundaryPropertyValuesAndConditions

template bool readBoxBoundaryPropertyValuesAndConditions<1U>(Model<1U>&,std::ifstream&,char*,size_t,bool);
template bool readBoxBoundaryPropertyValuesAndConditions<2U>(Model<2U>&,std::ifstream&,char*,size_t,bool);
template bool readBoxBoundaryPropertyValuesAndConditions<3U>(Model<3U>&,std::ifstream&,char*,size_t,bool);

/** ReadBoundaryPropertyValuesAndConditions

Example
=======
Sets variable at the boundary of interest to the values and conditions imposed in the textline.
@code
# a fluid pressure of 200 bar is assigned to the interior of a boundary named
# 'BOUNDARY3'
BOUNDARY3	tab interior tab DIRICH tab fluid pressure	tab   2.0E+7
#
# TOP is assigned values for the vector variable 'displacement
# Dirichlet names the condition that is applied; however only to the Y-component
TOP	complete	Dirichlet	displacement	Plain	Dirichlet	Plain	0.	-4.356	0.
@endcode

@attention Note the duplicate use of flags in the case of vector or tensor
properties. The first flag specifies which flag the method shall not overwrite in case
a value is already set to this flag at the boundary.
This controls the situation, for instance, where nodes that are shared among multiple
boundaries have already been set a Dirichlet constraint for one of the vector components
while the other is still up for assignment.
Setting the first flag to Dirichlet allows this because it signals to the method
not to overwrite preexisting Dirichlet constraints (anywhere on the boundary).

@attention This applies only to unique boundaries.

@attention In case of vector and tensor variables, flags are assigned to each component.
This means whatever was on the boundary before is overwritten.

*/
template<size_t dim>
bool readBoundaryPropertyValuesAndConditions( Model<dim>& sg,
                                              std::ifstream& ifs, char* text_line, size_t line_length, bool verbose )
{
    char*        token(0);
    const char*  delims =":,\t,\n,\r";
    string       assignment_spec("not defined"), prop_name("not specified"), prop_flag("not specified");

    cout <<"\n\nAssigning essential property values and flags to Boundaries..."<< endl;

    do {
        if ( !isCommentLine(text_line) )
        {
            // 1. name of boundary & ref to it
            string  boundaryName = strtok( text_line, delims );
            Boundary<dim>& boundary( sg.Boundary( boundaryName ) );

            // 2a. assignment whereto
            token = strtok( NULL, delims );
            if ( token == NULL )
            {
                cout <<"\n"<< text_line << endl;
                throw csmp::Exception(ERROR,
                                      "InputDataManager<dim>::ReadBoundaryPropertyValuesAndConditions",
                                      "Property assignment specifier missing");
                return false;
            }
            else
            {
                assignment_spec = token;
                std::transform(assignment_spec.begin(),assignment_spec.end(),assignment_spec.begin(),::toupper);
            }

            // 2b. flag
            token = strtok( NULL, delims );
            if ( token == NULL )
            {
                cout <<"\n"<< text_line << endl;
                throw csmp::Exception(ERROR,
                                      "InputDataManager<dim>::ReadBoundaryPropertyValuesAndConditions",
                                      "Condition type specifier missing (Dirichlet, Neumann...");
                return false;
            }
            else prop_flag = token;

            // 3a. property name
            token = strtok( NULL, delims );
            if ( token == NULL )
            {
                cout <<"\n"<< text_line << endl;
                throw csmp::Exception(ERROR,
                                      "InputDataManager<dim>::ReadBoundaryPropertyValuesAndConditions",
                                      "Property name could not be read");
                return false;
            }
            else
                prop_name = token;
            string  property(" '"); property += prop_name; property +="' ";

            // 3b. property name and type

            // 3b. property value
            VARIABLE_TYPE  prop_type = sg.Database().Type( prop_name.c_str() );
            size_t  length = sg.Database().Components( prop_name.c_str() );
            string  unit(" ["); unit += sg.Database().Unit( prop_name.c_str() ); unit +="] ";

            // overwrite protection flag
            VARIABLE_FLAG do_not_overwrite = parseStatus( prop_flag.c_str() );

            if ( prop_type == SCALAR )
            {
                ScalarVariable sc;
                sc.Flag() = parseCondition(prop_flag);
                readPropertyValue( sc );
                sg.Database().CheckRange( prop_name.c_str(), sc() );
                if      ( strcmp( "INTERIOR", assignment_spec.c_str() ) == 0 )
                    boundary.InputPropertyValue( prop_name.c_str(), sc, do_not_overwrite, INTERIOR );
                else if ( strcmp( "BOUNDARY", assignment_spec.c_str() ) == 0 )
                    boundary.InputPropertyValue( prop_name.c_str(), sc, do_not_overwrite, PERIMETER );
                else if ( strcmp( "COMPLETE", assignment_spec.c_str() ) == 0 )
                    boundary.InputPropertyValue( prop_name.c_str(), sc, do_not_overwrite, COMPLETE );
                else
                {
                    throw csmp::Exception( ERROR,
                                           "InputDataManager<dim>::ReadBoundaryPropertyValuesAndConditions",
                                           "Missing specifier (interior, boundary or complete) for", prop_name.c_str() );
                    return false;
                }

                string group(" '"); group += boundaryName; group +="' ";
                cout <<"\tinitialized"<< setw(35) << property << setw(15) << left << unit <<"in ";
                cout << setw(25) << right << group <<"to  "<< sc() << endl;
            }
            else if ( prop_type == VECTOR )
            {
                VectorVariable<dim> vc;
                for ( size_t i = 0; i < dim; ++i ) {
                    prop_flag = strtok( NULL, delims );
                    vc.Flag(i) = parseCondition(prop_flag);
                }
                readPropertyValue( vc );
                double64  val = vc.Length();
                sg.Database().CheckRange( prop_name.c_str(), val );
                if      ( strcmp( "INTERIOR", assignment_spec.c_str() ) == 0 )
                    boundary.InputPropertyValue( prop_name.c_str(), vc, do_not_overwrite, INTERIOR );
                else if ( strcmp( "BOUNDARY", assignment_spec.c_str() ) == 0 )
                    boundary.InputPropertyValue( prop_name.c_str(), vc, do_not_overwrite, PERIMETER );
                else if ( strcmp( "COMPLETE", assignment_spec.c_str() ) == 0 )
                    boundary.InputPropertyValue( prop_name.c_str(), vc, do_not_overwrite, COMPLETE );
                else
                {
                    throw csmp::Exception( ERROR,
                                           "InputDataManager<dim>::ReadBoundaryPropertyValuesAndConditions",
                                           "Missing specifier (interior, boundary or complete) for", prop_name.c_str() );
                    return false;
                }

                string group(" '"); group += boundaryName; group +="' ";
                cout <<"\tinitialized"<< setw(35) << property << setw(15) << left << unit <<"in ";
                cout << setw(25) << right << group <<"to  ";
                printPropertyValue( vc ); cout << endl;
            }
            else if ( prop_type == TENSOR )
            {
                TensorVariable<dim> ts;
                for( size_t i = 0; i < dim; ++i )
                    ts.Flag(i) = parseCondition(prop_flag);
                readPropertyValue( ts );
                double64  val = ts.MinElement();
                sg.Database().CheckRange( prop_name.c_str(), val );
                val = ts.MaxElement();
                sg.Database().CheckRange( prop_name.c_str(), val );
                if      ( strcmp( "INTERIOR", assignment_spec.c_str() ) == 0 )
                    boundary.InputPropertyValue( prop_name.c_str(), ts, do_not_overwrite, INTERIOR );
                else if ( strcmp( "BOUNDARY", assignment_spec.c_str() ) == 0 )
                    boundary.InputPropertyValue( prop_name.c_str(), ts, do_not_overwrite, PERIMETER );
                else if ( strcmp( "COMPLETE", assignment_spec.c_str() ) == 0 ) {
                    boundary.InputPropertyValue( prop_name.c_str(), ts, do_not_overwrite, COMPLETE );
                }
                else
                {
                    throw csmp::Exception( ERROR,
                                           "InputDataManager<dim>::ReadBoundaryPropertyValuesAndConditions",
                                           "Missing specifier (interior, boundary or complete) for", prop_name.c_str() );
                    return false;
                }

                string group(" '"); group += boundaryName; group +="' ";
                cout <<"\tinitialized"<< setw(35) << property << setw(15) << left << unit <<"in ";
                cout << setw(20) << right << group <<"to  ";
                printPropertyValue( ts ); cout << endl;
            }
            else if ( prop_type == ARRAY )
            {
                ArrayVariable av;
                av.Resize( length );
                av.Flag() = parseCondition(prop_flag);
                readPropertyValue( av );
                for( size_t i = 0; i< av.Size(); i++){
                    double64  val = av(i);
                    sg.Database().CheckRange( prop_name.c_str(), val );
                }
                if      ( strcmp( "INTERIOR", assignment_spec.c_str() ) == 0 )
                    boundary.InputPropertyValue( prop_name.c_str(), av, do_not_overwrite, INTERIOR );
                else if ( strcmp( "BOUNDARY", assignment_spec.c_str() ) == 0 )
                    boundary.InputPropertyValue( prop_name.c_str(), av, do_not_overwrite, PERIMETER );
                else if ( strcmp( "COMPLETE", assignment_spec.c_str() ) == 0 )
                    boundary.InputPropertyValue( prop_name.c_str(), av, do_not_overwrite, COMPLETE );
                else
                {
                    throw csmp::Exception( ERROR,
                                           "InputDataManager<dim>::ReadBoundaryPropertyValuesAndConditions",
                                           "Missing specifier (interior, boundary or complete) for", prop_name.c_str() );
                    return false;
                }

                string group(" '"); group += boundaryName; group +="' ";
                cout <<"\tinitialized"<< setw(35) << property << setw(15) << left << unit <<"in ";
                cout << setw(25) << right << group <<"to  ";
                printPropertyValue( av ); cout << endl;
            }
            else if ( prop_type == FLAGGEDARRAY )
            {
                FlaggedArrayVariable fv;
                fv.Resize( length );
                for ( size_t i = 0; i < fv.Size(); ++i ) {
                    prop_flag = strtok( NULL, delims );
                    fv.Flag(i) = parseCondition(prop_flag);
                }
                readPropertyValue( fv );
                for( size_t i = 0; i< fv.Size(); i++){
                    double64  val = fv(i);
                    sg.Database().CheckRange( prop_name.c_str(), val );
                }
                if      ( strcmp( "INTERIOR", assignment_spec.c_str() ) == 0 )
                    boundary.InputPropertyValue( prop_name.c_str(), fv, do_not_overwrite, INTERIOR );
                else if ( strcmp( "BOUNDARY", assignment_spec.c_str() ) == 0 )
                    boundary.InputPropertyValue( prop_name.c_str(), fv, do_not_overwrite, PERIMETER );
                else if ( strcmp( "COMPLETE", assignment_spec.c_str() ) == 0 )
                    boundary.InputPropertyValue( prop_name.c_str(), fv, do_not_overwrite, COMPLETE );
                else
                {
                    throw csmp::Exception( ERROR,
                                           "InputDataManager<dim>::ReadBoundaryPropertyValuesAndConditions",
                                           "Missing specifier (interior, boundary or complete) for", prop_name.c_str() );
                    return false;
                }

                string group(" '"); group += boundaryName; group +="' ";
                cout <<"\tinitialized"<< setw(35) << property << setw(15) << left << unit <<"in ";
                cout << setw(25) << right << group <<"to  ";
                printPropertyValue( fv ); cout << endl;
            }
        }
        ifs.getline( text_line, line_length );
    }
    while ( !isBlankLine(text_line) && !ifs.eof() );

    return true;

} // end readBoundaryPropertyValuesAndConditions



template bool readBoundaryPropertyValuesAndConditions<1U>(Model<1U>&,std::ifstream&,char*,size_t,bool);
template bool readBoundaryPropertyValuesAndConditions<2U>(Model<2U>&,std::ifstream&,char*,size_t,bool);
template bool readBoundaryPropertyValuesAndConditions<3U>(Model<3U>&,std::ifstream&,char*,size_t,bool);






















// READ BOUNDARY CONDITIONS


/**

reads the 'property flags assigned to boundary' (fifth config file data block)
from file and assigns the resulting property flags in the target region.
The format of any entry in the data block should have the form.

@code
# regionname	 	qualifier			input property			target flag
BOUNDARY3 	(tab)		interior	(tab) 	fluid pressure 	(tab) 	Dirichlet
@endcode

The target flag is the one that will be assigned to variables of the
target type in the boundary identified by name.

@section arguments Input Arguments

A reference to the current Model object and the line of text
that contains the boundary property flag information.

@section messages Messages

The method reports the data that are read from file to stdout and
reports any possible parsing errors.
*/
template<size_t dim>
bool readBoundaryPropertyConditions( Model<dim>& model,
                                     std::ifstream& ifs, char* text_line, size_t line_length, bool verbose )
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    char*          token(0);
    const char*    delims =":,\t,\n,\r";
    string         boundary_name, assignment_spec, prop_name;
    VARIABLE_TYPE  prop_type;
    vector<VARIABLE_FLAG>  flags;
    size_t         length;
    size_t         position;
    bool           digit( false );

    if( verbose )
        cout <<"\n\nAssigning essential property flags to Boundaries..."<< endl;

    do {
        if ( !isCommentLine(text_line) )
        {
            // 1. name of region
            boundary_name = strtok( text_line, delims );
            Boundary<dim>& boundary( model.Boundary( boundary_name ) );

            // 2. assignment whereto
            token = strtok( NULL, delims );
            if ( token == NULL )
            {
                std::cerr <<"\nFor region '"<< boundary_name <<"' text line: "<< text_line << endl;
                csmp_error.notice( ERROR,
                                  "readBoundaryPropertyConditions:",
                                  "Qualifier could not be parsed (interior, boundary or complete");
                return false;
            }
            else
            {
                assignment_spec = token;
                std::transform(assignment_spec.begin(),assignment_spec.end(),assignment_spec.begin(),::toupper);
            }

            // 3. property name
            token = strtok( NULL, delims );
            if ( token == NULL )
            {
                std::cerr <<"\nFor region '"<< boundary_name <<"' text_line: "<< text_line << endl;
                csmp_error.notice( ERROR,
                                  "readBoundaryPropertyConditions:",
                                  "Property name could not be read");
                return false;
            }
            else
                prop_name = token;

            // 3. property flag(s)
            prop_type = model.Database().Type( prop_name.c_str() );
            length    = model.Database().Components( prop_name.c_str() );

            if ( prop_type == SCALAR )
                readPropertyStatusOfScalar( 1U, digit, position, flags );
            else if ( prop_type == VECTOR )
                readPropertyStatusOfVector( dim, digit, position, flags );
            else if ( prop_type == TENSOR )
                readPropertyStatusOfTensor( dim, digit, position, flags );
            else if( prop_type == ARRAY )
                readPropertyStatusOfArray( 1U, digit, position, flags );
            else if( prop_type == FLAGGEDARRAY )
                readPropertyStatusOfFlaggedArray( length, digit, position, flags );

            if ( !digit )
            {
                if ( strcmp( "INTERIOR", assignment_spec.c_str() ) == 0 )
                    boundary.ChangePropertyStatus( prop_name.c_str(), flags, INTERIOR );

                else if ( strcmp( "BOUNDARY", assignment_spec.c_str() ) == 0 )
                    boundary.ChangePropertyStatus( prop_name.c_str(), flags, PERIMETER );

                else if ( strcmp( "COMPLETE", assignment_spec.c_str() ) == 0 )
                    boundary.ChangePropertyStatus( prop_name.c_str(), flags, COMPLETE );
                if( verbose )
                {
                    cout <<"\tassigned boundary condition to '"<< prop_name;
                    cout <<"' in '"<< assignment_spec <<"' boundary '"<< boundary_name <<"'"<< endl;
                }
            }
            else
            {
                if ( strcmp( "INTERIOR", assignment_spec.c_str() ) == 0 )
                    boundary.ChangePropertyStatus( prop_name.c_str(), position, flags[0], INTERIOR );

                else if ( strcmp( "BOUNDARY", assignment_spec.c_str() ) == 0 )
                    boundary.ChangePropertyStatus( prop_name.c_str(), position, flags[0], PERIMETER );

                else if ( strcmp( "COMPLETE", assignment_spec.c_str() ) == 0 )
                    boundary.ChangePropertyStatus( prop_name.c_str(), position, flags[0], COMPLETE );
                if( verbose )
                {
                    cout <<"\tassigned boundary condition to '"<< prop_name<<"' at position = "<<position;
                    cout <<" in '"<< assignment_spec <<"' boundary '"<< boundary_name <<"'"<< endl;
                }
            }
        }
        ifs.getline( text_line, line_length );
    }
    while ( !isBlankLine(text_line) && !ifs.eof() );

    return true;

} // end readBoundaryPropertyConditions

template bool readBoundaryPropertyConditions<1U>(Model<1U>&,std::ifstream&,char*,size_t,bool);
template bool readBoundaryPropertyConditions<2U>(Model<2U>&,std::ifstream&,char*,size_t,bool);
template bool readBoundaryPropertyConditions<3U>(Model<3U>&,std::ifstream&,char*,size_t,bool);

/**

reads the 'property flags assigned to region' (fifth config file data block)
from file and assigns the resulting property flags in the target region.
The format of any entry in the data block should have the form.

@code
# regionname	 	qualifier			input property			target flag
well 	(tab)		interior	(tab) 	fluid pressure 	(tab) 	Dirichlet
@endcode

The target flag is the one that will be assigned to variables of the
target type in the region identified by name.

@section arguments Input Arguments

A reference to the current Model object and the line of text
that contains the region property flag information.

@section messages Messages

The method reports the data that are read from file to stdout and
reports any possible parsing errors.
*/


template<size_t dim>
bool readRegionPropertyConditions( Model<dim>& model,
                                   std::ifstream& ifs, char* text_line, size_t line_length, bool verbose )
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    char*          token(0);
    const char*    delims =":,\t,\n,\r";
    string         group_name, assignment_spec, prop_name;
    VARIABLE_TYPE  prop_type;
    vector<VARIABLE_FLAG>  flags;
    size_t         length;
    size_t         position;
    bool           digit( false );

    if( verbose )
        cout <<"\n\nAssigning essential property flags to Regions..."<< endl;

    do {
        if ( !isCommentLine(text_line) )
        {
            // 1. name of region
            group_name = strtok( text_line, delims );

            // 2. assignment whereto
            token = strtok( NULL, delims );
            if ( token == NULL )
            {
                std::cerr <<"\nFor region '"<< group_name <<"' text line: "<< text_line << endl;
                csmp_error.notice( ERROR,
                                   "readRegionPropertyConditions:",
                                  "Qualifier could not be parsed (interior, boundary or complete");
                return false;
            }
            else
            {
                assignment_spec = token;
                std::transform(assignment_spec.begin(),assignment_spec.end(),assignment_spec.begin(),::toupper);
            }

            // 3. property name
            token = strtok( NULL, delims );
            if ( token == NULL )
            {
                std::cerr <<"\nFor region '"<< group_name <<"' text_line: "<< text_line << endl;
                csmp_error.notice( ERROR,
                                   "readRegionPropertyConditions:",
                                  "Property name could not be read");
                return false;
            }
            else prop_name = token;

            // 3. property flag(s)
            prop_type = model.Database().Type( prop_name.c_str() );
            length    = model.Database().Components( prop_name.c_str() );

            if ( prop_type == SCALAR )
                readPropertyStatusOfScalar( 1U, digit, position, flags );
            else if ( prop_type == VECTOR )
                readPropertyStatusOfVector( dim, digit, position, flags );
            else if ( prop_type == TENSOR )
                readPropertyStatusOfTensor( dim, digit, position, flags );
            else if( prop_type == ARRAY )
                readPropertyStatusOfArray( 1U, digit, position, flags );
            else if( prop_type == FLAGGEDARRAY )
                readPropertyStatusOfFlaggedArray( length, digit, position, flags );

            if ( !digit )
            {
                if ( strcmp( "INTERIOR", assignment_spec.c_str() ) == 0 )
                    model.Region( group_name.c_str() ).ChangePropertyStatus( prop_name.c_str(), flags, INTERIOR );

                else if ( strcmp( "BOUNDARY", assignment_spec.c_str() ) == 0 )
                    model.Region( group_name.c_str() ).ChangePropertyStatus( prop_name.c_str(), flags, PERIMETER );

                else if ( strcmp( "COMPLETE", assignment_spec.c_str() ) == 0 )
                    model.Region( group_name.c_str() ).ChangePropertyStatus( prop_name.c_str(), flags, COMPLETE );

                if( verbose )
                {
                    if ( prop_type == SCALAR ) cout <<"\tassigned "<< parseStatus(flags[0]) <<" condition to '"<< prop_name;
                    else cout <<"\tassigned VARIABLE_FLAG condition(s) to '"<< prop_name;
                    cout <<"' in '"<< assignment_spec <<"' region '"<< group_name <<"'"<< endl;
                }
            }
            else
            {
                if ( strcmp( "INTERIOR", assignment_spec.c_str() ) == 0 )
                    model.Region( group_name.c_str() ).ChangePropertyStatus( prop_name.c_str(), position, flags[0], INTERIOR );

                else if ( strcmp( "BOUNDARY", assignment_spec.c_str() ) == 0 )
                    model.Region( group_name.c_str() ).ChangePropertyStatus( prop_name.c_str(), position, flags[0], PERIMETER );

                else if ( strcmp( "COMPLETE", assignment_spec.c_str() ) == 0 )
                    model.Region( group_name.c_str() ).ChangePropertyStatus( prop_name.c_str(), position, flags[0], COMPLETE );

                if( verbose )
                {
                    if ( prop_type == SCALAR ) cout <<"\tassigned "<< parseStatus(flags[0]) <<" condition to '"<< prop_name<<"' at position = "<<position;
                    else cout <<"\tassigned VARIABLE_FLAG  condition(s) to '"<< prop_name<<"' at position = "<<position;
                    cout <<" in '"<< assignment_spec <<"' region '"<< group_name <<"'"<< endl;
                }
            }
        }
        ifs.getline( text_line, line_length );
    }
    while ( !isBlankLine(text_line) && !ifs.eof() );

    return true;

} // end readRegionPropertyConditions


template bool readRegionPropertyConditions<1U>(Model<1U>&,std::ifstream&,char*,size_t,bool);
template bool readRegionPropertyConditions<2U>(Model<2U>&,std::ifstream&,char*,size_t,bool);
template bool readRegionPropertyConditions<3U>(Model<3U>&,std::ifstream&,char*,size_t,bool);








// READ PROPERTIES



/** readBoundaryPropertyValues

Example
@code
# a fluid pressure of 200 bar is assigned to the interior of a boundary named
# 'BOUNDARY3'
BOUNDARY3	(tab)	interior	(tab)	fluid pressure	(tab)	2.0E+7
@endcode

@attention This applies only to unique boundary names.
*/
template<size_t dim>
bool readBoundaryPropertyValues( Model<dim>& model,
                                 std::ifstream& ifs, char* text_line, size_t line_length, bool verbose )
{
    char*        token(0);
    const char*  delims =":,\t,\n,\r";
    string       assignment_spec("not defined"), prop_name("not specified"), prop_flag("not specified");

    if( verbose )
        cout <<"\n\nAssigning essential property values to Boundaries..."<< endl;

    do {
        if ( !isCommentLine(text_line) )
        {
            // 1. name of boundary & ref to it
            string  boundaryName = strtok( text_line, delims );
            Boundary<dim>& boundary( model.Boundary( boundaryName ) );

            // 2a. assignment whereto
            token = strtok( NULL, delims );
            if ( token == NULL )
            {
                std::cerr <<"\n"<< text_line << endl;
                throw csmp::Exception( ERROR,
                                       "readBoundaryPropertyValues",
                                       "Property assignment specifier missing");
                return false;
            }
            else
            {
                assignment_spec = token;
                std::transform(assignment_spec.begin(),assignment_spec.end(),assignment_spec.begin(),::toupper);
            }

            // 3a. property name
            token = strtok( NULL, delims );
            if ( token == NULL )
            {
                std::cerr <<"\n"<< text_line << endl;
                throw csmp::Exception( ERROR,
                                       "readBoundaryPropertyValues",
                                       "Property name could not be read");
                return false;
            }
            else
                prop_name = token;
            string  property(" '"); property += prop_name; property +="' ";

            // 3b. property value
            VARIABLE_TYPE  prop_type = model.Database().Type( prop_name.c_str() );
            size_t  length = model.Database().Components( prop_name.c_str() );
            string  unit(" ["); unit += model.Database().Unit( prop_name.c_str() ); unit +="] ";

            if ( prop_type == SCALAR )
            {
                ScalarVariable sc;
                readPropertyValue( sc );
                model.Database().CheckRange( prop_name.c_str(), sc() );
                if      ( strcmp( "INTERIOR", assignment_spec.c_str() ) == 0 )
                    boundary.InputPropertyValue( prop_name.c_str(), sc, INTERIOR );
                else if ( strcmp( "BOUNDARY", assignment_spec.c_str() ) == 0 )
                    boundary.InputPropertyValue( prop_name.c_str(), sc, PERIMETER );
                else if ( strcmp( "COMPLETE", assignment_spec.c_str() ) == 0 )
                    boundary.InputPropertyValue( prop_name.c_str(), sc, COMPLETE );
                else
                {
                    throw csmp::Exception( ERROR,
                                           "readBoundaryPropertyValues",
                                           "Missing specifier (interior, boundary or complete) for", prop_name.c_str() );
                    return false;
                }
                if( verbose )
                {
                    string group(" '"); group += boundaryName; group +="' ";
                    cout <<"\tinitialized"<< setw(35) << property << setw(15) << left << unit <<"in ";
                    cout << setw(25) << right << group <<"to  "<< sc() << endl;
                }
            }
            else if ( prop_type == VECTOR )
            {
                VectorVariable<dim> vc;
                readPropertyValue( vc );
                double64  val = vc.Length();
                model.Database().CheckRange( prop_name.c_str(), val );
                if      ( strcmp( "INTERIOR", assignment_spec.c_str() ) == 0 )
                    boundary.InputPropertyValue( prop_name.c_str(), vc, INTERIOR );
                else if ( strcmp( "BOUNDARY", assignment_spec.c_str() ) == 0 )
                    boundary.InputPropertyValue( prop_name.c_str(), vc, PERIMETER );
                else if ( strcmp( "COMPLETE", assignment_spec.c_str() ) == 0 )
                    boundary.InputPropertyValue( prop_name.c_str(), vc, COMPLETE );
                else
                {
                    throw csmp::Exception( ERROR,
                                           "readBoundaryPropertyValues",
                                           "Missing specifier (interior, boundary or complete) for", prop_name.c_str() );
                    return false;
                }
                if( verbose )
                {
                    string group(" '"); group += boundaryName; group +="' ";
                    cout <<"\tinitialized"<< setw(35) << property << setw(15) << left << unit <<"in ";
                    cout << setw(25) << right << group <<"to  ";
                    printPropertyValue( vc ); cout << endl;
                }
            }
            else if ( prop_type == TENSOR )
            {
                TensorVariable<dim> ts;
                readPropertyValue( ts );
                double64  val = ts.MinElement();
                model.Database().CheckRange( prop_name.c_str(), val );
                val = ts.MaxElement();
                model.Database().CheckRange( prop_name.c_str(), val );
                if      ( strcmp( "INTERIOR", assignment_spec.c_str() ) == 0 )
                    boundary.InputPropertyValue( prop_name.c_str(), ts, INTERIOR );
                else if ( strcmp( "BOUNDARY", assignment_spec.c_str() ) == 0 )
                    boundary.InputPropertyValue( prop_name.c_str(), ts, PERIMETER );
                else if ( strcmp( "COMPLETE", assignment_spec.c_str() ) == 0 ) {
                    boundary.InputPropertyValue( prop_name.c_str(), ts, COMPLETE );
                }
                else
                {
                    throw csmp::Exception( ERROR,
                                           "readBoundaryPropertyValues",
                                           "Missing specifier (interior, boundary or complete) for", prop_name.c_str() );
                    return false;
                }
                if( verbose )
                {
                    string group(" '"); group += boundaryName; group +="' ";
                    cout <<"\tinitialized"<< setw(35) << property << setw(15) << left << unit <<"in ";
                    cout << setw(20) << right << group <<"to  ";
                    printPropertyValue( ts ); cout << endl;
                }
            }
            else if ( prop_type == ARRAY )
            {
                ArrayVariable av;
                av.Resize( length );
                readPropertyValue( av );
                for( size_t i = 0; i< av.Size(); i++){
                    double64  val = av(i);
                    model.Database().CheckRange( prop_name.c_str(), val );
                }
                if      ( strcmp( "INTERIOR", assignment_spec.c_str() ) == 0 )
                    boundary.InputPropertyValue( prop_name.c_str(), av, INTERIOR );
                else if ( strcmp( "BOUNDARY", assignment_spec.c_str() ) == 0 )
                    boundary.InputPropertyValue( prop_name.c_str(), av, PERIMETER );
                else if ( strcmp( "COMPLETE", assignment_spec.c_str() ) == 0 )
                    boundary.InputPropertyValue( prop_name.c_str(), av, COMPLETE );
                else
                {
                    throw csmp::Exception( ERROR,
                                           "readBoundaryPropertyValues",
                                           "Missing specifier (interior, boundary or complete) for", prop_name.c_str() );
                    return false;
                }
                if( verbose )
                {
                    string group(" '"); group += boundaryName; group +="' ";
                    cout <<"\tinitialized"<< setw(35) << property << setw(15) << left << unit <<"in ";
                    cout << setw(25) << right << group <<"to  ";
                    printPropertyValue( av ); cout << endl;
                }
            }
            else if ( prop_type == FLAGGEDARRAY )
            {
                FlaggedArrayVariable fv;
                fv.Resize( length );
                readPropertyValue( fv );
                for( size_t i = 0; i< fv.Size(); i++){
                    double64  val = fv(i);
                    model.Database().CheckRange( prop_name.c_str(), val );
                }
                if      ( strcmp( "INTERIOR", assignment_spec.c_str() ) == 0 )
                    boundary.InputPropertyValue( prop_name.c_str(), fv, INTERIOR );
                else if ( strcmp( "BOUNDARY", assignment_spec.c_str() ) == 0 )
                    boundary.InputPropertyValue( prop_name.c_str(), fv, PERIMETER );
                else if ( strcmp( "COMPLETE", assignment_spec.c_str() ) == 0 )
                    boundary.InputPropertyValue( prop_name.c_str(), fv, COMPLETE );
                else
                {
                    throw csmp::Exception( ERROR,
                                           "readBoundaryPropertyValues",
                                           "Missing specifier (interior, boundary or complete) for", prop_name.c_str() );
                    return false;
                }
                if( verbose )
                {
                    string group(" '"); group += boundaryName; group +="' ";
                    cout <<"\tinitialized"<< setw(35) << property << setw(15) << left << unit <<"in ";
                    cout << setw(25) << right << group <<"to  ";
                    printPropertyValue( fv ); cout << endl;
                }
            }
        }
        ifs.getline( text_line, line_length );
    }
    while ( !isBlankLine(text_line) && !ifs.eof() );

    return true;

} // end readBoundaryPropertyValues

template bool readBoundaryPropertyValues<1U>(Model<1U>&,std::ifstream&,char*,size_t,bool);
template bool readBoundaryPropertyValues<2U>(Model<2U>&,std::ifstream&,char*,size_t,bool);
template bool readBoundaryPropertyValues<3U>(Model<3U>&,std::ifstream&,char*,size_t,bool);

/**

This method reads data of the type Block 3, see class documentation below.
It uses the input string to infer from it which variable shall
be assigned the user-specified value.

@code
# a porosity of 25% is assigned to the interior of a region named
# 'fault zone'
fault zone	tab  	interior 	tab		porosity	tab   0.25
@endcode

Dependent on whether the property is a scalar, vector, or tensor variable
a different number of values are expected as input. This number also
depends on the spatial dimension of the model (a vector has 2 components
in 2D, but 3 in 3D, a 2D tensor is 2 x 2 and a 3D one 3 x 3 in size).

@section arguments Input Arguments

The method takes a reference to the Model object which contains the
region in which the variable values shall be modified.

The regionname, type of assignment (interior, boundary, or complete),
variable name and property value(s) are extracted from the supplied
character string.

@section messages Messages

The method attempts to correctly identify parsing errors that may
occur when the textline is interpreted.
*/
template<size_t dim>
bool readRegionPropertyValues( Model<dim>& model,
                               std::ifstream& ifs, char* text_line, size_t line_length, bool verbose )
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    char*        token(0);
    const char*  delims =":,\t,\n,\r";
    string   assignment_spec("not defined"), prop_name("not specified");

    if( verbose )
        cout <<"\n\nAssigning initial property values to Regions..."<< endl << endl;

    do {
        if ( !isCommentLine(text_line) )
        {
            // 1. name of region
            string  group_name = strtok( text_line, delims );

            // 2. assignment whereto
            token = strtok( NULL, delims );
            if ( token == NULL )
            {
                std::cerr <<"\n"<< text_line << endl;
                throw csmp::Exception( ERROR,
                                       "readRegionPropertyValues",
                                       "Property assignment specifier missing");
                return false;
            }
            else
            {
                assignment_spec = token;
                std::transform(assignment_spec.begin(),assignment_spec.end(),assignment_spec.begin(),::toupper);
            }

            // 3. property name
            token = strtok( NULL, delims );
            if ( token == NULL )
            {
                std::cerr <<"\n"<< text_line << endl;
                throw csmp::Exception( ERROR,
                                       "readRegionPropertyValues",
                                       "Property name could not be read");
                return false;
            }
            else
                prop_name = token;
            string  property(" '"); property += prop_name; property +="' ";

            // 3. property value
            VARIABLE_TYPE  prop_type = model.Database().Type( prop_name.c_str() );
            size_t  length = model.Database().Components( prop_name.c_str() );
            string  unit(" ["); unit += model.Database().Unit( prop_name.c_str() ); unit +="] ";

            if ( prop_type == SCALAR ) {
                ScalarVariable sc;
                readPropertyValue( sc );
                model.Database().CheckRange( prop_name.c_str(), sc() );
                if      ( strcmp( "INTERIOR", assignment_spec.c_str() ) == 0 )
                    model.Region(group_name.c_str()).InputPropertyValue( prop_name.c_str(), sc, INTERIOR );
                else if ( strcmp( "BOUNDARY", assignment_spec.c_str() ) == 0 )
                    model.Region(group_name.c_str()).InputPropertyValue( prop_name.c_str(), sc, PERIMETER );
                else if ( strcmp( "COMPLETE", assignment_spec.c_str() ) == 0 ) {
                    model.Region(group_name.c_str()).InputPropertyValue( prop_name.c_str(), sc, COMPLETE );
                }
                else
                {
                    throw csmp::Exception( ERROR,
                                           "readRegionPropertyValues",
                                           "Missing specifier (interior, boundary or complete) for", prop_name.c_str() );
                    return false;
                }
                if( verbose )
                {
                    string group(" '"); group += group_name; group +="' ";
                    cout <<"\tinitialized"<< setw(35) << property << setw(15) << left << unit <<"in ";
                    cout << setw(25) << right << group <<"to  "<< sc() << endl;
                }
            }
            else if ( prop_type == VECTOR ) {
                VectorVariable<dim> vc;
                readPropertyValue( vc );
                double64  val = vc.Length();
                model.Database().CheckRange( prop_name.c_str(), val );
                if      ( strcmp( "INTERIOR", assignment_spec.c_str() ) == 0 )
                    model.Region(group_name.c_str()).InputPropertyValue( prop_name.c_str(), vc, INTERIOR );
                else if ( strcmp( "BOUNDARY", assignment_spec.c_str() ) == 0 )
                    model.Region(group_name.c_str()).InputPropertyValue( prop_name.c_str(), vc, PERIMETER );
                else if ( strcmp( "COMPLETE", assignment_spec.c_str() ) == 0 ) {
                    model.Region(group_name.c_str()).InputPropertyValue( prop_name.c_str(), vc, COMPLETE );
                }
                else
                {
                    throw csmp::Exception( ERROR,
                                           "readRegionPropertyValues",
                                           "Missing specifier (interior, boundary or complete) for", prop_name.c_str() );
                    return false;
                }
                if( verbose )
                {
                    string group(" '"); group += group_name; group +="' ";
                    cout <<"\tinitialized"<< setw(35) << property << setw(15) << left << unit <<"in ";
                    cout << setw(25) << right << group <<"to  ";
                    printPropertyValue( vc ); cout << endl;
                }
            }
            else if ( prop_type == TENSOR ) {
              TensorVariable<dim> ts;
              readPropertyValue( ts );

              VectorVariable<dim> eigVals;
              TensorVariable<dim> eigVecs;

              // Check that minimum and maximum eigenvalues are in range
              if (!ts.EigenNonSymmetric( eigVals, eigVecs )) {
                csmp_error.notice(FATAL_ERROR,"readRegionPropertyValues",
                        "Cannot eigendecompose the tensor for property", prop_name.c_str());
              }
              for ( size_t i = 0; i < dim; ++i ) {
                double64  val = eigVals(i);
                model.Database().CheckRange( prop_name.c_str(), val );
              }
              if      ( strcmp( "INTERIOR", assignment_spec.c_str() ) == 0 )
                  model.Region(group_name.c_str()).InputPropertyValue( prop_name.c_str(), ts, INTERIOR );
              else if ( strcmp( "BOUNDARY", assignment_spec.c_str() ) == 0 )
                  model.Region(group_name.c_str()).InputPropertyValue( prop_name.c_str(), ts, PERIMETER );
              else if ( strcmp( "COMPLETE", assignment_spec.c_str() ) == 0 ) {
                  model.Region(group_name.c_str()).InputPropertyValue( prop_name.c_str(), ts, COMPLETE );
              }
              else
              {
                  throw csmp::Exception( ERROR,
                                         "readRegionPropertyValues",
                                         "Missing specifier (interior, boundary or complete) for", prop_name.c_str() );
                  return false;
              }
              if( verbose )
              {
                  string group(" '"); group += group_name; group +="' ";
                  cout <<"\tinitialized"<< setw(35) << property << setw(15) << left << unit <<"in ";
                  cout << setw(20) << right << group <<"to  ";
                  printPropertyValue( ts ); cout << endl;
              }
            }
            else if ( prop_type == ARRAY ) {
                ArrayVariable av;
                av.Resize( length );
                readPropertyValue( av );
                for( size_t i = 0; i< av.Size(); i++){
                    double64  val = av(i);
                    model.Database().CheckRange( prop_name.c_str(), val );
                }
                if      ( strcmp( "INTERIOR", assignment_spec.c_str() ) == 0 )
                    model.Region(group_name.c_str()).InputPropertyValue( prop_name.c_str(), av, INTERIOR );
                else if ( strcmp( "BOUNDARY", assignment_spec.c_str() ) == 0 )
                    model.Region(group_name.c_str()).InputPropertyValue( prop_name.c_str(), av, PERIMETER );
                else if ( strcmp( "COMPLETE", assignment_spec.c_str() ) == 0 ) {
                    model.Region(group_name.c_str()).InputPropertyValue( prop_name.c_str(), av, COMPLETE );
                }
                else
                {
                    throw csmp::Exception( ERROR,
                                           "readRegionPropertyValues",
                                           "Missing specifier (interior, boundary or complete) for", prop_name.c_str() );
                    return false;
                }
                if( verbose )
                {
                    string group(" '"); group += group_name; group +="' ";
                    cout <<"\tinitialized"<< setw(35) << property << setw(15) << left << unit <<"in ";
                    cout << setw(25) << right << group <<"to  ";
                    printPropertyValue( av ); cout << endl;
                }
            }
            else if ( prop_type == FLAGGEDARRAY ) {
                FlaggedArrayVariable fv;
                fv.Resize( length );
                readPropertyValue( fv );
                for( size_t i = 0; i< fv.Size(); i++){
                    double64  val = fv(i);
                    model.Database().CheckRange( prop_name.c_str(), val );
                }
                if      ( strcmp( "INTERIOR", assignment_spec.c_str() ) == 0 )
                    model.Region(group_name.c_str()).InputPropertyValue( prop_name.c_str(), fv, INTERIOR );
                else if ( strcmp( "BOUNDARY", assignment_spec.c_str() ) == 0 )
                    model.Region(group_name.c_str()).InputPropertyValue( prop_name.c_str(), fv, PERIMETER );
                else if ( strcmp( "COMPLETE", assignment_spec.c_str() ) == 0 ) {
                    model.Region(group_name.c_str()).InputPropertyValue( prop_name.c_str(), fv, COMPLETE );
                }
                else
                {
                    throw csmp::Exception( ERROR,
                                           "readRegionPropertyValues",
                                           "Missing specifier (interior, boundary or complete) for", prop_name.c_str() );
                    return false;
                }
                if( verbose )
                {
                    string group(" '"); group += group_name; group +="' ";
                    cout <<"\tinitialized"<< setw(35) << property << setw(15) << left << unit <<"in ";
                    cout << setw(25) << right << group <<"to  ";
                    printPropertyValue( fv ); cout << endl;
                }
            }
        }
        ifs.getline( text_line, line_length );
    }
    while ( !isBlankLine(text_line) && !ifs.eof() );

    return true;

} // end readRegionPropertyValues


template bool readRegionPropertyValues<1U>(Model<1U>&,std::ifstream&,char*,size_t,bool);
template bool readRegionPropertyValues<2U>(Model<2U>&,std::ifstream&,char*,size_t,bool);
template bool readRegionPropertyValues<3U>(Model<3U>&,std::ifstream&,char*,size_t,bool);






/**

This method reads data of the type Block 2, see class documentation below.
It assigns default values to the entire model. These values are interpreted
from the supplied input string that should have the form:

@code
property name  tab  value(s)
@endcode

The input property value is compared with the physically meaningful
range that the user specified in the '*-variables.txt' file. If this
range is exceeded, the value will be corrected to the nearest of the
two extrema.

Dependent on whether the property is a scalar, vector, or tensor variable
a different number of values are expected as input. This number also
depends on the spatial dimension of the model (a vector has 2 components
in 2D, but 3 in 3D, a 2D tensor is 2 x 2 and a 3D one 3 x 3 in size).

@section arguments Input Arguments

The method assigns the correctly parsed property value to its first
argument - a Model object. The second argument is character string
that should contain the tab-separated variable name and value(s).

@section messages Messages

The method reports potential parsing errors and echoes the values that
it reads to the screen.
 */
template<size_t dim>
bool readDefaultPropertyValues( Model<dim>& model,
                                std::ifstream& ifs, char* text_line, size_t line_length, bool verbose )
{
    ErrorHandler& error_handler ( ErrorHandler::Instance() );

    const char* delims =":,\t,\n,\r";

    if( verbose )
        cout <<"\n\nAssigning initial property values..."<< endl << endl;

    do {
        if ( !isCommentLine(text_line) )
        {
            string         prop_name = strtok( text_line, delims );
            string         property(" '"); property += prop_name; property +="' ";
            VARIABLE_TYPE  prop_type = model.Database().Type( prop_name.c_str() );
            size_t         length = model.Database().Components( prop_name.c_str() );
            string         unit(" ["); unit += model.Database().Unit( prop_name.c_str() ); unit +="] ";

            cout.setf( ios::scientific, ios::floatfield );
            cout.precision(3);

            if ( prop_type == SCALAR ) {
                ScalarVariable sc;
                readPropertyValue( sc );
                model.Database().CheckRange( prop_name.c_str(), sc() );
                model.InputPropertyValue( prop_name.c_str(), sc );
                if( verbose )
                    cout <<"\tinitialized"<< setw(35) << property << setw(15) << left << unit <<"to  "<< right << sc() << endl;
            }
            else if ( prop_type == VECTOR ) {
                VectorVariable<dim> vc;
                readPropertyValue( vc );
                double64  val = vc.Length();
                model.Database().CheckRange( prop_name.c_str(), val );
                model.InputPropertyValue( prop_name.c_str(), vc );
                if( verbose )
                {
                    cout <<"\tinitialized"<< setw(35) << property << setw(15) << left << unit <<"to  "<< right;
                    printPropertyValue( vc ); cout << endl;
                }
            }
            else if ( prop_type == TENSOR ) {
              TensorVariable<dim> ts;
              readPropertyValue( ts );
              
              VectorVariable<dim> eigVals;
              TensorVariable<dim> eigVecs;

              // Check that minimum and maximum eigenvalues are in range
              if (!ts.EigenNonSymmetric( eigVals, eigVecs )) {
                error_handler.notice(FATAL_ERROR,"readDefaultPropertyValues",
                        "Cannot eigendecompose the tensor for property", prop_name.c_str());
              }
              for ( size_t i = 0; i < dim; ++i ) {
                double64  val = eigVals(i);
                model.Database().CheckRange( prop_name.c_str(), val );
              }
              model.InputPropertyValue( prop_name.c_str(), ts );
              if( verbose )
              {
                cout <<"\tinitialized"<< setw(35) << property << setw(15) << left << unit <<"to  "<< right << endl;
                printPropertyValue( ts );
              }
            }
            else if ( prop_type == ARRAY ) {
                ArrayVariable av;
                av.Resize( length );
                readPropertyValue( av );
                for( size_t i = 0; i< av.Size(); i++){
                    double64  val = av(i);
                    model.Database().CheckRange( prop_name.c_str(), val );
                }
                model.InputPropertyValue( prop_name.c_str(), av );
                if( verbose )
                {
                    cout <<"\tinitialized"<< setw(35) << property << setw(15) << left << unit <<"to  "<< right;
                    printPropertyValue( av ); cout << endl;
                }
            }
            else if ( prop_type == FLAGGEDARRAY ) {
                FlaggedArrayVariable fv;
                fv.Resize( length );
                readPropertyValue( fv );
                for( size_t i = 0; i< fv.Size(); i++){
                    double64  val = fv(i);
                    model.Database().CheckRange( prop_name.c_str(), val );
                }
                model.InputPropertyValue( prop_name.c_str(), fv );
                if( verbose )
                {
                    cout <<"\tinitialized"<< setw(35) << property << setw(15) << left << unit <<"to  "<< right;
                    printPropertyValue( fv ); cout << endl;
                }
            }
        }
        ifs.getline( text_line, line_length );
    }
    while ( !isBlankLine(text_line) && !ifs.eof() );

    return true;

} // end readDefaultPropertyValues



template bool readDefaultPropertyValues<1U>(Model<1U>&,std::ifstream&,char*,size_t,bool);
template bool readDefaultPropertyValues<2U>(Model<2U>&,std::ifstream&,char*,size_t,bool);
template bool readDefaultPropertyValues<3U>(Model<3U>&,std::ifstream&,char*,size_t,bool);








// READ COMPUTATIONAL SETTINGS


/**

This method permits to initialize a computational settings object from
a file stream. This can be done by the data input manager so that
the runtime configuration can be read from the same input file as
the material properties.

This type of data effectively constitutes an additional data block
(Block 6).

@section arguments Input Arguments

The first argument must be a file stream that points to the beginning
of a computational settings data block of the form:

 @code
# (6) computational settings
time stepping	pedantic
duration		86400.
output time		10000.
output time		43200.
 @endcode

The third method argument represents the first line of the computational
settings datablock.

@return The method returns true if the data was read successfully.

@section messages Messages

Any reading problems are reported.
*/
template<size_t dim>
bool readComputationalSettings( ComputationalSettings& settings,
                                std::ifstream& ifs, char* text_line, size_t line_length, bool verbose )
{
    char*       token(0);
    const char* delims =":,\t,\n,\r";

    if( verbose )
        cout <<"\n\nreading computational settings..."<< endl << endl;

    double64 time_unit        ( 1.0 );
    double64 output_time_unit ( 1.0 );
    double64 monitor_time_unit( 1.0 );

    do {
        if ( !isCommentLine(text_line) )
        {
            token = strtok( text_line, delims );

            // time stepping strategy
            if ( strcmp( "time stepping", token ) == 0 ) {
                token = strtok( NULL, delims );
                settings.TimeSteppingApproach( parseTimeStrategy(token) );
            }

            // time unit
            else if ( strcmp( "time unit", token ) == 0 ) {
                token = strtok( NULL, delims );
                time_unit = atof(token);
            }

            // time increment
            else if ( strcmp( "time increment", token ) == 0 ) {
                token = strtok( NULL, delims );
                settings.TimeIncrement( atof(token) * time_unit );
            }

            // duration
            else if ( strcmp( "duration", token ) == 0 ) {
                token = strtok( NULL, delims );
                settings.Duration( atof(token) * time_unit );
            }

            // output time unit
            else if ( strcmp( "output time unit", token ) == 0 ) {
                token = strtok( NULL, delims );
                output_time_unit = atof(token);
            }

            // output times
            else if ( strcmp( "output time", token ) == 0 ) {
                token = strtok( NULL, delims );
                settings.AddOutputTime( atof(token) * output_time_unit );
            }

            // monitor time unit
            else if ( strcmp( "monitor time unit", token ) == 0 ) {
                token = strtok( NULL, delims );
                monitor_time_unit =  atof(token);
            }

            // monitor times
            else if ( strcmp( "monitor time", token ) == 0 ) {
                token = strtok( NULL, delims );
                settings.AddMonitorTime( atof(token) * monitor_time_unit );
            }
            else
                throw csmp::Exception( ERROR, "readComputationalSettings",
                                       "Unable to parse settings qualifier.");
        }
        ifs.getline( text_line, line_length );
    }
    while ( !isBlankLine(text_line) && !ifs.eof() );

    if( verbose )
        settings.Out();

    return true;

} // end readComputationalSettings

template bool readComputationalSettings<1U>(ComputationalSettings&,std::ifstream&,char*,size_t,bool);
template bool readComputationalSettings<2U>(ComputationalSettings&,std::ifstream&,char*,size_t,bool);
template bool readComputationalSettings<3U>(ComputationalSettings&,std::ifstream&,char*,size_t,bool);


/// skip unknown or invalid block of data
bool readUnknownOrInvalidBlock( const std::string& keyword,
                                std::ifstream& ifs, char* text_line, size_t line_length, bool verbose )
{
    if( verbose )
        cout <<"\nreadUnknownBlock: reading "<< keyword << " block of data...";

    do{
        ifs.getline( text_line, line_length );
    }
    while ( !isBlankLine( text_line ) && !ifs.eof() );

    if( verbose )
        cout <<"\nreadUnknownBlock: block "<< keyword << " have been successfully skiped!\n";
    return true;
}


} // end csmp

