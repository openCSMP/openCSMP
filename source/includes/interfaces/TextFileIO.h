#ifndef TEXT_FILE_IO_H
#define TEXT_FILE_IO_H

#include "Exception.h"
#include "ErrorHandler.h"

#include "Node.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"

#include "ScalarVariable.h"
#include "VectorVariable.h"
#include "TensorVariable.h"
#include "ArrayVariable.h"
#include "FlaggedArrayVariable.h"

#include "Model.h"

#include "ComputationalSettings.h"

namespace csmp {

/**

Collection of highlevel functions for seeking and reading in textfiles.

Support for comments, which must be preceded by the hash tag.

@note single-hash tags will not be parsed correctly.

@author S.K. Matthai
@author S. Geiger
@author S.G. Roberts
@author R. Manasipov
@date 2001,2012,2014

*/


// Common text file interface functionality

// Check existance of file
bool doesFileExist( const std::string& );
bool checkExistance( const char* filename );

// Open file
template<class FileStream>
bool openFile( FileStream& fs, const std::string& fname );
template<class FileStream>
bool openFile( FileStream& fs, const std::string& fname, const std::string& file_extension );
template<class FileStream>
bool openFile( FileStream& fs, const std::string& fname, const std::vector<std::string>& file_extensions );

// Reading specific parts of file
bool readFileHeader( std::ifstream& ifs, std::string& header, bool verbose = false );
bool readLineTellIfBlank( std::ifstream& ifs, char* text_line, size_t line_length );
bool readNonBlankLine( std::ifstream& ifs, char* text_line, size_t line_length );
int  readFirstLineInBlock( std::ifstream& ifs, char* text_line, size_t line_length );
int  readIncludeFileName( std::ifstream& ifs, char* text_line, size_t line_length,
                          std::string& fname );
bool readKeyword( std::ifstream& ifs, char* text_line, size_t line_length,
                  std::string& keyword, std::vector<std::string>& keyword_parameters );
bool readKeywordsAndParameters( std::ifstream& ifs, char* text_line, size_t line_length, bool verbose,
                                std::vector<std::pair<std::string, std::vector<std::string> > >& );
bool isBlankLine( const char* str );
/// detects lines beginning with #, %, --; char* not constant because first comment char gets replaced by '\0'
bool isCommentLine( char* str );
void advancePastCommentLine( std::ifstream& ifs );
int  yes_or_no( const std::string& );
bool isYES( const std::string& );
bool isNO( const std::string& );
bool isRealNumber( const std::string&, int len=2 );
bool isRealNumber( const char*, int len=2 );
bool isIntegerNumber( const std::string&, int len=2 );
bool isIntegerNumber( const char*, int len=2 );

/// tokenising a string with a single chosen control character; default = blank
std::vector<std::string> split( const char* , char c = ' ' );

/// tokenising a string with a set of chosen control characters using C++ regular expressions
std::vector<std::string> tokenise( std::string, const std::string regular_expression="[#%^\t\r\n]" );

/// removes whitespace from string, returning the remaining character sequence
std::string  withoutSpaces( std::string );


// Reading property values

/// parses scalars including "nodata" and NaN values
double parseDataValue( const std::string& value );

void readPropertyValue( ScalarVariable& sc );

template<size_t dim>
void readPropertyValue( VectorVariable<dim>& vc );

template<size_t dim>
void readPropertyValue( TensorVariable<dim>& ts );

void readPropertyValue( ArrayVariable& av );

void readPropertyValue( FlaggedArrayVariable& fv );

// Reading property status

void readPropertyStatusOfScalar( size_t depth, bool& digit, size_t& position, std::vector<VARIABLE_FLAG>& flags );

void readPropertyStatusOfVector( size_t depth, bool& digit, size_t& position, std::vector<VARIABLE_FLAG>& flags );

void readPropertyStatusOfTensor( size_t depth, bool& digit, size_t& position, std::vector<VARIABLE_FLAG>& flags );

void readPropertyStatusOfArray( size_t depth, bool& digit, size_t& position, std::vector<VARIABLE_FLAG>& flags );

void readPropertyStatusOfFlaggedArray( size_t depth, bool& digit, size_t& position, std::vector<VARIABLE_FLAG>& flags );


// Printing property values

void printPropertyValue( const ScalarVariable& sc );

template<size_t dim>
void printPropertyValue( const VectorVariable<dim>& vc );

template<size_t dim>
void printPropertyValue( const TensorVariable<dim>& ts );

void printPropertyValue( const ArrayVariable& av );

void printPropertyValue( const FlaggedArrayVariable& fv );


// Reporting about read data

void reportAssignment( const std::string& prop, const std::string& unit,
                       const std::string& bound,
                       ScalarVariable& sc1, ScalarVariable& sc2 );

void reportAssignment( const std::string& prop, const std::string& unit,
                       const std::string& bound,
                       ScalarVariable& sc1, ScalarVariable& sc2,
                       ScalarVariable& sc3, ScalarVariable& sc4 );

template<size_t dim>
void reportAssignment( const std::string& prop, const std::string& unit,
                       const std::string& bound,
                       VectorVariable<dim>& vc1, VectorVariable<dim>& vc2 );

template<size_t dim>
void reportAssignment( const std::string& prop, const std::string& unit,
                       const std::string& bound,
                       VectorVariable<dim>& vc1, VectorVariable<dim>& vc2,
                       VectorVariable<dim>& vc3, VectorVariable<dim>& vc4 );

template<size_t dim>
void reportAssignment( const std::string& prop, const std::string& unit,
                       const std::string& bound,
                       TensorVariable<dim>& ts1, TensorVariable<dim>& ts2 );

template<size_t dim>
void reportAssignment( const std::string& prop, const std::string& unit,
                       const std::string& bound,
                       TensorVariable<dim>& ts1, TensorVariable<dim>& ts2,
                       TensorVariable<dim>& ts3, TensorVariable<dim>& ts4 );

void reportAssignment( const std::string& prop, const std::string& unit,
                       const std::string& bound,
                       ArrayVariable& av1, ArrayVariable& av2 );

void reportAssignment( const std::string& prop, const std::string& unit,
                       const std::string& bound,
                       ArrayVariable& av1, ArrayVariable& av2,
                       ArrayVariable& av3, ArrayVariable& av4 );

void reportAssignment( const std::string& prop, const std::string& unit,
                       const std::string& bound,
                       FlaggedArrayVariable& fv1, FlaggedArrayVariable& fv2 );

void reportAssignment( const std::string& prop, const std::string& unit,
                       const std::string& bound,
                       FlaggedArrayVariable& fv1, FlaggedArrayVariable& fv2,
                       FlaggedArrayVariable& fv3, FlaggedArrayVariable& fv4 );


/// functionalty specificly related to the configuration file

bool doesConfigFileExist( const std::string& filename_prefix );
bool doesVariableFileExist( const std::string& filename_prefix );

/// read point data  such as well rates ( a.k.a block7 )
template<size_t dim>
bool readPointData( std::map<std::string,std::vector<double> >& well_data, int ndata,
                    std::ifstream& ifs, char* text_line, size_t line_length, bool verbose );

/// build Regions based on the range of particular property data ( a.k.a block1 )
template<size_t dim>
bool buildRegionsBasedOnPropertyRange( Model<dim>& model,
                                       std::set<std::string>& groups,
                                       std::ifstream& ifs, char* text_line, size_t line_length, bool verbose );

/// read box boundary values and conditions ( a.k.a. block4 )
template<size_t dim>
bool readBoxBoundaryPropertyValuesAndConditions( Model<dim>& model,
                                                 std::ifstream& ifs, char* text_line, size_t line_length, bool verbose );

/// read boundary values and conditions ( a.k.a. block6 )
template<size_t dim>
bool readBoundaryPropertyValuesAndConditions( Model<dim>& model,
                                              std::ifstream& ifs, char* text_line, size_t line_length, bool verbose );

/// read boundary conditions
template<size_t dim>
bool readBoundaryPropertyConditions( Model<dim>& model,
                                     std::ifstream& ifs, char* text_line, size_t line_length, bool verbose );

/// read region conditions ( a.k.a block5 )
template<size_t dim>
bool readRegionPropertyConditions( Model<dim>& model,
                                   std::ifstream& ifs, char* text_line, size_t line_length, bool verbose );

/// read boundary property values
template<size_t dim>
bool readBoundaryPropertyValues( Model<dim>& model,
                                 std::ifstream& ifs, char* text_line, size_t line_length, bool verbose );

/// read region property values ( a.k.a block3 )
template<size_t dim>
bool readRegionPropertyValues( Model<dim>& model,
                               std::ifstream& ifs, char* text_line, size_t line_length, bool verbose );

/// read default property values ( a.k.a block2 )
template<size_t dim>
bool readDefaultPropertyValues( Model<dim>& model,
                                std::ifstream& ifs, char* text_line, size_t line_length, bool verbose );

/// read computational settings such as duration time, output times, time stepping strategy ( a.k.a block8 )
template<size_t dim>
bool readComputationalSettings( ComputationalSettings& settings,
                                std::ifstream& ifs, char* text_line, size_t line_length, bool verbose );

/// skip unknown or invalid block of data
bool readUnknownOrInvalidBlock( const std::string& keyword,
                                std::ifstream& ifs, char* text_line, size_t line_length, bool verbose );

} // end csmp


#endif
