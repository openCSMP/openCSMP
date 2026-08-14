//
//  DataTable.cpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 20/12/21.
//  Copyright © 2021 Stephan Matthai. All rights reserved.
//

#include "DataTable.h"
#include "Exception.h"
#include <cassert>
#include <cstring>
#include <limits>
#include <fstream>
#include <ostream>

using namespace std;

namespace csmp {

DataTable::DataTable( string name, const vector<string>& col_titles, const vector<vector<double>>& rows_of_columns )
 : table_(rows_of_columns), col_titles_(col_titles), title_(name)
 {
 }


DataTable::DataTable( string name,
                      const vector<string>& col_titles,
                      const vector<string>& row_lables,
                      const vector<vector<double>>& rows_of_columns )
 : table_(rows_of_columns), col_titles_(col_titles), row_identifiers_(row_lables), title_(name)
 {
 }




/**
    If table does not already contain a data column with the indicated header it is appended on the right.
    
    @return true if a column was appended, false if the table was unmodified because it already contained
    a column with the desired header.
*/
bool DataTable::AppendColumn( std::string col_header )
 {
    if ( ContainsColumn(col_header) ) return false;
    
    col_titles_.push_back(col_header);
    const size_t n_columns_new{ col_titles_.size() };
    
    for ( auto& i : table_ )
      i.resize( n_columns_new, numeric_limits<double>::quiet_NaN() );
      
    return true;
 }





void DataTable::Out() const {
   cout <<"\nDataTable::Out: title '"<< title_ <<"':\n";
   for ( auto tit : col_titles_ ) cout << tit <<" ";
   cout <<"\n";
   size_t row{0};
   for ( auto i : table_ ) {
        if ( !row_identifiers_.empty() ) cout <<"\n"<< row_identifiers_[row] <<" ";
        else cout <<"\n ";
        for ( auto j : i ) cout << j <<" ";
     }
   cout << endl;
}


bool DataTable::ContainsColumn( std::string column_header ) const
 {
    for ( size_t idx{0}; idx<col_titles_.size(); ++idx )
      if ( col_titles_[idx] == column_header )
        return true;
        
    return false;
 }


    /// finds in which column the values of the target variable are stored; returns UINT_MAX if the variable is not contained
size_t DataTable::ColumnIndex( std::string variable ) const
 {
    for ( size_t idx{0}; idx<col_titles_.size(); ++idx )
      if ( col_titles_[idx] == variable )
        return idx;
        
    return std::numeric_limits<uint32_t>::max();
 }




    /// access of elements in the table by row and column indices
double DataTable::operator()( size_t row, size_t col ) const
 {
    assert( table_.size() >= 1 ); // at least one row
    assert( row < table_.size() );
    assert( col < table_[0].size() );
    
    return table_[row][col];
 }


double& DataTable::operator()( size_t row, size_t col )
 {
    assert( table_.size() >= 1 ); // at least one row
    assert( row < table_.size() );
    assert( col < table_[0].size() );
    
    return table_[row][col];
 }


    
    /// where empty sets are supplied as arguments, all column or rows will be extracted
DataTable  DataTable::ExtractSubTable( const set<string>& rows, const set<string>& cols ) const
 {
    assert( cols.size() <= col_titles_.size() );
    assert( rows.size() <= row_identifiers_.size() );
    assert( !cols.empty() );
    
    vector<vector<double> >  table( rows.size(), vector<double>( cols.size(), numeric_limits<double>::quiet_NaN() ) );
    vector<string>           col_titles( cols.size(), "column");
    vector<string>           row_labels( rows.size(), "row");
    
    // brute force copying
    bool   is_first_row{true};
    size_t col{0}, row{0};
    
    for ( auto i{0U}; i<table[0].size(); ++i )
      {
         if ( rows.find( row_identifiers_[i] ) != rows.end() ) {
             row_labels[row] = row_identifiers_[i];
             for ( size_t j{0}; j<table[i].size(); ++j ) {
                  if ( cols.find( col_titles_[j] ) != cols.end() ) {
                       if ( is_first_row ) {
                            col_titles[col] = col_titles_[j];
                            is_first_row = false;
                         }
                       // value assignment
                       table[row][col] = table[i][j];
                       col++;
                    }
                }
              row++;
           }
      }
 
    // new title: appended x_ r#, c#
    string title( title_ );
    title +="_subset_r";
    title += to_string( rows.size() );
    title += "_c";
    title += to_string( cols.size() );
 
    return DataTable( title, col_titles, row_labels, table );
    
 } // end ExtractSubTable (strings)
 
 
 
/**
      Extract rows and columns by number
*/
DataTable  DataTable::ExtractSubTable( const vector<int>& cols, const vector<int>& rows ) const
 {
    assert( !rows.empty() );
    assert( !cols.empty() );
    assert( cols.size() <= col_titles_.size() );
    assert( rows.size() <= row_identifiers_.size() );
    
    for ( auto i : rows ) assert( i < table_.size() );
    for ( auto j : cols ) assert( j < table_[0].size() );
    
    vector<vector<double> >  table( rows.size(), vector<double>( cols.size(), numeric_limits<double>::quiet_NaN() ) );
    vector<string>           col_titles( cols.size(), "column");
    vector<string>           row_labels( rows.size(), "row");
    
    // selected_rows
    size_t row{0};
    for ( auto rit : rows )
      row_labels[row++] = row_identifiers_[rit];
      
    // selected columns
    size_t col{0};
    for ( auto cit : cols )
      col_titles[col++] = col_titles_[cit];
    
    // table data
    row = 0;
    for ( auto i{0U}; i<rows.size(); ++i ) {
         for ( auto j{0U}, k{0U}; j<cols.size(); ++j )
           table[row][k++] = table_[ rows[i] ][ cols[j] ];
         row++;
      }
        
    // new title: appended x_ r#, c#
    string title( title_ );
    title +="_subset_r";
    title += to_string( rows.size() );
    title += "_c";
    title += to_string( cols.size() );
    
    return DataTable( title, col_titles, row_labels, table );
    
 } // end ExtractSubTable
 
 
 
/// writes the data from the currently stored table into a comma delimited text file
void DataTable::Write_CSV_File( std::string filename ) const
 {
   ofstream ofs( filename + ".csv" );
   ofs <<"\nDataTable::Out: title '"<< title_ <<"': read from file: "<< filename <<"\n";
   for ( auto tit : col_titles_ ) ofs << tit <<" ";
   ofs <<"\n";
   size_t row{0};
   for ( auto i : table_ ) {
        if ( !row_identifiers_.empty() ) ofs <<"\n"<< row_identifiers_[row] <<" ";
        else ofs <<"\n ";
        for ( auto j : i ) cout << j <<" ";
     }
   ofs << endl;
   
 } // end Write_CSV_File

 
 
 
 /*
 DataTable initialiseFrom_CSV_File( std::string csv_file )
  {
  } // end initialiseFrom_CSV_File
 */


/**
    Reads comma delimited point-data from ASCII file (ext. .csv),  returning them into the supplied vector,
    where the first  column is the row label if any and the following three column represent the point coordinates x,y,z.
    Then the data follow as indicated by the column titles.
    
    The first line is a headline with column header names, the number of which determines the number of columns in the file.
*/
size_t  read_CSV_File( string filename,
                       vector<string>& row_labels, vector<string>& col_titles,
                       vector<vector<double>>& rows_of_columns )
 {
    ifstream  ifs( filename + ".csv" );

    if ( !ifs.is_open() )
      throw csmp::Exception( ERROR, "read_CSV_File:", "input file could not be opened:", filename );
      
    const long LMAX(1024);
    const char* const delims =",";
    char              text_line[LMAX];
    char*             token(0);

    // reading and parsing the file header: region name, coordinates (double x dim), scalar-variable names
    // ---------------------------------------------------------------------------------------------------
    ifs.getline( text_line, LMAX );
    cout <<"\nread_CSV_File: file '"<< filename;
    cout <<"' header:\n\n"<< text_line << endl;
    
    // extracting information
    string  data_name(strtok(text_line,delims));
    string  data_X(strtok(NULL,delims));
    string  data_Y(strtok(NULL,delims));
    string  data_Z(strtok(NULL,delims));

    // parsing the column header / variable names
    if ( !col_titles.empty() ) col_titles.clear();
    else col_titles.reserve( 6 ); // row label, three coordinates, a value and a region identifier
    while( (token = strtok(NULL,delims)) != NULL  )
      col_titles.push_back(token);
    
    cout <<"\nread_CSV_File: file contains the variables:\n";
    for ( auto lt : col_titles )
      cout <<"\t"<< lt << endl;
      
    // row titles
    if ( !row_labels.empty() ) row_labels.clear();
    bool with_row_identifier = ( (*col_titles.begin()) != "x" ) ? true : false;
    if ( with_row_identifier ) row_labels.reserve( 100 );
   
    const size_t n_properties(col_titles.size() - 3 ); // the coordinate values
   
   
    // reading the data records storing them in the respective regions
    // ---------------------------------------------------------------
    double xmin(1.0e30), ymin(1.0e30), zmin(1.0e30), xmax(-1.0e30), ymax(-1.0e30), zmax(-1.0e30);
    if ( !rows_of_columns.empty() ) rows_of_columns.clear();
    size_t row_count{0};

    while ( !ifs.eof() )
      {
         // if the end of file or another errror is encountered, the reading proces is interrupted
         if ( !ifs.getline( text_line, LMAX ) ) break;

         // line is read, starting with the row_label if any
         if ( with_row_identifier ) {
              data_name = strtok(NULL,delims);
              row_labels.push_back( data_name );
           }
         // reading the point coordinates
         double x = atof(strtok(NULL,delims));
         double y = atof(strtok(NULL,delims));
         double z = atof(strtok(NULL,delims));
         xmin = std::min(xmin,x);
         xmax = std::max(xmax,x);
         ymin = std::min(ymin,y);
         ymax = std::max(ymax,y);
         zmin = std::min(zmin,z);
         zmax = std::max(zmax,z);
         // assigning these
         rows_of_columns.emplace_back( vector<double>{ x, y, z } );

         // reading the property values
         for ( auto i{0U}; i<n_properties; i++ ) {
              double data_value = atof(strtok(NULL,delims));
              // TODO: perhaps add a value check against property database here
              rows_of_columns[row_count].push_back( data_value );
           }
           
         row_count++;
        
      } // end while !eof
   
    ifs.close(); // data text file
   
    cout <<"\n\nread_CSV_File: successfully finished reading '"<< filename <<"' file.\n";
    cout <<"\n\tpoint data fall into the bounding box:\n";
    cout <<"\t\tx: "<< xmin <<" - "<< xmax <<" m.\n";
    cout <<"\t\ty: "<< ymin <<" - "<< ymax <<" m.\n";
    cout <<"\t\tz: "<< zmin <<" - "<< zmax <<" m.\n";
    
    return row_count;
 
 } // end read_CSV_File



} // end csmp
