// Copyright © 2021 Stephan Matthai. All rights reserved.
// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_DATA_TABLE_H
#define CSMP_DATA_TABLE_H

#include <string>
#include <vector>
#include <set>
#include <iostream>

namespace csmp {

class DataTable {
  public:
    /// creates a table from values and column headers supplied
    DataTable( std::string name, const std::vector<std::string>& col_headers,
               const std::vector<std::vector<double>>& rows_of_columns );
               
    DataTable( std::string name,
               const std::vector<std::string>& col_headers,
               const std::vector<std::string>& row_lables,
               const std::vector<std::vector<double>>& rows_of_columns );

    ~DataTable() = default;
    
    /// if table does not already contain a data column with the indicated header it is appended on the right
    bool AppendColumn( std::string col_header );
    
    size_t Rows() const { return table_.size(); }
    size_t Columns() const { return (table_.empty()) ? 0 : table_[0].size(); }
    
    /// checks whether there is a column header that matches the supplied variable
    bool ContainsColumn( std::string column_header ) const;
    
    /// finds in which column the values of the target variable are stored; returns UINT_MAX if the variable is not contained
    size_t ColumnIndex( std::string variable ) const;
    
    /// read access of elements in the table by row and column indices
    double  operator()( size_t row, size_t col ) const;
    /// write access
    double& operator()( size_t row, size_t col );
    
    /// whre empty sets are supplied as arguments, all column or rows will be extracted
    DataTable  ExtractSubTable( const std::set<std::string>& rows, const std::set<std::string>& columns ) const;
    DataTable  ExtractSubTable( const std::vector<int>& rows, const std::vector<int>& columns ) const;
    
    /// writes the data from the currently stored table into a comma delimited text file
    void Write_CSV_File( std::string filename ) const;
    
    void Out() const;
    
  private:
     std::vector<std::vector<double> >  table_;           ///< body of the table
     std::vector<std::string>           col_titles_;      ///< headings indicating the variable names
     std::vector<std::string>           row_identifiers_; ///< optional, could be region names
     std::string                        title_;           ///<  the name of the table
};

  /// reads CSV file collecting all the data needed to initialise the DataTable
  size_t read_CSV_File( std::string filename,
                        std::vector<std::string>& row_labels, std::vector<std::string>& col_headers,
                        std::vector<std::vector<double>>& rows_of_columns );


} // end csmp

#endif /* CSMP_DATA_TABLE_H */
