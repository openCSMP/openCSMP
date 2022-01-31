#ifndef CSMP_TEXT_INTERFACE_H
#define CSMP_TEXT_INTERFACE_H

#include "CSMP_definitions.h"
#include "DenseMatrix.h"
#include "CoordinateTransformer.h"

namespace csmp {

class Matrix;
template<size_t>  class Model;

/**
     CSMP Output interface for writing model variables to text files.
     
     @author S. K. Matthai
     @date 14/5/1999
     @copyright SKM
*/
class TextInterface {
  public:
    TextInterface();
  
    /// single variable text column output
    template<size_t dim>
    void OutputDataAsTextColumns( const Model<dim>&,
                                  const char* file_name, const char* s ) const;
    
    /// single variable text column output for a user specified region
    template<size_t dim>
    void OutputDataAsTextColumns( const char* region, const Model<dim>&, 
                                  const char* file_name, const char* s ) const;

    /// output of a suite of variables with same placement for a specific region; user can transform coordinates
    template<size_t dim>
    void OutputDataAsTextColumns( const Model<dim>&,
                                  const char* file_name, const char* region,
                                  const CoordinateTransformer<dim>&,
                                  const std::list<std::string>& output_variables ) const;

    /// (row numbered) single variable output
    template<size_t dim>
    void OutputDataAsTextColumnsNumbered( const Model<dim>&,
                                          const char* file_name, const char* s ) const;

    /// (row numbered) single variable output for a user specifed region
    template<size_t dim>
    void OutputDataAsTextColumnsNumbered( const char* region, const Model<dim>&,
                                          const char* file_name, const char* s ) const;
  
    /// (row numbered) single variable output to files uniquely identified by timestep
    template<size_t dim>
    void OutputDataAsTextColumns( const Model<dim>&,
                                  const char* file_name, const char* s, long timestep, bool numbered=false ) const;

    template<size_t dim>
    void OutputRegionsToTextFiles( const Model<dim>& sg, const char* property ) const;
    
    // text output from models   
    void AppendDataToText( const char* fname, double time, long idx, double value );
    void WriteStringToTextFile( const char* fname, const std::string& str, bool overwrite );
    void WriteStringToTextFile( const char* fname, const char* s, bool overwrite );
    void WriteMatrixToTextfile( const char* fname, DenseMatrix<DM_MIN>& mtrx );

    // Regular-gridded image files
    void SizeofPixelTextImage256( const char* fname, size_t& m, size_t& n );
    void ReadPixelTextImage256(   const char* fname, Matrix& data );

  private:
      bool IsInNextLine( std::ifstream& ifs, const char* search_string ) const;
};

} // csmp

#endif
