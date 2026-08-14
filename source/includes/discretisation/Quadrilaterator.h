#ifndef QUADRILATERATOR_H
#define QUADRILATERATOR_H

#include "CSMP_definitions.h"
#include "Matrix.h"

namespace csmp {

template<uint32_t> class VSet;

/**
 
@brief Creates a quadrilateral finite-element mesh from the input bitmap.
The node and element numbering is from 0..n-1; negative numbers identify 
missing neighbor elements and corresponding model boundary identifiers.

@author Stephan K. Matthaei
@author Stephen G. Roberts
@date 1999

@section motivation Motivation
 
Have a simple procedure to go from a geometry drawn with a pixel-type
tool or defined by a processed image, to go to a simple triangular
finite-element mesh for CSMP computations.

*/
class Quadrilaterator {
  public:
    Quadrilaterator( bool harmonic_permeability_averaging=false ); 

    /// prompts user for color-coded permeability values (0..256) matrix in ascii text file and creates a regular quadrilateral mesh from this data
    void QuadrilateralsFromRegularGrid( VSet<2U>& vset, bool from_bitmap=true );
  
    /// reads 0..256 valued point matrix of color-coded permeability values from ascii text file and creates a regular, yet scaled quadrilateral mesh from it
    void QuadrilateralsFromRegularGrid( VSet<2U>& vset, const char* file_name,
                                        double x_extend, double y_extend, bool from_bitmap=true );
  
    /// creates a regular quadrilateral mesh with row_nodes-1 elements in the vertical, anr column_nodes-1 elements in the horizontal direction
    void QuadrilateralsFromRegularGrid( VSet<2U>& vset, double x_extend, double y_extend, size_t row_nodes, size_t column_nodes );

  private:
     double  HarmonicPermeabilityAverage( unsigned int m, unsigned int n, const csmp::Matrix& perm ) const;
     void      ReadPixelMatrix( bool from_bitmap );
     void      ReadPixelMatrix( const char* name, bool from_bitmap );
     void      ScaleModelRange( VSet<2U>& vset_scaled ) const;
     void      ScaleModelRange( VSet<2U>& vset_scaled, double x_dim, double y_dim ) const;
     size_t    MapVertex( unsigned int m, unsigned int n ) const;
     void      GenerateVSet( VSet<2U>& vset ) const;

     size_t  rows, cols, n_vertices, n_elements;
     bool    harmonic;
     Matrix  matrix_;
};  






} // end namespace csmp

#endif
