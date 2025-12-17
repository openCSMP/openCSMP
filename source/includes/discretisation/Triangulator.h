#ifndef CSMP_TRIANGULATOR_H
#define CSMP_TRIANGULATOR_H

#include "CSMP_definitions.h"

namespace csmp {

class Matrix;
template<uint32_t> class VSet;

/**
 
@brief Creates a triangular element mesh from the pixels of the regular grid.
Boundaries that are defined by changes in the pixel values are used to control
how the pixels are split such that a smooth boundary is created.

@author Stephan K. Matthai
@author Stephen G. Roberts
@date 1999

@section motivation Motivation
 
Have a simple procedure to go from a geometry drawn with a pixel-type
tool or defined by a processed image, to go to a simple triangular
finite-element mesh for CSP computations.

*/
class Triangulator {
  public:
    Triangulator(); 

    void TrianglesFromRegularGrid( const Matrix& grid, VSet<2U>& vs );

  private:
     size_t m_mtrx, n_mtrx;
     size_t n_vertices, n_elements;
     double elperm, elperm2;

     short TestOutline( unsigned int m, unsigned int n, const Matrix& perm, double& el_perm );
     unsigned long MapVertex( unsigned int m, unsigned int n );
};  

/// reads pixel-based text input file and use triangulator to generate VSet
VSet<2U> readTextPixelData();

#endif


} // end namespace csmp
