// Copyright © 2020 Stephan Matthai. All rights reserved.
// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef MODEL_PRE_PROCESSOR_2D_H
#define MODEL_PRE_PROCESSOR_2D_H

#include "CSMP_definitions.h"

namespace csmp {

template<uint32_t> class Model;


/**

Preparation of input models so that they have all the facilities (boundaries, splitboundaries, neighbor information etc.
to perform standard CSMP computations.

Output of CSMP native models into different formats for comparison runs with other software.

TODO: use VTK mesh and contouring of variables to perform a flux-adapted gridding.

*/
class ModelPreProcessor2D {
  public:
    /// reads CSMP native model from file to perform the preprocessing on
    explicit ModelPreProcessor2D( const std::string& csmp_model );
    ~ModelPreProcessor2D();
    
    /// writes log10 of permeability to regular grid; works only or quadrilateral (regular) meshes
    void ModelToMatrix( const std::string& output_file_name /* model.Name()-matrix".txt" */,
                        size_t rows, size_t columns ) const;
    
  // UTILITIES FOR PROPERTY MANIPULATION
  private:
    /// computes an average grid-spacing from the distance beween cell centers in a given coordinate direction
//    double AverageGridSpacing( char x_or_y ) const;
    
  private:
    Model<2U>*  model_ptr_ = 0;
};




} // end csmp

#endif /* MODEL_PRE_PROCESSOR_2D_H */
