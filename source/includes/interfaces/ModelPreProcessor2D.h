//
//  ModelPreProcessor2D.hpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 24/5/20.
//  Copyright © 2020 Stephan Matthai. All rights reserved.
//

#ifndef MODEL_PRE_PROCESSOR_2D_H
#define MODEL_PRE_PROCESSOR_2D_H

#include "CSMP_definitions.h"

namespace csmp {

template<size_t> class Model;


/**

Preparation of input models so that they have all the facilities (boundaries, splitboundaries, neighbor information etc.
to perform standard CSMP computations.

Output of CSMP native models into different formats for comparison runs with other software.

TODO: use VTK mesh and contouring of variables to perform a flux-adapted gridding.

*/
class ModelPreProcessor2D {
  public:
    explicit ModelPreProcessor2D( const std::string& csmp_model );
    ~ModelPreProcessor2D();
    
    /// writes log10 of permeability to regular grid; works only or quadrilateral (regular) meshes
    void ModelToMatrix( const std::string& output_file_name /* model.Name()-matrix".txt" */ ) const;
    
  // UTILITIES FOR PROPERTY MANIPULATION
    
  private:
    Model<2U>*  model_ptr_ = 0;
};




} // end csmp

#endif /* MODEL_PRE_PROCESSOR_2D_H */
