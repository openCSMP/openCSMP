//
//  ModelPreProcessor2D.hpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 24/5/20.
//  Copyright © 2020 Stephan Matthai. All rights reserved.
//

#ifndef MODEL_PRE_PROCESSOR_3D_H
#define MODEL_PRE_PROCESSOR_3D_H

#include "CSMP_definitions.h"

namespace csmp {

template<uint32_t> class Model;


/**

Preparation of input models so that they have all the facilities (boundaries, splitboundaries, neighbor information etc.
to perform standard CSMP computations.

Output of CSMP native models into different formats for comparison runs with other software.

TODO: use VTK mesh and contouring of variables to perform a flux-adapted gridding.

*/
class ModelPreProcessor3D {
  public:
    /// reads CSMP native model from file to perform the preprocessing on
    explicit ModelPreProcessor3D( const std::string& csmp_model );
    ~ModelPreProcessor3D();
    
    /// tests:
    void CheckIntegrity( const std::string& report_file /* default name: model.Name()-integrity_check".txt" */ ) const;
    
  // UTILITIES FOR PROPERTY MANIPULATION
  private:
    
  private:
    Model<3U>*  model_ptr_ = 0;
};

} // end csmp

#endif /* MODEL_PRE_PROCESSOR_3D_H */
