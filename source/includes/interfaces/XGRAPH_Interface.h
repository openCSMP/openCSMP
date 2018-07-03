#ifndef XGRAPH_INTERFACE_H
#define XGRAPH_INTERFACE_H

#include "CSMP_definitions.h"

namespace csmp {

template<size_t>  class Model;

/**
@file XGRAPH_Interface.h
*/

/**
@addtogroup CSMPglobalFunctions
*/


/**

Writes CSMP 1D model to file: user can choose whether to append to existing file or
whether to start a new set of files.
    
Results are written to a set of files with the extensions:
    
      _saturation-oil.txt
      _fluid-pressure.txt
      _absolute-fluid-pressure.txt
      _total-mobility.txt
      _volume-flux.txt
      
@attention If the corresponding variables are not present in the model, this function will fail.
 
*/
void outputToXGraph( const Model<1U>& model, const char* file_name, double64 model_time, bool append_to_files=false );

/**
  @}
  */

}

#endif
