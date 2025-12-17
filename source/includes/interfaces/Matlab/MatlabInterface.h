#ifndef CSP_MATLAB_INTERFACE_H
#define CSP_MATLAB_INTERFACE_H

#include "Model.h"

namespace csmp {

class MatlabInterface {
  public:
    MatlabInterface();
    ~MatlabInterface();
  
    // output methods for two-dimensional models
    
    void Write2DMatlabFile( Model<2U>&, const char* file_name, const char* variable_name, long step, const char* region_name="Model" );
    void ExtractAndWrite1DVariableProfileAlongX( Model<2U>&, double distance, const char* name, const char* variable, long step );
    void ExtractAndWrite1DVariableProfileAlongY( Model<2U>&, double distance, const char* name, const char* variable, long step );
    void WriteSavedTimeStepsFile( long step, const char* name="saved_time_steps.txt" );
};

} // end csmp namespace

#endif
