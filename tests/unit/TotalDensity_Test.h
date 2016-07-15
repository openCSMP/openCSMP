#ifndef __TotalDensity_Test_h__
#define __TotalDensity_Test_h__

#include "CSMP_definitions.h"
#include "Test.h"
#include "Model.h"
#include "VSet.h"
#include "TRIANGLE_Interface.h"
//#include "TotalDensity.h"

namespace csmp {

class TotalDensity_Test : public Test {
public:
    TotalDensity_Test( const std::string& mesh_name, 
                       const std::string& var_file );
    ~TotalDensity_Test();
    void run();
    void elementTest();
    void nodeTest();
/*
private:
   void setNodeVariable( std::vector<csp_float>& var, const char* var_name );

     double tol_;
     TRIANGLE_Interface* mi_;
     VSet<csp_float, 2>* vset_;
     SuperGroup<csp_float, 2>* sg_;
     */
};

} // end namespace csp

#endif
