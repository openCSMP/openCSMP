#ifndef CSMP_LOOPUNROLLER_H
#define CSMP_LOOPUNROLLER_H

#include "CSMP_definitions.h"

namespace csmp {

class LoopUnroller {
  public:
    LoopUnroller();
    ~LoopUnroller();
    
    // a[i] expr b[i];
    void Unroll( const char* file,
                 int i,
                 const char* var1, const char* var2,
                 const char* expression );

    // a[i][j]; expr b[i][j];
    void Unroll( const char* file,
                 int i, int j,
                 const char* var1, const char* var2,
                 const char* expression );

    void Unroll( const char* file,
                 const char* rows, const char* cols,
                 int i, int j,
                 const char* var1, 
                 const char* expression1, 
                 const char* var2 );

 };

} // csmp

#endif
