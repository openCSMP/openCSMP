#ifndef PROPERTYATPOINTVISITOR_TEST_H
#define PROPERTYATPOINTVISITOR_TEST_H

#include "Test.h"

#include "PropertyAtPointVisitor.h"

namespace csmp {

class PropertyAtPointVisitor_Test: public Test
   {
    public:
      PropertyAtPointVisitor_Test( bool verbose );
      ~PropertyAtPointVisitor_Test();

      virtual void run();

    private:
      void run1D_MeshTests();
      void run2D_MeshTests();
      void run3D_MeshTests();

      // 1D mesh test
      void IsoparametricLinear1DMesh_Test( const char* mesh_name, const char* regionfile_name, const char* varfile_name, const char* mesh_type, bool brute_force_search, double tolerance );
      // 2D mesh test
      void IsoparametricLinear2DMesh_Test( const char* mesh_name, const char* regionfile_name, const char* varfile_name, const char* mesh_type, bool brute_force_search, double tolerance );
      // 3D mesh test
      void IsoparametricLinear3DMesh_Test( const char* mesh_name, const char* regionfile_name, const char* varfile_name, const char* mesh_type, bool brute_force_search, double tolerance );

      void OutputElapsedTime( clock_t start, clock_t end );
     
      const bool verbose_;
};

}

#endif // PROPERTYATPOINTVISITOR_TEST_H
