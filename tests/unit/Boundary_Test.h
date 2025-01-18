#ifndef BOUNDARY_TEST_H
#define BOUNDARY_TEST_H

#include "Test.h"
#include "Model.h"
#include "VTU_Interface.h"

namespace csmp {

/** tests integrity of boundary also after retrieval from CSMP native binary model
 
        @author P. Lang
        @date 2011

        @author SKM (adding extra tests)
        @date 2024
 */
class Boundary_Test : public Test {
public:
  virtual void run();

private:
  const static bool verbose_ = false;

  void runLegacy(); // Ansys models TODO: refactor to use CSMP native models
  void runCurrent(); // Ansys models

  template<uint32_t dim>
  size_t InputElementAreaAsVolumeVariable( Model<dim>&, Boundary<dim>&, const char* variableName );

  template<uint32_t dim>
  void TestBoxBoundary( Model<dim>&, const std::string& boundary, VTU_Interface<dim>& );

  template<uint32_t dim>
  void CheckFaceNeighbors( const Boundary<dim>& );

  template<uint32_t dim>
  void CheckNodeFlags( const Boundary<dim>&, BOX_BOUNDARY flag, bool interiorOnly = false );

  template<uint32_t dim>
  void CheckFaceUnitNormalOrientation( const Boundary<dim>& );
  
  /// using prism_test because it has a host of element types
  void UnitNormalTest3D();

  template<uint32_t dim>
  void CheckNodeParents( const Boundary<dim>& );

  template<uint32_t dim>
  void ElementNodes( const csmp::Region<dim>& );

  template<uint32_t dim>
  void NoSurfaceElementsAsNodeParents( const Region<dim>& );

  /// checks that variable flags are set correctly for INTERIOR, PERIMETER and COMPLETE discriminators
  bool Test_ChangeBoundaryStatus();
};

} // csmp

#endif // BOUNDARY_TEST_H
