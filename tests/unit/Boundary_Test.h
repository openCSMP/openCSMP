#ifndef BOUNDARY_TEST_H
#define BOUNDARY_TEST_H

#include "Test.h"
#include "Model.h"
#include "VTU_Interface.h"

namespace csmp {

/** tests integrity of boundary also after retrieval from CSMP native binary model
 
        @author P. Lang
        @date 2011
 */
class Boundary_Test : public Test {
public:
  Boundary_Test( bool verbose=false ) : verbose_(verbose) {}

  virtual void run();

private:
  void runLegacy();
  void runCurrent();

  template<size_t dim>
  size_t InputElementAreaAsVolumeVariable( Model<dim>& model, Boundary<dim>& boundary, const char* variableName );

  template<size_t dim>
  void TestBoxBoundary( Model<dim>& model, const std::string& boundary, VTU_Interface<dim>& vtu );

  template<size_t dim>
  void CheckFaceNeighbors( const Boundary<dim>& boundary );

  template<size_t dim>
  void CheckNodeFlags( const Boundary<dim>& boundary, BOX_BOUNDARY flag, bool interiorOnly = false );

  template<size_t dim>
  void CheckFaceUnitNormalOrientation( const Boundary<dim>& boundary );
  
  /// using prism_test because it has a host of element types
  void UnitNormalTest3D();

  template<size_t dim>
  void CheckNodeParents( const Boundary<dim>& boundary );

  template<size_t dim>
  void ElementNodes( const csmp::Region<dim>& region );

  template<size_t dim>
  void NoSurfaceElementsAsNodeParents( const Region<dim>& region );
  
  bool verbose_;
};

} // csmp

#endif // BOUNDARY_TEST_H
