#ifndef THERMALCONVECTION_ETHZ_EXAMPLE_H
#define THERMALCONVECTION_ETHZ_EXAMPLE_H

#include "Example.h"

#include "CSMP_definitions.h"
#include "Model.h"
#include "H2OLookup.h"
#include "VTK_Interface.h"

namespace csmp {

    //template<size_t dim> class Model;
    //template<size_t dim> class Region;
    //class H2OLookup;
    //template<size_t dim> class VTK_Interface;
  
class ThermalConvectionETHZ_Example : public Example{
public:
  virtual void Run();
  virtual void Specifications();
    
private:
  /// sector volume, finite volume, FV pore volume
  void initialiseFiniteVolumeProperties( Model<3U>& );

  void initialFluidProperties( Model<3U>&, Region<3U>&, H2OLookup& );

  void transientFluidProperties( Model<3U>&, Region<3U>&, H2OLookup& );
  
  /// mobility, gravity term
  void steadyStatePressureOperands( Model<3U>&, Region<3U>&, H2OLookup& );

  /// system compressibility etc.
  void transientOperands( Model<3U>&, Region<3U>&, H2OLookup& );

  /// calculation of new density liquid and enthalpy in FV sectors
  void thermalEquilibrationAtSectorIntegrationPoints( Model<3U>&, Region<3U>&, H2OLookup& );

  /// nodal source-sink terms due to density amd energy content changes upon thermal equilibration
  void nodalSourceOperands( Model<3U>&, Region<3U>&, H2OLookup&, const double& time_increment );

  /// qm at integration point
  void massFlux( Model<3U>&, Region<3U>& );

  void outputToVTK( const Model<3U>&, VTK_Interface<3U>&, long time_step );
  
  void monitorVariableRange( Model<3U>& );
};

} // csmp

#endif // THERMALCONVECTIONETHZ_EXAMPLE_H
