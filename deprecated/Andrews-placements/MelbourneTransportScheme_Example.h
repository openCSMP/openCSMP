#ifndef MELBOURNE_TRANSPORT_SCHEME_EXAMPLE_H
#define MELBOURNE_TRANSPORT_SCHEME_EXAMPLE_H

#include "Example.h"

#include "CSMP_number_types.h"

namespace csmp {

  template<size_t dim> class Model;
  template<size_t dim> class VTK_Interface;
  template<size_t dim> class NodeCenteredFiniteVolumeTransport;

/// transport scheme illustrated with an example of the passive advection of a conservative tracer
class  MelbourneTransportScheme_Example : public Example {
  public:
    virtual void Run();
    virtual void Specifications();
    
  private:
    bool     NCFVT_methods( Model<3U>& , NodeCenteredFiniteVolumeTransport<3U>&, VTK_Interface<3U>& );
    
    double64 TestNodeCenteredFiniteVolumeTransport_PrescribedVelocity( Model<3U>& );
    void     TestNodeCenteredFiniteVolumeStencil( Model<3U>&, VTK_Interface<3U>& );
    void     TestNodeCenteredFiniteVolumeTransport( Model<3U>& );
    
    void AdvectVariableExplicit( Model<3U>&, bool second_order = true );
    void AdvectVariableExplicit( Model<3U>&, const char* group, bool second_order );
    void AdvectVariableFirstOrderImplicit( Model<3U>&, VTK_Interface<3U>& );
    void AdvectVariableFirstOrderImplicit( Model<3U>&, const char* group );
    void AdvectVariableSecondOrderImplicit( Model<3U>&, bool bijective_mapping );
    void AdvectVariableSecondOrderImplicitSecondOrderInTime( Model<3U>&, bool bijective_mapping );
};

} // csmp

#endif // PASSIVE_ADVECTION_OF_TRACER_EXAMPLE_H
