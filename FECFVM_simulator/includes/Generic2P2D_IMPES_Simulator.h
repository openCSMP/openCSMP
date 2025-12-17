#ifndef GENERIC_2P2D_IMPES_SIMULATOR_H
#define GENERIC_2P2D_IMPES_SIMULATOR_H

#include "CSMP_definitions.h"
#include "VTK_Interface.h"
#include "VTU_Interface.h"
#include "CSMP_definitions.h"
#include "IMPES_SimulatorKeys.h"
#include "ComputationalSettings.h"

namespace csmp {

template<uint32_t> class Model;
template<uint32_t> class Region;

template<uint32_t> class IMPES_Setup;
template<uint32_t> class TransportModel;
template<uint32_t> class ControlVolumeElement;
template<uint32_t> class TwoPhaseElementBasedTransport;
template<uint32_t> class TwoPhaseModel;
template<uint32_t> class PropertyHandle;

class Generic2P2D_IMPES_Simulator
{
public:
    Generic2P2D_IMPES_Simulator( const char* model_name,
                                 bool  ansys_model_true_or_csmp_binary_false,
                                 IMPES_Setup<2>& );
    ~Generic2P2D_IMPES_Simulator();
    void ConfigureModel();
    void InitializeProperties();
    void OutputToVTK(double time_unit = 60.);
    void OutputToVTU(double time_unit = 60.);
    void ComputeFlowProperties();
    void ComputeFlowPropertiesAllRegions();
    void ComputeFlowPropertiesActiveRegions();
    void ComputeTotalVelocity();
    void ComputeOilVelocity();
    void ComputeCapillaryPressure();
    // This method corrects pressure BC based on hydrostatic head. 
    // top_as_reference = true assigns pressure in configuration file
    // to the topest y-coordinates in the boundary and calculates 
    // the pressure at other points based on that reference.
    // top_as_reference = false does the calculation using the lowest 
    // point in the boundary as reference
    void CorrectPressureBoundaryConditionBasedonHydrostaticHead(bool top_as_reference);
    void WriteRestartFile();
    void ReadRestartFile();
    void ProcessPotentialSignals();
    void OutputFluidVolumes(double time, bool restart);
    void OutputFluidPressure(double time, bool restart) const;
    void OutputWellF_FluidPressure(double time, bool restart) const;
    void SolvePressure();
    void Run(bool restart = false);
    void CalculateEffectivePermeability();
    void OutputVelocityValues(double time);
    void OutputInflowOutflow(double time, bool restart);
    
    void ReadPermVolmod();
    void WritePermAndVolmodValues();
    void OutputCapillaryPressureGradient(double time);
    

    void FindContinuousLineRegionContainsCVE(ControlVolumeElement<2U>* elm,
                                             std::set<ControlVolumeElement<2U>* >& visited,
                                             std::vector<ControlVolumeElement<2U>* >& discovered );

    void SeparateContinuousLineElementsToDifferentRegionAndAssignProperty();

private:

    Generic2P2D_IMPES_Simulator();

    std::string                                     modelName_;
    Model<2U>*                                      reservoirModel_;
    VTK_Interface<2U>                               vtkOutput_;
    VTU_Interface<2U>*                              vtuOutput_;

    ComputationalSettings                           runSettings_;

    IMPES_Setup<2U>&                                setup_;
    IMPES_SimulatorKeys<2U>*                        propKeys_;

    TransportModel<2U>*                             transportModel_;
    TwoPhaseModel<2U>*                              relpermModel_;
    TwoPhaseElementBasedTransport<2U>*              advector_;
    PDE_Integrator<2U,Element>*                     pressureSolver_;
    DenseMatrix<DM_MIN>                             DERIV_;
};


} // end csmp

#endif
