// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef TWO_PHASE_ELEMENT_BASED_TRANSPORT_H
#define TWO_PHASE_ELEMENT_BASED_TRANSPORT_H

#include "CSMP_definitions.h"
#include "DenseMatrix.h"
#include "Element.h"

namespace csmp {

    template<uint32_t> class Model;
    template<uint32_t> class Region;
    template<uint32_t> class TwoPhaseModel;
    template<uint32_t> class TransportModel;
    template<uint32_t> class ElementFace;
    template<uint32_t> class ControlVolumeElement;
    template<uint32_t> class IMPES_Setup;

    template<uint32_t dim>
    class TwoPhaseElementBasedTransport
    {
    public:
        TwoPhaseElementBasedTransport(Model<dim>& model,
                                      TransportModel<dim>& transport_model,
                                      const char* region_name,
                                      const char* porosity,
                                      const char* permeability,
                                      const char* saturation_oil,
                                      const char* saturation_water, // it needs this key to update it
                                      const char* divergence,       // class uses this storage to store the divergence
                                      const char* total_velocity,
                                      const char* element_oil_injection_volume_rate, // only positive values are allowed
                                      const char* element_water_injection_volume_rate, // only positive values are allowed
                                      const char* element_total_production_volume_rate, // only negative values are allowed
                                      IMPES_Setup<dim>& setup);

        virtual ~TwoPhaseElementBasedTransport();

        void                                  AssignTransportBoundaryCondition(const char* boundary_node_oil_saturation);
        void                                  Divergence(double time);
        void                                  ConstructFluxes();
        virtual double                      Transport(double total_time_increment, TwoPhaseModel<dim>& flow_model);
        void                                  FindUpstreamProperties(const double& face_flux, const double& area, 
                                                                     const Point<dim>& face_normal, const double& perm,
                                                                     const double& lambda_w_out, const double& lambda_o_out,
                                                                     const double& lambda_w_in, const double& lambda_o_in,
                                                                     const double& ro_w_out, const double& ro_o_out,
                                                                     const double& ro_w_in, const double& ro_o_in,
                                                                     const double& normal_capillary_grad,
                                                                     double& lambda_w_upstream, double& lambda_o_upstream,
                                                                     double& ro_w_upstream, double& ro_o_upstream) const;

    
    protected:
        void                                  CalculateBarycentersDistanceDerivativeTerm();
        std::pair<double, double>         CFL_Timestep_MaxDivergence();

        double CalculateOilFluxInteriorLineFace(ElementFace<dim>& face, TwoPhaseModel<dim>& flow_model);
        double CalculateOilFluxBoundaryLineFace(ElementFace<dim>& face, TwoPhaseModel<dim>& flow_model);
        double CalculateOilFluxInteriorPointFace(ElementFace<dim>& face, TwoPhaseModel<dim>& flow_model);
        void TransportAcrossZeroDimensionalCVE_Faces(ControlVolumeElement<dim>& cve, TwoPhaseModel<dim>& flow_model);


    protected:
        Model<dim>&                                                       modelRef_;
        Region<dim>&                                                      regionRef_;
        TransportModel<dim>&                                              transportModelRef_;
        IMPES_Setup<dim>&                                                 setup_;
        DenseMatrix<DM_MIN>                                               DERIV_;
        csmp::Index                                                       totalVelocityKey_;
        csmp::Index                                                       permeabilityKey_;
        csmp::Index                                                       oilSaturationKey_;
        csmp::Index                                                       waterSaturationKey_;
        csmp::Index                                                       porosityKey_;
        csmp::Index                                                       divergenceKey_;
        csmp::Index                                                       totalProdRateKey_;
        csmp::Index                                                       oilInjRateKey_;
        csmp::Index                                                       waterInjRateKey_;
        const bool                                                        withGravity_;
        const bool                                                        withCapillary_;
        const double                                                    G_; // gravitational acceleration
        const double                                                    EPSILON_; // epsilon to avoid devision by zero

    };


} // end namespace csmp

#endif
