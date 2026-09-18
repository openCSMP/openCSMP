// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef VAN_GENUCHTEN_MODEL_H
#define VAN_GENUCHTEN_MODEL_H

#include "TwoPhaseModel.h"
#include "CubicSpline.h"

namespace csmp {

template<uint32_t dim>
class VanGenuchten : public TwoPhaseModel<dim> {
  public:
    VanGenuchten( const PropertyDatabase<dim>& database,
                  const char* model_parameter_n, // calculates m = 1 - 1/n
                  const char* model_parameter_alpha,
                  bool withFractureMatrixTransfer = false,
                  const bool sw_ro_mu_placement = true); // NODE=true ELEMENT=false


    /// constructor for constant capillary pressure equal to entry pressure
    VanGenuchten( const PropertyDatabase<dim>& database,
                 const char* permeability,
                 const char* viscosity_nw, const char* viscosity_w,
                 const char* density_nw, const char* density_w,
                 const char* model_parameter_n, // calculates m = 1 - 1/n
                 const char* model_parameter_alpha,
                 const char* sat_w, const char* res_sat_nw, const char* res_sat_w,
                 const char* permeability_matrix = "matrix permeability",
                 const char* res_sat_nw_matrix   = "residual saturation oil matrix",
                 const char* res_sat_w_matrix    = "residual saturation water matrix",
                 bool withFractureMatrixTransfer = false,
                 const bool sw_ro_mu_placement = true ); // NODE=true ELEMENT=fals


    virtual ~VanGenuchten();
    
    virtual void Initialize( const Element<dim>& e );
    void InitializeVirtualMatrix( const Element<dim>& e );

    // relative permeabilities
    virtual double krn_Phase() const;
    virtual double krw_Phase() const;

    // derivatives of relative permeabilities
    virtual double dkrnds_Phase() const;
    virtual double dkrwds_Phase() const;

    // capillary pressure
    virtual double pc_Phase( ) const;

    // capillary pressure derivatives
    virtual double dpcds_Phase( ) const;

    // inverse capillary pressure function
    virtual double Sw_Phase( double pc_Phase ) const;

    // inverse capillary pressure derivative
    virtual double dsdpc_Phase( double pc_Phase ) const;

    // maximum absolute value returned by dfdS
    virtual double MaxFractionalFlowDerivative() const;

    double PoreParametersCapillaryPressure( double pt, double pr, double tension ) const;

  protected:

    VanGenuchten();

  private:

    csmp::Index  n_key_,      // VG parameter (element property)
                 alpha_key_;  // alpha is the name of pd in this model

    csmp::Index  residualMatrixSaturationWKey_;
    csmp::Index  residualMatrixSaturationOKey_;
    csmp::Index  matrixPermeabilityKey_;
    bool         matrix_tensor_permeability_;
    mutable double  n_, alpha_;
    const   double  acc_gravity_,
                PC_LOW_SW_LIMIT_,
                PC_HIGH_SW_LIMIT_,
                KRW_HIGH_SW_LIMIT_,
                KRN_LOW_SW_LIMIT_;

   double m_from_n( double n_ ) const;

};

/**
@class VanGenuchten VanGenuchten "two_phase_flow/VanGenuchten.h"

@author S.K. Matthaei

Relative Permeability model according to original VanGenuchten formulation

*/

} // end namespace csmp

#endif
