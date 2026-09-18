// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef MOHR_COULOMB_FAILURE_CRITERION_H
#define MOHR_COULOMB_FAILURE_CRITERION_H

#include "Interrelation.h"

namespace csmp {

/**
       Using the current local stress state stored on the model in the variable "stress",
       this Interrelation subclass computes the failure criterion 'failure' performing a Mohr Coulomb analysis,
       but without consideration of the fluid pressure.
       
       To apply this interrelation the following discretised variables need to be defined:
       
       "stress", "mean stress", "cohesion", "failure".
*/
template<uint32_t dim>
class MohrCoulombFailure : public Interrelation<dim> {
  public:
    MohrCoulombFailure( const PropertyDatabase<dim>&, double friction_angle );
                        
    ~MohrCoulombFailure() = default;
    
    /// computes the local value of result property "failure" which gets written back to the model
    void Calculate() override final;

  protected:
    ///  average of the magnitude of the principal stress
    double    MeanStress( const TensorVariable<dim>& );
    
    /// computes the shear-stress related stress invariant,  returning into second argument; method difference between sigma1 and sigma3
    double    DeviatoricStress( const TensorVariable<dim>&, double& t_stress_invariant );
    
    /// computes stress state invariance theta (Smith & Griffith, 2014, p. 236
    double    Theta( const TensorVariable<dim>&, double t_stress_invariant );
        
    /// computes stress state invariance, see Zienkiewitz II, p. 89
    double    G_OfTheta( double theta );

  private:
    Operand<dim>&  MS;        ///< mean stress
    Operand<dim>&  STRESS;
    Operand<dim>&  CRIT;
    Operand<dim>&  COH;       ///< cohesion of material
    double         phi;       ///< angle of friction
    double         Ts;        ///< tensile strength

    TensorVariable<dim>  ts_;
    ScalarVariable       ch_;
};


} // csmp

#endif

