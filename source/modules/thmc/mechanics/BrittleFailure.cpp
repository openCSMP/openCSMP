//
//  BrittleFailure.cpp
//  CSMP_ReservoirSimulator
//
//  Created by Stephan Matthai on 10/25/14.
//  Copyright (c) 2014 Stephan K. Matthai. All rights reserved.
//

#include "BrittleFailure.h"
#include "CSMP_mathUtilities.h"

using namespace std;

namespace csmp {

std::string  parseFailure( FAILURE f )
 {
    if ( f ==  0 ) return "NONE";
    if ( f == -1 ) return "COMPRESSIVE";
    if ( f ==  1 ) return "FRICTIONAL_SLIDING";
    if ( f ==  2 ) return "SHEAR_FRACTURE";
    if ( f ==  3 ) return "TENSILE_OPENING";
    if ( f ==  4 ) return "MODE1_FRACTURE";
    return "failure mode undefined";
}

/// if the failure mode is not recognized, none is returned
FAILURE parseFailure( const std::string& failure )
 {
    if ( failure == "NONE" ) return NONE;
    if ( failure == "COMPRESSIVE" ) return COMPRESSIVE;
    if ( failure == "FRICTIONAL_SLIDING" ) return FRICTIONAL_SLIDING;
    if ( failure == "SHEAR_FRACTURE" ) return SHEAR_FRACTURE;
    if ( failure == "TENSILE_OPENING" ) return TENSILE_OPENING;
    if ( failure == "MODE1_FRACTURE" ) return MODE1_FRACTURE;
    return NONE;
}

/** 
     Evaluates failure criteria -1 to 4 without the explicit consideration
     of pore pressure (note that this works for fluid-saturated porous
     media as well, if the computed stresses are already effective stresses.
     
     @attention Since this approach is based on civil engineering concepts,
     tensile stresses are treated as positive and compressive ones as negative.
*/
FAILURE BrittleFailure::Evaluate( const MechanicalProperties& props,
                                  const csmp::StressInvariants& invars )
 {
    // 1. if maximum compressive stress is tensile, opening and tensile failure are considered first
    if ( invars.MaximumPrincipalStress1() > 0. ) {
          // evaluating tensile fracture
          if ( invars.LeastPrincipalStress3() >= props.TS_ ) return MODE1_FRACTURE;
          // just tensile opening
          return TENSILE_OPENING;
      }
   
    // establishing input parameters
    const double alpha(degreesToRadians(frictionAngle(props.mu_)));
    const double meanstress(invars.MeanStress());
    const double devstress(invars.DeviatoricStress());
   
    // 2. If the stress is essentially isostatic, failure may still occur
    // if the rock is porous and the compressive strength of the rock under isostatic
    // conditions is exceeded by the stress
    const double bar(1e5);
    if ( devstress < bar ) {
         // evaluation of the yield cap (Wong et al. 97, JGR) for isostatic state of stress
         // eqn. 2.29, Fjaer et al. 08', p. 68
         if ( (meanstress*meanstress + devstress*devstress) > (props.pstar_* props.pstar_) )
           return COMPRESSIVE;
      }
   
    // 3. if at least the maximum principal stress is compressive,
    // the likelyhood  of shear failure can be evaluated.
    // This is done here by calculating friction criterion Fmc
    // - material is described only in terms of its friction angle and cohesion.
    // (yielding occurs when Fmc >= 0)
    double FMC1 = meanstress * sin(alpha);
    FMC1         -= props.C_ * cos(alpha);
    FMC1         += devstress * invars.SmoothMohrCoulombYieldEnvelope( radiansToDegrees(alpha) );
    if ( FMC1 > 0. ) return SHEAR_FRACTURE;

    // 4. if shear failure does not occur, the likelyhood of frictional sliding on
    // pre-existing shear planes is evaluated.
    // This is done here by calculating friction criterion Fmc
    // - material is described only in terms of its friction angle and cohesion.
    // (yielding occurs when Fmc >= 0)
    // TODO: do this properly using the actual orientation of such a sliding plane
    double FMC2 = meanstress * sin(alpha);
    FMC2         += devstress * invars.SmoothMohrCoulombYieldEnvelope( radiansToDegrees(alpha) );
    if ( FMC2 > 0. ) return FRICTIONAL_SLIDING;
 
    return NONE;
 
 } // end Failure



/** 
     Evaluates failure criteria -1 to 4 considering pore pressure explicitly,
     establishing the effective stress:
     
     sigma_ij' = sigma_ij - pf * alpha * delta_ij
     
     where delta_ij is Kronecker delta.
 
     @attention Since this approach is based on civil engineering concepts,
     tensile stresses are treated as positive and compressive ones as negative.
     
     @attention The ordering of the principal stresses in engineering is such that
     sigma1 is the most negative = least tensile stress.
     
     @test looks OK, but not rigorously tested yet.
*/
FAILURE BrittleFailure::Evaluate( const MechanicalProperties& props,
                                  const csmp::StressInvariants& invars,
                                  double pf )
 {

    // 1. if even the maximum compressive stress + fluid pressure is tensile (=positive),
    //    a frictional / shear strength failure analysis does not apply and
    //    opening and tensile failure are considered first.
    if ( invars.MaximumPrincipalStress1() + props.alpha_ * pf > 0. ) {
          const double S3_eff(invars.LeastPrincipalStress3() + props.alpha_ * pf);
          // evaluating tensile fracture
          if ( S3_eff >= props.TS_ ) return MODE1_FRACTURE;
          // tensile opening
          return TENSILE_OPENING;
      }
   
    // 2. establishing input parameters for frictional analysis of failure
    const double alpha(degreesToRadians(frictionAngle(props.mu_)));
    // fluid pressure is added to meanstress because tension is positive
    const double meanstress(invars.MeanStress() + props.alpha_ * pf);
    const double devstress(invars.DeviatoricStress());
   
    // 3. If the stress is essentially isostatic, failure may still occur
    // if the rock is porous and the compressive strength of the rock under isostatic
    // conditions is exceeded by the stress
    const double bar(1e5);
    if ( devstress < bar ) {
         // evaluation of the yield cap (Wong et al. 97, JGR) for isostatic state of stress
         // eqn. 2.29, Fjaer et al. 08', p. 68
         if ( (meanstress*meanstress + devstress*devstress) > (props.pstar_* props.pstar_) )
           return COMPRESSIVE;
      }
   
    // 4. if the maximum principal stress is compressive, i.e. negative,
    // the likelyhood  of shear failure can be evaluated.
    // This is done here by calculating the friction criterion FMC1
    // only taking into account the friction angle and the cohesion.
    // (yielding occurs when Fmc >= 0)
    double FMC1 = meanstress * sin(alpha);
    FMC1         -= props.C_ * cos(alpha);
    FMC1         += devstress * invars.SmoothMohrCoulombYieldEnvelope( radiansToDegrees(alpha) );
    if ( FMC1 > 0. ) return SHEAR_FRACTURE;

    // 5. if shear failure does not occur, there is still the possibility of frictional sliding
    // on pre-existing shear planes. This is evaluated here by calculating a frictional sliding criterion FMC2.
    // TODO: do this properly using the actual orientation of such a sliding plane and a specific friction angle.
    double FMC2 = meanstress * sin(alpha);
    FMC2         += devstress * invars.SmoothMohrCoulombYieldEnvelope( radiansToDegrees(alpha) );
    if ( FMC2 > 0. ) return FRICTIONAL_SLIDING;
 
    return NONE;
 
 } // end Failure(with explicit fluid pressure effects)


/*  CHECKING
cerr <<"\nsigma 1,2,3: ";
cerr << invars.MaximumPrincipalStress1() <<" ";
cerr << invars.IntermediatePrincipalStress2() <<" ";
cerr << invars.LeastPrincipalStress3();
cerr <<" mean S: "<< invars.MeanStress();
cerr <<" dev S: "<< invars.DeviatoricStress();
cerr <<" pf: "<< pf;
props.Out();
*/





} // end csmp
