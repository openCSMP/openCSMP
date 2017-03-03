#include "MohrCoulombFailure_Visitor.h"
#include "Model.h"
#include "Node.h"
#include "Exception.h"
#include <cmath>


using namespace std;

namespace csmp {

/**
    @attention SKM 28/9/2014 - added documentation, fixed wrong output and formatted code
*/
template<size_t dim>
MohrCoulombFailure_Visitor<dim>::MohrCoulombFailure_Visitor( Model<dim>& model,
                                                             bool positive_compressive_stress_convention,
                                                             bool verbose )
    : Visitor<dim>( MODEL, ELEMENT ),
      stressref_(model.Database()),
      verbose_(verbose),
      // getting csmp::Index values for the involved Node and ELEMENT variables
      Stress_key_(stressref_.StorageKey("stress")),
      Cohesion_key_(stressref_.StorageKey("cohesion")),
      Friction_key_(stressref_.StorageKey("friction angle")),
      Failure_key_(stressref_.StorageKey("failure")),
      Failure01_key_(stressref_.StorageKey("tensile failure")),
      // converting the friction angle in 'degrees' to 'radians'
      degrees_to_radians_( (2*3.14159265358979)/360. ),
      sqrt3_(sqrt(3.)),
      sqrt32_(sqrt(3./2.)),
      fluid_pressure_(0.),
      sign_of_tensile_stress_( positive_compressive_stress_convention ? -1.0 : 1.0 )
{
   // fluid pressure is only considered if it os defined
   if ( model.Database().IsDefined("fluid pressure") )
       pf_key_ = model.Database().StorageKey("fluid pressure");

    cout <<"\nMohrCoulombFailure_Visitor: diagnostics: ";
    if ( (Cohesion_key_.place != ELEMENT_INTEGRATION_POINT and Cohesion_key_.place != ELEMENT) || Cohesion_key_.type != SCALAR )
        throw csmp::Exception( CSMP_ERROR, "MohrCoulombFailure_Visitor (constructor):",
                              "'cohesion' must be a scalar variable." );

    if ( (Failure_key_.place != ELEMENT_INTEGRATION_POINT and Failure_key_.place != ELEMENT) || Failure_key_.type != SCALAR )
        throw csmp::Exception( CSMP_ERROR, "MohrCoulombFailure_Visitor (constructor):",
                               "'failure' must be a scalar variable." );

    if ( (Failure01_key_.place != ELEMENT_INTEGRATION_POINT and Failure01_key_.place != ELEMENT) || Failure01_key_.type != SCALAR )
        throw csmp::Exception( CSMP_ERROR, "MohrCoulombFailure_Visitor (constructor)",
                               "'failure01' must be a scalar variable" );

    if ( (Stress_key_.place != ELEMENT_INTEGRATION_POINT and Stress_key_.place != ELEMENT) || Stress_key_.type != TENSOR )
        throw csmp::Exception( CSMP_ERROR, "MohrCoulombFailure_Visitor (constructor)",
                               "'stress' must be a tensor variable" );
}



template<size_t dim>
MohrCoulombFailure_Visitor<dim>::~MohrCoulombFailure_Visitor()
{
}




/**
    Loops over the element integration points and evaluates shear and tensile failure potential.
    Where failure occurred the variable 'failure' is set to 1.
    A distinction is made between tensile 'tensile failure' and shear failure 'failure' variables.
    
    @attention it is assumed that the tensile strength is about 1/10 of the cohesive strength of the rock
    
    @attention SKM 28/9/2014 - added case where stress is placed on the element.
*/
template<size_t dim>
void MohrCoulombFailure_Visitor<dim>::Visit( Element<dim>* e )
 {
    // cohesion, failure variables are in integration point
    if ( Stress_key_.place == ELEMENT_INTEGRATION_POINT ) {
        const size_t integration_points(e->FE()->IntegrationPoints());
        for ( size_t i=0U; i<integration_points; i++ )
          {
            // 1. reading input variables
            e->Read(i, Stress_key_, Cartesian_stress_ );
            double64 cohesion = e->Read(i,Cohesion_key_);
            double64 alpha    = e->Read(i,Friction_key_) * degrees_to_radians_;
            if ( pf_key_.place != UNDEFINED )
              fluid_pressure_ = e->PropertyValueAtIntegrationPoint( pf_key_, i );

            // six components of symmetric stress tensor
            double64 sigmaxx = Cartesian_stress_(0,0);
            double64 sigmayy = Cartesian_stress_(1,1);
            double64 sigmazz = Cartesian_stress_(2,2);
            double64 sigmaxy = Cartesian_stress_(0,1);
            double64 sigmaxz = Cartesian_stress_(0,2);
            double64 sigmayz = Cartesian_stress_(1,2);

            // Biot coefficient (alpha) should be included as
            // alpha = 1. - (K/K_{s}),
            // where K and K_{s} are effective and bulk moduli correspondingly
            // for the time being alpha assumed to be 1
            biot_coefficient_alpha_ = 1.0;

            // 2. computing parameter s
            //    Smith & Griffiths, eqn. 6.3, p. 227
            //    s ( distance from the origin to the pi-plane in which the stress point lies )
            const double64 s = (dim==2) ? (sigmaxx + sigmayy) / sqrt3_
                                        : (sigmaxx + sigmayy + sigmazz) / sqrt3_;

            // 3. computing parameter t
            //    Smith & Griffiths, eq. 6.3, p. 227
            //    t ( perpendicular distance of the stress point from the space diagonal: sigma1 = sigma2 = sigma3 )
            double64 t(0.);
            if ( dim == 2 ) {
                 t = ((sigmaxx - sigmayy)*(sigmaxx - sigmayy)) + 3. * sigmaxy * sigmaxy;
              }
            else if ( dim == 3 ) {
                 // assuming a symmetric stress tensor
                 t  = ((sigmaxx - sigmayy)*(sigmaxx - sigmayy));
                 t += ((sigmayy - sigmazz)*(sigmayy - sigmazz));
                 t += ((sigmazz - sigmaxx)*(sigmazz - sigmaxx));
                 t += 6. * (sigmaxy * sigmaxy) + 6. * (sigmayz * sigmayz) + 6. * (sigmaxz * sigmaxz);
              }
            t  = sqrt(t) / sqrt3_;

            // 4. calculate Lode angle(theta), in radians
            //    Smith & Griffiths, eq. 6.3. p. 227
            double64 sx, sy, sz, J3(std::numeric_limits<double64>::quiet_NaN());
            if ( dim == 2 ) {
                // 2D case not sure yet, search reference
                sx  = (2. * sigmaxx - sigmayy) / 2.;
                sy  = (2. * sigmayy - sigmaxx) / 2.;
                J3  = sx * sy + 2.* sigmaxy;
              }
            else if ( dim == 3 ) {
                sx  = (2. * sigmaxx - sigmayy - sigmazz) / 3.;
                sy  = (2. * sigmayy - sigmazz - sigmaxx) / 3.;
                sz  = (2. * sigmazz - sigmaxx - sigmayy) / 3.;
                J3  = sx * sy * sz;
                J3 -= sx * (sigmayz * sigmayz);
                J3 -= sy * (sigmaxz * sigmaxz);
                J3 -= sz * (sigmaxy * sigmaxy);
                J3 += 2. * sigmaxy * sigmaxz * sigmayz;
              }
            // the Lode angle theta
            const double64 theta = 1./3. * asin( (-3.* sqrt(6.) * J3) / (t * t * t) );

            // 5. computing mean stress,
            //    Smith & Griffiths, eq. 6.4, p. 228
            //    mean stress = sqrt(1/3)*s
            double64 meanstress = s/sqrt3_;

            // taking into account sign convention,
            // since current formulation is ment to be used for positive tensile stress convention,
            // the correction of sign must be perfomed in the mean stress
            // if one is applying positive compressive stress convention
            // ------------------------------
            meanstress *= sign_of_tensile_stress_;

            // taking into account pore pressure
            // note that plus is because of positive tensile stress convention
            // no conversion is needed here, because it is lined up with current formulation
            // ------------------------------
            meanstress += biot_coefficient_alpha_*fluid_pressure_;

            // 6. computing deviatoric stress,
            //    Smith & Griffiths, eq. 6.4, p. 228
            //    deviatoric stress = sqrt(3/2)*t
            const double64 devstress = t*sqrt32_;

            // 7. computing K_OfTheta ( original Mohr-Coulomb yield surface )
            //    Zienkiewitz, Finite element method for solid and structural mechanics,
            //    Chapter 4. Inelastic and non-linear materials
            //    4.5.1 Isotropic yield surfaces
            //double64 K_OfTheta = -sin(alpha);
            //K_OfTheta *= sin(theta);
            //K_OfTheta /= sqrt3_;
            //K_OfTheta += cos(theta);
            //K_OfTheta /= sqrt3_;

            // 7. computing K_OfTheta = 1/G_OfTheta ( modified Mohr-Coulomb envelope with smooth boundaries ):
            //    Zienkiewitz, Finite element method for solid and structural mechanics,
            //    Chapter 4. Inelastic and non-linear materials
            //    4.11. Non-uniqueness and localization in elasto-plastic deformations
            double64 K(sin(alpha));
            K = (3. - K) / (3. + K);
            const double64 G_OfTheta = (2.*K) / ((1+K) - sin(3. * theta) * (1-K));
            const double64 K_OfTheta = 1.0/G_OfTheta;

            // 8. calculating friction criterion Fmc
            //    Yielding can occur when Fmc >= 0,
            //    the material is described only in terms of its friction angle and cohesion.
            // ------------------------------
            // TODO: organise failure modes like in fault module
            // TODO: add failure if bulk modulus is exceeded
            // TODO: deal with pure tensile failure if normal stress is negative
            double64 Fmc = meanstress * sin(alpha) - cohesion * cos(alpha) + devstress * K_OfTheta;

            // tensile failure determination assuming that the tensile strength is about 0.1 of
            // the value of the cohesion
            const double64 tensile_strength(cohesion * 0.1);
            double64 F01 = ( (Fmc - tensile_strength ) <= 0. ) ? F01 = 1. : F01 = 0.;

            // 9. store variables
            e->Store(i, Failure_key_,   makeScalar(PLAIN,Fmc) );
            e->Store(i, Failure01_key_, makeScalar(PLAIN,F01) );

            // 10. output if desired
            if ( verbose_ ) {
                cout <<"\nshear failure criterion F "<< endl;
                cout.width(35); cout <<"at element integration point ";
                cout << i << ": "<< Fmc << endl;
                cout.flush();
            }
         }
       return;
     }
  
     
    // cohesion, failure variables are placed on the element
    if ( Stress_key_.place == ELEMENT )
      {
          // 1. reading input variables
          e->Read( Stress_key_, Cartesian_stress_ );
          const double64 alpha    = e->Read(Friction_key_) * degrees_to_radians_;
          const double64 cohesion = e->Read( Cohesion_key_ );
          if ( pf_key_.place != UNDEFINED ) {
               ScalarVariable sc;
               e->PropertyValueAtBaryCenter( pf_key_, sc );
               fluid_pressure_ = sc();
            }
          const double64 sigmaxx = Cartesian_stress_(0,0);
          const double64 sigmayy = Cartesian_stress_(1,1);
          const double64 sigmazz = Cartesian_stress_(2,2);
          const double64 sigmaxy = Cartesian_stress_(0,1);
          const double64 sigmaxz = Cartesian_stress_(0,2);
          const double64 sigmayz = Cartesian_stress_(1,2);

          // Biot coefficient (alpha) should be included as
          // alpha = 1. - (K/K_{s}),
          // where K and K_{s} are effective and bulk moduli correspondingly
          // for the time being alpha assumed to be 1
          biot_coefficient_alpha_ = 1.;

          // 2. computing parameter s
          //    Smith & Griffiths, eqn. 6.3, p. 227
          //    s ( distance from the origin to the pi-plane in which the stress point lies )
          const double64 s = (dim==2) ? (sigmaxx + sigmayy) / sqrt3_
                                      : (sigmaxx + sigmayy + sigmazz) / sqrt3_;

          // 3. computing parameter t
          //    Smith & Griffiths, eq. 6.3, p. 227
          //    t ( perpendicular distance of the stress point from the space diagonal: sigma1 = sigma2 = sigma3 )
          double64 t(0.);
          if ( dim == 2 ) {
               t = ((sigmaxx - sigmayy)*(sigmaxx - sigmayy)) + 3. * sigmaxy * sigmaxy;
            }
          else if ( dim == 3 ) {
               // assuming a symmetric stress tensor
               t  = ((sigmaxx - sigmayy)*(sigmaxx - sigmayy));
               t += ((sigmayy - sigmazz)*(sigmayy - sigmazz));
               t += ((sigmazz - sigmaxx)*(sigmazz - sigmaxx));
               t += 6. * (sigmaxy * sigmaxy) + 6. * (sigmayz * sigmayz) + 6. * (sigmaxz * sigmaxz);
            }
          t  = sqrt(t) / sqrt3_;

          // 4. calculate Lode angle(theta), in radians
          //    Smith & Griffiths, eq. 6.3. p. 227
          double64 sx, sy, sz, J3(std::numeric_limits<double64>::quiet_NaN());
          if ( dim == 2 ) {
              // 2D case not sure yet, search reference
              sx  = (2. * sigmaxx - sigmayy) / 2.;
              sy  = (2. * sigmayy - sigmaxx) / 2.;
              J3  = sx * sy + 2.* sigmaxy;
            }
          else if ( dim == 3 ) {
              sx  = (2. * sigmaxx - sigmayy - sigmazz) / 3.;
              sy  = (2. * sigmayy - sigmazz - sigmaxx) / 3.;
              sz  = (2. * sigmazz - sigmaxx - sigmayy) / 3.;
              J3  = sx * sy * sz;
              J3 -= sx * (sigmayz * sigmayz);
              J3 -= sy * (sigmaxz * sigmaxz);
              J3 -= sz * (sigmaxy * sigmaxy);
              J3 += 2. * sigmaxy * sigmaxz * sigmayz;
            }
          // the Lode angle theta
          const double64 theta = 1./3. * asin( (-3.* sqrt(6.) * J3) / (t * t * t) );

          // 5. computing mean stress,
          //    Smith & Griffiths, eq. 6.4, p. 228
          //    mean stress = sqrt(1/3)*s
          double64 meanstress = s/sqrt3_;

          // taking into account sign convention,
          // since current formulation is ment to be used for positive tensile stress convention,
          // the correction of sign must be perfomed in the mean stress
          // if one is applying positive compressive stress convention
          // ------------------------------
          meanstress *= sign_of_tensile_stress_;

          // taking into account pore pressure
          // note that plus is because of positive tensile stress convention
          // no conversion is needed here, because it is lined up with current formulation
          // ------------------------------
          meanstress += biot_coefficient_alpha_*fluid_pressure_;

          // 6. computing deviatoric stress,
          //    Smith & Griffiths, eq. 6.4, p. 228
          //    deviatoric stress = sqrt(3/2)*t
          const double64 devstress = t*sqrt32_;

          // 7. computing K_OfTheta ( original Mohr-Coulomb yield surface )
          //    Zienkiewitz, Finite element method for solid and structural mechanics,
          //    Chapter 4. Inelastic and non-linear materials
          //    4.5.1 Isotropic yield surfaces
          //double64 K_OfTheta = -sin(alpha);
          //K_OfTheta *= sin(theta);
          //K_OfTheta /= sqrt3_;
          //K_OfTheta += cos(theta);
          //K_OfTheta /= sqrt3_;

          // 7. computing K_OfTheta = 1/G_OfTheta ( modified Mohr-Coulomb envelope with smooth boundaries ):
          //    Zienkiewitz, Finite element method for solid and structural mechanics,
          //    Chapter 4. Inelastic and non-linear materials
          //    4.11. Non-uniqueness and localization in elasto-plastic deformations
          double64 K(sin(alpha));
          K = (3. - K) / (3. + K);
          const double64 G_OfTheta = (2.*K) / ((1+K) - sin(3. * theta) * (1-K));
          const double64 K_OfTheta = 1.0/G_OfTheta;

          // 8. calculating friction criterion Fmc
          //    Yielding can occur when Fmc >= 0,
          //    the material is described only in terms of its friction angle and cohesion.
          // ------------------------------
          // TODO: organise failure modes like in fault module
          // TODO: add failure if bulk modulus is exceeded
          // TODO: deal with pure tensile failure if normal stress is negative
          double64 Fmc = meanstress * sin(alpha) - cohesion * cos(alpha) + devstress * K_OfTheta;

          // tensile failure determination assuming that the tensile strength is about 0.1 of
          // the value of the cohesion
          const double64 tensile_strength(cohesion * 0.1);
          double64 F01 = ( (Fmc - tensile_strength ) <= 0. ) ? F01 = 1. : F01 = 0.;

          // 9. store variables
          e->Store( Failure_key_,   makeScalar(PLAIN,Fmc)   );
          e->Store( Failure01_key_, makeScalar(PLAIN,F01) );

          // 10. output if desired
          if ( verbose_ ) {
              cout <<"\nshear failure criterion F "<< endl;
              cout.width(35); cout <<" in element "<< e->Idx();
              cout <<": "<< Fmc << endl;
              cout.flush();
          }
      }
  
} // end Visit(Element)

template class MohrCoulombFailure_Visitor<1U>;
template class MohrCoulombFailure_Visitor<2U>;
template class MohrCoulombFailure_Visitor<3U>;

} // end csmp

