#include "BrittleFailureVisitor.h"
#include "BrittleFailure.h"
#include "MechanicalProperties.h"
#include "Model.h"
#include "Node.h"
#include "Exception.h"
#include "CSMP_mathUtilities.h"
#include <cmath>


using namespace std;

namespace csmp {

/**
    @attention SKM 28/9/2014 - added documentation, fixed wrong output and formatted code
*/
template<size_t dim>
BrittleFailureVisitor<dim>::BrittleFailureVisitor( Model<dim>& model,
                                                   bool verbose )
    : Visitor<dim>( MODEL, ELEMENT ),
      verbose_(verbose),
      // getting csmp::Index values for the involved Node and ELEMENT variables
      Youngs_key_(model.Database().StorageKey("Youngs modulus")),
      Poissons_key_(model.Database().StorageKey("Poissons ratio")),
      Stress_key_(model.Database().StorageKey("stress")),
      Cohesion_key_(model.Database().StorageKey("cohesion")),
      Friction_key_(model.Database().StorageKey("friction angle")),
      Failure_key_(model.Database().StorageKey("failure")),
      // converting the friction angle in 'degrees' to 'radians'
      fluid_pressure_(PLAIN,0.)
{
   // fluid pressure is only considered if it is defined
   if ( model.Database().IsDefined("fluid pressure") )
       Pressure_key_ = model.Database().StorageKey("fluid pressure");
   // Biot coefficient can only be considered if the fluid pressure is defined
   if ( model.Database().IsDefined("Biot alpha") )
       Biot_key_ = model.Database().StorageKey("Biot alpha");

    cout <<"\nBrittleFailureVisitor: diagnostics: ";
    if ( (Cohesion_key_.place != ELEMENT_INTEGRATION_POINT and Cohesion_key_.place != ELEMENT) || Cohesion_key_.type != SCALAR )
        throw csmp::Exception( CSMP_ERROR, "BrittleFailureVisitor (constructor):",
                              "'cohesion' must be a scalar variable." );

    if ( (Failure_key_.place != ELEMENT_INTEGRATION_POINT and Failure_key_.place != ELEMENT) || Failure_key_.type != SCALAR )
        throw csmp::Exception( CSMP_ERROR, "BrittleFailureVisitor (constructor):",
                               "'failure' must be a scalar variable." );

    if ( (Stress_key_.place != ELEMENT_INTEGRATION_POINT and Stress_key_.place != ELEMENT) || Stress_key_.type != TENSOR )
        throw csmp::Exception( CSMP_ERROR, "BrittleFailureVisitor (constructor)",
                               "'stress' must be a tensor variable" );
}



template<size_t dim>
BrittleFailureVisitor<dim>::~BrittleFailureVisitor()
{
}



/** 
    Reads all relevant rock and fluid properties from the element.
    Missing parameters are computed from correlations:
    
    - UCS is calculated from E and the Hoek-Brown failure envelope.
 
    - pstar = compactive strength is estimated from UCS
 
    - G and B are calculated from E and nu assuming perfect linear elasticity
*/
template<size_t dim>
void BrittleFailureVisitor<dim>::InitializeInputProperties( Element<dim>* eptr )
  {
     if ( Stress_key_.place == ELEMENT )
       eptr->Read( Stress_key_, Cartesian_stress_ );

     mprops_.E  = eptr->Read( Youngs_key_ );          // Young's modulus
     mprops_.nu = eptr->Read( Poissons_key_ );        // Poisson's ratio
     // friction coefficient computed from friction angle
     mprops_.mu = tan(degreesToRadians(eptr->Read( Friction_key_)));
     if ( Biot_key_.place != UNDEFINED ) mprops_.alpha = eptr->Read( Biot_key_ );
     else mprops_.alpha = (0.6);                      // Biot coefficient phi <= alpha <= 1
     mprops_.UCS = 0.1 * mprops_.E;                   // Hoek-Brown criterion inferred value
     mprops_.TS = tensileStrengthFromUCS_Griffith( mprops_.UCS, mprops_.mu );
     mprops_.pstar = 0.6 * mprops_.E;                 // about 5-7 times of the UCS of the rock
     mprops_.C = eptr->Read( Cohesion_key_ );         // cohesion = inherent shear strength

     // calculated using E and nu from relations between elastic moduli
     mprops_.G = (0.5 * mprops_.E)/(1 + mprops_.nu);
     mprops_.B = (2./3.) * mprops_.G * (1. + mprops_.nu)/(1. - 2. * mprops_.nu); // K_dry

     if ( Pressure_key_.place != UNDEFINED )
       eptr->PropertyValueAtBaryCenter( Pressure_key_, fluid_pressure_ );
   
  } // end InitializeMechanicalProperties




/**
    Loops over the element integration points and evaluates shear and tensile failure potential.
    Where failure occurred the variable 'failure' is set to 1.
    A distinction is made between tensile 'tensile failure' and shear failure 'failure' variables.
    
    @attention it is assumed that the tensile strength is about 1/10 of the cohesive strength of the rock
    
    @attention SKM 28/9/2014 - added case where stress is placed on the element.
*/
template<size_t dim>
void BrittleFailureVisitor<dim>::Visit( Element<dim>* e )
{
   InitializeInputProperties( e );
  
   // reading extra variables at integration points and calculating criteria
   if ( Stress_key_.place == ELEMENT_INTEGRATION_POINT ) {
        const size_t integration_points(e->FE()->IntegrationPoints());
        for ( size_t i=0U; i<integration_points; i++ )
          {
            // 1. reading extra input variables
            if ( Pressure_key_.place != UNDEFINED )
              fluid_pressure_ = e->PropertyValueAtIntegrationPoint( Pressure_key_, i );

            // 2. stress tensor and its invariants
            e->Read(i, Stress_key_, Cartesian_stress_ );
            StressInvariants  invars( Cartesian_stress_ );
            
            // 3. failure analysis
            FAILURE failure = BrittleFailure::Evaluate( mprops_, invars, fluid_pressure_() );

            // 4. store variables
            e->Store(i, Failure_key_, makeScalar(PLAIN,failure) );

            // 5. output if desired
            if ( verbose_ ) {
                cout <<"\nshear failure criterion F "<< endl;
                cout.width(35); cout <<"at element integration point ";
                cout << i << ": "<< failure << endl;
                cout.flush();
            }
         }
       return;
     }
  
     
  // cohesion, failure variables are placed on the element
  if ( Stress_key_.place == ELEMENT )
    {
        // 1. stress tensor, its invariants, and fluid pressure
        e->Read( Stress_key_, Cartesian_stress_ );
        StressInvariants  invars( Cartesian_stress_ );
        // NB: pore pressure has already been read by InitializeInputProperties

        // 2. failure analysis
        FAILURE failure = BrittleFailure::Evaluate( mprops_, invars, fluid_pressure_() );

        // 3. store variables
        e->Store( Failure_key_, makeScalar(PLAIN,failure) );

        // 10. output if desired
        if ( verbose_ ) {
            cout <<"\nshear failure criterion F "<< endl;
            cout.width(35); cout <<" in element "<< e->Idx();
            cout <<": "<< failure << endl;
            cout.flush();
        }
    }
  
} // end Visit(Element)

template class BrittleFailureVisitor<2U>;
template class BrittleFailureVisitor<3U>;

} // end csmp

