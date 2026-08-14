#include "MohrCoulombFailure.h"

using namespace std;

namespace csmp {

template<uint32_t dim>
MohrCoulombFailure<dim>::MohrCoulombFailure( const PropertyDatabase<dim>& p,
                                             double friction_angle ) 
      : Interrelation<dim>(p),
        STRESS( Interrelation<dim>::GlobalProperty("stress") ),
        MS( Interrelation<dim>::GlobalProperty("mean stress") ),
        COH( Interrelation<dim>::GlobalProperty("cohesion") ),
        CRIT( Interrelation<dim>::GlobalProperty("failure") ),
        Ts(1.0e+7) // 10 MPa
 {
    Interrelation<dim>::Name("MohrCoulombFailure");
    Interrelation<dim>::OutputCondition( CRIT, PLAIN );
    Interrelation<dim>::ResultProperty("failure");
    
    // converting the friction angle in 'degrees' to 'radians'
    double  degrees_to_radians( (2*3.14159265358979)/360. );
    
    phi = friction_angle * degrees_to_radians;
 }




/** Zienkiewicz II, p. 89, modified Mohr-Coulomb envelope with smooth boundaries.

Yielding can occur when F >= 0, the material is described only in terms of 
its friction angle and cohesion.
*/
template<uint32_t dim>
void MohrCoulombFailure<dim>::Calculate()
 {
    STRESS.AssignTo( ts_ );
    COH.AssignTo( ch_ );
    
    double  sm   = MeanStress( ts_ );
    double  t;
    double  sd   = DeviatoricStress( ts_, t );
    double  ta   = Theta( ts_, t );
    double  gt   = G_OfTheta( ta );
    
    double  F = sm * sin(phi) - ch_() * cos(phi) + sd / gt;
    
    CRIT = 1.0;
    
    cout <<"\nMohrCoulombFailure::Calculate: Friction Criterion F: "<< F << endl;
    cout <<"\nParameters mean-, deviatoric-stress, lode angle: ";
    cout << sm <<", "<< sd <<", "<< ta << endl;
    
    MS.AssignTo( ch_ );
    cout <<"\nactual mean stress: "<< ch_ << endl;

 } // end Calculate




/// calculate lode angle, Smith & Griffiths, p. 233, computed in radians
template<uint32_t dim>
double MohrCoulombFailure<dim>::Theta( const TensorVariable<dim>& ts, double t )
 {
    double J3 = std::numeric_limits<double>::quiet_NaN();
  
    if constexpr ( dim == 2 )
      {
         // 2D case not sure yet, search reference
         double sx  = (2*ts(0,0) - ts(1,1)) / 2.;
         double sy  = (2*ts(1,1) - ts(0,0)) / 2.;
         J3  = sx * sy + 2.* ts(0,1);
      }

    if constexpr ( dim == 3 )
      {
         double sx  = (2*ts(0,0) - ts(1,1) - ts(2,2)) / 3.;
         double sy  = (2*ts(1,1) - ts(2,2) - ts(0,0)) / 3.;
         double sz  = (2*ts(2,2) - ts(0,0) - ts(1,1)) / 3.;
         J3  = sx * sy * sz;
         J3 -= sx * (ts(1,2)*ts(1,2));
         J3 -= sy * (ts(0,2)*ts(0,2));
         J3 -= sz * (ts(0,1)*ts(0,1));
         J3 += 2. * ts(0,1) * ts(0,2) * ts(1,2);
      }
      
    return 1./3. * std::asin( (-3.* std::sqrt(6.) * J3)/(t*t*t) );
 }


/// Zienkiewitz II, p. 89
template<uint32_t dim>
double  MohrCoulombFailure<dim>::G_OfTheta( double theta ) 
 {
    double K(std::sin(phi));
    
    K = (3. - K) / (3. + K);

    return (2.*K) / ((1+K) - std::sin(3.*theta)*(1-K));
    
 } // end G_OfTheta




/// Smith & Griffiths, p. 233
template<uint32_t dim>
double  MohrCoulombFailure<dim>::MeanStress( const TensorVariable<dim>& ts ) 
 {
    // Smith & Griffiths, p. 233
    if ( dim == 2 ) return (ts(0,0) + ts(1,1)) / std::sqrt(3.);
    return (ts(0,0) + ts(1,1) + ts(2,2)) / std::sqrt(3.);
    
 } // end



/// Zienkiewitz II, p. 61 bottom, p. 62 top, sigma-dash = deviatoric stress
template<uint32_t dim>
double  MohrCoulombFailure<dim>::DeviatoricStress( const TensorVariable<dim>& ts, double& t ) 
 {
    // Smith & Griffiths, p. 236 (t is one of the stress invariants)
    if constexpr ( dim == 2U )
      {
         t = ((ts(0,0)-ts(1,1))*(ts(0,0)-ts(1,1))) + 3.*ts(0,1)*ts(0,1);
         t = std::sqrt(t) / std::sqrt(3.);
      }
    else
      {
         // assuming symmetric stress tensor
         t  = ((ts(0,0)-ts(1,1))*(ts(0,0)-ts(1,1)));
         t += ((ts(1,1)-ts(2,2))*(ts(1,1)-ts(2,2)));
         t += ((ts(2,2)-ts(0,0))*(ts(2,2)-ts(0,0)));
         t += 6.*ts(0,1)*ts(0,1) + 6.*ts(1,2)*ts(1,2) + 6.*ts(0,2)*ts(0,2);
         t  = std::sqrt(t) / std::sqrt(3.);
      }
    
    static_assert( dim != 1U, "MohrCoulombFailure<dim>::DeviatoricStress: does not work in 1D.");
    
    return t * std::sqrt(3./2.);
    
 } // end DeviatoricStress


//template class MohrCoulombFailure<1U>;
template class MohrCoulombFailure<2U>;
template class MohrCoulombFailure<3U>;

} // csmp


