#include "QuartzPrecipitation.h"

using namespace std;

namespace csmp {

template<size_t dim>
QuartzPrecipitation<dim>::QuartzPrecipitation( const PropertyDatabase<dim>& p, double64 time_increment ) 
      : Interrelation<dim>(p), dt(time_increment),
        F( Interrelation<dim>::GlobalProperty("velocity") ),
        DS( Interrelation<dim>::GlobalProperty("quartz solubility gradient") ),
        FD( Interrelation<dim>::GlobalProperty("fluid density") ),
        Q( Interrelation<dim>::GlobalProperty("quartz") )
 {
    Interrelation<dim>::Name("QuartzPrecipitation");
    Interrelation<dim>::OutputCondition( Q, PLAIN );
    Interrelation<dim>::ResultProperty("quartz");
 }


/** Calculating how much quartz gets precipitated

1. calculate quartz solubility at the nodes (done in QuartzSolubility.h)
2. get vector properties flux and solubility gradient in
   the element
3. project solubility gradient on fluid flux
4. compute product of mass flux, projected solubility gradient, and
   time. This product gives the amount of quartz that is dissolved
   or precipitated
   */
template<size_t dim>
void QuartzPrecipitation<dim>::Calculate()
 {
     // 1. getting input data for the calculation 
     // -----------------------------------------
     F.AssignTo(  flux );
     DS.AssignTo( dSdz );
     FD.AssignTo( rho );
     
     // 2. calculating the pressure and temperature gradients in the
     //    direction of fluid flow a_b = |a| . cos( angle between a & b) 
     // ----------------------------------------------------------------
     // quartz solubility

     if ( flux(0) != 0.0 || flux(1) != 0.0 )
       {
          // 1. time-integrated mass flux
          // ----------------------------
          mass_flux = flux.Length() * rho() * dt;

          // 2. gradient of quartz solubility in the direction of fluid flow
          // ---------------------------------------------------------------
          dSpr  = dSdz.ProjectOnto( flux );
          if ( dSpr(0) < 0.0 || dSpr(1) < 0.0 ) sign = -1.0;
          else                                  sign =  1.0;
          grad_S = dSpr.Length() * sign;

          // 3. compute quartz precipitation/dissolution in kg
          // -------------------------------------------------
          Q += mass_flux * -grad_S;
       }
 } 

template class QuartzPrecipitation<1U>;
template class QuartzPrecipitation<2U>;
template class QuartzPrecipitation<3U>;

} // csmp
