#include "TransportStepSize.h"

using namespace std;

namespace csmp {

template<uint32_t dim>
TransportStepSize<dim>::TransportStepSize( const PropertyDatabase<dim>& p, 
                                              const char* n_velo,
                                              const char* e_velo, 
                                              const char* i_radius ) 
      : Interrelation<dim>(p),
        NV( Interrelation<dim>::GlobalProperty(n_velo) ),
        EV( Interrelation<dim>::GlobalProperty(e_velo) ),
        IR( Interrelation<dim>::GlobalProperty(i_radius) ),
        advection_increment(0.0),
        called(false)
 {
    Interrelation<dim>::Name("TransportStepSize");
    Interrelation<dim>::OutputCondition( IR, PLAIN );
    // dummy, since nothings output
    Interrelation<dim>::ResultProperty(i_radius);
 }


template<uint32_t dim>
TransportStepSize<dim>::~TransportStepSize() {}
  
  
    
template<uint32_t dim>
double TransportStepSize<dim>::AdvectionTimeIncrement() const
 {
    if ( advection_increment == 0.0 )
      {
          cout <<"\nTransportStepSize::AdvectionTimeIncrement: ";
          cout <<" Error: Current increment is zero, ";
          cout <<"interrelation may not have been applied yet."<< endl;
      }
    return advection_increment;
 }


template<uint32_t dim>
void TransportStepSize<dim>::Calculate()
 {
    // 1. Calculate peak velocity in the element
    // -----------------------------------------
    NV.AssignTo( nv );
    EV.AssignTo( ev );
    n_velo = nv.Length();
    e_velo = ev.Length();
    if ( n_velo >= e_velo ) scalar_velocity = n_velo;
    else                    scalar_velocity = e_velo;
    
    // 2. Calculate maximum permitted transport increment in element
    //    (twice the inner radius is the diameter of the inner circle)
    //    - scalar velocity should be less than 100 m s-1
    // ---------------------------------------------------------------
    assert( scalar_velocity >= 0.0 );
    IR.AssignTo( ir );
    if (scalar_velocity == 0.0 )
      e_increment = (ir()*2.0) / 1.0e-30;
    else
      e_increment = (ir()*2.0) / scalar_velocity;    
//    e_increment = ir() / scalar_velocity;    
    
    // 3. Compare against thus far permitted increment, if smaller 
    //    the size of the permitted increment is decreased.
    // ----------------------------------------------------
    if ( !called )
      {
         advection_increment = e_increment;
         called              = true;
      }
    else if ( e_increment < advection_increment ) advection_increment = e_increment;
 } 

template class TransportStepSize<1U>;
template class TransportStepSize<2U>;
template class TransportStepSize<3U>;

} // csmp
