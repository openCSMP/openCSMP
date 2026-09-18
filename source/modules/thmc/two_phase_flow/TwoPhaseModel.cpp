// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "TwoPhaseModel.h"
#include "CSMP_physical_constants.h"

using namespace std;

namespace csmp {


template<uint32_t dim>
TwoPhaseModel<dim>::TwoPhaseModel()
 : mun_(2.0e-3), muw_(1.6e-3),
   rhn_(800.), rhw_(1000.),
   sat_(std::numeric_limits<double>::quiet_NaN()),
   seff_(std::numeric_limits<double>::quiet_NaN()),
   swr_(0.), snr_(0.),
   k_(std::numeric_limits<double>::quiet_NaN()),
   K_(PLAIN,std::numeric_limits<double>::quiet_NaN()),
   ift_(0.05), // N m-1 Danesh (2003), p. 292
   acc_gravity_(ACC_GRAVITY),
   tolerance_(numeric_limits<double>::epsilon()),
   MAX_CAPILLARY_PRESSURE_(4e7),       // maximum tensile strength of a rock
   MAX_CAPILLARY_PRESSURE_SLOPE_(1e7), // maximum slope of the capillary pressure curve
   interpolate_fluid_properties_(false),
   sw_ro_mu_placement_(true)
 {
 }
 
 
template<uint32_t dim>
TwoPhaseModel<dim>::TwoPhaseModel( const PropertyDatabase<dim>& database,
                                      double viscosity_nw, double viscosity_w,
                                      double density_nw, double density_w,
                                      const char* kkk,
                                      const char* sat,
                                      const char* snr, 
                                      const char* swr,
                                      const bool sw_ro_mu_placement )
 : perm_key_(database.StorageKey(kkk)),
   sat_key_(database.StorageKey(sat)),
   snr_key_(database.StorageKey(snr)),
   swr_key_(database.StorageKey(swr)),
   mun_(viscosity_nw), muw_(viscosity_w),
   rhn_(density_nw), rhw_(density_w),
   lt_key_(database.StorageKey("total mobility")),
   sat_(std::numeric_limits<double>::quiet_NaN()),
   seff_(std::numeric_limits<double>::quiet_NaN()),
   swr_(0.), snr_(0.),
   k_(std::numeric_limits<double>::quiet_NaN()),
   K_(PLAIN,std::numeric_limits<double>::quiet_NaN()),
   ift_(0.05), // N m-1 Danesh (2003), p. 292
   acc_gravity_(ACC_GRAVITY),
   tolerance_(1.0e-17),
   MAX_CAPILLARY_PRESSURE_(4e7), // maximum strength of a rock
   MAX_CAPILLARY_PRESSURE_SLOPE_(1e7),
   interpolate_fluid_properties_(false),
   sw_ro_mu_placement_(sw_ro_mu_placement),
   tensor_permeability_(false)
 {

    if ( perm_key_.place != ELEMENT || (perm_key_.type != SCALAR && perm_key_.type != TENSOR) )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseModel(constructor)",
                      kkk, " variable must be scalar or tensor element property" );

    if( perm_key_.type == TENSOR)
        tensor_permeability_ = true;

    if ( snr_key_.place != ELEMENT || snr_key_.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseModel(constructor)",
                      snr, " variable must be scalar element property" );

    if ( swr_key_.place != ELEMENT || swr_key_.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseModel(constructor)",
                      swr, " variable must be scalar element property" );

    if (sw_ro_mu_placement_)
    {
        if ( sat_key_.place != NODE || sat_key_.type != SCALAR )
          throw csmp::Exception( FATAL_ERROR, "TwoPhaseModel(constructor)",
                          sat, " variable must be scalar node property" );

    }
    else
    {
        if ( sat_key_.place != ELEMENT || sat_key_.type != SCALAR )
          throw csmp::Exception( FATAL_ERROR, "TwoPhaseModel(constructor)",
                          sat, " variable must be scalar element property" );

    }

}



template<uint32_t dim>
TwoPhaseModel<dim>::TwoPhaseModel( const PropertyDatabase<dim>& database,
                                  const char* kkk,
                                  const char* mun,
                                  const char* muw,
                                  const char* rhn,
                                  const char* rhw,
                                  const char* sat,
                                  const char* snr, 
                                  const char* swr,
                                  const bool sw_ro_mu_placement )
 : perm_key_(database.StorageKey(kkk)),
   mun_key_(database.StorageKey(mun)),
   muw_key_(database.StorageKey(muw)),
   rhn_key_(database.StorageKey(rhn)),
   rhw_key_(database.StorageKey(rhw)),
   lt_key_(database.StorageKey("total mobility")),
   sat_key_(database.StorageKey(sat)),
   snr_key_(database.StorageKey(snr)),
   swr_key_(database.StorageKey(swr)),
   sat_(std::numeric_limits<double>::quiet_NaN()),
   seff_(std::numeric_limits<double>::quiet_NaN()),
   swr_(0.), snr_(0.),
   k_(std::numeric_limits<double>::quiet_NaN()),
   ift_(0.05), // N m-1 Danesh (2003), p. 292
   mun_(std::numeric_limits<double>::quiet_NaN()),
   muw_(std::numeric_limits<double>::quiet_NaN()),
   rhn_(std::numeric_limits<double>::quiet_NaN()),
   rhw_(std::numeric_limits<double>::quiet_NaN()),
   acc_gravity_(ACC_GRAVITY),
   tolerance_(1.0e-17),
   MAX_CAPILLARY_PRESSURE_(4e7), // maximum strength of a rock
   MAX_CAPILLARY_PRESSURE_SLOPE_(1e7),
   interpolate_fluid_properties_(true),
   sw_ro_mu_placement_(sw_ro_mu_placement),
   tensor_permeability_(false)
 {
    if ( perm_key_.place != ELEMENT || (perm_key_.type != SCALAR && perm_key_.type != TENSOR) )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseModel(constructor)", 
                      kkk, " variable must be scalar or tensor element property" );

    if( perm_key_.type == TENSOR)
        tensor_permeability_ = true;
                      
    if ( snr_key_.place != ELEMENT || snr_key_.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseModel(constructor)", 
                      snr, " variable must be scalar element property" ); 
                      
    if ( swr_key_.place != ELEMENT || swr_key_.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TwoPhaseModel(constructor)", 
                      swr, " variable must be scalar element property" );
    
    if (sw_ro_mu_placement_)
    {
        if ( mun_key_.place != NODE || mun_key_.type != SCALAR )
          throw csmp::Exception( FATAL_ERROR, "TwoPhaseModel(constructor)", 
                          mun, " variable must be scalar node property" ); 
                          
        if ( muw_key_.place != NODE || muw_key_.type != SCALAR )
          throw csmp::Exception( FATAL_ERROR, "TwoPhaseModel(constructor)", 
                          muw, " variable must be scalar node property" ); 
                          
        if ( rhn_key_.place != NODE || rhn_key_.type != SCALAR )
          throw csmp::Exception( FATAL_ERROR, "TwoPhaseModel(constructor)", 
                          rhn, " variable must be scalar node property" ); 
                          
        if ( rhw_key_.place != NODE || rhw_key_.type != SCALAR )
          throw csmp::Exception( FATAL_ERROR, "TwoPhaseModel(constructor)", 
                          rhw, " variable must be scalar node property" ); 
                          
        if ( sat_key_.place != NODE || sat_key_.type != SCALAR )
          throw csmp::Exception( FATAL_ERROR, "TwoPhaseModel(constructor)", 
                          sat, " variable must be scalar node property" ); 
        
    }
    else
    {
        if ( mun_key_.place != ELEMENT || mun_key_.type != SCALAR )
          throw csmp::Exception( FATAL_ERROR, "TwoPhaseModel(constructor)", 
                          mun, " variable must be scalar element property" ); 
                          
        if ( muw_key_.place != ELEMENT || muw_key_.type != SCALAR )
          throw csmp::Exception( FATAL_ERROR, "TwoPhaseModel(constructor)", 
                          muw, " variable must be scalar element property" ); 
                          
        if ( rhn_key_.place != ELEMENT || rhn_key_.type != SCALAR )
          throw csmp::Exception( FATAL_ERROR, "TwoPhaseModel(constructor)", 
                          rhn, " variable must be scalar element property" ); 
                          
        if ( rhw_key_.place != ELEMENT || rhw_key_.type != SCALAR )
          throw csmp::Exception( FATAL_ERROR, "TwoPhaseModel(constructor)", 
                          rhw, " variable must be scalar element property" ); 
                          
        if ( sat_key_.place != ELEMENT || sat_key_.type != SCALAR )
          throw csmp::Exception( FATAL_ERROR, "TwoPhaseModel(constructor)", 
                          sat, " variable must be scalar element property" ); 
    
    }

 }
 
 
 
template<uint32_t dim>
TwoPhaseModel<dim>::~TwoPhaseModel()
 {
 }



template<uint32_t dim>
double csmp::TwoPhaseModel<dim>::Swr() const
  {
    return swr_;
  }

template<uint32_t dim>
double csmp::TwoPhaseModel<dim>::Snr() const
  {
    return snr_;
  }

template<uint32_t dim>
void csmp::TwoPhaseModel<dim>::Swr( double wettingResidual )
  {
     swr_ = wettingResidual;
  }

template<uint32_t dim>
void csmp::TwoPhaseModel<dim>::Snr( double wettingNonResidual )
  {
    snr_ = wettingNonResidual;
  }


template<uint32_t dim>
void TwoPhaseModel<dim>::SaturationWettingPhase( double s_wetting )
{ sat_ = s_wetting; }

template<uint32_t dim>
void TwoPhaseModel<dim>::ViscosityNonWettingPhase( double visc )
{ mun_ = visc; }

template<uint32_t dim>
void TwoPhaseModel<dim>::ViscosityWettingPhase( double visc )
{ muw_ = visc; }

template<uint32_t dim>
void TwoPhaseModel<dim>::DensityNonWettingPhase( double dens )
{ rhn_ = dens; }

template<uint32_t dim>
void TwoPhaseModel<dim>::DensityWettingPhase( double dens )
{ rhw_ = dens; }


template<uint32_t dim>
double TwoPhaseModel<dim>::Permeability() const
{ return k_; }

template<uint32_t dim>
void TwoPhaseModel<dim>::Permeability( double permeability )
  { k_ = permeability; }

template<uint32_t dim>
TensorVariable<dim> TwoPhaseModel<dim>::TensorPermeability() const
  { return K_; }

template<uint32_t dim>
void TwoPhaseModel<dim>::TensorPermeability( TensorVariable<dim> permeability )
  { K_ = permeability; }

template<uint32_t dim>
double TwoPhaseModel<dim>::ViscosityNonWettingPhase() const
{ return mun_; }

template<uint32_t dim>
double TwoPhaseModel<dim>::ViscosityWettingPhase() const
{ return muw_; }

template<uint32_t dim>
double TwoPhaseModel<dim>::DensityNonWettingPhase() const
{ return rhn_; }

template<uint32_t dim>
double TwoPhaseModel<dim>::DensityWettingPhase() const
{ return rhw_; }

template<uint32_t dim>
double TwoPhaseModel<dim>::MaxCapillaryPressure( uint32_t phase ) const
{ 
   assert( phase == 1U or phase == 2U );
   if ( phase == 2U ) return MAX_CAPILLARY_PRESSURE_;
   return static_cast<double>(0.); 
}


template<uint32_t dim>
double TwoPhaseModel<dim>::Saturation( uint32_t phase ) const
{ 
   assert( phase == 1U or phase == 2U );
   if ( phase == 2U ) return static_cast<double>(1.) - sat_;
   return sat_;
}


template<uint32_t dim>
double TwoPhaseModel<dim>::MobilityPhase( uint32_t phase ) const
 {
    assert( phase == 1U or phase == 2U );
    if ( phase == 1U ) return krw_Phase() / muw_;
    return krn_Phase() / mun_;
 }


template<uint32_t dim>
double TwoPhaseModel<dim>::TotalMobilityMultiplier() const 
 {
    return krn_Phase() / mun_ + krw_Phase() / muw_;
 }

template<uint32_t dim>
double TwoPhaseModel<dim>::TotalMobility() const 
 {
    return k_ * TotalMobilityMultiplier();
 }

template<uint32_t dim>
double TwoPhaseModel<dim>::ViscosityRatio() const
{
  return mun_ / muw_;
}

/** 
    calculates lambda_t (total mobility) = sum of phase mobilities

    @todo SKM averages tensor properties, but should create tensor mobilities in stead
*/
template<uint32_t dim>
void TwoPhaseModel<dim>::Initialize( const Element<dim>& e )
 {
    swr_ = e.Read( swr_key_ );
    snr_ = e.Read( snr_key_ );

    if( tensor_permeability_){
        e.Read( perm_key_, K_);
        // TODO: Skm: fix these unwanted averages of the tensor k
        k_ = K_.Trace()/static_cast<double>(dim);
    }else{
        k_ = e.Read( perm_key_ );
        K_.operator=( VectorVariable<dim>(PLAIN, k_ ) );
    }
    
    // if properties are discretized on element
    if (!sw_ro_mu_placement_) {
        sat_ = e.Read( sat_key_ );
        mun_ = e.Read( mun_key_ );
        muw_ = e.Read( muw_key_ );
        rhn_ = e.Read( rhn_key_ );
        rhw_ = e.Read( rhw_key_ );
    }
    
 }
 
/**
    non const version which allows to set parameters on element 
*/
template<uint32_t dim>
void TwoPhaseModel<dim>::InitializeAndStore( Element<dim>& e )
 {
    swr_ = e.Read( swr_key_ );
    snr_ = e.Read( snr_key_ );

    if( tensor_permeability_){
        e.Read( perm_key_, K_);
        // TODO: Skm: fix these unwanted averages of the tensor k
        k_ = K_.Trace()/static_cast<double>(dim);
    }else{
        k_ = e.Read( perm_key_ );
        K_.operator=( VectorVariable<dim>(PLAIN, k_ ) );
    }
   
    // if properties are discretized on element
    if (!sw_ro_mu_placement_) {
        sat_ = e.Read( sat_key_ );
        mun_ = e.Read( mun_key_ );
        muw_ = e.Read( muw_key_ );
        rhn_ = e.Read( rhn_key_ );
        rhw_ = e.Read( rhw_key_ );
    }
   
 }

 
 
 
template<uint32_t dim>
void TwoPhaseModel<dim>::InitializeForBaryCenter( const Element<dim>& e )
 {
    if (sw_ro_mu_placement_) {
        e.N_AtBaryCenter( e.FE()->NRST );
        InterpolateNodeProperties( e );
    }
 }  


template<uint32_t dim>
void TwoPhaseModel<dim>::InitializeForNode( const Element<dim>& e,
                                                   uint32_t fem_node )
 {
    if (sw_ro_mu_placement_) {
        sat_ = e.N(fem_node)->Read( sat_key_ );
        mun_ = e.N(fem_node)->Read( mun_key_ );
        muw_ = e.N(fem_node)->Read( muw_key_ );
        rhn_ = e.N(fem_node)->Read( rhn_key_ );
        rhw_ = e.N(fem_node)->Read( rhw_key_ );
    }

 }  


template<uint32_t dim>
void TwoPhaseModel<dim>::InitializeForIntegrationPoint( uint32_t ip,
                                                               const Element<dim>& e )
 {
    if (sw_ro_mu_placement_) {
        e.N_AtIntegrationPoint( ip, e.FE()->NRST );
        InterpolateNodeProperties( e );
    }
    
 } // end 
 
 
 
 
template<uint32_t dim>
void TwoPhaseModel<dim>::InitializeForFacetIntegrationPoint( uint32_t facet, uint32_t ip,
                                                                    const Element<dim>& e )
 {
    if (sw_ro_mu_placement_) {
        (e).N_AtFacetIntegrationPoint( facet, ip );  
        InterpolateNodeProperties( e );
    }

 }  



template<uint32_t dim>
void TwoPhaseModel<dim>::InitializeForSectorIntegrationPoint( uint32_t sector, uint32_t ip,
                                                                     const Element<dim>& e )
 {
    if (sw_ro_mu_placement_) {
        (e).N_AtSectorIntegrationPoint( sector, ip );  
        InterpolateNodeProperties( e );
    }

 } 





 

/**
 
The effective saturation is always computed w.r.t. the wetting phase. 
(compare with Helmig 1997, p. 71 and note that phase 1 is the wetting 
phase; note also that Initialize() must be called first to get input 
parameters like residual saturations).  
*/
template<uint32_t dim>
double TwoPhaseModel<dim>::EffectiveSaturation() const 
 {
    return seff_ = std::min( std::max( (sat_ - swr_) / (1. - swr_ - snr_), 0. ), 1. );
    
 } // end EffectiveSaturation


/// return wetting-phase saturation based on effective saturation
template<uint32_t dim>
double TwoPhaseModel<dim>::SeffToSw() const
 {
    return sat_ = seff_*(1. - swr_ - snr_) + swr_;

 }

/// return wetting-phase saturation based on effective saturation ( useful for numerical calculations )
template<uint32_t dim>
double TwoPhaseModel<dim>::SeffToSw( double seff ) const
 {
    return sat_ = seff*(1. - swr_ - snr_) + swr_;

 }

/**
 
Computes the fractional flow of the wetting (phase=1) and non-wetting
(phase=2) phases using the relative k's. and viscosities. Note that 
Initialize() must be called first.  
*/
template<uint32_t dim>
double TwoPhaseModel<dim>::f_Phase( uint32_t phase ) const
 {
    assert( phase == 1U or phase == 2U );
      if ( phase == 1U )
        return (krw_Phase() / muw_) / TotalMobilityMultiplier();

      return   (krn_Phase() / mun_) / TotalMobilityMultiplier();
 }

      
                                  
/**
Computes G = lamdba_w * lambda_n / (lambda_w + lambda_n), cf., van Duijn 
and de Neef (1998). Note that Initialize() must be called first.  
 */
template<uint32_t dim>
double TwoPhaseModel<dim>::G() const
 {
    const double lambda_w(krw_Phase() / muw_),
             lambda_n(krn_Phase() / mun_);
             
	return (lambda_w * lambda_n) / (lambda_w + lambda_n);
 }
 



/**
Returns the diffusion coefficient for the phase of interest. If not 
overloaeded, the hydraulic conductivity is returned.  
*/
template<uint32_t dim>
double TwoPhaseModel<dim>::DiffusionMultiplier( size_t phase ) const
{
     assert( phase == 1U or phase == 2U );
     return k_ / ( (phase==1u) ? muw_ : mun_ );
} 
 


/**
See Helmig, 1997, p. 108, eqn. 3.74, term 1. This takes into account the
permeability in direction of flow  x  lambda_overbar  x pc-gradient.  
*/
template<uint32_t dim>
double TwoPhaseModel<dim>::CapillaryDiffusionMultiplier( ) const
{
   return k_ * G() * dpcds_Phase( );
} 



/**
 
If not overloaded, this returns the derivative of the fractional flow
function at the current saturation of the wetting phase (see Helmig, 1997, 
p. 108, eqn. 3.74, term 2 (first part).  
*/
template<uint32_t dim>
double TwoPhaseModel<dim>::AdvectionMultiplier( ) const
 {
    return dfds();
 } 

template<uint32_t dim>
double TwoPhaseModel<dim>::GravityTerm() const
{
  // note that the projected gravity acts opposite the y-axis
  const double k_g_drho = k_ * -acc_gravity_ * (rhw_ - rhn_);

    // economizing the calculation
    if ( std::fabs(k_g_drho) < tolerance_ ) return static_cast<double>(0.);

    // else compute result using G saturation derivative
    return k_g_drho;
}

/** See Sebastian Geiger's thesis (2004), closed form.
*/
template<uint32_t dim>
double TwoPhaseModel<dim>::GravityMultiplier_G( ) const
{
  // note that the projected gravity acts opposite the y-axis
  const double k_g_drho = k_ * -acc_gravity_ * (rhw_ - rhn_);
	
	// economizing the calculation
    if ( std::fabs(k_g_drho) < tolerance_ ) return static_cast<double>(0.);

	// else compute result using G saturation derivative
	return k_g_drho * G();	
} 


 
/**
 
Computes multiplier for advection gravity coefficient. The divergence 
lamda_ div k g (rhw-rhn) must be dealt with separately, see Helmig, 1997, 
p. 108, eqn. 3.74, term 2 (second part).  
*/
template<uint32_t dim>
double TwoPhaseModel<dim>::GravityMultiplier_dGds( ) const
{
    // SKM flow equations worked out with Adrian
    const double k_g_drho = k_ * -acc_gravity_ * (rhw_ - rhn_);
	
	// economizing the calculation
    if ( std::fabs(k_g_drho) < tolerance_ ) return static_cast<double>(0.);

	// else compute result using G saturation derivative
    return k_g_drho * dGds(  );
} 

template<uint32_t dim>
csmp::Index TwoPhaseModel<dim>::WettingPhaseSaturationKey() const
{
  return sat_key_;
}

template<uint32_t dim>
csmp::Index TwoPhaseModel<dim>::NonWettingPhaseSaturationKey() const
{
    throw csmp::Exception( ERROR, "TwoPhaseModel<dim>::NonWettingPhaseSaturationKey()",
                           "It appears someone decided this should be the irreducible non-wet. phase saturation. Wrong." );
    return snr_key_;
}

template<uint32_t dim>
csmp::Index TwoPhaseModel<dim>::TotalMobilityKey() const
{
  return lt_key_;
}



template<uint32_t dim>
void TwoPhaseModel<dim>::InterpolateNodeProperties( const Element<dim>& e )
 {
    if ( interpolate_fluid_properties_ ) {
         mun_ = muw_ = rhn_ = rhw_ = sat_ = static_cast<double>(0.);
         for ( auto i{0U}; i<e.Nodes(); i++ ) {
              sat_ += e.FE()->NRST[i] * e.N(i)->Read( sat_key_ );
              mun_ += e.FE()->NRST[i] * e.N(i)->Read( mun_key_ );
              muw_ += e.FE()->NRST[i] * e.N(i)->Read( muw_key_ );
              rhn_ += e.FE()->NRST[i] * e.N(i)->Read( rhn_key_ );
              rhw_ += e.FE()->NRST[i] * e.N(i)->Read( rhw_key_ );
           }
      }
    else {
         sat_ = static_cast<double>(0.);
         for ( auto i{0U}; i<e.Nodes(); i++ )
           sat_ += e.FE()->NRST[i] * e.N(i)->Read( sat_key_ );
      }
 }


/// relative permeabilities
template<uint32_t dim>
double TwoPhaseModel<dim>::krw_Phase() const
 {
    cout <<"\nTwoPhaseModel<dim>::krw_Phase (base class): ";
    cout <<"This method needs to be defined in this subclass to achieve desired functionality."<< endl;
    throw logic_error("TwoPhaseModel<dim>::krw_Phase: Method not defined in subclass");
    return std::numeric_limits<double>::signaling_NaN();
 }



template<uint32_t dim>
double TwoPhaseModel<dim>::krn_Phase() const
 {
    cout <<"\nTwoPhaseModel<"<< dim <<">::krn_Phase (base class): ";
    cout <<"This method needs to be defined in this subclass to achieve desired functionality."<< endl;
    throw logic_error("TwoPhaseModel<dim>::krn_Phase: Method not defined in subclass");
    return std::numeric_limits<double>::signaling_NaN();
 }

/// derivatives of relative permeabilities
template<uint32_t dim>
double TwoPhaseModel<dim>::dkrwds_Phase() const
 {
    return dkrwds_numerical();
 }



template<uint32_t dim>
double TwoPhaseModel<dim>::dkrnds_Phase() const
 {
    return dkrnds_numerical();
 }

/// capillary pressure
template<uint32_t dim>
double TwoPhaseModel<dim>::pc_Phase() const
 {
    cout <<"\nTwoPhaseModel<"<< dim <<">::pc_Phase (base class): ";
    cout <<"This method needs to be defined in this subclass to achieve desired functionality."<< endl;
    throw logic_error("TwoPhaseModel<dim>::pc_Phase: Method not defined in subclass");
    return std::numeric_limits<double>::signaling_NaN();
 }

/// inverse capillary pressure function
template<uint32_t dim>
double TwoPhaseModel<dim>::Sw_Phase(double ) const
 {
    cout <<"\nTwoPhaseModel<"<< dim <<">::SwFromPc(base class): ";
    cout <<"This method needs to be defined in this subclass to achieve desired functionality."<< endl;
    throw logic_error("TwoPhaseModel<dim>::SwFromPc: Method not defined in subclass");
    return std::numeric_limits<double>::signaling_NaN();
 }

    /// capillary pressure derivatives
template<uint32_t dim>
double TwoPhaseModel<dim>::dpcds_Phase( ) const
 {
    return dpcds_numerical();
 }

/// derivatives of inverse capillary pressure function
template<uint32_t dim>
double TwoPhaseModel<dim>::dsdpc_Phase( double ) const
 {
    cout <<"\nTwoPhaseModel<"<<  dim <<">::dpcdsw_Phase (base class): ";
    cout <<"This method needs to be defined in this subclass to achieve desired functionality."<< endl;
    throw logic_error("TwoPhaseModel<dim>::dpcdsw_Phase: Method not defined in subclass");
    return std::numeric_limits<double>::signaling_NaN();
 }

// linearized diffusion multiplier for large-timestep calculations
template<uint32_t dim>
double TwoPhaseModel<dim>::DiffusionCharacteristic( size_t ) const
{
cout <<"\nTwoPhaseModel<"<<  dim <<">::DiffusionCharacteristic (base class): ";
cout <<"This method needs to be defined in this subclass to achieve desired functionality."<< endl;
return std::numeric_limits<double>::signaling_NaN();
}


/// derivative of wetting phase mobility
template<uint32_t dim>
double TwoPhaseModel<dim>::dlwds() const
 {
    return dkrwds_Phase() / muw_;

 }

/// derivative of non-wetting phase mobility
template<uint32_t dim>
double TwoPhaseModel<dim>::dlnds() const
 {
    return dkrnds_Phase() / mun_;

 }

/// derivative of fractional flow function (advection multipliers)
template<uint32_t dim>
double TwoPhaseModel<dim>::dfds() const
 {

    //if ( ( seff_ < static_cast<double>(0.) ) || ( seff_ > static_cast<double>(1.) ) )
    //    return static_cast<double>(0.);

    const double lw  = krw_Phase() / muw_;
    const double ln  = krn_Phase() / mun_;
    const double lt  = lw + ln;
    const double lt2 = lt*lt;

    const double dlwds = dkrwds_Phase() / muw_;
    const double dlnds = dkrnds_Phase() / mun_;

    return ( dlwds*ln - dlnds*lw )/lt2;

 }
       

     
/// derivatives of gravitational flow
template<uint32_t dim>
double TwoPhaseModel<dim>::dGds( ) const
 {

    //if ( ( seff_ < static_cast<double>(0.) ) || ( seff_ > static_cast<double>(1.) ) )
    //    return static_cast<double>(0.);

    const double lw  = krw_Phase() / muw_;
    const double ln  = krn_Phase() / mun_;
    const double lt  = lw + ln;
    const double lt2 = lt*lt;
    const double ln2 = ln*ln;
    const double lw2 = lw*lw;


    const double dlwds = dkrwds_Phase() / muw_;
    const double dlnds = dkrnds_Phase() / mun_;

    return ( dlwds*ln2 + dlnds*lw2 )/lt2;
 }



/// Numerical derivatives

/// derivative of wetting phase mobility
template<uint32_t dim>
double TwoPhaseModel<dim>::dlwds_numerical( double h ) const
 {
    return dkrwds_numerical( h ) / muw_;

 }

/// derivative of non-wetting phase mobility
template<uint32_t dim>
double TwoPhaseModel<dim>::dlnds_numerical( double h) const
 {
    return dkrnds_numerical( h ) / mun_;

 }

template<uint32_t dim>
double TwoPhaseModel<dim>::dkrwds_numerical( double h ) const
{

  //if ( seff_ < 0.0 || seff_ > 1.0 )
  //      return static_cast<double>(0.);

  const double dSedSw( 1.0/ (1.0 - swr_ - snr_ ) );

  if( seff_ < 0.+h )
      return (krw_at( seff_ + h ) - krw_at( seff_ ) ) / h * dSedSw;
  if( seff_ > 1.-h )
      return ( krw_at( seff_ ) - krw_at( seff_ - h ) ) / h * dSedSw;

  return ( krw_at( seff_+h ) - krw_at( seff_-h ) )/ (2.0*h) * dSedSw;

}

template<uint32_t dim>
double TwoPhaseModel<dim>::dkrnds_numerical( double h ) const
{
  //if ( seff_ < 0.0 || seff_ > 1.0 )
  //      return static_cast<double>(0.);

  const double dSedSw( 1.0/ (1.0 - swr_ - snr_ ) );

  if( seff_ < 0.+h )
      return ( krn_at( seff_ + h ) - krn_at( seff_ ) ) / h * dSedSw;
  if( seff_ > 1.-h )
      return ( krn_at( seff_ ) - krn_at( seff_ - h ) ) / h * dSedSw;

  return ( krn_at( seff_+h ) - krn_at( seff_-h ) )/ (2.0*h) * dSedSw;
}



template<uint32_t dim>
double TwoPhaseModel<dim>::dfds_numerical( double h) const
{

  //if ( seff_ < 0.0 || seff_ > 1.0 )
  //    return static_cast<double>(0.);

  /*
  // first version: direct differentiation
  const double dSedSw( 1.0/ (1.0 - swr_ - snr_ ) );

  if( seff_ < 0.+h )
      return ( fw_at( seff_ + h ) - fw_at( seff_ ) ) / h * dSedSw;
  if( seff_ > 1.-h )
      return ( fw_at( seff_ ) - fw_at( seff_ - h ) ) / h * dSedSw;

  return ( fw_at( seff_+h ) - fw_at( seff_-h ) )/ (2.0*h)* dSedSw;

  */

  ///*
  // second version: mixed analytical and numerical differentiation
  const double lw  = krw_Phase() / muw_;
  const double ln  = krn_Phase() / mun_;
  const double lt  = lw + ln;
  const double lt2 = lt*lt;

  const double dlwds = dkrwds_numerical( h )/muw_;
  const double dlnds = dkrnds_numerical( h )/mun_;

  return ( dlwds*ln - dlnds*lw )/lt2;
  //*/

}

template<uint32_t dim>
double TwoPhaseModel<dim>::dGds_numerical( double h ) const
{
  //if ( seff_ < 0.0 || seff_ > 1.0 )
  //      return static_cast<double>(0.);

  /*
  // first version: direct differentiation
  const double dSedSw( 1.0/ (1.0 - swr_ - snr_ ) );
  h = 0.0000001;
  if( seff_ < 0.+h )
      return ( G_at( seff_ + h ) - G_at( seff_ ) ) / h * dSedSw;
  if( seff_ > 1.-h )
      return ( G_at( seff_ ) - G_at( seff_ - h ) ) / h * dSedSw;


  return ( G_at( seff_+h ) - G_at( seff_-h ) )/ (2.0*h) * dSedSw;
  */

  // second version: mixed analytical and numerical differentiation
  const double lw  = krw_Phase() / muw_;
  const double ln  = krn_Phase() / mun_;
  const double lt  = lw + ln;
  const double lt2 = lt*lt;
  const double ln2 = ln*ln;
  const double lw2 = lw*lw;


  const double dlwds = dkrwds_numerical( h )/muw_;
  const double dlnds = dkrnds_numerical( h )/mun_;

  return ( dlwds*ln2 + dlnds*lw2 )/lt2;

}

template<uint32_t dim>
double TwoPhaseModel<dim>::dpcds_numerical( double h) const
{

  //if ( seff_ < 0.0 || seff_ > 1.0 )
  //    return static_cast<double>(0.);

  const double dSedSw( 1.0/ (1.0 - swr_ - snr_ ) );
  double dpcds;

  if( seff_ < 0.+h )
      dpcds = ( pc_at( seff_ + h ) - pc_at( seff_ ) ) / h * dSedSw;
  else if( seff_ > 1.-h )
      dpcds = ( pc_at( seff_ ) - pc_at( seff_ - h ) ) / h * dSedSw;
  else
      dpcds = ( pc_at( seff_+h ) - pc_at( seff_-h ) )/ (2.0*h) * dSedSw;

  //if( dpcds < -MAX_CAPILLARY_PRESSURE_SLOPE_ )
  //    return -MAX_CAPILLARY_PRESSURE_SLOPE_;

  return dpcds;

}


template<uint32_t dim>
double TwoPhaseModel<dim>::MaxFractionalFlowDerivative() const
 {
    double speed, height;
    ShockSpeedHeight( speed, height );
    return speed;
 }

    /// linearized fractional flow derivative
template<uint32_t dim>
double TwoPhaseModel<dim>::ShockSpeed() const
 {
    double speed, height;
    ShockSpeedHeight( speed, height );
    return speed;
 }


template<uint32_t dim>
double TwoPhaseModel<dim>::ShockHeight() const
 {
    double speed, height;
    ShockSpeedHeight( speed, height );
    return height;
 }

template<uint32_t dim>
void TwoPhaseModel<dim>::ShockSpeedHeight( double& speed, double& height)const
{
  double se(0.);
  speed = 0;

  for( auto i = 0; i<=100; ++i){
    if( speed < dfds_at( se ) )
    {
      speed  = dfds_at( se );
      height = dfds_at( se ) * se;
    }
    se += 0.01;
  }

}

/// General accessory functions

template<uint32_t dim>
double TwoPhaseModel<dim>::mobility_w_at( double se ) const
{
  return krw_at( se ) / muw_;
}

template<uint32_t dim>
double TwoPhaseModel<dim>::mobility_n_at( double se ) const
{
  return krn_at( se ) / mun_;
}


template<uint32_t dim>
double TwoPhaseModel<dim>::fw_at( double se ) const
{
  return 1.0 / ( 1.0 + mobility_n_at( se ) / mobility_w_at( se ) );
}

template<uint32_t dim>
double TwoPhaseModel<dim>::fn_at( double se ) const
{
  return 1.0 / ( 1.0 + mobility_w_at( se ) / mobility_n_at( se ) );
}

template<uint32_t dim>
double TwoPhaseModel<dim>::krw_at( double se) const
{

  const double cache = seff_;
  seff_ = se;
  const double kr( krw_Phase() );
  seff_ = cache;

  return kr;
}

template<uint32_t dim>
double TwoPhaseModel<dim>::krn_at( double se) const
{

  const double cache = seff_;
  seff_ = se;
  const double kr( krn_Phase() );
  seff_ = cache;

  return kr;
}

template<uint32_t dim>
double TwoPhaseModel<dim>::dkrwds_at( double se) const
{
  const double cache = seff_;
  seff_ = se;
  const double dkrwds( dkrwds_Phase() );
  seff_ = cache;

  return dkrwds;
}

template<uint32_t dim>
double TwoPhaseModel<dim>::dkrnds_at( double se) const
{

  const double cache = seff_;
  seff_ = se;
  const double dkrnds( dkrnds_Phase() );
  seff_ = cache;

  return dkrnds;
}


template<uint32_t dim>
double TwoPhaseModel<dim>::dfds_at( double se ) const
{

  const double cache = seff_;
  seff_ = se;
  const double df = dfds();
  seff_ = cache;

  return df;
}

template<uint32_t dim>
double TwoPhaseModel<dim>::G_at( double se ) const
{

  const double cache = seff_;
  seff_ = se;
  const double g = G();
  seff_ = cache;

  return g;
}

template<uint32_t dim>
double TwoPhaseModel<dim>::pc_at( double se ) const
{

  const double cache = seff_;
  seff_ = se;
  const double pc = pc_Phase( );
  seff_ = cache;

  return pc;
}

template<uint32_t dim>
double TwoPhaseModel<dim>::dpcds_at( double se) const
{
    double cache = seff_;
    seff_ = se;
    double dpcds( dpcds_Phase( ) );
    seff_ = cache;

    return dpcds;
}

template<uint32_t dim>
double TwoPhaseModel<dim>::Sw_at( double pc ) const
{
  const double cache = seff_;
  const double Sw = Sw_Phase( pc );
  seff_ = cache;

  return Sw;
}

template<uint32_t dim>
double TwoPhaseModel<dim>::dsdpc_at( double pc ) const
{
    const double cache = seff_;
    const double dsdpc = dsdpc_Phase( pc );
    seff_ = cache;

    return dsdpc;
}



template<uint32_t dim>
void TwoPhaseModel<dim>::Out( uint32_t phase ) const
 { 
    assert( phase == 1U or phase == 2U );
    cout <<"\nTwoPhaseModel<" << dim <<">::Out: ";
    if ( phase == 1U ) cout <<" data for the wetting phase: ";
    else               cout <<" data for the non-wetting phase: ";
    if ( !interpolate_fluid_properties_ )
      cout <<"fluid properties are fixed. "<< endl;
    else 
      cout <<"fluid properties are interpolated. "<< endl;
      
    //perm_key,    // permeability = scalar element property
    //sat_key,     // saturation of the wetting phase (nodal property)
    //snr_key,     // irreducible saturation of non-wetting phase (element property)
    //swr_key,     // irreducible saturation of wetting phase (element property)
    //mtrl_key,    // Brooks-Corey or vanGenuchten Parameter
    //ift_key,     // interfacial tension key
    //mun_key,     // viscosity of non-wetting phase (nodal property)
    //muw_key,     // viscosity of wetting phase (nodal property)
    //rhw_key,     // density of non-wetting phase (nodal property)
    //rhn_key;     // density of non-wetting phase (nodal property)
    cout <<"\nelement properties:";
    cout <<"\n       irreducible saturation of wetting phase: "<< swr_;
    cout <<"\n   irreducible saturation of non-wetting phase: "<< snr_;
    cout <<"\n                                  permeability: "<< k_;
    cout <<"\nnode properties:";
    cout <<"\n                    viscosity of wetting phase: "<< muw_;
    cout <<"\n                viscosity of non-wetting phase: "<< mun_;
    cout <<"\n                      density of wetting phase: "<< rhw_;
    cout <<"\n                  density of non-wetting phase: "<< rhn_;
    cout <<"\ninterfacial tension (surface tension of fluid): "<< ift_;
    cout <<"\n                      saturation wetting phase: "<< sat_;
    cout <<"\n            effective saturation wetting phase: "<< seff_;
    cout <<"\n                               fractional flow: "<< f_Phase(phase);
    cout <<"\n                                         df/dS: "<< dfds();
    cout <<"\n                                   df/dSn(max): "<< MaxFractionalFlowDerivative();
    cout <<"\n                                             G: "<< G();
    cout <<"\n                                         dG/dS: "<< dGds();
    cout <<"\n                                         pc(2): "<< pc_Phase();
    cout <<"\n                                      dpcdS(2): "<< dpcds_Phase();
    cout << endl << endl;
}



template class TwoPhaseModel<1U>;
template class TwoPhaseModel<2U>;
template class TwoPhaseModel<3U>;


} // end namespace csp





