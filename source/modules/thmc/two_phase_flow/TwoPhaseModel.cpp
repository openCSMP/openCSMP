#include "TwoPhaseModel.h"

using namespace std;

namespace csmp {


template<size_t dim>
TwoPhaseModel<dim>::TwoPhaseModel()
 : mun_(2.0e-3), muw_(1.6e-3),
   rhn_(800.), rhw_(1000.),
   sat_(std::numeric_limits<double64>::quiet_NaN()),
   seff_(std::numeric_limits<double64>::quiet_NaN()),
   swr_(0.), snr_(0.),
   k_(std::numeric_limits<double64>::quiet_NaN()),
   K_(PLAIN,std::numeric_limits<double64>::quiet_NaN()),
   ift_(0.05), // N m-1 Danesh (2003), p. 292
   acc_gravity_(9.8066),
   tolerance_(numeric_limits<double64>::epsilon()),
   MAX_CAPILLARY_PRESSURE_(4e7), // maximum strength of a rock
   MAX_CAPILLARY_PRESSURE_SLOPE_(1e7), // maximum strength of a rock
   interpolate_fluid_properties_(false),
   sw_ro_mu_placement_(true)
 {
 }
 
 
template<size_t dim>
TwoPhaseModel<dim>::TwoPhaseModel( const PropertyDatabase<dim>& database,
                                      double64 viscosity_nw, double64 viscosity_w,
                                      double64 density_nw, double64 density_w,
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
   sat_(std::numeric_limits<double64>::quiet_NaN()),
   seff_(std::numeric_limits<double64>::quiet_NaN()),
   swr_(0.), snr_(0.),
   k_(std::numeric_limits<double64>::quiet_NaN()),
   K_(PLAIN,std::numeric_limits<double64>::quiet_NaN()),
   ift_(0.05), // N m-1 Danesh (2003), p. 292
   acc_gravity_(9.8066),
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



template<size_t dim>
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
   sat_(std::numeric_limits<double64>::quiet_NaN()),
   seff_(std::numeric_limits<double64>::quiet_NaN()),
   swr_(0.), snr_(0.),
   k_(std::numeric_limits<double64>::quiet_NaN()),
   ift_(0.05), // N m-1 Danesh (2003), p. 292
   mun_(std::numeric_limits<double64>::quiet_NaN()),
   muw_(std::numeric_limits<double64>::quiet_NaN()),
   rhn_(std::numeric_limits<double64>::quiet_NaN()),
   rhw_(std::numeric_limits<double64>::quiet_NaN()),
   acc_gravity_(9.8066),
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
 
 
 
template<size_t dim>
TwoPhaseModel<dim>::~TwoPhaseModel()
 {
 }



template<size_t dim>
void TwoPhaseModel<dim>::InterpolateNodeProperties( const Element<dim>& e )
 {
    if ( interpolate_fluid_properties_ ) {
         mun_ = muw_ = rhn_ = rhw_ = sat_ = static_cast<double64>(0.);
         for ( size_t i=0U; i<e.Nodes(); i++ ) {
              sat_ += e.FE()->NRST[i] * e.N(i)->Read( sat_key_ );
              mun_ += e.FE()->NRST[i] * e.N(i)->Read( mun_key_ );
              muw_ += e.FE()->NRST[i] * e.N(i)->Read( muw_key_ );
              rhn_ += e.FE()->NRST[i] * e.N(i)->Read( rhn_key_ );
              rhw_ += e.FE()->NRST[i] * e.N(i)->Read( rhw_key_ );
           }
      }
    else {
         sat_ = static_cast<double64>(0.);
         for ( size_t i=0U; i<e.Nodes(); i++ )
           sat_ += e.FE()->NRST[i] * e.N(i)->Read( sat_key_ );
      }
 }


/// relative permeabilities
template<size_t dim>
double64 TwoPhaseModel<dim>::krw_Phase() const
 {
    cout <<"\nTwoPhaseModel<dim>::krw_Phase (base class): ";
    cout <<"This method needs to be defined in this subclass to achieve desired functionality."<< endl;
    throw logic_error("TwoPhaseModel<dim>::krw_Phase: Method not defined in subclass");
    return std::numeric_limits<double64>::quiet_NaN();
 }



template<size_t dim>
double64 TwoPhaseModel<dim>::krn_Phase() const
 {
    cout <<"\nTwoPhaseModel<"<< dim <<">::krn_Phase (base class): ";
    cout <<"This method needs to be defined in this subclass to achieve desired functionality."<< endl;
    throw logic_error("TwoPhaseModel<dim>::krn_Phase: Method not defined in subclass");
    return std::numeric_limits<double64>::quiet_NaN();
 }

/// derivatives of relative permeabilities
template<size_t dim>
double64 TwoPhaseModel<dim>::dkrwds_Phase() const
 {
    return dkrwds_numerical();
 }



template<size_t dim>
double64 TwoPhaseModel<dim>::dkrnds_Phase() const
 {
    return dkrnds_numerical();
 }

/// capillary pressure
template<size_t dim>
double64 TwoPhaseModel<dim>::pc_Phase() const
 {
    cout <<"\nTwoPhaseModel<"<< dim <<">::pc_Phase (base class): ";
    cout <<"This method needs to be defined in this subclass to achieve desired functionality."<< endl;
    throw logic_error("TwoPhaseModel<dim>::pc_Phase: Method not defined in subclass");
    return std::numeric_limits<double64>::quiet_NaN();
 }

/// inverse capillary pressure function
template<size_t dim>
double64 TwoPhaseModel<dim>::Sw_Phase(double64 ) const
 {
    cout <<"\nTwoPhaseModel<"<< dim <<">::SwFromPc(base class): ";
    cout <<"This method needs to be defined in this subclass to achieve desired functionality."<< endl;
    throw logic_error("TwoPhaseModel<dim>::SwFromPc: Method not defined in subclass");
    return std::numeric_limits<double64>::quiet_NaN();
 }

    /// capillary pressure derivatives
template<size_t dim>
double64 TwoPhaseModel<dim>::dpcds_Phase( ) const
 {
    return dpcds_numerical();
 }

/// derivatives of inverse capillary pressure function
template<size_t dim>
double64 TwoPhaseModel<dim>::dsdpc_Phase( double64 ) const
 {
    cout <<"\nTwoPhaseModel<"<<  dim <<">::dpcdsw_Phase (base class): ";
    cout <<"This method needs to be defined in this subclass to achieve desired functionality."<< endl;
    throw logic_error("TwoPhaseModel<dim>::dpcdsw_Phase: Method not defined in subclass");
    return std::numeric_limits<double64>::quiet_NaN();
 }

// linearized diffusion multiplier for large-timestep calculations
template<size_t dim>
double64 TwoPhaseModel<dim>::DiffusionCharacteristic( size_t ) const
{
cout <<"\nTwoPhaseModel<"<<  dim <<">::DiffusionCharacteristic (base class): ";
cout <<"This method needs to be defined in this subclass to achieve desired functionality."<< endl;
return std::numeric_limits<double64>::quiet_NaN();
}


/// derivative of wetting phase mobility
template<size_t dim>
double64 TwoPhaseModel<dim>::dlwds() const
 {
    return dkrwds_Phase() / muw_;

 }

/// derivative of non-wetting phase mobility
template<size_t dim>
double64 TwoPhaseModel<dim>::dlnds() const
 {
    return dkrnds_Phase() / mun_;

 }

/// derivative of fractional flow function (advection multipliers)
template<size_t dim>
double64 TwoPhaseModel<dim>::dfds() const
 {

    //if ( ( seff_ < static_cast<double64>(0.) ) || ( seff_ > static_cast<double64>(1.) ) )
    //    return static_cast<double64>(0.);

    const double64 lw  = krw_Phase() / muw_;
    const double64 ln  = krn_Phase() / mun_;
    const double64 lt  = lw + ln;
    const double64 lt2 = lt*lt;

    const double64 dlwds = dkrwds_Phase() / muw_;
    const double64 dlnds = dkrnds_Phase() / mun_;

    return ( dlwds*ln - dlnds*lw )/lt2;

 }
       

     
/// derivatives of gravitational flow
template<size_t dim>
double64 TwoPhaseModel<dim>::dGds( ) const
 {

    //if ( ( seff_ < static_cast<double64>(0.) ) || ( seff_ > static_cast<double64>(1.) ) )
    //    return static_cast<double64>(0.);

    const double64 lw  = krw_Phase() / muw_;
    const double64 ln  = krn_Phase() / mun_;
    const double64 lt  = lw + ln;
    const double64 lt2 = lt*lt;
    const double64 ln2 = ln*ln;
    const double64 lw2 = lw*lw;


    const double64 dlwds = dkrwds_Phase() / muw_;
    const double64 dlnds = dkrnds_Phase() / mun_;

    return ( dlwds*ln2 + dlnds*lw2 )/lt2;
 }



/// Numerical derivatives

/// derivative of wetting phase mobility
template<size_t dim>
double64 TwoPhaseModel<dim>::dlwds_numerical( double64 h ) const
 {
    return dkrwds_numerical( h ) / muw_;

 }

/// derivative of non-wetting phase mobility
template<size_t dim>
double64 TwoPhaseModel<dim>::dlnds_numerical( double64 h) const
 {
    return dkrnds_numerical( h ) / mun_;

 }

template<size_t dim>
double64 TwoPhaseModel<dim>::dkrwds_numerical( double64 h ) const
{

  //if ( seff_ < 0.0 || seff_ > 1.0 )
  //      return static_cast<double64>(0.);

  const double64 dSedSw( 1.0/ (1.0 - swr_ - snr_ ) );

  if( seff_ < 0.+h )
      return (krw_at( seff_ + h ) - krw_at( seff_ ) ) / h * dSedSw;
  if( seff_ > 1.-h )
      return ( krw_at( seff_ ) - krw_at( seff_ - h ) ) / h * dSedSw;

  return ( krw_at( seff_+h ) - krw_at( seff_-h ) )/ (2.0*h) * dSedSw;

}

template<size_t dim>
double64 TwoPhaseModel<dim>::dkrnds_numerical( double64 h ) const
{
  //if ( seff_ < 0.0 || seff_ > 1.0 )
  //      return static_cast<double64>(0.);

  const double64 dSedSw( 1.0/ (1.0 - swr_ - snr_ ) );

  if( seff_ < 0.+h )
      return ( krn_at( seff_ + h ) - krn_at( seff_ ) ) / h * dSedSw;
  if( seff_ > 1.-h )
      return ( krn_at( seff_ ) - krn_at( seff_ - h ) ) / h * dSedSw;

  return ( krn_at( seff_+h ) - krn_at( seff_-h ) )/ (2.0*h) * dSedSw;

}

template<size_t dim>
double64 TwoPhaseModel<dim>::dfds_numerical( double64 h) const
{

  //if ( seff_ < 0.0 || seff_ > 1.0 )
  //    return static_cast<double64>(0.);

  /*
  // first version: direct differentiation
  const double64 dSedSw( 1.0/ (1.0 - swr_ - snr_ ) );

  if( seff_ < 0.+h )
      return ( fw_at( seff_ + h ) - fw_at( seff_ ) ) / h * dSedSw;
  if( seff_ > 1.-h )
      return ( fw_at( seff_ ) - fw_at( seff_ - h ) ) / h * dSedSw;

  return ( fw_at( seff_+h ) - fw_at( seff_-h ) )/ (2.0*h)* dSedSw;

  */

  ///*
  // second version: mixed analytical and numerical differentiation
  const double64 lw  = krw_Phase() / muw_;
  const double64 ln  = krn_Phase() / mun_;
  const double64 lt  = lw + ln;
  const double64 lt2 = lt*lt;

  const double64 dlwds = dkrwds_numerical( h )/muw_;
  const double64 dlnds = dkrnds_numerical( h )/mun_;

  return ( dlwds*ln - dlnds*lw )/lt2;
  //*/

}

template<size_t dim>
double64 TwoPhaseModel<dim>::dGds_numerical( double64 h ) const
{
  //if ( seff_ < 0.0 || seff_ > 1.0 )
  //      return static_cast<double64>(0.);

  /*
  // first version: direct differentiation
  const double64 dSedSw( 1.0/ (1.0 - swr_ - snr_ ) );
  h = 0.0000001;
  if( seff_ < 0.+h )
      return ( G_at( seff_ + h ) - G_at( seff_ ) ) / h * dSedSw;
  if( seff_ > 1.-h )
      return ( G_at( seff_ ) - G_at( seff_ - h ) ) / h * dSedSw;


  return ( G_at( seff_+h ) - G_at( seff_-h ) )/ (2.0*h) * dSedSw;
  */

  // second version: mixed analytical and numerical differentiation
  const double64 lw  = krw_Phase() / muw_;
  const double64 ln  = krn_Phase() / mun_;
  const double64 lt  = lw + ln;
  const double64 lt2 = lt*lt;
  const double64 ln2 = ln*ln;
  const double64 lw2 = lw*lw;


  const double64 dlwds = dkrwds_numerical( h )/muw_;
  const double64 dlnds = dkrnds_numerical( h )/mun_;

  return ( dlwds*ln2 + dlnds*lw2 )/lt2;

}

template<size_t dim>
double64 TwoPhaseModel<dim>::dpcds_numerical( double64 h) const
{

  //if ( seff_ < 0.0 || seff_ > 1.0 )
  //    return static_cast<double64>(0.);

  const double64 dSedSw( 1.0/ (1.0 - swr_ - snr_ ) );
  double64 dpcds;

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


template<size_t dim>
double64 TwoPhaseModel<dim>::MaxFractionalFlowDerivative() const
 {
    double64 speed, height;
    ShockSpeedHeight( speed, height );
    return speed;
 }

    /// linearized fractional flow derivative
template<size_t dim>
double64 TwoPhaseModel<dim>::ShockSpeed() const
 {
    double64 speed, height;
    ShockSpeedHeight( speed, height );
    return speed;
 }


template<size_t dim>
double64 TwoPhaseModel<dim>::ShockHeight() const
 {
    double64 speed, height;
    ShockSpeedHeight( speed, height );
    return height;
 }

template<size_t dim>
void TwoPhaseModel<dim>::ShockSpeedHeight( double64& speed, double64& height)const
{
  double64 se(0.);
  speed = 0;

  for( size_t i = 0; i<=100; ++i){
    if( speed < dfds_at( se ) )
    {
      speed  = dfds_at( se );
      height = dfds_at( se ) * se;
    }
    se += 0.01;
  }

}

/// General accessory functions

template<size_t dim>
inline double64 TwoPhaseModel<dim>::mobility_w_at( double64 se ) const
{
  return krw_at( se ) / muw_;
}

template<size_t dim>
inline double64 TwoPhaseModel<dim>::mobility_n_at( double64 se ) const
{
  return krn_at( se ) / mun_;
}


template<size_t dim>
inline double64 TwoPhaseModel<dim>::fw_at( double64 se ) const
{
  return 1.0 / ( 1.0 + mobility_n_at( se ) / mobility_w_at( se ) );
}

template<size_t dim>
inline double64 TwoPhaseModel<dim>::fn_at( double64 se ) const
{
  return 1.0 / ( 1.0 + mobility_w_at( se ) / mobility_n_at( se ) );
}

template<size_t dim>
inline double64 TwoPhaseModel<dim>::krw_at( double64 se) const
{

  const double64 cache = seff_;
  seff_ = se;
  const double64 kr( krw_Phase() );
  seff_ = cache;

  return kr;
}

template<size_t dim>
inline double64 TwoPhaseModel<dim>::krn_at( double64 se) const
{

  const double64 cache = seff_;
  seff_ = se;
  const double64 kr( krn_Phase() );
  seff_ = cache;

  return kr;
}

template<size_t dim>
inline double64 TwoPhaseModel<dim>::dkrwds_at( double64 se) const
{
  const double64 cache = seff_;
  seff_ = se;
  const double64 dkrwds( dkrwds_Phase() );
  seff_ = cache;

  return dkrwds;
}

template<size_t dim>
inline double64 TwoPhaseModel<dim>::dkrnds_at( double64 se) const
{

  const double64 cache = seff_;
  seff_ = se;
  const double64 dkrnds( dkrnds_Phase() );
  seff_ = cache;

  return dkrnds;
}


template<size_t dim>
inline double64 TwoPhaseModel<dim>::dfds_at( double64 se ) const
{

  const double64 cache = seff_;
  seff_ = se;
  const double64 df = dfds();
  seff_ = cache;

  return df;
}

template<size_t dim>
inline double64 TwoPhaseModel<dim>::G_at( double64 se ) const
{

  const double64 cache = seff_;
  seff_ = se;
  const double64 g = G();
  seff_ = cache;

  return g;
}

template<size_t dim>
inline double64 TwoPhaseModel<dim>::pc_at( double64 se ) const
{

  const double64 cache = seff_;
  seff_ = se;
  const double64 pc = pc_Phase( );
  seff_ = cache;

  return pc;
}

template<size_t dim>
inline double64 TwoPhaseModel<dim>::dpcds_at( double64 se) const
{
    double64 cache = seff_;
    seff_ = se;
    double64 dpcds( dpcds_Phase( ) );
    seff_ = cache;

    return dpcds;
}

template<size_t dim>
inline double64 TwoPhaseModel<dim>::Sw_at( double64 pc ) const
{
  const double64 cache = seff_;
  const double64 Sw = Sw_Phase( pc );
  seff_ = cache;

  return Sw;
}

template<size_t dim>
inline double64 TwoPhaseModel<dim>::dsdpc_at( double64 pc ) const
{
    const double64 cache = seff_;
    const double64 dsdpc = dsdpc_Phase( pc );
    seff_ = cache;

    return dsdpc;
}


/// Interpolations
template<size_t dim>
double64 TwoPhaseModel<dim>::spline_value( double64 x, double64 x1, double64 x2, double64 y1, double64 y2, double64 k1, double64 k2) const
{
    const double64 a =  k1*( x2-x1 ) - ( y2 - y1 );
    const double64 b = -k2*( x2-x1 ) + ( y2 - y1 );
    const double64 t = ( x - x1) / ( x2 - x1 );

    return (1. - t)*y1 + t*y2 + t*(1.-t)*( a*(1.-t) + b*t);

}

template<size_t dim>
double64 TwoPhaseModel<dim>::spline_derivative( double64 x, double64 x1, double64 x2, double64 y1, double64 y2, double64 k1, double64 k2) const
{
    const double64 a =  k1*( x2-x1 ) - ( y2 - y1 );
    const double64 b = -k2*( x2-x1 ) + ( y2 - y1 );
    const double64 t = ( x - x1) / ( x2 - x1 );

    return (y2-y1)/( x2-x1 ) + (1.-2.*t)*( a*(1.-t)+b*t)/(x2-x1) + t*(1.-t)*(b-a)/(x2-x1);

}

template<size_t dim>
double64 TwoPhaseModel<dim>::spline_second_derivative( double64 x, double64 x1, double64 x2, double64 y1, double64 y2, double64 k1, double64 k2) const
{
    const double64 a =  k1*( x2-x1 ) - ( y2 - y1 );
    const double64 b = -k2*( x2-x1 ) + ( y2 - y1 );
    const double64 t = ( x - x1) / ( x2 - x1 );

    return 2.*( b-2.*a +(a-b)*3.*t)/(x2-x1)/(x2-x1);

}

template<size_t dim>
void TwoPhaseModel<dim>::Out( std::ostream& os, size_t phase ) const
 { 
    assert( phase == 1U or phase == 2U );
    os <<"\nTwoPhaseModel<" << dim <<">::Out: ";
    if ( phase == 1U ) os <<" data for the wetting phase: ";
    else               os <<" data for the non-wetting phase: ";
    if ( !interpolate_fluid_properties_ )
      os <<"fluid properties are fixed. "<< endl;
    else 
      os <<"fluid properties are interpolated. "<< endl;
      
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
    os <<"\nelement properties:";
    os <<"\n       irreducible saturation of wetting phase: "<< swr_;
    os <<"\n   irreducible saturation of non-wetting phase: "<< snr_;
    os <<"\n                                  permeability: "<< k_;
    os <<"\nnode properties:";
    os <<"\n                    viscosity of wetting phase: "<< muw_;
    os <<"\n                viscosity of non-wetting phase: "<< mun_;
    os <<"\n                      density of wetting phase: "<< rhw_;
    os <<"\n                  density of non-wetting phase: "<< rhn_;
    os <<"\ninterfacial tension (surface tension of fluid): "<< ift_;
    os <<"\n                      saturation wetting phase: "<< sat_;
    os <<"\n            effective saturation wetting phase: "<< seff_;
    os <<"\n                               fractional flow: "<< f_Phase(phase);
    os <<"\n                                         df/dS: "<< dfds();
    os <<"\n                                   df/dSn(max): "<< MaxFractionalFlowDerivative();
    os <<"\n                                             G: "<< G();
    os <<"\n                                         dG/dS: "<< dGds();
    os <<"\n                                         pc(2): "<< pc_Phase();
    os <<"\n                                      dpcdS(2): "<< dpcds_Phase();
    os << endl << endl;
}



template class TwoPhaseModel<1U>;
template class TwoPhaseModel<2U>;
template class TwoPhaseModel<3U>;


} // end namespace csp





