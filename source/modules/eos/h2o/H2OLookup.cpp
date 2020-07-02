#include "H2OLookup.h"
#include "ConvertConcentrationUnitsNaCl.h"
#include "binaryReadWrite.h"

#include <string>

using namespace std;

namespace csmp {


H2OLookup::H2OLookup()
{
  // 1. Determine t_dim and p_dim
  //    These determinations upon construction are done to allow changes in the temperature
  //    and pressure resolution scheme of the lookup tables without affecting the rest of the code
  //    in this file.
  //    Known issues:
  //    - if your really change these resolutions, code still tries to read in existing lookup tables 
  //      without range check! So, delete old lookup-tables, they will be re-computed
  it = GetTemperatureIndex( 2000.0001 );
  t_dim = it + 1;
  ip = GetPressureIndex( 10000.0001*1.0e5 );
  p_dim = ip + 1;
  cout << "H2OLookup t_dim, p_dim = " << t_dim << "\t" << p_dim << endl;

  // 2. Parameters for viscosity computations. Ideally to be made static consts
  ak[0] = 0.0181583, ak[1] = 0.0177624, ak[2] = 0.0105287, ak[3] = -0.0036744;

  bij[0][0] = 0.501938, bij[1][0] = 0.162888, bij[2][0] = -0.130356, bij[3][0] = 0.907919,
    bij[4][0] = -0.551119, bij[5][0] = 0.146543, bij[0][1] = 0.235622, bij[1][1] = 0.789393,
    bij[2][1] = 0.673665, bij[3][1] = 1.207552, bij[4][1] = 0.0670665, bij[5][1] = -0.084337,
    bij[0][2] = -0.274637, bij[1][2] = -0.743539, bij[2][2] = -0.959456, bij[3][2] = -0.687343,
    bij[4][2] = -0.497089, bij[5][2] = 0.195286, bij[0][3] = 0.145831, bij[1][3] = 0.263129,
    bij[2][3] = 0.347247, bij[3][3] = 0.213486, bij[4][3] = 0.100754, bij[5][3] = -0.032932,
    bij[0][4] = -0.0270448, bij[1][4] = -0.0253093, bij[2][4] = -0.0267758, bij[3][4] = -0.0822904,
    bij[4][4] = 0.0602253, bij[5][4] = -0.0202595;

  // 3. Set up table parameters
  n_tables = 3;
  strcpy( filename[0], "H2OPropertiesLookupTableSaturatedVapor.bin" );
  strcpy( filename[1], "H2OPropertiesLookupTableSaturatedLiquid.bin" );
  strcpy( filename[2], "H2OPropertiesLookupTableSinglePhase.bin" );

  t_dim_table.resize( n_tables );
  t_dim_table[0] = GetTemperatureIndex( 374.00001 ) + 1;
  t_dim_table[1] = t_dim_table[0];
  t_dim_table[2] = t_dim;
  cout << "H2OLookup found t_dim 0 and 1 = " << t_dim_table[0] << ", t_dim 2 = " << t_dim_table[2] << endl;
  satvap.resize( t_dim_table[0] * max_index );
  satliq.resize( t_dim_table[1] * max_index );
  singlephase.resize( t_dim_table[2] * p_dim*max_index );

  // 4. Try to locate and read existing lookup tables files declared above.
  //    If a file doesn't exist, compute lookup vectors and store in file.
  char vpname[200], lpname[200];
  strcpy( vpname, filename[0] );
  strcpy( lpname, filename[1] );

  fstream vp( vpname, ios::in | ios::binary );
  if ( !vp.is_open() )
  {
    cout << "\nH2OLookup::H2OLookup() - File: " << vpname;
    cout << " could not be opened, computing ..." << endl;
    BuildTable0And1();
    // ... 
    cout << "now reading newly built file " << vpname << " ... ";
    vp.open( vpname, ios::in | ios::binary );
    binaryFileRead( vp, satvap );
    vp.close();
    cout << "done!\n";
  }
  else
  {
    cout << "reading file " << vpname << " ... ";
    binaryFileRead( vp, satvap );
    vp.close();
    cout << "done!\n";
  }

  fstream lp( lpname, ios::in | ios::binary );
  if ( !lp.is_open() )
  {
    cout << "\nH2OLookup::H2OLookup() - File: " << lpname;
    cout << " could not be opened, computing ..." << endl;
    BuildTable0And1();
    // ...
    cout << "now reading newly built file " << lpname << " ... ";
    lp.open( lpname, ios::in | ios::binary );
    binaryFileRead( lp, satliq );
    lp.close();
    cout << "done!\n";
  }
  else
  {
    cout << "reading file " << lpname << " ... ";
    binaryFileRead( lp, satliq );
    lp.close();
    cout << "done !\n";
  }

  // Single phase stuff
  char spname[200];
  strcpy( spname, filename[2] );
  fstream sp( spname, ios::in | ios::binary );
  if ( !sp.is_open() )
  {
    cout << "\nH2OLookup::H2OLookup() - File: " << spname;
    cout << " could not be opened, computing ..." << endl;
    BuildTable2();
    // ...
    cout << "now reading newly built file " << spname << " ... ";
    sp.open( spname, ios::in | ios::binary );
    binaryFileRead( sp, singlephase );
    sp.close();
    cout << "done!\n";
  }
  else
  {
    cout << "reading file " << spname << " ... ";
    binaryFileRead( sp, singlephase );
    sp.close();
    cout << "done!\n";
  }

  cout << "H2OLookup, leaving constructor ...\n\n";

  return;
}


H2OLookup::~H2OLookup()
{
}


// functions that mimic the interface of Sebstian's former H2OPropertiesLookupTable

double64 H2OLookup::SaturationTemperatureFromP( double64 p )
{
  //	cout << "entering H2OLookup::SaturationTemperatureFromP\n"; 
  myt = TfromP( p );
  return myt;
  //	cout << "leaving H2OLookup::SaturationTemperatureFromP\n"; 
}
double64 H2OLookup::LiquidDensityFromP( double64 p ) { myt = TfromP( p ); return LiquidProperty( myt, density_index ); }
double64 H2OLookup::LiquidEnthalpyFromP( double64 p ) { myt = TfromP( p ); return LiquidProperty( myt, enthalpy_index ); }
double64 H2OLookup::LiquidHeatCapacityFromP( double64 p ) { myt = TfromP( p ); return LiquidProperty( myt, heatcapacity_index ); }
double64 H2OLookup::LiquidViscosityFromP( double64 p ) { myt = TfromP( p ); return LiquidProperty( myt, viscosity_index ); }
double64 H2OLookup::LiquidCompressibilityFromP( double64 p ) { myt = TfromP( p ); return LiquidProperty( myt, compressibility_index ); }
double64 H2OLookup::VaporDensityFromP( double64 p ) { myt = TfromP( p ); return VaporProperty( myt, density_index ); }
double64 H2OLookup::VaporEnthalpyFromP( double64 p ) { myt = TfromP( p ); return VaporProperty( myt, enthalpy_index ); }
double64 H2OLookup::VaporHeatCapacityFromP( double64 p ) { myt = TfromP( p ); return VaporProperty( myt, heatcapacity_index ); }
double64 H2OLookup::VaporViscosityFromP( double64 p ) { myt = TfromP( p ); return VaporProperty( myt, viscosity_index ); }
double64 H2OLookup::VaporCompressibilityFromP( double64 p ) { myt = TfromP( p ); return VaporProperty( myt, compressibility_index ); }
double64 H2OLookup::SaturationPressureFromT( double64 t ) { return LiquidProperty( t, pressure_index ); }
double64 H2OLookup::DSaturationPressureFromTDT( double64 t ) {
  if ( t <= 5.0 ) return (SaturationPressureFromT( t + 0.01 ) - SaturationPressureFromT( t )) / 0.01;
  else return (SaturationPressureFromT( t + 0.01 ) - SaturationPressureFromT( t )) / 0.01;
}
double64 H2OLookup::LiquidDensityFromT( double64 t ) { return LiquidProperty( t, density_index ); }
double64 H2OLookup::LiquidEnthalpyFromT( double64 t ) { return LiquidProperty( t, enthalpy_index ); }
double64 H2OLookup::LiquidHeatCapacityFromT( double64 t ) { return LiquidProperty( t, heatcapacity_index ); }
double64 H2OLookup::LiquidViscosityFromT( double64 t ) { return LiquidProperty( t, viscosity_index ); }
double64 H2OLookup::LiquidCompressibilityFromT( double64 t ) { return LiquidProperty( t, compressibility_index ); }
double64 H2OLookup::VaporDensityFromT( double64 t ) { return VaporProperty( t, density_index ); }
double64 H2OLookup::VaporEnthalpyFromT( double64 t ) { return VaporProperty( t, enthalpy_index ); }
double64 H2OLookup::VaporHeatCapacityFromT( double64 t ) { return VaporProperty( t, heatcapacity_index ); }
double64 H2OLookup::VaporViscosityFromT( double64 t ) { return VaporProperty( t, viscosity_index ); }
double64 H2OLookup::VaporCompressibilityFromT( double64 t ) { return VaporProperty( t, compressibility_index ); }
double64 H2OLookup::Density( double64 t, double64 p ) { return SinglePhaseProperty( t, p, density_index ); }
double64 H2OLookup::HeatCapacity( double64 t, double64 p ) { return SinglePhaseProperty( t, p, heatcapacity_index ); }
double64 H2OLookup::Enthalpy( double64 t, double64 p ) { return SinglePhaseProperty( t, p, enthalpy_index ); }
double64 H2OLookup::Viscosity( double64 t, double64 p ) { return SinglePhaseProperty( t, p, viscosity_index ); }
double64 H2OLookup::Compressibility( double64 t, double64 p ) { return SinglePhaseProperty( t, p, compressibility_index ); }
double64 H2OLookup::MolarVolume( double64 t, double64 p ) { return DensityAndX2MolarVolume( SinglePhaseProperty( t, p, density_index ), 0.0 ); }
double64 H2OLookup::LiquidMolarVolumeFromP( double64 p ) { myt = TfromP( p ); return DensityAndX2MolarVolume( LiquidProperty( myt, density_index ), 0.0 ); }
double64 H2OLookup::VaporMolarVolumeFromP( double64 p ) { myt = TfromP( p ); return DensityAndX2MolarVolume( VaporProperty( myt, density_index ), 0.0 ); }
double64 H2OLookup::Dvdt( double64 t, double64 p )
{
  if ( t > cp_h2o.Temperature() ) dvdt = (MolarVolume( t + 0.01, p ) - MolarVolume( t, p )) / 0.01;
  else {
    if ( p > SaturationPressureFromT( t ) )
    {
      if ( t < 1.0e0 ) dvdt = (MolarVolume( t + 0.01, p ) - MolarVolume( t, p )) / 0.01;
      else          dvdt = (MolarVolume( t, p ) - MolarVolume( t - 0.01, p )) / 0.01;
    }
    else dvdt = (MolarVolume( t + 0.01, p ) - MolarVolume( t, p )) / 0.01;
  }
  return dvdt;
}
double64 H2OLookup::Dvdpbar( double64 t, double64 p ) { return -18015.0e0*Compressibility( t, p ) / Density( t, p )*1.0e5; }


// Functions that retrieve data from the lookup table
// No range check done, it's YOUR responsibility !!!
//
// May 2009:
// - simple range check implemented in SinglePhaseProperty(const double64& t, const double64& p, const int& property_index).
//   Any suggestions for better and/or more inexpensive and/or more intelligent range checks are welcome, 
//   send to thomas.driesner@erdw.ethz.ch

// These three were rather for debugging. You don't want to use these, rather got to the three following after these
double64 H2OLookup::LiquidProperty( const long& it, const int& property_index ) { return satliq[t_dim_table[0] * property_index + it]; }

double64 H2OLookup::VaporProperty( const long& it, const int& property_index ) { return satvap[t_dim_table[1] * property_index + it]; }

double64 H2OLookup::SinglePhaseProperty( const long& it, const long& ip, const int& property_index ) {
  return singlephase[t_dim_table[2] * p_dim * property_index + p_dim*it + ip];
}

// The real stuff
double64 H2OLookup::LiquidProperty( const double64& t, const int& property_index )
{
  // Given t and property, return value on liquid side of boiling curve
  // Caution: range ends at tcrit
  // 	cout << "H2OLookup::LiquidProperty(const double64& t, const int& property_index) : computing for property_index = " << property_index << endl;
  // 	cout << "                                                                           t was " << t << endl;
  double64 tnorm;
  it = GetTemperatureIndex( t );
  tnorm = t - satliq[t_dim_table[0] * temperature_index + it];
  tnorm /= satliq[t_dim_table[0] * temperature_index + it + 1] - satliq[t_dim_table[0] * temperature_index + it];
  return satliq[t_dim_table[0] * property_index + it] + tnorm*(satliq[t_dim_table[0] * property_index + it + 1] - satliq[t_dim_table[0] * property_index + it]);
}

double64 H2OLookup::VaporProperty( const double64& t, const int& property_index )
{
  // Given t and property, return value on vapor side of boiling curve
  // Caution: range ends at tcrit
  double64 tnorm;
  it = GetTemperatureIndex( t );
  tnorm = t - satvap[t_dim_table[0] * temperature_index + it];
  tnorm /= satvap[t_dim_table[0] * temperature_index + it + 1] - satvap[t_dim_table[0] * temperature_index + it];
  return satvap[t_dim_table[0] * property_index + it] + tnorm*(satvap[t_dim_table[0] * property_index + it + 1] - satvap[t_dim_table[0] * property_index + it]);
}

double64 H2OLookup::SinglePhaseProperty( const double64& t, const double64& p, const int& property_index )
{
  ErrorHandler& csmp_err( ErrorHandler::Instance() );
  // Given t,p and property, return value for single phase
  int yesno;
  // first, check for extreme temperature
  if ( t >= 2000.0e0 )
  {
    if ( p >= 10000.0e5 )
    {
      // top back corner of lookup table
      csmp_err.notice( WARNING,
                       "H2OLookup::SinglePhaseProperty(const double64& t, const double64& p, const int& property_index) -",
                       "temperature > 2000C and pressure > 1000MPa, i.e., out of valid range. Returning index for 2000C, 1000MPa properties" );
      cout << t << "\t" << p << endl;
      cin >> yesno;
      return singlephase[t_dim_table[2] * p_dim * (property_index + 1) - 1];
    }
    else
    {
      // t extreme but p ok
      // p-interpolation is still missing but unlikely to make any difference
      csmp_err.notice( WARNING,
                       "H2OLookup::SinglePhaseProperty(const double64& t, const double64& p, const int& property_index) -",
                       "temperature > 2000C, i.e., out of valid range. Returning index for 2000C properties" );
      cout << p << endl;
      it = GetTemperatureIndex( t );
      ip = GetPressureIndex( p );
      cin >> yesno;
      return singlephase[t_dim_table[2] * p_dim * property_index + p_dim*it + ip];
    }
  }

  // ok, temperature test passed, now check for remaining extreme pressure
  if ( p == 10000.0e5 )
  {
    // t ok but p extreme
    // t-interpolation is still missing but unlikely to make any real difference
    it = GetTemperatureIndex( t );
    ip = GetPressureIndex( p );
    return singlephase[t_dim_table[2] * p_dim * property_index + p_dim*it + ip];
  }
  if ( p > 10000.0e5 )
  {
    // t ok but p extreme
    // t-interpolation is still missing but unlikely to make any real difference
    csmp_err.notice( WARNING,
                     "H2OLookup::SinglePhaseProperty(const double64& t, const double64& p, const int& property_index) -",
                     "pressure > 1000 MPa, i.e., out of valid range. Returning index for 1000 MPa properties" );
    cout << p << endl;
    it = GetTemperatureIndex( t );
    ip = GetPressureIndex( p );
    // this is based on the assumption that ip = pdim-1 is determined correctly
    cin >> yesno;
    return singlephase[t_dim_table[2] * p_dim * property_index + p_dim*it + ip];
  }

  // both temperature and pressure checks passed, now do normal interpolation
  GetIndex_iA( t, p, property_index );

  if ( t >= 374.0e0 ) return NormalInterpolation( t, p, property_index );
  else if ( (t >= 373.9e0) && (t < 374.0e0) && (p >= 220.5e5) && (p < 220.6e5) )
  {
    return NearCriticalInterpolation( t, p, property_index );
  }
  else
  {
    // perform phase boundary check
    if ( BoilingCurveInInterpolationCell() == true )
    {
      pboil = VaporProperty( t, pressure_index );
      if ( p > pboil )      return LiquidInterpolationNearBoilingCurve( t, p, property_index );
      else if ( p < pboil ) return VaporInterpolationNearBoilingCurve( t, p, property_index );
      else               return AccidentalBoilingCurveEncounter( t, p, property_index );
    }
    else return NormalInterpolation( t, p, property_index );
  }

  return singlephase[t_dim_table[2] * p_dim * property_index + p_dim*it + ip];
}




// some auxilliary functions
double64 H2OLookup::AccidentalBoilingCurveEncounter( const double64& t, const double64& p, const int& property_index )
{
  cout << "******************* accidentally hit boiling curve !!!! ***********************\n";
  cout << "\n have no idea, what to return as only property of one phase was requested!!! \n";
  cout << "\n my advice: best stop run and get back to Thomas Driesner for debugging      \n";
  cout << " otherwise I return bogus value of 9.e99 and csmp will stop anyhow.            \n";
  cout << "\nEnter any key if you want to continue anyhow ...:";
  cout << "t was " << t << endl;
  cout << "p was " << p << endl;

  long myinput;
  cin >> myinput;
  return 9.e99;
}

bool H2OLookup::BoilingCurveInInterpolationCell()
{
  p_boil_low = satvap[t_dim_table[0] * pressure_index + it];
  p_boil_high = satvap[t_dim_table[0] * pressure_index + it + 1];
  //d	cout << "boiling ps = " << p_boil_low << "\t" << p_boil_high << endl;
  if ( (p_boil_low <= p_iD) && (p_boil_high >= p_iB) ) return true;
  else return false;
}

void H2OLookup::GetIndex_iA( const double64& tcurrent, const double64& pcurrent, const int& property_index )
{
  // May 2009: CAUTION - this scheme may fail in case of tcurrent > tmax and pcurrent > pmax
  //           to do: improve range control!

  it = GetTemperatureIndex( tcurrent );
  ip = GetPressureIndex( pcurrent );
  iA = t_dim_table[2] * p_dim*property_index + it*p_dim + ip;
  iB = iA + p_dim;
  iC = iB + 1;
  iD = iA + 1;

  i_dummy = t_dim_table[2] * p_dim*temperature_index + it*p_dim + ip;
  t_iA = singlephase[i_dummy];
  // the range problem is here:
  t_iB = singlephase[i_dummy + p_dim];
  t_iC = singlephase[i_dummy + 1 + p_dim];
  t_iD = singlephase[i_dummy + 1];
  // end range problem for extreme t

  i_dummy = t_dim_table[2] * p_dim*pressure_index + it*p_dim + ip;
  p_iA = singlephase[i_dummy];
  // and here:
  p_iB = singlephase[i_dummy + p_dim];
  p_iD = singlephase[i_dummy + 1];
  p_iC = p_iD;
  // end range problem for extreme p

  //d 	cout << "H2OLookup::GetIndex_iA values are:\n";
  //d 	cout << "t_iA = " << t_iA << endl;
  //d 	cout << "t_iB = " << t_iB << endl;
  //d 	cout << "t_iC = " << t_iC << endl;
  //d 	cout << "t_iD = " << t_iD << endl;
  //d 	cout << "p_iA = " << p_iA << endl;
  //d 	cout << "p_iB = " << p_iB << endl;
  //d 	cout << "p_iC = " << p_iC << endl;
  //d 	cout << "p_iD = " << p_iD << endl;

  return;
}


double64 H2OLookup::DynamicViscosity( const double64& T, const double64& rho )  const
{
  const double64 tstar( 647.27e0 ), rhostar( 317.763e0 );
  double64 trat, trat1, rhorat, rhorat1, n0, n;
  int i, j, k;

  trat = T / tstar;
  trat1 = 1.0e0 / trat - 1.0e0;
  rhorat = rho / rhostar;
  rhorat1 = rhorat - 1.0e0;

  n0 = 0;
  for ( k = 0; k<4; ++k )
  {
    n0 += (ak[k] * pow( (1.0e0 / trat), k ));
  }
  n0 = 1.0e0 / n0;
  n0 *= (sqrt( trat ));

  n = 0;
  for ( i = 0; i<6; ++i )
  {
    for ( j = 0; j<5; ++j )
    {
      n += (bij[i][j] * pow( trat1, i )*pow( rhorat1, j ));
    }
  }
  n = n0*exp( rhorat*n );

  return n*1.0e-06; // converts to kg/s/m
}




// the interpolation routines

double64 H2OLookup::TfromP( const double64& press )
{
  double64 pnorm;
  i_max = t_dim_table[0] - 1;
  i_guess = i_max / 2;
  i_min = 0;
  int        myi = 0;
  while ( i_max - i_guess > 1 )
  {
    if ( press == satliq[i_guess + t_dim_table[0] * pressure_index] )
    {
      return satliq[i_guess]; // if p is exactly correct at i_guess
    }
    else if ( press < satliq[i_guess + t_dim_table[0] * pressure_index] ) { i_max = i_guess; i_guess = (i_guess + i_min) / 2; }
    else { i_min = i_guess; i_guess = (i_max + i_guess) / 2; }
    myi++;
    if ( myi > 50 )
    {
      cout << "not converged ..." << endl;
      return 0.0e0;
    }
  }
  if ( press < satliq[i_guess + t_dim_table[0] * pressure_index] ) i_guess -= 1;
  pnorm = (press - satliq[i_guess + t_dim_table[0] * pressure_index]) / (satliq[i_guess + 1 + t_dim_table[0] * pressure_index] - satliq[i_guess + t_dim_table[0] * pressure_index]);
  return satliq[i_guess] + pnorm*(satliq[i_guess + 1] - satliq[i_guess]);
}


double64 H2OLookup::NormalInterpolation( const double64& tcurrent, const double64& pcurrent, const int& property_index )
{
  //d cout << "Using H2OLookup::NormalInterpolation( const int& property_index ) ...\n";
  // iA etc. must be known before this is called!
  //   P
  //   ^
  //   |         |
  //  -D---------C-
  //   |  x      |
  //   |         |
  //  -A---------B--> T
  //
  double64 tnorm, pnorm;
  tnorm = (tcurrent - t_iA) / (t_iB - t_iA);
  pnorm = (pcurrent - p_iA) / (p_iD - p_iA);
  v_bottom = singlephase[iA] + tnorm*(singlephase[iB] - singlephase[iA]);
  v_top = singlephase[iD] + tnorm*(singlephase[iC] - singlephase[iD]);
  v_interpolated = v_bottom + pnorm*(v_top - v_bottom);
  return v_interpolated;
}



double64 H2OLookup::NearCriticalInterpolation( const double64& tcurrent, const double64& pcurrent, const int& property_index )
{
  //d cout << "using H2OLookup::NearCriticalInterpolation ...\n";
  it_A = GetTemperatureIndex( cp_h2o.Temperature() );
  it_B = it_A + 1;
  ip_A = GetPressureIndex( cp_h2o.Pressure() );
  ip_D = ip_A + 1;

  double64 tnorm, pnorm;
  // topology is as follows:
  //
  //    D|_________|C    
  //     |         |
  //     |         |
  //     |     x   |
  //     |    /    |
  //    _|___/_____|_
  //    A|  /      |B
  //
  // p_D - p_A = 0.1 bar
  // t_B - t_A = 0.1 C

  if ( pcurrent >= cp_h2o.Pressure() )
  {
    if ( tcurrent >= cp_h2o.Temperature() )
    {
      pnorm = (pcurrent - cp_h2o.Pressure()) / (220.6e5 - cp_h2o.Pressure());
      tnorm = (cp_h2o.Temperature() - 373.9e0) / 0.1e0;
      v_top = SinglePhaseProperty( it_A, ip_D, property_index )
        + tnorm*(SinglePhaseProperty( it_B, ip_D, property_index ) - SinglePhaseProperty( it_A, ip_D, property_index ));
      v_behind = SinglePhaseProperty( it_B, ip_A, property_index ) + (pcurrent - 220.5e5) / 0.1e0*
        (SinglePhaseProperty( it_B, ip_D, property_index ) - SinglePhaseProperty( it_B, ip_A, property_index ));
      it = GetTemperatureIndex( 374.001 );
      v_interpolated = LiquidProperty( it, property_index ) + pnorm*(v_top - LiquidProperty( it, property_index ));
      tnorm = (tcurrent - cp_h2o.Temperature()) / (374.0e0 - cp_h2o.Temperature());
      v_interpolated += tnorm * (v_behind - v_interpolated);
    }
    else
    {
      pnorm = (pcurrent - cp_h2o.Pressure()) / (220.6e5 - cp_h2o.Pressure());
      tnorm = (cp_h2o.Temperature() - 373.9e0) / 0.1;
      v_top = SinglePhaseProperty( it_A, ip_D, property_index )
        + tnorm*(SinglePhaseProperty( it_B, ip_D, property_index ) - SinglePhaseProperty( it_A, ip_D, property_index ));
      v_before = SinglePhaseProperty( it_A, ip_A, property_index ) + (pcurrent - 220.5e5) / 0.1e0*
        (SinglePhaseProperty( it_A, ip_D, property_index ) - SinglePhaseProperty( it_A, ip_A, property_index ));
      it = GetTemperatureIndex( 374.001 );
      v_behind = LiquidProperty( it, property_index ) + pnorm*(v_top - LiquidProperty( it, property_index ));
      tnorm = (tcurrent - 373.9e0) / (cp_h2o.Temperature() - 373.9e0);
      v_interpolated = v_before + tnorm * (v_behind - v_before);
    }
  }
  else
  {
    if ( tcurrent >= cp_h2o.Temperature() )
    {
      tboil = TfromP( 220.5e5 );
      pnorm = (pcurrent - 220.5e5) / (cp_h2o.Pressure() - 220.5e5);
      tnorm = (cp_h2o.Temperature() - tboil) / (374.0e0 - tboil);
      v_bottom = VaporProperty( tboil, property_index )
        + tnorm*(SinglePhaseProperty( it_B, ip_A, property_index ) - VaporProperty( tboil, property_index ));
      v_behind = SinglePhaseProperty( it_B, ip_A, property_index ) + (pcurrent - 220.5e5) / 0.1e0*
        (SinglePhaseProperty( it_B, ip_D, property_index ) - SinglePhaseProperty( it_B, ip_A, property_index ));
      it = GetTemperatureIndex( 374.001 );
      v_interpolated = v_bottom + pnorm*(LiquidProperty( it, property_index ) - v_bottom);
      tnorm = (tcurrent - cp_h2o.Temperature()) / (374.0e0 - cp_h2o.Temperature());
      v_interpolated += tnorm * (v_behind - v_interpolated);
    }
    else
    {
      if ( pcurrent > LiquidProperty( tcurrent, pressure_index ) )
      {
        tboil = TfromP( pcurrent );
        tnorm = (tcurrent - 373.9e0) / (tboil - 373.9e0);
        v_before = SinglePhaseProperty( it_A, ip_A, property_index ) + (pcurrent - 220.5e5) / 0.1e0*
          (SinglePhaseProperty( it_A, ip_D, property_index ) - SinglePhaseProperty( it_A, ip_A, property_index ));
        v_boil = LiquidProperty( tboil, property_index );
        v_interpolated = v_before + tnorm*(v_boil - v_before);
      }
      else if ( pcurrent < LiquidProperty( tcurrent, pressure_index ) )
      {
        tboil = TfromP( 220.5e5 );
        pnorm = (pcurrent - 220.5e5) / (cp_h2o.Pressure() - 220.5e5);
        tnorm = (tcurrent - tboil) / (cp_h2o.Temperature() - tboil);
        v_bottom = VaporProperty( tboil, property_index )
          + tnorm*(SinglePhaseProperty( it_B, ip_A, property_index ) - VaporProperty( tboil, property_index ));
        it = GetTemperatureIndex( 374.001 );
        v_behind = v_bottom + pnorm*(LiquidProperty( it, property_index ) - v_bottom);
        tboil = TfromP( pcurrent );
        tnorm = (tcurrent - tboil) / (cp_h2o.Temperature() - tboil);
        v_boil = VaporProperty( tboil, property_index );
        v_interpolated = v_boil + tnorm*(v_behind - v_boil);
      }
      else return AccidentalBoilingCurveEncounter( tcurrent, pcurrent, property_index );
    }
  }
  //d cout << "leaving near-critical interpolation\n";
  return v_interpolated;
}


double64 H2OLookup::LiquidInterpolationNearBoilingCurve( const double64& tcurrent, const double64& pcurrent, const int& property_index )
{
  //d	cout << "Using H2OLookup::LiquidInterpolationNearBoilingCurve(t,p,property_index) ...\n";
  double64 tnorm, pnorm;

  tboil = TfromP( pcurrent );

  if ( tboil <= t_iC )
  {

    if ( satliq[t_dim_table[0] * pressure_index + it] <= p_iA )
    {
      // case (1) and (2): Interpolation to be done between A-D and boil_curve_at_pcurrent
      //
      //  (1)                                (2)                          
      //                                   
      //   P                                  P                             
      //   ^            boilcurve              ^         boilcurve          
      //   |         | /                      |        /                    
      //  -D---------C/-                     -D-------/-C-              
      //   |         |                        |      /  |           
      //   |        /|                        |     /   |        
      // --|--x----/-|------ pcurrent       --|-x--/----|------ pcurrent  
      //   |      /  |                        |   /     |         
      //  -A-----/---B--> T                  -A--/------B--> T      
      //   |    /    |                        |         |   
      //

      v_iA = singlephase[iA];
      v_iD = singlephase[iD];
      pnorm = (pcurrent - p_iD) / (p_iA - p_iD);
      v_before = singlephase[iD] + pnorm*(singlephase[iA] - singlephase[iD]);
      v_boil = LiquidProperty( tboil, property_index );
      tnorm = (tcurrent - t_iA) / (tboil - t_iA);
      v_interpolated = v_before + tnorm*(v_boil - v_before);
    }
    else
    {
      // case (3) and (4): Interpolation to be done between boil_at_t_iD-D-segment 
      // and boil_curve_at_pcurrent
      //
      //      (4)                              (2)
      //                                                               
      //   P                                P
      //   ^    boilcurve                    ^
      //   |   /                            |         |   boilcurve
      //  -D--/------C-                    -D---------C__/
      // --|x/-------|------ pcurrent       |     ___/| 
      //   |/        |                    --|-x__/----|------ pcurrent
      //   /         |                    __|/        |
      //  /|         |                   /  |         |
      //  -A---------B--> T                -A---------B--> T
      //   |         |                      |         |  
      //

      // notice: it is still t = t_iD;
      pnorm = (pcurrent - p_iD) / (LiquidProperty( t_iD, pressure_index ) - p_iD);
      v_iD = singlephase[iD];
      v_boil_before = LiquidProperty( t_iD, property_index );
      v_boil = LiquidProperty( tboil, property_index );
      v_before = v_iD + pnorm*(v_boil_before - v_iD);
      tnorm = (tcurrent - t_iD) / (tboil - t_iD);
      v_interpolated = v_before + tnorm*(v_boil - v_before);
    }
  }
  else
  {
    // tboil is > t_iB 
    if ( LiquidProperty( t_iD, pressure_index ) >= p_iA )
    {
      // case (5): Interpolation to be done between D-boil_at_t_iD-segment and C-boil_at_t_iC-segment
      //
      //  (5)                      
      //
      //   P                        
      //   ^                        
      //   |                        
      //  -D---------C ___  boilcurve  
      // --|-----x--_|/-----pcurrent           
      //   | _____/  |              
      // __|/        |              
      //   |         |     
      //  -A---------B--> T         
      //   |         |              
      //
      //
      pnorm = (pcurrent - p_iD) / (LiquidProperty( t_iD, pressure_index ) - p_iD);
      v_iC = singlephase[iC];
      v_iD = singlephase[iD];
      v_boil_before = LiquidProperty( t_iD, property_index ); // (between A and D)
      v_boil_behind = LiquidProperty( t_iB, property_index );
      v_before = v_iD + pnorm*(v_boil_before - v_iD);
      pnorm = (pcurrent - p_iC) / (LiquidProperty( t_iC, pressure_index ) - p_iC);
      v_behind = v_iC + pnorm*(v_boil_behind - v_iC);
      tnorm = (tcurrent - t_iD) / (t_iB - t_iD);
      v_interpolated = v_before + tnorm*(v_behind - v_before);
    }
    else
    {
      //  (6)
      //
      //   P
      //   ^            boilcurve
      //   |           /
      //  -D--------C-/
      //  -|--x-----|/--- pcurrent
      //   |        /
      //   |       /|
      //   |      / |
      //  -A-----/--B--> T
      //   |    /   |        
      //
      pnorm = (pcurrent - p_iD) / (p_iA - p_iD);
      v_iD = singlephase[iD];
      v_iC = singlephase[iC];
      v_iA = singlephase[iA];
      v_boil_behind = LiquidProperty( t_iC, property_index );
      v_before = v_iD + pnorm*(v_iA - v_iD);
      pnorm = (pcurrent - p_iC) / (LiquidProperty( t_iC, pressure_index ) - p_iC);
      v_behind = v_iC + pnorm*(v_boil_behind - v_iC);
      tnorm = (tcurrent - t_iD) / (t_iB - t_iD);
      v_interpolated = v_before + tnorm*(v_behind - v_before);
    }
  }
  return v_interpolated;
}


double64 H2OLookup::VaporInterpolationNearBoilingCurve( const double64& tcurrent, const double64& pcurrent, const int& property_index )
{
  //d cout << "Using H2OLookup::VaporInterpolationNearBoilingCurve( const int& property_index ) ...\n";
  double64 tnorm, pnorm;

  // a few data needed to determine topology
  tboil = TfromP( pcurrent );
  //d cout << "tboil is " << tboil << endl;
  pboil_iB = VaporProperty( t_iB, pressure_index );

  // There is a number of topologies possible (x is the value to be interpolated):

  if ( tboil >= t_iA )
  {
    if ( pboil_iB <= p_iC )
    {
      //d 		cout << "H2OLookup::VaporInterpolationNearBoilingCurve for t,p = " << tcurrent << ", " << pcurrent << ": case 1/2 for property_index = " << property_index << endl;
      //d 		cout << "t_iA, t_iB, t_iC, t_iD = " << t_iA << ", " << t_iB << ", " << t_iC << ", " << t_iD << endl;
      //d 		cout << "p_iA, p_iB, p_iC, p_iD = " << p_iA << ", " << p_iB << ", " << p_iC << ", " << p_iD << endl;
      // case (1) and (2): Interpolation to be done between boil_v_at_pcurrent and B-boil_v_at_t_iB-segment
      //
      //  (1)                               (2)
      //
      //   P                                 P
      //   ^            boil_v                ^
      //   |         | /                     |         |   boil_v
      //  -D---------C/-                    -D---------C__/
      //   |         |                       |     ___/| 
      //   |        /|                     --|--__/--x-|------ pcurrent
      // --|-------/x|------ pcurrent      __|/        |
      //   |      /  |                    /  |         |
      //  -A-----/---B--> T                 -A---------B--> T
      //   |    /    |                       |         |
      //
      pnorm = (pcurrent - p_iB) / (pboil_iB - p_iB);
      v_iB = singlephase[iB];
      v_boil_behind = VaporProperty( t_iB, property_index );
      v_boil = VaporProperty( tboil, property_index );
      v_behind = v_iB + pnorm*(v_boil_behind - v_iB);
      tnorm = (tcurrent - tboil) / (t_iB - tboil);
      v_interpolated = v_boil + tnorm*(v_behind - v_boil);
    }
    else
    {
      //d 		cout << "H2OLookup::VaporInterpolationNearBoilingCurve for t,p = " << tcurrent << ", " << pcurrent << ": case 3/4 for property_index = " << property_index << endl;
      //d 		cout << "t_iA, t_iB, t_iC, t_iD = " << t_iA << ", " << t_iB << ", " << t_iC << ", " << t_iD << endl;
      //d 		cout << "p_iA, p_iB, p_iC, p_iD = " << p_iA << ", " << p_iB << ", " << p_iC << ", " << p_iD << endl;
      // case (3) and (4): Interpolation to be done between boil_v_at_pcurrent and B-C
      //
      //  (3)                               (4)
      //
      //   P                                 P
      //   ^         boil_v               ^    boil_v
      //   |        /                        |   /        
      //  -D-------/-C-                     -D--/------C-
      //   |      /  |                       | /       | 
      //   |     /   |                     --|/------x-|------ pcurrent
      // --|----/---x|------ pcurrent        /         |
      //   |   /     |                      /|         |
      //  -A--/------B--> T                 -A---------B--> T
      //   |         |                       |         |
      //
      //
      //
      //
      pnorm = (pcurrent - p_iB) / (p_iC - p_iB);
      v_behind = singlephase[iB] + pnorm*(singlephase[iC] - singlephase[iB]);
      tnorm = (tcurrent - tboil) / (t_iB - tboil);
      v_boil = VaporProperty( tboil, property_index );
      v_interpolated = v_boil + tnorm*(v_behind - v_boil);
    }
  }
  else
  { // tboil_pcurrent is < T-iA
    if ( pboil_iB <= p_iC )
    {
      //d 	    cout << "H2OLookup::VaporInterpolationNearBoilingCurve for t,p = " << tcurrent << ", " << pcurrent << ": case 5 for property_index = " << property_index << endl;
      //d 	    cout << "t_iA, t_iB, t_iC, t_iD = " << t_iA << ", " << t_iB << ", " << t_iC << ", " << t_iD << endl;
      //d 	    cout << "p_iA, p_iB, p_iC, p_iD = " << p_iA << ", " << p_iB << ", " << p_iC << ", " << p_iD << endl;
      // case (5): Interpolation to be done between A-boil_v_at_iA-segment and B-boil_v_at_t_iB-segment
      //
      //  (5)                      
      //
      //   P                        
      //   ^                        
      //   |                        
      //  -D---------C   boil_v  
      //   |       __|__/           
      //   | _____/  |              
      // __|/        |              
      //   |------x--|-pcurrent     
      //  -A---------B--> T         
      //   |         |              
      //
      //
      pnorm = (pcurrent - p_iB) / (pboil_iB - p_iB);
      v_iB = singlephase[iB];
      v_boil_behind = VaporProperty( t_iB, property_index );
      v_behind = v_iB + pnorm*(v_boil_behind - v_iB);
      pboil_iA = VaporProperty( t_iA, pressure_index );
      pnorm = (pcurrent - p_iA) / (pboil_iA - p_iA);
      v_iA = singlephase[iA];
      v_boil = VaporProperty( t_iA, property_index );
      v_before = v_iA + pnorm*(v_boil - v_iA);
      tnorm = (tcurrent - t_iA) / (t_iB - t_iA);
      v_interpolated = v_before + tnorm*(v_behind - v_before);
    }
    else
    {
      //d 		cout << "H2OLookup::VaporInterpolationNearBoilingCurve for t,p = " << tcurrent << ", " << pcurrent << ": case 6 for property_index = " << property_index << endl;
      //d 		cout << "t_iA, t_iB, t_iC, t_iD = " << t_iA << ", " << t_iB << ", " << t_iC << ", " << t_iD << endl;
      //d 		cout << "p_iA, p_iB, p_iC, p_iD = " << p_iA << ", " << p_iB << ", " << p_iC << ", " << p_iD << endl;
      // case (6): Interpolation to be done between A-boil_v_at_iA-segment and B-C 
      //
      //  (6)
      //
      //   P
      //   ^    boil_v
      //   |   /        
      //  -D--/------C-
      //   | /       | 
      //   |/        |
      //   /         |
      // -/|------x--|-pcurrent
      //  -A---------B--> T
      //   |         |        
      //
      pnorm = (pcurrent - p_iB) / (p_iC - p_iB);
      v_behind = singlephase[iB] + pnorm*(singlephase[iC] - singlephase[iB]);
      tnorm = (tcurrent - t_iA) / (t_iB - t_iA);
      pboil_iA = VaporProperty( t_iA, pressure_index );
      pnorm = (pcurrent - p_iA) / (pboil_iA - p_iA);
      v_iA = singlephase[iA];
      v_boil = VaporProperty( t_iA, property_index );
      v_before = v_iA + pnorm*(v_boil - v_iA);
      tnorm = (tcurrent - t_iA) / (t_iB - t_iA);
      v_interpolated = v_before + tnorm*(v_behind - v_before);
    }
  }
  return v_interpolated;
}


// functions that build the lookup tables

void H2OLookup::BuildTable0And1()
{
  Prop *vapprops, *liqprops, *properties;
  vapprops = newProp( 'T', 'p', 1 );
  liqprops = newProp( 'T', 'p', 1 );
  properties = newProp( 'T', 'p', 1 );

  std::vector<double64> vap_vector, liq_vector;
  vap_vector.resize( t_dim_table[0] * max_index );
  liq_vector.resize( t_dim_table[0] * max_index );

  for ( tcurrent = 0.0e0; tcurrent < cp_h2o.Temperature(); tcurrent += t_res )
  {
    if ( tcurrent > cp_h2o.Temperature() ) break;
    t_res = GetTemperatureResolution( tcurrent + 1.0e-3 );
    it = GetTemperatureIndex( tcurrent + 1.0e-3 );
    T = tcurrent + 273.15;
    sat_t( T, liqprops, vapprops );
    if ( tcurrent == 0.0e0 ) sat_t( T + 0.011, liqprops, vapprops ); // insurance against non-convergence for derived properties below triple point (0.01 C)

    pcurrent = 0.5*(liqprops->p + vapprops->p); // *** PROST CAUTION liqprops->p != vapprops->p !!! ***

    vap_vector[t_dim_table[0] * temperature_index + it] = tcurrent;
    vap_vector[t_dim_table[0] * pressure_index + it] = pcurrent;
    vap_vector[t_dim_table[0] * composition_index + it] = 0.0e0;
    vap_vector[t_dim_table[0] * density_index + it] = vapprops->d;
    vap_vector[t_dim_table[0] * enthalpy_index + it] = vapprops->h;
    vap_vector[t_dim_table[0] * heatcapacity_index + it] = vapprops->cp;
    vap_vector[t_dim_table[0] * compressibility_index + it] = 1.0e0 / vapprops->d / (vapprops->dp->d_CT); //? 
    vap_vector[t_dim_table[0] * viscosity_index + it] = DynamicViscosity( T, vapprops->d );

    // this is an insurance against problems of PROST to report correct derivatives for saturated vapor
    // very close to the critical point (is also related to the pboil-problem seen a few lines above)
    if ( tcurrent > 373.85 )
    {
      water_td( T, vapprops->d - 2.0, properties );
      vap_vector[t_dim_table[0] * heatcapacity_index + it] = properties->cp;
      vap_vector[t_dim_table[0] * compressibility_index + it] = 1.0e0 / properties->d / (properties->dp->d_CT); //? 
    }

    liq_vector[t_dim_table[0] * temperature_index + it] = tcurrent;
    liq_vector[t_dim_table[0] * pressure_index + it] = pcurrent;
    liq_vector[t_dim_table[0] * composition_index + it] = 0.0e0;
    liq_vector[t_dim_table[0] * density_index + it] = liqprops->d;
    liq_vector[t_dim_table[0] * enthalpy_index + it] = liqprops->h;
    liq_vector[t_dim_table[0] * heatcapacity_index + it] = liqprops->cp;
    liq_vector[t_dim_table[0] * compressibility_index + it] = 1.0e0 / liqprops->d / (liqprops->dp->d_CT); //? 
    liq_vector[t_dim_table[0] * viscosity_index + it] = DynamicViscosity( T, liqprops->d );
  }

  // critical properties
  t_res = GetTemperatureResolution( 374.001 );
  it = GetTemperatureIndex( 374.001 );
  vap_vector[t_dim_table[0] * temperature_index + it] = cp_h2o.Temperature();
  vap_vector[t_dim_table[0] * pressure_index + it] = cp_h2o.Pressure();
  vap_vector[t_dim_table[0] * composition_index + it] = 0.0e0;
  vap_vector[t_dim_table[0] * density_index + it] = cp_h2o.Density();
  vap_vector[t_dim_table[0] * enthalpy_index + it] = cp_h2o.Enthalpy()*1.0e3;
  vap_vector[t_dim_table[0] * heatcapacity_index + it] = 1.0e10; // since infinity is not wanted 
  vap_vector[t_dim_table[0] * compressibility_index + it] = 1.0e10;
  vap_vector[t_dim_table[0] * viscosity_index + it] = DynamicViscosity( cp_h2o.Temperature() + 273.15, cp_h2o.Density() );

  liq_vector[t_dim_table[0] * temperature_index + it] = cp_h2o.Temperature();
  liq_vector[t_dim_table[0] * pressure_index + it] = cp_h2o.Pressure();
  liq_vector[t_dim_table[0] * composition_index + it] = 0.0e0;
  liq_vector[t_dim_table[0] * density_index + it] = cp_h2o.Density();
  liq_vector[t_dim_table[0] * enthalpy_index + it] = cp_h2o.Enthalpy()*1.0e3;
  liq_vector[t_dim_table[0] * heatcapacity_index + it] = 1.0e10;
  liq_vector[t_dim_table[0] * compressibility_index + it] = 1.0e10;
  liq_vector[t_dim_table[0] * viscosity_index + it] = DynamicViscosity( cp_h2o.Temperature() + 273.15, cp_h2o.Density() );

  freeProp( vapprops );
  freeProp( liqprops );
  freeProp( properties );

  char vname[200];
  strcpy( vname, filename[0] );
  fstream vp( vname, ios::out | ios::binary );
  if ( !vp.is_open() )
    cout << "\nH2OLookup::BinaryOut: File: " << filename[0] << " could not be opened" << endl;
  cout << "writing file " << vname << endl;
  binaryFileWrite( vp, vap_vector );
  vp.close();

  char lname[200];
  strcpy( lname, filename[1] );
  fstream lp( lname, ios::out | ios::binary );
  if ( !lp.is_open() )
    cout << "\nH2OLookup::BinaryOut: File: " << filename[1] << " could not be opened" << endl;
  cout << "writing file " << lname << endl;
  binaryFileWrite( lp, liq_vector );
  lp.close();

}


void H2OLookup::BuildTable2()
{
  Prop *properties, *liqprops, *vapprops;
  double64 psat;

  properties = newProp( 'T', 'p', 1 );
  liqprops = newProp( 'T', 'p', 1 );
  vapprops = newProp( 'T', 'p', 1 );

  std::vector<double64> singlephase_vector;
  singlephase_vector.resize( t_dim_table[2] * p_dim*max_index );

  // new variables for region check
  // 9.9.2010 JPW
  double64 dl, dv;
  int    reg;
  S_mliq MLiq;
  S_mpro MPro;

  for ( tcurrent = 0.0e0; tcurrent <= 2000.1e0; tcurrent += t_res )
  {
    t_res = GetTemperatureResolution( tcurrent + 1.0e-3 );
    it = GetTemperatureIndex( tcurrent + 1.0e-3 );
    if ( it > t_dim_table[2] - 1 ) break;
    T = tcurrent + 273.15e0;
    if ( tcurrent < cp_h2o.Temperature() )
    {
      sat_t( T, liqprops, vapprops );
      psat = 0.5*(liqprops->p + vapprops->p);
    }

    for ( pcurrent = 0.0e0; pcurrent <= 10000.1e5; pcurrent += p_res )
    {
      p_res = GetPressureResolution( pcurrent + 1.0e-3 );
      ip = GetPressureIndex( pcurrent + 1.0e-3 );
      if ( ip > p_dim - 1 ) break;

      // these 2 lines create the basis for a pseudo-zero-bar baseline in the lookup-table
      pdummy = pcurrent;
      if ( pcurrent < 0.1e5 ) pcurrent = 5.0e2;

      /*  August 21, 2009, Thomas Driesner : originally, I had the following line here:

      water_tp(T, pcurrent, 500.0e0, 1.0e-8, properties);

      However, at near-critical pressures below 22.1 MPa, this gives problems at 373.9
      They were not visible if PROST was run stand-alone, though!

      The problem seems to be related to convergence of PROST's water_tp function.
      When digging deep into the PROST iaps.c source code, one reads that for correct
      determinations of properties, the suggested density (500.0 in above code line)in
      the vapor region should be below critical density (321.89...) and above critical
      density in the liquid region.

      Therefore, I introduced a psat determination and added the following if's:
      */

      // check for saturated region
      // 9.9.2010 JPW

      singlephase_vector[t_dim_table[2] * p_dim * temperature_index + p_dim*it + ip] = tcurrent;
      singlephase_vector[t_dim_table[2] * p_dim * pressure_index + p_dim*it + ip] = pcurrent;
      singlephase_vector[t_dim_table[2] * p_dim * composition_index + p_dim*it + ip] = 0.0e0;

      if ( tcurrent < cp_h2o.Temperature() )
      {
        reg = region_tp( T, pcurrent*0.1, &dl, &dv, &MLiq, &MPro );
        if ( reg == 4 )
        {
          if ( pcurrent < psat )
          {
            singlephase_vector[t_dim_table[2] * p_dim * density_index + p_dim*it + ip] = vapprops->d;
            singlephase_vector[t_dim_table[2] * p_dim * enthalpy_index + p_dim*it + ip] = vapprops->h;
            singlephase_vector[t_dim_table[2] * p_dim * heatcapacity_index + p_dim*it + ip] = vapprops->cp;
            singlephase_vector[t_dim_table[2] * p_dim * compressibility_index + p_dim*it + ip] = 1.0e0 / vapprops->d / (vapprops->dp->d_CT);
            //		       cout << "vapor compressibility: " << 1.0e0/vapprops->d/(vapprops->dp->d_CT) << endl;
            //		       cout << "vapprops->d: " << vapprops->d << endl;
            //		       cout << "vapprops->dp->d_CT: " << vapprops->dp->d_CT << endl;		       
            singlephase_vector[t_dim_table[2] * p_dim * viscosity_index + p_dim*it + ip] = DynamicViscosity( T, vapprops->d );
          }
          else
          {
            singlephase_vector[t_dim_table[2] * p_dim * density_index + p_dim*it + ip] = liqprops->d;
            singlephase_vector[t_dim_table[2] * p_dim * enthalpy_index + p_dim*it + ip] = liqprops->h;
            singlephase_vector[t_dim_table[2] * p_dim * heatcapacity_index + p_dim*it + ip] = liqprops->cp;
            singlephase_vector[t_dim_table[2] * p_dim * compressibility_index + p_dim*it + ip] = 1.0e0 / liqprops->d / (liqprops->dp->d_CT);
            //		       cout << "liquid compressibility: " << 1.0e0/liqprops->d/(liqprops->dp->d_CT) << endl;
            singlephase_vector[t_dim_table[2] * p_dim * viscosity_index + p_dim*it + ip] = DynamicViscosity( T, liqprops->d );
          }
        }
        else
        {
          if ( pcurrent < psat ) water_tp( T, pcurrent, 200.0e0, 1.0e-8, properties );
          else                water_tp( T, pcurrent, 500.0e0, 1.0e-8, properties );
          singlephase_vector[t_dim_table[2] * p_dim * density_index + p_dim*it + ip] = properties->d;
          singlephase_vector[t_dim_table[2] * p_dim * enthalpy_index + p_dim*it + ip] = properties->h;
          singlephase_vector[t_dim_table[2] * p_dim * heatcapacity_index + p_dim*it + ip] = properties->cp;
          singlephase_vector[t_dim_table[2] * p_dim * compressibility_index + p_dim*it + ip] = 1.0e0 / properties->d / (properties->dp->d_CT);
          //		      cout << "compressibility: " << 1.0e0/properties->d/(properties->dp->d_CT) << endl;
          singlephase_vector[t_dim_table[2] * p_dim * viscosity_index + p_dim*it + ip] = DynamicViscosity( T, properties->d );
        }
      }
      else
      {
        water_tp( T, pcurrent, 500.0e0, 1.0e-8, properties );
        singlephase_vector[t_dim_table[2] * p_dim * density_index + p_dim*it + ip] = properties->d;
        singlephase_vector[t_dim_table[2] * p_dim * enthalpy_index + p_dim*it + ip] = properties->h;
        singlephase_vector[t_dim_table[2] * p_dim * heatcapacity_index + p_dim*it + ip] = properties->cp;
        singlephase_vector[t_dim_table[2] * p_dim * compressibility_index + p_dim*it + ip] = 1.0e0 / properties->d / (properties->dp->d_CT);
        //		       cout << "compressibility: " << 1.0e0/properties->d/(properties->dp->d_CT) << endl;
        singlephase_vector[t_dim_table[2] * p_dim * viscosity_index + p_dim*it + ip] = DynamicViscosity( T, properties->d );
      }

      pcurrent = pdummy;
    }
  }

  freeProp( properties );
  freeProp( vapprops );
  freeProp( liqprops );

  char sname[200];
  strcpy( sname, filename[2] );
  fstream sp( sname, ios::out | ios::binary );
  if ( sp.is_open() )
    cout << "\nH2OLookup::BinaryOut: File: " << filename[2] << " could not be opened" << endl;
  cout << "writing file " << sname << endl;
  binaryFileWrite( sp, singlephase_vector );
  sp.close();

}

long H2OLookup::GetTemperatureIndex( const double64& t )
{
  double64 t_res;
  long      it;
  if ( t <    0.0e0 )
  {
    cout << "H2OLookup::GetTemperatureIndex(const double64& t) : t < 0 (t = " << t << "), better terminate ...\n";
    cout << "or to continue, enter any key : ";
    char yesno;
    cin >> yesno;
  }
  else if ( t <= 250.0e0 ) { t_res = 5.0; it = static_cast<long>((t - 0.0) / t_res); }
  else if ( t <= 350.0e0 ) { t_res = 2.0; it = 50 + static_cast<long>((t - 250.0) / t_res); }
  else if ( t <= 360.0e0 ) { t_res = 1.0; it = 100 + static_cast<long>((t - 350.0) / t_res); }
  else if ( t <= 370.0e0 ) { t_res = 0.5; it = 110 + static_cast<long>((t - 360.0) / t_res); }
  else if ( t <= 372.0e0 ) { t_res = 0.2; it = 130 + static_cast<long>((t - 370.0) / t_res); }
  else if ( t <= 378.0e0 ) { t_res = 0.1; it = 140 + static_cast<long>((t - 372.0) / t_res); }
  else if ( t <= 380.0e0 ) { t_res = 0.2; it = 200 + static_cast<long>((t - 378.0) / t_res); }
  else if ( t <= 390.0e0 ) { t_res = 0.5; it = 210 + static_cast<long>((t - 380.0) / t_res); }
  else if ( t <= 400.0e0 ) { t_res = 1.0; it = 230 + static_cast<long>((t - 390.0) / t_res); }
  else if ( t <= 450.0e0 ) { t_res = 2.0; it = 240 + static_cast<long>((t - 400.0) / t_res); }
  else if ( t <= 550.0e0 ) { t_res = 5.0; it = 265 + static_cast<long>((t - 450.0) / t_res); }
  else if ( t <= 1000.0e0 ) { t_res = 10.0; it = 285 + static_cast<long>((t - 550.0) / t_res); }
  else if ( t <= 2000.0e0 ) { t_res = 50.0; it = 330 + static_cast<long>((t - 1000.0) / t_res); }
  else { it = 350; }; //                  t_res = 10.0; it = 285+static_cast<long>( (t-550.0)/t_res ); /* throw out of range ? */ }
                      // t_dim is therefore (2000-1000)/50+330 + 1 = 351
  return it;
}

long H2OLookup::GetPressureIndex( const double64& p )
{
  double64 p_res;
  long      ip;
  /*!
  New version, August 24, 2009, Thomas Driesner

  In order to resolve low-Temperature vapor reasonably well, the resolution
  in that range could either be increased or a fake "0 bar" value be invented.
  I decided to try the latter, i.e., i will use - in Table2 - a value at 0.005 bar,
  which is always below psat. This increases p_dim by 1 compared to earlier versions.
  This should probably also be invented in the H2O-NaCl Tables ...
  */
  if ( p <= 20.0e5 ) { p_res = 0.5e5; ip = static_cast<long>((p) / p_res); }
  else if ( p <= 210.0e5 ) { p_res = 1.0e5; ip = 40 + static_cast<long>((p - 20.0e5) / p_res); }
  else if ( p <= 215.0e5 ) { p_res = 0.5e5; ip = 230 + static_cast<long>((p - 210.0e5) / p_res); }
  else if ( p <= 225.0e5 ) { p_res = 0.1e5; ip = 240 + static_cast<long>((p - 215.0e5) / p_res); }
  else if ( p <= 230.0e5 ) { p_res = 0.5e5; ip = 340 + static_cast<long>((p - 225.0e5) / p_res); }
  else if ( p <= 250.0e5 ) { p_res = 1.0e5; ip = 350 + static_cast<long>((p - 230.0e5) / p_res); }
  else if ( p <= 300.0e5 ) { p_res = 2.0e5; ip = 370 + static_cast<long>((p - 250.0e5) / p_res); }
  else if ( p <= 400.0e5 ) { p_res = 5.0e5; ip = 395 + static_cast<long>((p - 300.0e5) / p_res); }
  else if ( p <= 700.0e5 ) { p_res = 10.0e5; ip = 415 + static_cast<long>((p - 400.0e5) / p_res); }
  else if ( p <= 1000.0e5 ) { p_res = 25.0e5; ip = 445 + static_cast<long>((p - 700.0e5) / p_res); }
  else if ( p <= 2000.0e5 ) { p_res = 50.0e5; ip = 457 + static_cast<long>((p - 1000.0e5) / p_res); }
  else if ( p <= 5000.0e5 ) { p_res = 100.0e5; ip = 477 + static_cast<long>((p - 2000.0e5) / p_res); }
  else if ( p <= 10000.0e5 ) { p_res = 500.0e5; ip = 507 + static_cast<long>((p - 5000.0e5) / p_res); }
  else { ip = 517; }//                  p_res = 100.0e5; ip = 476+static_cast<long>( (p-2000.0e5)/p_res );/* throw out of range ? */ }
                    // p_dim is therefore (10000e5-5000e5)/500e5+506 + 1 = 518
  return ip;
}

} // end namespace csmp
