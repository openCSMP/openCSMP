#include <limits>

#include "NaClSaturatedVaporLookup.h"
#include "binaryReadWrite.h"
#include "Brine.h"
#include "NaClSaturatedVapor.h"
#include "LookupPropertyIndex.h"


using namespace std;

// 16.04.2008 known issues: interpolations near VLH boundary do not catch some topologies.
//                          however, unlike in the case of the twophase surfaces near the critical curve,
//                          this should be no problem as no sharp gradients appear. Nevertheless, the 
//                          safe versions should be implemented

// Feb 2011: probably resolved the above stuff, further testing required, TD

namespace csmp
{
  NaClSaturatedVaporLookup::NaClSaturatedVaporLookup(const double64& externaltemperature, 
                                                     const double64& externalpressure)
    : temperature(externaltemperature), 
      pressure(externalpressure), 
      tcurrent(0.0),
      pcurrent(0.0),
      xcurrent(0.0),
      tdummy(0.0),
      t_res(0.0), 
      p_res(0.0), 
      tnorm(0.0), 
      pnorm(0.0),
      t_iA(0.0),
      t_iB(0.0), 
      t_iC(0.0), 
      t_iD(0.0), 
      p_iA(0.0), 
      p_iB(0.0), 
      p_iC(0.0), 
      p_iD(0.0), 
    pvlh_iB(0.0), 
    pvlh_iA(0.0),
    value_bottom(0.0), 
    value_top(0.0), 
    value_interpolated(0.0), 
    value_iA(0.0), 
    value_iB(0.0), 
    value_iC(0.0), 
    value_iD(0.0), 
    value_before(0.0), 
    value_behind(0.0), 
    value_vlh(0.0), 
    value_vlh_behind(0.0), 
    value_vlh_before(0.0),
    tvlh(0.0), 
    x_low(0.0), 
    x_high(0.0),
    pmax_low(0.0),
    pmax_high(0.0),
    it(0), 
    ip(0), 
    t_dim(0), 
    p_dim(0), 
    iA(0), 
    iB(0), 
    iC(0), 
    iD(0), 
    i_dummy(0),
    it_p_max(0),
    ip_p_max_low(0),
    ip_p_max_high(0),
    state(none), 
    state_iA(none), 
    state_iB(none), 
    state_iC(none), 
    state_iD(none),
    vlh_vapor(tdummy), 
    water(),
    csmp_error( ErrorHandler::Instance() )
  {
    GetTemperatureIndex(vlh_vapor.Tmax());
    it_p_max      = it;
    tdummy        = 590.; // *** dangerous if resolution of lookup table changes; this is intended to be the temperature at it 
    pmax_low      = vlh_vapor.Pressure();
    GetPressureIndex(pmax_low);
    ip_p_max_low  = ip;
    tdummy        = 600.; // *** dito but at it+1
    pmax_high     = vlh_vapor.Pressure();
    GetPressureIndex(pmax_high);
    ip_p_max_high = ip;
    
    if( ip_p_max_low != ip_p_max_high )
      {
        csmp_error.notice( FATAL_ERROR, 
                           "NaClSaturatedVaporLookup constructor -",
                           "ip_p_max_low != ip_p_max_high.\nThis probably implies that the lookup table resolution was changed.\nIf that heppened, lookup near the pressure maximum of the vlh curve is likely to fail.\nReport issue to Thomas Driesner, thomas.driesner@erdw.ethz.ch"); 
      }

    GetTemperatureIndex(1000.0001);
    t_dim = it+1;
    GetPressureIndex(5000.0e5+0.0001);
    p_dim = ip+1;
    //d cout << "NaClSaturatedVaporLookup t_dim and p_dim = " << t_dim << "\t" << p_dim << endl;
    
    storage_vector.resize(t_dim*p_dim*max_index);
    state_vector.resize(t_dim*p_dim);
    
    char filename[60], statefilename[60];
    strcpy( filename, "NaClSaturatedVaporPropertiesLookupTable.bin" );
    strcpy( statefilename, "NaClSaturatedVaporStateLookupTable.bin" );
	fstream infile1(filename, ios::in | ios::binary);
	fstream infile2(statefilename, ios::in | ios::binary);

	if (!infile1.is_open() || !infile2.is_open())
      {
        //d cout << "NaClSaturatedVaporLookup : at least one lookup file missing, computing ...\n\n";
        it = 0;
      
        state    = V; // for brine only!
        NaClSaturatedVapor                  naclsatvap(tcurrent,pcurrent);
        Brine                               brine(tcurrent,pcurrent,xcurrent);
      
        cout << "NaClSaturatedVaporLookup::NaClSaturatedVaporLookup(...): computing subcritical data ...\n";
        for(tcurrent = 0.0e0; tcurrent < 1000.1e0; tcurrent += t_res)
          {
            GetTemperatureIndex(tcurrent+1.0e-3);
            cout << "computing for t = " << tcurrent << ", data set " << it << endl << endl;
            if(it > t_dim-1)
              {
                cout << "too high t-index ...\n";
                break;
              }
	
            tdummy = tcurrent;
	
            //@todo: don't have to go to 5000, vlh_vapor.Pmax is sufficient!!!
            for(pcurrent = 1.0e5; pcurrent <= 5000.0e5; pcurrent += p_res)
              {
                GetPressureIndex(pcurrent+1.0e-3);
                if(ip > p_dim-1)
                  {
                    cout << "too high p-index ...\n";
                    break;
                  }
                // The following three cases are distinguished to make sure that lookup can properly
                // be done even if one to three corners of the cell are NOT in the VH field
                // Assignment of states in the state_vector is purely formal and for identification
                // of such cases, it does NOT imply that the respective state (except VH) is
                // indeed existing at that T-P coordinate at xsat
                if( (tcurrent <= 800.7) && (pcurrent < vlh_vapor.Pressure()) )
                  {
                    xcurrent = naclsatvap.MassFractionNaCl();
                    storage_vector[t_dim*p_dim*temperature_index     + it*p_dim + ip] = tcurrent;
                    storage_vector[t_dim*p_dim*pressure_index        + it*p_dim + ip] = pcurrent;
                    storage_vector[t_dim*p_dim*composition_index     + it*p_dim + ip] = xcurrent;
                    storage_vector[t_dim*p_dim*density_index         + it*p_dim + ip] = brine.Density();
                    storage_vector[t_dim*p_dim*enthalpy_index        + it*p_dim + ip] = brine.Enthalpy();
                    storage_vector[t_dim*p_dim*heatcapacity_index    + it*p_dim + ip] = brine.HeatCapacity();
                    storage_vector[t_dim*p_dim*compressibility_index + it*p_dim + ip] = brine.Compressibility();
                    storage_vector[t_dim*p_dim*viscosity_index       + it*p_dim + ip] = brine.Viscosity();
                    state_vector[it*p_dim + ip ] = VH;
                  }
                else if( (tcurrent < 374.0e0) && (pcurrent > water.SaturationPressureFromT(tcurrent)) )
                  { // t> 374 not necessary becuase of smaller discetization
                    storage_vector[t_dim*p_dim*temperature_index     + it*p_dim + ip] = tcurrent;
                    storage_vector[t_dim*p_dim*pressure_index        + it*p_dim + ip] = pcurrent;
                    storage_vector[t_dim*p_dim*composition_index     + it*p_dim + ip] = 0.0e0;
                    storage_vector[t_dim*p_dim*density_index         + it*p_dim + ip] = 0.0e0;
                    storage_vector[t_dim*p_dim*enthalpy_index        + it*p_dim + ip] = 0.0e0;
                    storage_vector[t_dim*p_dim*heatcapacity_index    + it*p_dim + ip] = 0.0e0;
                    storage_vector[t_dim*p_dim*compressibility_index + it*p_dim + ip] = 0.0e0;
                    storage_vector[t_dim*p_dim*viscosity_index       + it*p_dim + ip] = 0.0e0;
                    state_vector[it*p_dim + ip ] = F;
                  }
                else
                  { // not comprehensively correct
                    storage_vector[t_dim*p_dim*temperature_index     + it*p_dim + ip] = tcurrent;
                    storage_vector[t_dim*p_dim*pressure_index        + it*p_dim + ip] = pcurrent;
                    storage_vector[t_dim*p_dim*composition_index     + it*p_dim + ip] = 0.0e0;
                    storage_vector[t_dim*p_dim*density_index         + it*p_dim + ip] = 0.0e0;
                    storage_vector[t_dim*p_dim*enthalpy_index        + it*p_dim + ip] = 0.0e0;
                    storage_vector[t_dim*p_dim*heatcapacity_index    + it*p_dim + ip] = 0.0e0;
                    storage_vector[t_dim*p_dim*compressibility_index + it*p_dim + ip] = 0.0e0;
                    storage_vector[t_dim*p_dim*viscosity_index       + it*p_dim + ip] = 0.0e0;
                    state_vector[it*p_dim + ip ] = VL;
                  }
              }
          }
        
		fstream outfile1(filename, ios::out | ios::binary);
		if (!outfile1.is_open())
          {
            cout << "could not even it for writing, please stop program and debug !!!!\n";
            char yesno;
            cout << "or enter any key to continue (simulation likely to crash or give wrong results!) :";
            cin  >> yesno;
          }
        else
          {
            cout << "writing file " << filename << " ... ";
            skm_C_fwrite( outfile1, storage_vector );
            outfile1.close();
            cout << "done!\n";
          }

        fstream outfile2(statefilename, ios::out | ios::binary);
		if (!outfile2.is_open())
          {
            cout << "could not even it for writing, please stop program and debug !!!!\n";
            char yesno;
            cout << "or enter any key to continue (simulation likely to crash or give wrong results!) :";
            cin  >> yesno;
          }
        else
          {
            cout << "writing file " << statefilename << " ... ";
            skm_C_fwrite( outfile2, state_vector );
			outfile2.close();
            cout << "done!\n";
          }

      }
	
    else
      {
        cout << "reading file " << filename << " ... ";
        skm_C_fread( infile1, storage_vector );
		infile1.close();		
        cout << "done!\n";

        cout << "reading file " << statefilename << " ... ";
        skm_C_fread( infile2, state_vector );
		infile2.close();
        cout << "done!\n";
	    
      }
    cout << "NaClSaturatedVaporLookup, leaving constructor ...\n\n";
  }


  NaClSaturatedVaporLookup::~NaClSaturatedVaporLookup()
  {
  }



  // The data interpolation routines
  double64 NaClSaturatedVaporLookup::Temperature(){     SetTemperatureAndPressure(); return ValueOf(temperature_index);    }
  double64 NaClSaturatedVaporLookup::Pressure(){        SetTemperatureAndPressure(); return ValueOf(pressure_index);       }
  double64 NaClSaturatedVaporLookup::MassFractionNaCl(){     SetTemperatureAndPressure(); return ValueOf(composition_index);    }
  double64 NaClSaturatedVaporLookup::Density(){         SetTemperatureAndPressure(); return ValueOf(density_index);        }
  double64 NaClSaturatedVaporLookup::Enthalpy(){        SetTemperatureAndPressure(); return ValueOf(enthalpy_index);       }
  double64 NaClSaturatedVaporLookup::HeatCapacity(){    SetTemperatureAndPressure(); return ValueOf(heatcapacity_index);   }
  double64 NaClSaturatedVaporLookup::Compressibility(){ SetTemperatureAndPressure(); return ValueOf(compressibility_index);}
  double64 NaClSaturatedVaporLookup::Viscosity(){       SetTemperatureAndPressure(); return ValueOf(viscosity_index);      }

  double64 NaClSaturatedVaporLookup::ReportMassFractionNaCl(){ return ValueOf(composition_index); }
  double64 NaClSaturatedVaporLookup::ReportEnthalpy(){    return ValueOf(enthalpy_index); }


  void NaClSaturatedVaporLookup::SetTemperatureAndPressure()
  {
    tcurrent = temperature;
    pcurrent = pressure;
    tdummy   = tcurrent;
    return;
  }


  double64 NaClSaturatedVaporLookup::DEnthalpyDT()
  {
    SetTemperatureAndPressure(); 
    tdummy   = tcurrent;

    if(tcurrent <= vlh_vapor.Tmax())
      {
        tcurrent += 0.1;
        x_high   = ReportEnthalpy();
        tcurrent = tdummy;
        return (x_high-ReportEnthalpy())/0.1;
      }
    else
      {
        tvlh = vlh_vapor.TfromP(pcurrent,700.0);
        if(tcurrent >= tvlh-0.101)
          {
            tcurrent -= 0.1;
            x_low     = ReportEnthalpy();
            tcurrent  = tdummy;
            return (ReportEnthalpy()-x_low)/0.1;
          }
        else
          {
            tcurrent += 0.1;
            x_high    = ReportEnthalpy();
            tcurrent  = tdummy;
            return (x_high-ReportEnthalpy())/0.1;
          }
      }
  } 


  double64 NaClSaturatedVaporLookup::DCompositionDT()
  {
    SetTemperatureAndPressure(); 
    tdummy   = tcurrent;

    if(tcurrent <= vlh_vapor.Tmax())
      {
        tcurrent += 0.1;
        x_high   = MassFractionNaCl();
        tcurrent = tdummy;
        return (x_high-MassFractionNaCl())/0.1;
      }
    else
      {
        tvlh = vlh_vapor.TfromP(pcurrent,700.0);
        if(tcurrent >= tvlh-0.101)
          {
            tcurrent -= 0.1;
            x_low   = MassFractionNaCl();
            tcurrent = tdummy;
            return (MassFractionNaCl()-x_low)/0.1;
          }
        else
          {
            tcurrent += 0.1;
            x_high    = MassFractionNaCl();
            tcurrent = tdummy;
            return (x_high-MassFractionNaCl())/0.1;
          }
      }
  } 


  double64 NaClSaturatedVaporLookup::DMassFractionNaClDT()
  {
    SetTemperatureAndPressure(); 
    tdummy   = tcurrent;

    if(tcurrent <= vlh_vapor.Tmax())
      {
        tcurrent += 0.1;
        x_high   = XNaCl2Massfraction(MassFractionNaCl());
        tcurrent = tdummy;
        return (x_high-XNaCl2Massfraction(MassFractionNaCl()))/0.1;
      }
    else
      {
        tvlh = vlh_vapor.TfromP(pcurrent,700.0);
        if(tcurrent >= tvlh-0.101)
          {
            tcurrent -= 0.1;
            x_low   = XNaCl2Massfraction(MassFractionNaCl());
            tcurrent = tdummy;
            return (XNaCl2Massfraction(MassFractionNaCl())-x_low)/0.1;
          }
        else
          {
            tcurrent += 0.1;
            x_high    = XNaCl2Massfraction(MassFractionNaCl());
            tcurrent  = tdummy;
            return (x_high-XNaCl2Massfraction(MassFractionNaCl()))/0.1;
          }
      }
  }

  void NaClSaturatedVaporLookup::GetIndex_iA(const int& property_index)
  {
    // CAUTION: tcurrent and pcurrent must be known before this function is called !!!
    GetTemperatureIndex(tcurrent);
    GetPressureIndex(pcurrent);
    iA  = t_dim*p_dim*property_index + it*p_dim + ip; //it*p_dim*(property_index+1) + ip;
    iB  = iA + p_dim;
    iC  = iB + 1;
    iD  = iA + 1;

    return;
  }


  double64 NaClSaturatedVaporLookup::ValueOf(const int& property_index)
  {
    // do NOT set tcurrent and pcurrent here, do it outside this function too keep it versatile !!!
   
    GetIndex_iA(property_index);
    //    cout << "it = " << it << ", ip = " << ip << endl;

    //CheckForOutOfRange();
    state_iA = state_vector[ it*p_dim + ip ];
    state_iB = state_vector[ it*p_dim + ip + p_dim ];
    state_iC = state_vector[ it*p_dim + ip + p_dim +1 ];
    state_iD = state_vector[ it*p_dim + ip + 1 ];

    i_dummy  = t_dim*p_dim*temperature_index + it*p_dim + ip;
    //    cout << "i_dummy t = " << i_dummy << endl;
    t_iA     = storage_vector[i_dummy];
    t_iB     = storage_vector[i_dummy + p_dim];
    t_iC     = storage_vector[i_dummy +1 + p_dim];
    t_iD     = storage_vector[i_dummy + 1];
 
    i_dummy  = t_dim*p_dim*pressure_index + it*p_dim + ip;
    //    cout << "i_dummy p = " << i_dummy << endl;
    p_iA     = storage_vector[i_dummy];
    p_iB     = storage_vector[i_dummy + p_dim];
    p_iD     = storage_vector[i_dummy + 1];
    p_iC     = p_iD;
 
    //d     cout << "Temperatures and Pressures iA-iD:\n";
    //d     cout << t_iA << "\t" << t_iB << "\t" << t_iC << "\t" << t_iD << endl;
    //d     cout << p_iA << "\t" << p_iB << "\t" << p_iC << "\t" << p_iD << endl;

    // The following if-statements try to detect the position of the vlh curve
    // relative to the interpolation cell, naming of points is as follows:
    //
    //   P
    //   ^            vlhcurve
    //   |         | /
    //  -D---------C/-
    //   |         |
    //   |        /|
    //   |       / |
    //   |      /x |
    //  -A-----/---B--> T
    //        / 
    //

    // The simplest case is when the vlh curve is outside the cell:
    if( (state_iA == VH) && (state_iB == VH) && (state_iC == VH) && (state_iD == VH) )
      {
        return NormalInterpolation(property_index); 
      }
    else
      {
        // *** WARNING: it and ip must be correct for tcurrent and pcurrent here!
        //              otherwise if-statements will very likely fail
        //	cout << state_iA << "\t" << state_iB << "\t" << state_iC << "\t" << state_iD << endl;
        if( it == it_p_max && (ip == ip_p_max_low || ip == ip_p_max_low+1) ) 
          // we are in the region that contains the pressure maximum of the vlh curve
          // see also comment at begin of constructor body; 
          //ip_p_max_low must be ip_p_max_high for this to be strict!
          return NearVLHMaxInterpolation( property_index );
        else if( state_iB == VH && state_iD == VL && it < it_p_max ) 
          return NearVLHInterpolationLowT( property_index);
        else if( state_iA == VH && state_iC == VL && it > it_p_max ) 
          return NearVLHInterpolationHighT(property_index);
        else if( state_iD == F  && state_iB == VH && it < it_p_max ) 
          return NearVLHInterpolationLowT( property_index);
        else
          {
            cerr << "\nNaClSaturatedVaporLookup::missed to find interpolation type in case 2, returning bogus value ...\n";
            cerr << "tcurrent was " << tcurrent << endl;
            cerr << "pcurrent was " << pcurrent << endl;
            cerr << "tmax, pmax   " << vlh_vapor.Tmax() << ", " << vlh_vapor.Pmax() << endl;
            cerr << "state_iA = " << state_iA << endl;
            cerr << "state_iB = " << state_iB << endl;
            cerr << "state_iC = " << state_iC << endl;
            cerr << "state_iD = " << state_iD << endl;
            tdummy = t_iA;
            cerr << "t_iA     = " << t_iA << endl;
            cerr << "pvlh iA  = " << vlh_vapor.Pressure() << endl;
            tdummy = t_iB;
            cerr << "t_iB     = " << t_iB << endl;
            cerr << "pvlh iB  = " << vlh_vapor.Pressure() << endl;
            return 9.9e99;
          }
      }
  }
  
  double64 NaClSaturatedVaporLookup::NormalInterpolation( const int& property_index )
  {
    //d cout << "Using NaClSaturatedVaporLookup::NormalInterpolation( const int& property_index ) ...\n";

    //   P
    //   ^
    //   |         |
    //  -D---------C-
    //   |  x      |
    //   |         |
    //  -A---------B--> T
    //

    // ****** now compute iA  ******

    tnorm          = (tcurrent - t_iA ) / ( t_iB - t_iA );
    pnorm          = (pcurrent - p_iA ) / ( p_iD - p_iA );
    value_bottom       =  storage_vector[iA] + tnorm*(storage_vector[iB] - storage_vector[iA]);
    value_top          =  storage_vector[iD] + tnorm*(storage_vector[iC] - storage_vector[iD]);
    value_interpolated =  value_bottom + pnorm*(value_top-value_bottom);
    return value_interpolated;
  }

  double64 NaClSaturatedVaporLookup::NearVLHMaxInterpolation( const int& property_index )
  {
    tdummy = tcurrent; // to make sure that vlh_vapor returns what we need;
    // unfortunately, the topology is awkward, hope to catch everything:
    //
    //   P
    //   ^   |         |
    //   |  -------------
    //   |   |         |  + is Tmax, Pmax // this function applies only here
    //   |   |   _+_   |
    //   |  ------------- ip_p_max_low/high+1
    //   |   |_/     \_|
    //   | _/|  x      |\_
    //   |  ------------- ip_p_max_low/high
    //   |   |         |
    //    -------------> T
    //
    // *** CAUTION : usage of this function strictly requires that pcurrent < p_vlh at tcurrent
    if( ip == ip_p_max_low ) // lower cell of the two shown above
      {
        tnorm              = (tcurrent - t_iA) / ( t_iB - t_iA );
        value_bottom       =  storage_vector[iA] + tnorm * ( storage_vector[iB] - storage_vector[iA] );
        pnorm              = (pcurrent - p_iA) / ( vlh_vapor.Pressure() - p_iA );
        value_interpolated =  value_bottom + pnorm * ( vlh_vapor.ValueOf( property_index ) - value_bottom );
        return value_interpolated;
      }
    else if( ip == ip_p_max_low+1 )
      {
        // 1. re-caculate iA, iB because we need those from the lower of the two cells
        iA  = t_dim*p_dim*property_index + it*p_dim + ip_p_max_low; 
        iB  = iA + p_dim;
        iC  = iB + 1;
        iD  = iA + 1;

        // 2. use above scheme
        tnorm              = (tcurrent - t_iA) / ( t_iB - t_iA );
        value_bottom       = storage_vector[iA] + tnorm * ( storage_vector[iB] - storage_vector[iA] );
        pnorm              = (pcurrent - p_iA) / ( vlh_vapor.Pressure() - p_iA );
        value_interpolated = value_bottom + pnorm * ( vlh_vapor.ValueOf( property_index ) - value_bottom );
        return value_interpolated;
      }
    else
      {
        csmp_error.notice( FATAL_ERROR, 
                           "NaClSaturatedVaporLookup::NearVLHMaxInterpolation( const int& property_index ) -",
                           "ip_p_max_low != ip_p_max_high.missed both if-statements.\nReport issue to Thomas Driesner, thomas.driesner@erdw.ethz.ch"); 


        return std::numeric_limits<double>::signaling_NaN();
      }
  }

  double64 NaClSaturatedVaporLookup::NearVLHInterpolationLowT( const int& property_index )
  {
    //d    cout << "Using NaClSaturatedVaporLookup::NearVLHInterpolationLowT( const int& property_index ) ...\n";

    tdummy = tcurrent;
    //    cout << "VLH Pressure at T of interest is " << vlh_vapor.Pressure() << endl;
    if(pcurrent > vlh_vapor.Pressure())
      {
        // this next statement should be fine as this object is queried only if an external object (FluidLookup)
        // determined that we are in the VH field. The problem is - as often - the numerical precision for > and <
        // comparisons.
        if( (vlh_vapor.Pressure()-pcurrent) < 0.0e0) pcurrent -= 1.0e-3;
        else
          {
            cout << "NaClSaturatedVaporLookup::NearVLHInterpolationLowT detected that T-P are rather in the VL field ...!\n";
            cout.setf(ios::scientific);
            cout << "with a pressure difference of pcurrent - vlh_vapor.Pressure() = " << pcurrent - vlh_vapor.Pressure() << endl;
            cout << "returning bogus value ...\n";
            return 9.9e99;
          }
      }
    
    // a few data needed to determine topology
    tvlh     = vlh_vapor.TfromP(pcurrent,200.0);
    tdummy    = t_iB;
    pvlh_iB  = vlh_vapor.Pressure();
    tdummy    = tcurrent;
    
    // There is a number of topologies possible (x is the value to be interpolated):
    
    if(tvlh >= t_iA)
      {
	
        if(pvlh_iB <= p_iC)
          { 
            // case (1) and (2): Interpolation to be done between vlh_vapor_at_pcurrent and B-vlh_vapor_at_t_iB-segment
            //
            //  (1)                               (2)
            //
            //   P                                 P
            //   ^            vlh_v                ^
            //   |         | /                     |         |   vlh_v
            //  -D---------C/-                    -D---------C__/
            //   |         |                       |     ___/| 
            //   |        /|                     --|--__/--x-|------ pcurrent
            // --|-------/x|------ pcurrent      __|/        |
            //   |      /  |                    /  |         |
            //  -A-----/---B--> T                 -A---------B--> T
            //   |    /    |                       |         |
            //
            pnorm             = (pcurrent-p_iB) / (pvlh_iB-p_iB);
            value_iB              = storage_vector[iB];
	    
            tdummy       = t_iB;
            value_vlh_behind = vlh_vapor.ValueOf(property_index);
            tdummy       = tvlh;
            value_vlh        = vlh_vapor.ValueOf(property_index);
	    
            value_behind          =  value_iB + pnorm*( value_vlh_behind - value_iB );
            tnorm             = (tcurrent-tvlh) / (t_iB-tvlh);
            value_interpolated    =  value_vlh + tnorm*(value_behind - value_vlh);
          }
        else
          {
	    
            // case (3) and (4): Interpolation to be done between vlh_vapor_at_pcurrent and B-C
            //
            //  (3)                               (4)
            //
            //   P                                 P
            //   ^         vlh_v               ^    vlh_v
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
            pnorm          = (pcurrent-p_iB) / (p_iC-p_iB);
            value_behind       =  storage_vector[iB] + pnorm*( storage_vector[iC] - storage_vector[iB] );
            tnorm          = (tcurrent-tvlh) / (t_iB-tvlh);
            tdummy         = tvlh;
            value_vlh          = vlh_vapor.ValueOf(property_index);
            value_interpolated =  value_vlh + tnorm*(value_behind - value_vlh);    
          }
      }
    else
      { // tvlh_pcurrent is < T-iA
	
        if(pvlh_iB <= p_iC)
          {
            // case (5): Interpolation to be done between A-vlh_vapor_at_iA-segment and B-vlh_vapor_at_t_iB-segment
            //
            //  (5)                      
            //
            //   P                        
            //   ^                        
            //   |                        
            //  -D---------C   vlh_v  
            //   |       __|__/           
            //   | _____/  |              
            // __|/        |              
            //   |------x--|-pcurrent     
            //  -A---------B--> T         
            //   |         |              
            //
            //
            pnorm        = (pcurrent-p_iB) / (pvlh_iB-p_iB);
            value_iB         = storage_vector[iB];
            tdummy       = t_iB;            
            value_vlh_behind = vlh_vapor.ValueOf(property_index);
            value_behind     =  value_iB + pnorm*( value_vlh_behind - value_iB );
            // needs update
            // ******************* why? *********************
	    
            tdummy       = t_iA;
            pvlh_iA      = vlh_vapor.Pressure();
            pnorm        = (pcurrent-p_iA) / (pvlh_iA-p_iA);
            value_iA         = storage_vector[iA];
            value_vlh        = vlh_vapor.ValueOf(property_index);
	    
            value_before           =  value_iA + pnorm*( value_vlh - value_iA );
            tnorm             = (tcurrent-t_iA) / (t_iB-t_iA);
            value_interpolated    =  value_before + tnorm*(value_behind - value_before);
          }
        else
          {
            // case (6): Interpolation to be done between A-vlh_vapor_at_iA-segment and B-C 
            //
            //  (6)
            //
            //   P
            //   ^    vlh_v
            //   |   /        
            //  -D--/------C-
            //   | /       | 
            //   |/        |
            //   /         |
            // -/|------x--|-pcurrent
            //  -A---------B--> T
            //   |         |        
            //
            pnorm          = (pcurrent-p_iB) / (p_iC-p_iB);
            value_behind       =  storage_vector[iB] + pnorm*( storage_vector[iC] - storage_vector[iB] );
            tnorm          = (tcurrent-t_iA) / (t_iB-t_iA);

            tdummy    = t_iA;
            pvlh_iA   = vlh_vapor.Pressure();
            pnorm     = (pcurrent-p_iA) / (pvlh_iA-p_iA);
            value_iA      = storage_vector[iA];

            value_vlh     = vlh_vapor.ValueOf(property_index);

            value_before           =  value_iA + pnorm*( value_vlh - value_iA );
            tnorm             = (tcurrent-t_iA) / (t_iB-t_iA);
            value_interpolated    =  value_before + tnorm*(value_behind - value_before);
          }
      }

    //    cout << "value_behind = " << value_behind << ", value_vlh = " << value_vlh << ", value_interpolated = " << value_interpolated << endl;

    return value_interpolated;

  }




  double64 NaClSaturatedVaporLookup::NearVLHInterpolationHighT( const int& property_index )
  {
    //d    cout << "Using NaClSaturatedVaporLookup::NearVLHInterpolationHighT( const int& property_index ) ...\n";

    tdummy = tcurrent;
    //    cout << "VLH Pressure at T of interest is " << vlh_vapor.Pressure() << endl;
    if(pcurrent > vlh_vapor.Pressure())
      {
        // this next statement should be fine as this object is queried only if an external object (FluidLookup)
        // determined that we are in the VH field. The problem is - as often - the numerical precision for > and <
        // comparisons.
        if( (vlh_vapor.Pressure() - pcurrent)< 0.0e0) pcurrent -= 1.0e-3;
        else
          {
            cout << "NaClSaturatedVaporLookup::NearVLHInterpolationHighT detected that T-P are rather in the VL field ...!\n";
            cout.setf(ios::scientific);
            cout << "with a pressure difference of pcurrent - vlh_vapor.Pressure() = " << pcurrent - vlh_vapor.Pressure() << endl;
            cout << "returning bogus value ...\n";
            return 9.9e99;
          }
      }
    tvlh      = vlh_vapor.TfromP(pcurrent,700.0e0);

    if(tvlh  <= t_iB)
      {

        tdummy    = t_iD;
        if(vlh_vapor.Pressure() >= p_iD)
          {
	    
            //  (1)                     //  (2)                    
            //                          //                         
            //   P                      //   P                     
            //   ^ vlh                  //   ^ vlh                 
            //   | \       |            //   |    \    |           
            //  -D--\------C--          //  -D-----\---C--         
            //   |   \     |            //  -|-x----\--|-pcurrent  
            //  -|-x--\----|-pcurrent   //   |       \ |           
            //   |     \   |            //   |        \|           
            //   |      \  |            //   |         \
            //  -A-------\-B--> T       //  -A---------B--> T      
            //   |         |            //   |         |           
            //                                 

            pnorm          = (pcurrent-p_iA) / (p_iD-p_iA);
            value_before       =  storage_vector[iA] + pnorm*( storage_vector[iD] - storage_vector[iA] );
            tnorm          = (tcurrent-t_iA) / (tvlh-t_iA);
            tdummy         =  tvlh;
            value_vlh          = vlh_vapor.ValueOf(property_index);
            value_interpolated =  value_before + tnorm*(value_vlh - value_before);   
          }

        else
          {
            // vlh_vapor.Pressure(t_iD) < p_iD
            // tdummy is still t_iD, make this sure if you modify the code

            //    (3)                       // (4)
            //                              //
            //     P                        //    P
            //     ^                        //    ^
            //     |          |             //  \ |          |
            //    -D----------C-            //   \D----------C-
            //     |          |             //    \          |
            //   \_|__        |             //    |\         |
            //   --|x-\-_-_---|--pcurrent   //    | \        |
            //     |        \_|__           //   -|x-\-------|--pcurrent
            //    -A----------B- \__vlh     //   -A---\------B-
            //     |          |             //    |    \     |

            pnorm             = (pcurrent-p_iA) / (vlh_vapor.Pressure()-p_iA);
            value_iA              = storage_vector[iA];
            // tdummy    = t_iD !
            value_vlh_before      = vlh_vapor.ValueOf(property_index);
            tdummy            = tvlh;
            value_vlh             = vlh_vapor.ValueOf(property_index);
	    
            value_before          =  value_iA + pnorm*( value_vlh_before - value_iA );
            tnorm             = (tcurrent-t_iA) / (tvlh-t_iA);
            value_interpolated    =  value_before + tnorm*(value_vlh- value_before);
          }
      }
    else
      {
        // tvlh > t_iB
        tdummy = t_iA;
        if(vlh_vapor.Pressure() <= p_iD)
          {

            //  (5)     
            //  
            //   P    
            //   ^          
            //   |          |    
            //  -D----------C-     
            // \_|___       |         
            //   |   \____  |        
            //   |        \_|__vlh         
            // --|--------x-|--pcurrent 
            //  -A----------B-     
            //   |          |       


            pnorm             = (pcurrent-p_iA) / (vlh_vapor.Pressure()-p_iA);
            value_iA              = storage_vector[iA];
            // tdummy    = t_iA !
            value_vlh_before      = vlh_vapor.ValueOf(property_index);
            value_before          =  value_iA + pnorm*( value_vlh_before - value_iA );

            tdummy            = t_iB;
            pnorm             = (pcurrent-p_iB)/(vlh_vapor.Pressure()-p_iB);
            value_iB              = storage_vector[iB];
            value_vlh_behind      = vlh_vapor.ValueOf(property_index);
            value_behind       = value_iB + pnorm*(value_vlh_behind-value_iB);
            tnorm          = (tcurrent-t_iA) / (t_iB-t_iA);
            value_interpolated = value_before + tnorm*(value_behind-value_before);
          }
        else
          {
	    
            //  (6)
            //
            //   P
            //   w^
            //   |        \ |
            //  -D---------\C-
            //   |          \
            //   |          |\
            // --|--------x-|--pcurrent
            //   |          |
            //  -A----------B-
            //   |          |
            //
     
            // tdummy    = t_iA !

            pnorm        = (pcurrent-p_iA) / (p_iD-p_iA);
            value_before     =  storage_vector[iA] + pnorm*( storage_vector[iD] - storage_vector[iA] );

            tdummy       = t_iB;
            pnorm        = (pcurrent-p_iB)/(vlh_vapor.Pressure()-p_iB);
            value_iB         = storage_vector[iB];
            value_vlh_behind = vlh_vapor.ValueOf(property_index);
            value_behind       = value_iB + pnorm*(value_vlh_behind-value_iB);
            tnorm          = (tcurrent-t_iA) / (t_iB-t_iA);
            value_interpolated = value_before + tnorm*(value_behind-value_before);
          }
      }

    return value_interpolated;
  }

  void NaClSaturatedVaporLookup::GetTemperatureIndex(const double64& t)
  {
    // new version
    if(     t <    0.0e0){  t_res =  5.0; it = 0; }
    else if(t <= 250.0e0){  t_res =  5.0; it =     static_cast<long>( (t-  0.0)/t_res ); }
    else if(t <= 350.0e0){  t_res =  2.0; it =  50+static_cast<long>( (t-250.0)/t_res ); }
    else if(t <= 360.0e0){  t_res =  1.0; it = 100+static_cast<long>( (t-350.0)/t_res ); }
    else if(t <= 370.0e0){  t_res =  0.5; it = 110+static_cast<long>( (t-360.0)/t_res ); }
    else if(t <= 372.0e0){  t_res =  0.2; it = 130+static_cast<long>( (t-370.0)/t_res ); }
    else if(t <= 378.0e0){  t_res =  0.1; it = 140+static_cast<long>( (t-372.0)/t_res ); }
    else if(t <= 380.0e0){  t_res =  0.2; it = 200+static_cast<long>( (t-378.0)/t_res ); }
    else if(t <= 390.0e0){  t_res =  0.5; it = 210+static_cast<long>( (t-380.0)/t_res ); }
    else if(t <= 400.0e0){  t_res =  1.0; it = 230+static_cast<long>( (t-390.0)/t_res ); }
    else if(t <= 450.0e0){  t_res =  2.0; it = 240+static_cast<long>( (t-400.0)/t_res ); }
    else if(t <= 550.0e0){  t_res =  5.0; it = 265+static_cast<long>( (t-450.0)/t_res ); }
    else if(t <= 1000.0e0){ t_res = 10.0; it = 285+static_cast<long>( (t-550.0)/t_res ); }
    else{ it = 330; }//                  t_res = 10.0; it = 285+static_cast<long>( (t-550.0)/t_res ); /* throw out of range ? */ }
    // t_dim is therefore (1000-550)/10+285 + 1 = 331
  }
    
  void NaClSaturatedVaporLookup::GetPressureIndex(const double64& p)
  {
    // new version
    if(     p <=   20.0e5){ p_res =   0.5e5; ip =     static_cast<long>( (p-   0.5e5)/p_res ); }    
    else if(p <=  210.0e5){ p_res =   1.0e5; ip =  39+static_cast<long>( (p-  20.0e5)/p_res ); }
    else if(p <=  215.0e5){ p_res =   0.5e5; ip = 229+static_cast<long>( (p- 210.0e5)/p_res ); }
    else if(p <=  225.0e5){ p_res =   0.1e5; ip = 239+static_cast<long>( (p- 215.0e5)/p_res ); }
    else if(p <=  230.0e5){ p_res =   0.5e5; ip = 339+static_cast<long>( (p- 225.0e5)/p_res ); }
    else if(p <=  250.0e5){ p_res =   1.0e5; ip = 349+static_cast<long>( (p- 230.0e5)/p_res ); }
    else if(p <=  300.0e5){ p_res =   2.0e5; ip = 369+static_cast<long>( (p- 250.0e5)/p_res ); }
    else if(p <=  400.0e5){ p_res =   5.0e5; ip = 394+static_cast<long>( (p- 300.0e5)/p_res ); }
    else if(p <=  700.0e5){ p_res =  10.0e5; ip = 414+static_cast<long>( (p- 400.0e5)/p_res ); }
    else if(p <= 1000.0e5){ p_res =  25.0e5; ip = 444+static_cast<long>( (p- 700.0e5)/p_res ); }
    else if(p <= 2000.0e5){ p_res =  50.0e5; ip = 456+static_cast<long>( (p-1000.0e5)/p_res ); }
    else if(p <= 5000.0e5){ p_res = 100.0e5; ip = 476+static_cast<long>( (p-2000.0e5)/p_res ); }
    else{ ip = 506; }//                  p_res = 100.0; ip = 476+static_cast<long>( (p-2000.0)/p_res );/* throw out of range ? */ }
    // p_dim is therefore (5000-2000)/100+476 + 1 = 507
  }

} // namespace csmp
