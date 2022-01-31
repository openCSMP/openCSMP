#include <limits>

#include "HaliteLiquidusLookup.h"

#include "binaryReadWrite.h"

#include "Brine.h"
#include "HaliteLiquidus.h"
#include "LookupPropertyIndex.h"

#include <cmath>

using namespace std;

namespace csmp
{
  HaliteLiquidusLookup::HaliteLiquidusLookup(const double& externaltemperature, 
                                             const double& externalpressure)
    : temperature(externaltemperature), 
      pressure(externalpressure), 
      tcurrent(0.0),
      pcurrent(0.0),
      xcurrent(0.0),
      tdummy(0.0),
      pdummy(0.0),
      t_res(0.0),
      p_res(0.0),
      tnorm(0.0),
      pnorm(0.0),
      pvlh(0.0),
      t_iA(0.0),
      t_iB(0.0),
      t_iC(0.0),
      t_iD(0.0),
      p_iA(0.0),
      p_iB(0.0),
      p_iC(0.0),
    p_iD(0.0),
    value_iA(0.0),
    value_iB(0.0),
    value_iC(0.0),
    value_iD(0.0),
    value_before(0.0),
    value_behind(0.0),
    value_bottom(0.0),
    value_top(0.0),
    value_interpolated(0.0),
    value_vlh(0.0),
    value_vlh_behind(0.0),
    value_vlh_before(0.0),
    tcrit(0.0),
    tvlh(0.0),
    x_high(0.0),
    vlh_pmax(vlh_liquid.Pmax()),
    vlh_tmax(vlh_liquid.Tmax()),
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
    ip_p_max(0),
    state(none),
    state_iA(none),
    state_iB(none),
    state_iC(none),
    state_iD(none),
    vlh_liquid(tdummy), 
    naclmelt_liquid(tdummy,pdummy),
    csmp_error( ErrorHandler::Instance() )
  {
    // Prepare values that are necessary to identify cases where we are very close to Tmax,Pmax of the VLH surface
    GetTemperatureIndex(vlh_tmax);
    it_p_max          = it;
    GetPressureIndex(vlh_pmax);
    ip_p_max          = ip;
      
    // Setup table dimensions and rescale stroage vector sizes
    GetTemperatureIndex(1000.0001);
    t_dim = it+1;
    GetPressureIndex(5000.0e5+0.001);
    p_dim = ip+1;
    //d  cout << "HaliteLiquidusLookup t_dim and p_dim = " << t_dim << "\t" << p_dim << endl;
	
    storage_vector.resize(t_dim*p_dim*max_index);
    state_vector.resize(t_dim*p_dim);

    char filename[60], statefilename[60];
    strcpy( filename, "HaliteLiquidusPropertiesLookupTable.bin" );
    strcpy( statefilename, "HaliteLiquidusStateLookupTable.bin" );
    
	fstream infile1(filename, ios::in | ios::binary);
	fstream infile2(statefilename, ios::in | ios::binary);

    if( !infile1.is_open() || !infile2.is_open() )
      {
        cout << "HaliteLiquidusLookup : at least one lookup file missing, computing ...\n\n";
        it = 0;
        state    = L; // relevant for viscosity computations, don't change to VL
        HaliteLiquidus                       liquidus(tcurrent,pcurrent);
        Brine                                brine(tcurrent,pcurrent,xcurrent);

        cout << "HaliteLiquidusLookup::HaliteLiquidusLookup(...): computing data ...\n";
        for(tcurrent = 0.0e0; tcurrent < 1000.1e0; tcurrent += t_res)
          {
            GetTemperatureIndex(tcurrent+1.0e-3);
            cout << "HaliteLiquidusLookup computing for t = " << tcurrent << ", data set " << it << endl;
            if(it > t_dim-1){
              cout << "too high t-index ...\n";
              break;
            }

            tdummy = tcurrent;

            for(pcurrent = 0.5e5; pcurrent <= 5000.0e5; pcurrent += p_res)
              {
                pdummy = pcurrent;
                GetPressureIndex(pcurrent+1.0e-3);
                if(ip > p_dim-1)
                  {
                    cout << "too high p-index ...\n";
                    break;
                  }

                if( (tcurrent <= 800.7) && (pcurrent < vlh_liquid.Pressure() ) )
                  {
                    xcurrent = liquidus.MassFractionNaCl();
                    storage_vector[t_dim*p_dim*temperature_index     + it*p_dim + ip] = tcurrent;
                    storage_vector[t_dim*p_dim*pressure_index        + it*p_dim + ip] = pcurrent;
                    storage_vector[t_dim*p_dim*composition_index     + it*p_dim + ip] = xcurrent;
                    storage_vector[t_dim*p_dim*density_index         + it*p_dim + ip] = brine.Density();
                    storage_vector[t_dim*p_dim*enthalpy_index        + it*p_dim + ip] = brine.Enthalpy();
                    storage_vector[t_dim*p_dim*heatcapacity_index    + it*p_dim + ip] = brine.HeatCapacity();
                    storage_vector[t_dim*p_dim*compressibility_index + it*p_dim + ip] = brine.Compressibility();
                    storage_vector[t_dim*p_dim*viscosity_index       + it*p_dim + ip] = brine.Viscosity();
                    state_vector[it*p_dim + ip ] = V;
                  }
                else if( tcurrent > naclmelt_liquid.TmeltFromP() )
                  {
                    storage_vector[t_dim*p_dim*temperature_index     + it*p_dim + ip] = tcurrent;
                    storage_vector[t_dim*p_dim*pressure_index        + it*p_dim + ip] = pcurrent;
                    storage_vector[t_dim*p_dim*composition_index     + it*p_dim + ip] = 0.0e0;
                    storage_vector[t_dim*p_dim*density_index         + it*p_dim + ip] = 0.0e0;
                    storage_vector[t_dim*p_dim*enthalpy_index        + it*p_dim + ip] = 0.0e0;
                    storage_vector[t_dim*p_dim*heatcapacity_index    + it*p_dim + ip] = 0.0e0;
                    storage_vector[t_dim*p_dim*compressibility_index + it*p_dim + ip] = 0.0e0;
                    storage_vector[t_dim*p_dim*viscosity_index       + it*p_dim + ip] = 0.0e0;
                    state_vector[it*p_dim + ip ] = L;
                  }
                else
                  {
                    xcurrent = liquidus.MassFractionNaCl();
                    state    = L;
                    storage_vector[t_dim*p_dim*temperature_index     + it*p_dim + ip] = tcurrent;
                    storage_vector[t_dim*p_dim*pressure_index        + it*p_dim + ip] = pcurrent;
                    storage_vector[t_dim*p_dim*composition_index     + it*p_dim + ip] = xcurrent;
                    storage_vector[t_dim*p_dim*density_index         + it*p_dim + ip] = brine.Density();
                    storage_vector[t_dim*p_dim*enthalpy_index        + it*p_dim + ip] = brine.Enthalpy();
                    storage_vector[t_dim*p_dim*heatcapacity_index    + it*p_dim + ip] = brine.HeatCapacity();
                    storage_vector[t_dim*p_dim*compressibility_index + it*p_dim + ip] = brine.Compressibility();
                    storage_vector[t_dim*p_dim*viscosity_index       + it*p_dim + ip] = brine.Viscosity();
                    state_vector[it*p_dim + ip ] = LH;
                  }
              }
          }
	
		fstream outfile1(filename, ios::out | ios::binary);
        if( !outfile1.is_open() )
          {
            cout << "could not even it for writing, please stop program and debug !!!!\n";
            char yesno;
            cout << "or enter any key to continue (simulation likely to crash or give wrong results!) :";
            cin  >> yesno;
          }
        else
          {
            cout << "writing file " << filename << " ... ";
            binaryFileWrite( outfile1, storage_vector );
			      outfile1.close();
            cout << "done!\n";
          }

        fstream outfile2(statefilename, ios::out | ios::binary);
        if( !outfile2.is_open() )
          {
            cout << "could not even it for writing, please stop program and debug !!!!\n";
            char yesno;
            cout << "or enter any key to continue (simulation likely to crash or give wrong results!) :";
            cin  >> yesno;
          }
        else
          {
            cout << "writing file " << statefilename << " ... ";
            binaryFileWrite( outfile2, state_vector );
			      outfile2.close();
            cout << "done!\n";
          }

      }
	
    else
      {
        cout << "reading file " << filename << " ... ";
        binaryFileRead( infile1, storage_vector );
		    infile1.close();
        cout << "done!\n";

        cout << "reading file " << statefilename << " ... ";
        binaryFileRead( infile2, state_vector );
		    infile2.close();
        cout << "done!\n";
	    
      }
    cout << "HaliteLiquidusLookup, leaving constructor ...\n\n";
  }


  HaliteLiquidusLookup::~HaliteLiquidusLookup()
  {
  }

  void HaliteLiquidusLookup::SetTemperatureAndPressure()
  {
    tcurrent = temperature;
    pcurrent = pressure;
    // ********** this is my personal p_min **********
    // if you want to change it, you must probably adapt p_dim and GetPressureIndex 
    // in all Lookup files
    if(pcurrent < 0.5e5) pcurrent = 0.5e5;
    tdummy   = tcurrent;
    pdummy   = pcurrent;
    return;
  }



  void HaliteLiquidusLookup::GetIndex_iA(const int& property_index)
  {
    // CAUTION: tcurrent and pcurrent must be known before this function is called !!!
    GetTemperatureIndex(tcurrent);
    GetPressureIndex(pcurrent);
    iA  = t_dim*p_dim*property_index + it*p_dim + ip;
    iB  = iA + p_dim;
    iC  = iB + 1;
    iD  = iA + 1;

    return;
  }


  // The data interpolation routines
  double HaliteLiquidusLookup::Temperature(){     SetTemperatureAndPressure(); return ValueOf(temperature_index);    }
  double HaliteLiquidusLookup::Pressure(){        SetTemperatureAndPressure(); return ValueOf(pressure_index);       }
  double HaliteLiquidusLookup::MassFractionNaCl(){     SetTemperatureAndPressure(); return ValueOf(composition_index);    }
  double HaliteLiquidusLookup::Density(){         SetTemperatureAndPressure(); return ValueOf(density_index);        }
  double HaliteLiquidusLookup::Enthalpy(){        SetTemperatureAndPressure(); return ValueOf(enthalpy_index);       }
  double HaliteLiquidusLookup::HeatCapacity(){    SetTemperatureAndPressure(); return ValueOf(heatcapacity_index);   }
  double HaliteLiquidusLookup::Compressibility(){ SetTemperatureAndPressure(); return ValueOf(compressibility_index);}
  double HaliteLiquidusLookup::Viscosity(){       SetTemperatureAndPressure(); return ValueOf(viscosity_index);      }

  double HaliteLiquidusLookup::ReportComposition(){ return ValueOf(composition_index); }
  double HaliteLiquidusLookup::ReportEnthalpy(){    return ValueOf(enthalpy_index); }

  // *******************************************************************************
  // for the following three ones, the case of being very close to the melting curve 
  // probably needs special attention
  // *******************************************************************************

  double HaliteLiquidusLookup::DCompositionDT()
  {
    SetTemperatureAndPressure(); 
    tdummy    = tcurrent;
    tcurrent += 0.1;
    x_high    = ReportComposition();
    tcurrent  = tdummy;
    return (x_high-ReportComposition())/0.1;
  }

  double HaliteLiquidusLookup::DEnthalpyDT()
  {
    SetTemperatureAndPressure(); 
    tdummy    = tcurrent;
    tcurrent += 0.1;
    x_high    = ReportEnthalpy();
    tcurrent  = tdummy;
    return (x_high-ReportEnthalpy())/0.1;
  }

  double HaliteLiquidusLookup::DSaltMassFractionDT()
  {
    SetTemperatureAndPressure(); 
    tdummy    = tcurrent;
    tcurrent += 0.1;
    x_high    = XNaCl2Massfraction(ReportComposition());
    tcurrent  = tdummy;
    return (x_high- XNaCl2Massfraction(ReportComposition()))/0.1;
  }



  double HaliteLiquidusLookup::ValueOf(const int& property_index)
  {
    SetTemperatureAndPressure();

    GetIndex_iA(property_index);
    //d cout << "\n\nHaliteLiquidusLookup::ValueOf for property_index " << property_index << endl;
    //d     cout << "it = " << it << ", ip = " << ip << endl;

    if(ip >= p_dim-1)
      {
        if(it >= t_dim-1)
          {
            // if p >= 5000 bar and t >=  1000C, return value at 5000 bar, 1000C
            return storage_vector[t_dim*p_dim*property_index + (t_dim-1)*p_dim + (p_dim-1)];
          } 
        else
          {
            // if p >= 5000 bar and t < 1000C, return value at 5000 bar and t
            return storage_vector[t_dim*p_dim*property_index + it*p_dim + (p_dim-1)];
          }
      }
    else if(it >= t_dim-1)
      {
        // if p < 5000 bar and t >= 1000C, return value at p and 1000C
        return storage_vector[t_dim*p_dim*property_index + (t_dim-1)*p_dim + ip];
      }
    else if(ip == 0 )
      {
        return storage_vector[t_dim*p_dim*property_index + t_dim*ip];
      }
    else
      {
        state_iA = state_vector[ it*p_dim + ip ];
        state_iB = state_vector[ it*p_dim + ip + p_dim ];
        state_iC = state_vector[ it*p_dim + ip + p_dim +1 ];
        state_iD = state_vector[ it*p_dim + ip + 1 ];
	    
        i_dummy  = t_dim*p_dim*temperature_index + it*p_dim + ip;
        //d     cout << "i_dummy t = " << i_dummy << endl;
        t_iA     = storage_vector[i_dummy];
        t_iB     = storage_vector[i_dummy + p_dim];
        t_iC     = storage_vector[i_dummy +1 + p_dim];
        t_iD     = storage_vector[i_dummy + 1];
	    
        i_dummy  = t_dim*p_dim*pressure_index + it*p_dim + ip;
        //d     cout << "i_dummy p = " << i_dummy << endl;
        p_iA     = storage_vector[i_dummy];
        p_iB     = storage_vector[i_dummy + p_dim];
        p_iD     = storage_vector[i_dummy + 1];
        p_iC     = p_iD;
	    
        //d    cout << "Temperatures and Pressures iA-iD:\n";
        //d     cout << t_iA << "\t" << t_iB << "\t" << t_iC << "\t" << t_iD << endl;
        //d     cout << p_iA << "\t" << p_iB << "\t" << p_iC << "\t" << p_iD << endl;
	    
        if( it == it_p_max && (ip == ip_p_max || ip == ip_p_max-1/*+1*/) )
          {
            return NearVLHMaxInterpolation(property_index);
          }

	
        if( (state_iA == LH) && (state_iB == LH) && (state_iC == LH) && (state_iD == LH) )
          {
            //d 	cout << state_iA << "\t" << state_iB << "\t" << state_iC << "\t" << state_iD << endl;
            return NormalInterpolation(property_index); 
          }
        //   P
        //   ^            vlh
        //   |         | /
        //  -D---------C/-
        //   |         |
        //   |        /|
        //   |       / |
        //   |      /x |
        //  -A-----/---B--> T
        //        / 
	    
        else
          {
            //d 	cout << state_iA << "\t" << state_iB << "\t" << state_iC << "\t" << state_iD << endl;
		
            if( (tcurrent <= vlh_liquid.Tmax()) && state_iB == V && state_iD == LH ) return NearVLHInterpolationLowT(property_index);
            else if( (tcurrent > vlh_liquid.Tmax()) && state_iA == V && state_iC == LH ) return NearVLHInterpolationHighT(property_index);
            else if( state_iD == LH && state_iB == L ) return NearNaClMeltInterpolation(property_index);
            else
              {
                char yesno;
                cout << "HaliteLiquidusLookup::missed to find interpolation type in case 2, returning bogus value ...\n";
                cout << state_iA << "\t" << state_iB << "\t" << state_iC << "\t" << state_iD << endl;
                cout << "conditions were : tcurrent = " << tcurrent << ", pcurrent = " << pcurrent << endl; 
                cin >> yesno;
                return 9.9e99;
              }
          }
      }
  }

	
  double HaliteLiquidusLookup::NormalInterpolation( const int& property_index )
  {
    //d     cout << "Using HaliteLiquidusLookup::NormalInterpolation( const int& property_index ) ...\n";

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
    //d     cout << "property index = " << property_index << endl;
    //d     cout << "tnorm    = " << tnorm << ", pnorm = " << pnorm << endl;
    //d     cout << "tcurrent = " << tcurrent << ", t_iA = " << t_iA << endl;
    value_bottom       =  storage_vector[iA] + tnorm*(storage_vector[iB] - storage_vector[iA]);
    value_top          =  storage_vector[iD] + tnorm*(storage_vector[iC] - storage_vector[iD]);
    value_interpolated =  value_bottom + pnorm*(value_top-value_bottom);
    return value_interpolated;
  }

  double HaliteLiquidusLookup::NearVLHMaxInterpolation(const int& property_index)
  {
    // Feb. 2011
    // The following treatment strongl;y relies on a correct detrmination that we are
    // indeed in the T-P range for VL coexistence; if this is not the case (e.g., if
    // T and P were in the VH field, this appraoch is very likely to FAIL!
    //
    // Geometry is as follows (CAN CHANGE IF RESOLUTION OF LOOKUP TABLES IS CHANGED!!!):
    //
    //     |                 |
    //   --D-----------------C-- 
    //     |         VL      |
    //     |                 |
    //     |   tmax,pmax     |
    //     |       *         |
    //     |      / \        |
    //   -D/A---/----\------C/B- ip_p_max // wrong: +1
    //     | _/       V\     |
    //     |/           L\   |
    //   _/|              H\ |  
    //  /  |       VH          \ 
    //     |                 | \
    //   --A-----------------B-- ip_p_max // needs: -1
    //     |                 |
    //   it_p_max
    //     

    // iA etc. have already been determined when this function is called
    if( ip == ip_p_max /*+1*/ ) // in this case, ABDC in the above figure are moved up by one cell
      {
        if( pcurrent >= vlh_pmax ) 
          {
            // this may seem overly complicated but is required to maintain consistency with
            // VLH_LiquidLookup, which provides exact values tmax, pmax
            if( tcurrent <= vlh_tmax )
              {
                // first, value_before
                pnorm              = (pcurrent - p_iA) / (p_iD - p_iA);
                value_before       =  storage_vector[iA] + pnorm*(storage_vector[iD]-storage_vector[iA]);
                // then value_behind = @vlh_tmax
                tnorm              = (vlh_tmax - t_iD) / (t_iC - t_iD);
                value_top          =  storage_vector[iD] + tnorm*(storage_vector[iC] - storage_vector[iD]);
                pnorm              = (pcurrent - vlh_pmax) / (p_iD - vlh_pmax);
                tdummy             =  vlh_tmax; // needed for vlh_liquid in following line
                value_behind       =  vlh_liquid.ValueOf(property_index) + pnorm*(value_top - vlh_liquid.ValueOf(property_index));
                // then interpolate
                tnorm              = (tcurrent - t_iD) / (vlh_tmax-t_iD);
                value_interpolated =  value_before + tnorm*(value_behind-value_before);
                return                value_interpolated;
              }
            else
              {
                // first, value_before
                tnorm              = (vlh_tmax - t_iD) / (t_iC - t_iD);
                value_top          =  storage_vector[iD] + tnorm*(storage_vector[iC] - storage_vector[iD]);
                pnorm              = (pcurrent - vlh_pmax) / (p_iD - vlh_pmax);
                tdummy             =  vlh_tmax; // needed for vlh_liquid in following line
                value_before       =  vlh_liquid.ValueOf(property_index) + pnorm*(value_top - vlh_liquid.ValueOf(property_index));
                // then value_behind = @vlh_tmax
                pnorm              = (pcurrent - p_iB) / (p_iC - p_iB);
                value_behind       =  storage_vector[iB] + pnorm*(storage_vector[iC]-storage_vector[iB]);
                // then interpolate
                tnorm              = (tcurrent - vlh_tmax) / (t_iC - vlh_tmax);
                value_interpolated =  value_before + tnorm*(value_behind-value_before);
                return                value_interpolated;
              }
          }
        else
          {
            if(tcurrent < vlh_tmax )
              {
                // first, value_before
                pnorm              = (pcurrent - p_iA) / (p_iD - p_iA);
                value_before       =  storage_vector[iA] + pnorm*(storage_vector[iD]-storage_vector[iA]);
                // then value_behind @vlh curve at p
                tdummy             =  vlh_liquid.TfromP(pcurrent,300.); // needed for vlh_liquid in following line
                value_behind       =  vlh_liquid.ValueOf(property_index);
                // then interpolate
                tnorm              = (tcurrent - t_iA) / (tdummy - t_iA);
                value_interpolated =  value_before + tnorm*(value_behind-value_before);
                return                value_interpolated;		
              }
            else if( tcurrent > vlh_tmax )
              {
                // first, value_before @vlh curve at p
                tdummy             =  vlh_liquid.TfromP(pcurrent,700.); // needed for vlh_liquid in following line
                value_before       =  vlh_liquid.ValueOf(property_index);
                // then value_behind 
                pnorm              = (pcurrent - p_iB) / (p_iC - p_iB);
                value_behind       =  storage_vector[iB] + pnorm*(storage_vector[iC]-storage_vector[iB]);
                // then interpolate
                tnorm              = (tcurrent - tdummy) / (t_iB - tdummy);
                value_interpolated =  value_before + tnorm*(value_behind-value_before);
                return                value_interpolated;		
              }
            else
              {
                csmp_error.notice( FATAL_ERROR, 
                                 "HaliteLiquidusLookup::NearVLHMaxInterpolation(const int& property_index) - ",
                                 "p < vlh_pmax (in cell ip_p_max+1) but neither tcurrent<vlh_tmax nor tcurrent>vlh_tmax condition worked!\nReport issue to Thomas Driesner, thomas.driesner@erdw.ethz.ch"); 
                return std::numeric_limits<double>::signaling_NaN();
              }
          }
      }
    else if( ip == ip_p_max-1 )
      {
        //d cout << "HaliteLiquidusLookup::NearVLHMaxInterpolation(const int& property_index) is at ip_p_max-1\n";
        //d cout << "and interpolating for property_index " << property_index << endl;
        if(tcurrent < vlh_tmax )
          {
            //d cout << "tcurrent < t_vlh_max\n";
            // first, value_before
            tdummy             =  t_iA;	    
            pnorm              = (pcurrent - vlh_liquid.Pressure()) / (p_iD -  vlh_liquid.Pressure());
            value_before       =  vlh_liquid.ValueOf(property_index) + pnorm*(storage_vector[iD]-vlh_liquid.ValueOf(property_index));
            // then value_behind @vlh curve at p
            tdummy             =  vlh_liquid.TfromP(pcurrent,300.); // tdummy id the one that vlh_liquid uses, see constructor!
            value_behind       =  vlh_liquid.ValueOf(property_index);
            // then interpolate
            tnorm              = (tcurrent - t_iA) / (tdummy - t_iA);
            value_interpolated =  value_before + tnorm*(value_behind-value_before);
            return                value_interpolated;		
          }
        else if( tcurrent > vlh_tmax )
          {
            tdummy             =  vlh_liquid.TfromP(pcurrent,700.); // tdummy id the one that vlh_liquid uses, see constructor!
            tvlh               = tdummy; 
            value_before       =  vlh_liquid.ValueOf(property_index);
            tdummy             =  t_iB;	    
            pnorm              = (pcurrent - vlh_liquid.Pressure()) / (p_iC - vlh_liquid.Pressure());
            value_behind       =  vlh_liquid.ValueOf(property_index) + pnorm*(storage_vector[iC]-vlh_liquid.ValueOf(property_index));
            tnorm              = (tcurrent - tvlh) / (t_iB - tvlh);
            value_interpolated =  value_before + tnorm*(value_behind-value_before);
            return                value_interpolated;		
          }
        else
          {
            csmp_error.notice( FATAL_ERROR, 
                             "TwophaseLiquidLookup::NearVLHMaxInterpolation(const int& property_index) - ",
                             "p < vlh_pmax (in cell ip_p_max) but neither tcurrent<vlh_tmax nor tcurrent>vlh_tmax condition worked!\nReport issue to Thomas Driesner, thomas.driesner@erdw.ethz.ch"); 
            return std::numeric_limits<double>::signaling_NaN();
          }
      }
    else
      {
        csmp_error.notice( FATAL_ERROR, 
                         "TwophaseLiquidLookup::NearVLHMaxInterpolation(const int& property_index) - ",
                         "Missed ALL if-statements!\nReport issue to Thomas Driesner, thomas.driesner@erdw.ethz.ch"); 

        return std::numeric_limits<double>::signaling_NaN();
      }
  }

  double HaliteLiquidusLookup::NearVLHInterpolationLowT( const int& property_index )
  {
    //d 	cout << "Using HaliteLiquidusLookup::NearVLHInterpolationLowT( const int& property_index ) ...\n";
    tvlh      = vlh_liquid.TfromP(pcurrent,300.0e0);
    if(tvlh  <= t_iC)
      {
        tdummy    = t_iD;
        if(vlh_liquid.Pressure() <= p_iA)
          {
            // case (1) and (2): Interpolation to be done between A-D and vlh_curve_at_pcurrent
            //
            //  (1)                                (2)                          
            //                                   
            //   P                                  P                             
            //   ^            vlhcurve              ^         vlhcurve          
            //   |         | /                      |        /                    
            //  -D---------C/-                     -D-------/-C-              
            //   |         |                        |      /  |           
            //   |        /|                        |     /   |        
            // --|--x----/-|------ pcurrent       --|-x--/----|------ pcurrent  
            //   |      /  |                        |   /     |         
            //  -A-----/---B--> T                  -A--/------B--> T      
            //   |    /    |                        |         |   
            //

            value_iA              = storage_vector[iA];
            value_iD              = storage_vector[iD];
	    
            pnorm             = (pcurrent-p_iD) / (p_iA-p_iD);
            value_before          = storage_vector[iD] + pnorm*( storage_vector[iA] - storage_vector[iD] );
            tdummy            = tvlh;
            value_vlh             = vlh_liquid.ValueOf(property_index);
            tnorm             = (tcurrent-t_iA) / (tvlh-t_iA);
            value_interpolated    =  value_before + tnorm*(value_vlh-value_before);
            //d  cout << "computed via if-2 as value_before = " << value_before << ", " << "value_vlh = " 
            //d  << value_vlh << ", " << value_interpolated << endl;
          }
        else
          {
            // case (3) and (4): Interpolation to be done between vlh_at_t_iD-D-segment 
            // and vlh_curve_at_pcurrent
            //
            //      (4)                              (2)
            //                                                               
            //   P                                P
            //   ^    vlhcurve                    ^
            //   |   /                            |         |   vlhcurve
            //  -D--/------C-                    -D---------C__/
            // --|x/-------|------ pcurrent       |     ___/| 
            //   |/        |                    --|-x__/----|------ pcurrent
            //   /         |                    __|/        |
            //  /|         |                   /  |         |
            //  -A---------B--> T                -A---------B--> T
            //   |         |                      |         |  
            //
	    
            // notice: it is still tdummy = t_iD;
            pnorm             = (pcurrent-p_iD) / (vlh_liquid.Pressure()-p_iD);
            value_iD              = storage_vector[iD];
            value_vlh_before      = vlh_liquid.ValueOf(property_index); // (between A and D)
            tdummy            = tvlh;
            value_vlh             = vlh_liquid.ValueOf(property_index);
            value_before          =  value_iD + pnorm*( value_vlh_before - value_iD );
            tnorm             = (tcurrent-t_iD) / (tvlh-t_iD);
            value_interpolated    =  value_before + tnorm*(value_vlh-value_before);
            //d  cout << "value_iD = " << value_iD << ", tnorm = " << tnorm << endl;
            //d  cout << "computed via if-1 as value_before = " << value_before << ", " 
            //d       << "value_vlh = " << value_vlh << ", " << value_interpolated << endl;
          }
      }
    else
      {
        // tvlh is > t_iB 
        tdummy    = t_iD;
        if(vlh_liquid.Pressure() >= p_iA)
          {
            // case (5): Interpolation to be done between D-vlh_at_t_iD-segment and C-vlh_at_t_iC-segment
            //
            //  (5)                      
            //
            //   P                        
            //   ^                        
            //   |                        
            //  -D---------C ___  vlhcurve  
            // --|-----x--_|/-----pcurrent           
            //   | _____/  |              
            // __|/        |              
            //   |         |     
            //  -A---------B--> T         
            //   |         |              
            //
            //
            pnorm             = (pcurrent-p_iD) / (vlh_liquid.Pressure()-p_iD);
            value_iC              = storage_vector[iC];
            value_iD              = storage_vector[iD];
            value_vlh_before      = vlh_liquid.ValueOf(property_index); // (between A and D)
            tdummy            = t_iB;
            value_vlh_behind      = vlh_liquid.ValueOf(property_index);

            value_before          =  value_iD + pnorm*( value_vlh_before - value_iD );
            tdummy            =  t_iC;
            pnorm             = (pcurrent-p_iC) / (vlh_liquid.Pressure()-p_iC);
            value_behind          =  value_iC + pnorm*( value_vlh_behind - value_iC); 
            tnorm             = (tcurrent-t_iD) / (t_iB-t_iD);
            value_interpolated    =  value_before + tnorm*(value_behind-value_before);
          }
        else
          {
	
            //  (6)
            //
            //   P
            //   ^            vlhcurve
            //   |           /
            //  -D--------C-/
            //  -|--x-----|/--- pcurrent
            //   |        /
            //   |       /|
            //   |      / |
            //  -A-----/--B--> T
            //   |    /   |        
            //
            pnorm             = (pcurrent-p_iD) / (p_iA-p_iD);
            value_iD              = storage_vector[iD];
            value_iC              = storage_vector[iC];
            value_iA              = storage_vector[iA];

            tdummy            = t_iC;
            value_vlh_behind      = vlh_liquid.ValueOf(property_index);
	    
            value_before          =  value_iD + pnorm*( value_iA - value_iD );
            tdummy            =  t_iC;
            pnorm             = (pcurrent-p_iC) / (vlh_liquid.Pressure()-p_iC);
            value_behind          =  value_iC + pnorm*( value_vlh_behind - value_iC); 
            tnorm             = (tcurrent-t_iD) / (t_iB-t_iD);
            value_interpolated    =  value_before + tnorm*(value_behind-value_before);
          }
      }
    return value_interpolated;


  }




  double HaliteLiquidusLookup::NearVLHInterpolationHighT( const int& property_index )
  {
    //d  cout << "Using HaliteLiquidusLookup::NearVLHInterpolationHighT( const int& property_index ) ...\n";
    tvlh      = vlh_liquid.TfromP(pcurrent,700.0e0);
    //d  cout << "tvlh found as " << tvlh << endl;
    if(tvlh  >= t_iA)
      {
        tdummy    = t_iB;
        if(vlh_liquid.Pressure() <= p_iB)
          {
            //d  cout << "case 1 or 2 \n";
            //  (1)                        (2)
            //
            //   P                          P
            //   ^ vlh                      ^
            //   | \       |                |          |
            //  -D--\------C--             -D----------C-
            //   |   \     |              \ |          |
            //  -|----\--x-|-pcurrent      \|          |
            //   |     \   |                \          |
            //   |      \  |               -|\--x------|--pcurrent
            //  -A-------\-B--> T          -A-\--------B-
            //   |         |                |  \       |
            //                                  vlh
            pnorm          = (pcurrent-p_iB) / (p_iC-p_iB);
            value_behind       =  storage_vector[iB] + pnorm*( storage_vector[iC] - storage_vector[iB] );
            tnorm          = (tcurrent-tvlh) / (t_iB-tvlh);
            tdummy         =  tvlh;
            value_vlh          = vlh_liquid.ValueOf(property_index);
	    
            value_interpolated =  value_vlh + tnorm*(value_behind - value_vlh);   
          }
        else
          {
            //d  cout << "case 3 or 4 \n";
            // vlh_liquid.Pressure() > p_iB

            //  (3)                        (4)
            //
            //   P                          P
            //   ^ vlh                      ^
            //   |    \    |                |          |
            //  -D-----\---C--             -D----------C-
            //  -|------\x-|-pcurrent       |          |
            //   |       \ |              \_|__        |
            //   |        \|              --|--\-_-_--x|--pcurrent
            //   |         \                |        \_|__
            //  -A---------B--> T          -A----------B- \__vlh
            //   |         |                |          |
            pnorm             = (pcurrent-p_iC) / (vlh_liquid.Pressure()-p_iC);
            value_iC              = storage_vector[iC];
            // tdummy    = t_iB !
            value_vlh_behind      = vlh_liquid.ValueOf(property_index);
            tdummy            = tvlh;
            value_vlh             = vlh_liquid.ValueOf(property_index);
	    
            value_behind          =  value_iC + pnorm*( value_vlh_behind - value_iC );
            tnorm             = (tcurrent-tvlh) / (t_iB-tvlh);
            value_interpolated    =  value_vlh + tnorm*(value_behind - value_vlh);
          }
      }
    else
      {
        //d  cout << "case 5 \n";
        // tvlh < t_iA
        tdummy = t_iB;
        if(vlh_liquid.Pressure() >= p_iB)
          {

            //  (4)     
            //  
            //   P    
            //   ^          
            //   |          |    
            //  -D----------C-     
            // --|--------x-|--pcurrent 
            // \_|___       |         
            //   |   \____  |        
            //   |        \_|__         
            //  -A----------B- \__vlh    
            //   |          |       


            pnorm             = (pcurrent-p_iC) / (vlh_liquid.Pressure()-p_iC);
            value_iC              = storage_vector[iC];
            // tdummy    = t_iB !
            value_vlh_behind      = vlh_liquid.ValueOf(property_index);
            tdummy            = tvlh;
            value_vlh             = vlh_liquid.ValueOf(property_index);
            value_behind          =  value_iC + pnorm*( value_vlh_behind - value_iC );

            tdummy       = t_iD;
            pnorm        = (pcurrent-p_iD)/(vlh_liquid.Pressure()-p_iD);
            value_iD         = storage_vector[iD];
            value_vlh_before  = vlh_liquid.ValueOf(property_index);

            value_before       = value_iD + pnorm*(value_vlh_before-value_iD);
            tnorm          = (tcurrent-t_iD) / (t_iB-t_iD);
            value_interpolated = value_before + tnorm*(value_behind-value_before);
          }
        else
          {
            //d  cout << "case 6 \n";
	    
            //  (6)
            //
            //   P
            //   w^
            //   |          |
            //  -D----------C-
            // --|--------x-|--pcurrent
            // \ |          |
            //  \|          |
            //   \          |
            //  -A\---------B-
            //   | \        |
            //      vlh

            pnorm             = (pcurrent-p_iC) / (p_iB-p_iC);
            value_behind          =  storage_vector[iC] + pnorm*( storage_vector[iB] - storage_vector[iC] );

            tdummy       = t_iD;
            pnorm        = (pcurrent-p_iD)/(vlh_liquid.Pressure()-p_iD);
            value_iD         = storage_vector[iD];
            value_vlh_before  = vlh_liquid.ValueOf(property_index);
            value_before       = value_iD + pnorm*(value_vlh_before-value_iD);
            tnorm          = (tcurrent-t_iD) / (t_iB-t_iD);
            value_interpolated = value_before + tnorm*(value_behind-value_before);
          }
      }

    return value_interpolated;
    
  }



  double HaliteLiquidusLookup::NearNaClMeltInterpolation( const int& property_index )
  {
    //    cout << "Using HaliteLiquidusLookup::NearNaClMeltInterpolation( const int& property_index ) ...\n";

    pdummy = pcurrent;
    tvlh   = naclmelt_liquid.TmeltFromP();
    if(tvlh  <= t_iC)
      {

        tdummy    = t_iD;
        if(naclmelt_liquid.PmeltFromT() <= p_iA)
          {
            // case (1) and (2): Interpolation to be done between A-D and vlh_curve_at_pcurrent
            //
            //  (1)                                (2)                          
            //                                   
            //   P                                  P                             
            //   ^            meltcurve             ^         meltcurve          
            //   |         | /                      |        /                    
            //  -D---------C/-                     -D-------/-C-              
            //   |         |                        |      /  |           
            //   |        /|                        |     /   |        
            // --|--x----/-|------ pcurrent       --|-x--/----|------ pcurrent  
            //   |      /  |                        |   /     |         
            //  -A-----/---B--> T                  -A--/------B--> T      
            //   |    /    |                        |         |   
            //

            value_iA              = storage_vector[iA];
            value_iD              = storage_vector[iD];
	    
            pnorm             = (pcurrent-p_iD) / (p_iA-p_iD);
            value_before          = storage_vector[iD] + pnorm*( storage_vector[iA] - storage_vector[iD] );
            tdummy            = tvlh;
            value_vlh             = naclmelt_liquid.ValueOf(property_index);
            tnorm             = (tcurrent-t_iA) / (tvlh-t_iA);
            value_interpolated    =  value_before + tnorm*(value_vlh-value_before);
            //d  cout << "computed via if-2 as value_before = " << value_before << ", " << "value_vlh = " 
            //d  << value_vlh << ", " << value_interpolated << endl;
          }
        else
          {
            // case (3) and (4): Interpolation to be done between vlh_at_t_iD-D-segment 
            // and vlh_curve_at_pcurrent
            //
            //      (4)                              (2)
            //                                                               
            //   P                                P
            //   ^    meltcurve                   ^
            //   |   /                            |         |   meltcurve
            //  -D--/------C-                    -D---------C__/
            // --|x/-------|------ pcurrent       |     ___/| 
            //   |/        |                    --|-x__/----|------ pcurrent
            //   /         |                    __|/        |
            //  /|         |                   /  |         |
            //  -A---------B--> T                -A---------B--> T
            //   |         |                      |         |  
            //
	    
            // notice: it is still tdummy = t_iD;
            pnorm             = (pcurrent-p_iD) / (naclmelt_liquid.PmeltFromT()-p_iD);
            value_iD              = storage_vector[iD];
            value_vlh_before      = naclmelt_liquid.ValueOf(property_index); // (between A and D)
            tdummy            = tvlh;
            value_vlh             = naclmelt_liquid.ValueOf(property_index);
	    
            value_before          =  value_iD + pnorm*( value_vlh_before - value_iD );
            tnorm             = (tcurrent-t_iD) / (tvlh-t_iD);
            value_interpolated    =  value_before + tnorm*(value_vlh-value_before);
            //d  cout << "value_iD = " << value_iD << ", tnorm = " << tnorm << endl;
            //d  cout << "computed via if-1 as value_before = " << value_before << ", " 
            //d       << "value_vlh = " << value_vlh << ", " << value_interpolated << endl;
          }
      }
    else
      {
        // tvlh is > t_iB 
        tdummy    = t_iD;
        if(naclmelt_liquid.PmeltFromT() >= p_iA)
          {
            // case (5): Interpolation to be done between D-vlh_at_t_iD-segment and C-vlh_at_t_iC-segment
            //
            //  (5)                      
            //
            //   P                        
            //   ^                        
            //   |                        
            //  -D---------C ___  meltcurve  
            // --|-----x--_|/-----pcurrent           
            //   | _____/  |              
            // __|/        |              
            //   |         |     
            //  -A---------B--> T         
            //   |         |              
            //
            //
            pnorm             = (pcurrent-p_iD) / (naclmelt_liquid.PmeltFromT()-p_iD);
            value_iC              = storage_vector[iC];
            value_iD              = storage_vector[iD];
            value_vlh_before      = naclmelt_liquid.ValueOf(property_index); // (between A and D)
            tdummy            = t_iB;
            value_vlh_behind      = naclmelt_liquid.ValueOf(property_index);
	    
            value_before          =  value_iD + pnorm*( value_vlh_before - value_iD );
            tdummy            =  t_iC;
            pnorm             = (pcurrent-p_iC) / (naclmelt_liquid.PmeltFromT()-p_iC);
            value_behind          =  value_iC + pnorm*( value_vlh_behind - value_iC); 
            tnorm             = (tcurrent-t_iD) / (t_iB-t_iD);
            value_interpolated    =  value_before + tnorm*(value_behind-value_before);
          }
        else
          {
	
            //  (6)
            //
            //   P
            //   ^            meltcurve
            //   |           /
            //  -D--------C-/
            //  -|--x-----|/--- pcurrent
            //   |        /
            //   |       /|
            //   |      / |
            //  -A-----/--B--> T
            //   |    /   |        
            //
            pnorm             = (pcurrent-p_iD) / (p_iA-p_iD);
            value_iD              = storage_vector[iD];
            value_iC              = storage_vector[iC];
            value_iA              = storage_vector[iA];

            tdummy            = t_iC;
            value_vlh_behind      = naclmelt_liquid.ValueOf(property_index);

            value_before          =  value_iD + pnorm*( value_iA - value_iD );
            tdummy            =  t_iC;
            pnorm             = (pcurrent-p_iC) / (naclmelt_liquid.PmeltFromT()-p_iC);
            value_behind          =  value_iC + pnorm*( value_vlh_behind - value_iC); 
            tnorm             = (tcurrent-t_iD) / (t_iB-t_iD);
            value_interpolated    =  value_before + tnorm*(value_behind-value_before);
          }
      }
    return value_interpolated;
  }


  void HaliteLiquidusLookup::GetTemperatureIndex(const double& t)
  {
    // new version
    if(     t <= 250.0e0){  t_res =  5.0; it =     static_cast<long>( (t-  0.0)/t_res ); }
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
    
  void HaliteLiquidusLookup::GetPressureIndex(const double& p)
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
