#include <limits>

#include "TwophaseLiquidLookup.h"

#include "Brine.h"
#include "TwophaseLiquid.h"
#include "LookupPropertyIndex.h"

#include "binaryReadWrite.h"

#include <cmath>

using namespace std;

namespace csmp
{
  TwophaseLiquidLookup::TwophaseLiquidLookup(const double64& externaltemperature, 
                                             const double64& externalpressure)
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
      ph2o(0.0),
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
    value_bottom(0.0),
    value_top(0.0),
    value_before(0.0),
    value_behind(0.0),
    value_interpolated(0.0),
    value_crit(0.0),
    value_crit_behind(0.0),
    value_vlh(0.0),
    value_vlh_behind(0.0),
    value_vlh_before(0.0),
    value_boil(0.0),
    value_boil_behind(0.0),
    pboil_iA(0.0),
    pboil_iB(0.0),
    tboil(0.0),
    tcrit(0.0),
    tvlh(0.0),
    x_low(0.0),
    x_high(0.0),
    pcrit_iB(0.0),
    pcrit_iA(0.0),
    it(0),
    ip(0),
    t_dim(0),
    p_dim(0),
    iA(0),
    iB(0),
    iC(0),
    iD(0),
    i_dummy(0),
    eqtype(0),
    it_p_max(0),
    ip_p_max(0),
    state(none),
    state_iA(none),
    state_iB(none),
    state_iC(none),
    state_iD(none),
    critcurve(tdummy), 
    vlh_liquid(tdummy),
    water(),
    vlh_pmax(vlh_liquid.Pmax()),
    vlh_tmax(vlh_liquid.Tmax()),
    csmp_error( ErrorHandler::Instance())
  {

    // Prepare values that are necessary to identify cases where we are very close to Tmax,Pmax of the VLH surface
    GetTemperatureIndex(vlh_tmax);
    it_p_max          = it;
    GetPressureIndex(vlh_pmax);
    ip_p_max          = ip;
      
    // Setup table dimensions and rescale stroage vector sizes
    GetTemperatureIndex(1000.0001);
    t_dim             = it+1;
    GetPressureIndex(5000.0e5+0.001);
    p_dim             = ip+1;
    //d cout << "TwophaseLiquidLookup t_dim and p_dim = " << t_dim << "\t" << p_dim << endl;
	
    storage_vector.resize(t_dim*p_dim*max_index);
    state_vector.resize(  t_dim*p_dim);

    char                filename[60], statefilename[60];
    strcpy( filename, "TwophaseLiquidPropertiesLookupTable.bin" );
    strcpy( statefilename, "TwophaseLiquidStateLookupTable.bin" );
    FILE*               infile1;
    FILE*               infile2;
    infile1           = fopen( filename, "rb" );
    infile2           = fopen( statefilename, "rb" );
	
    if( (infile1 == NULL) || (infile2 == NULL) )
      {
        cout << "TwophaseLiquidLookup : at least one lookup file missing, computing ...\n\n";
        it = 0;
	    
        state    = L; // infor for brine only, don't change!!!
	    
        TwophaseLiquid                       twophase_l(tcurrent,pcurrent,ph2o);
        Brine                                brine(     tcurrent,pcurrent,xcurrent);

        //d cout << "TwophaseLiquidLookup::TwophaseLiquidLookup(...): computing subcritical data ...\n";
        for(tcurrent = 0.0e0; tcurrent < cp_h2o.Temperature(); tcurrent += t_res)
          {
            GetTemperatureIndex(tcurrent+1.0e-3);
            cout << "computing for t = " << tcurrent << ", data set " << it << endl << endl;

            ph2o    = water.SaturationPressureFromT(tcurrent);
            //d cout << "done ph2o stuff\n"; 
            if(it > t_dim-1)
              {
                cout << "too high t-index ...\n";
                break;
              }

            tdummy = tcurrent;
		
            for(pcurrent = 0.5e5; pcurrent <= 5000.0e5; pcurrent += p_res)
              {
                GetPressureIndex(pcurrent+1.0e-3);
                if(ip > p_dim-1)
                  {
                    cout << "too high p-index ...\n";
                    break;
                  }
		    
                if( pcurrent < vlh_liquid.Pressure() )
                  {
                    storage_vector[t_dim*p_dim*temperature_index     + it*p_dim + ip] = tcurrent;
                    storage_vector[t_dim*p_dim*pressure_index        + it*p_dim + ip] = pcurrent;
                    storage_vector[t_dim*p_dim*composition_index     + it*p_dim + ip] = 0.0e0;
                    storage_vector[t_dim*p_dim*density_index         + it*p_dim + ip] = 0.0e0;
                    storage_vector[t_dim*p_dim*enthalpy_index        + it*p_dim + ip] = 0.0e0;
                    storage_vector[t_dim*p_dim*heatcapacity_index    + it*p_dim + ip] = 0.0e0;
                    storage_vector[t_dim*p_dim*compressibility_index + it*p_dim + ip] = 0.0e0;
                    storage_vector[t_dim*p_dim*viscosity_index       + it*p_dim + ip] = 0.0e0;
                    state_vector[it*p_dim + ip ] = V;
                  }
                else if(pcurrent > water.SaturationPressureFromT(tcurrent))
                  {
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
                  {
                    xcurrent = twophase_l.MassFractionNaCl();
                    storage_vector[t_dim*p_dim*temperature_index     + it*p_dim + ip] = tcurrent;
                    storage_vector[t_dim*p_dim*pressure_index        + it*p_dim + ip] = pcurrent;
                    storage_vector[t_dim*p_dim*composition_index     + it*p_dim + ip] = xcurrent;
                    storage_vector[t_dim*p_dim*density_index         + it*p_dim + ip] = brine.Density();
                    storage_vector[t_dim*p_dim*enthalpy_index        + it*p_dim + ip] = brine.Enthalpy();
                    storage_vector[t_dim*p_dim*heatcapacity_index    + it*p_dim + ip] = brine.HeatCapacity();
                    storage_vector[t_dim*p_dim*compressibility_index + it*p_dim + ip] = brine.Compressibility();
                    storage_vector[t_dim*p_dim*viscosity_index       + it*p_dim + ip] = brine.Viscosity();
                    state_vector[it*p_dim + ip ] = VL;
                  }
              }
          }
        //d cout << "TwophaseLiquidLookup::TwophaseLiquidLookup(...): computing supercritical data ...\n";
	    
        // ***** check if the usage of 374.0 causes problems
        for(tcurrent = 374.0e0; tcurrent < 1000.1e0; tcurrent += t_res)
          {
            GetTemperatureIndex(tcurrent+1.0e-3);
            cout << "computing for t = " << tcurrent << ", data set " << it << endl << endl;
            if(it > t_dim-1)
              {
                cout << "too high t-index ...\n";
                break;
              }
            tdummy = tcurrent; // for critcurve and vlh_liquid_lookup
		
            for(pcurrent = 1.0e5; pcurrent <= 5000.0e5; pcurrent += p_res)
              {
                GetPressureIndex(pcurrent+1.0e-3);
                if(ip > p_dim-1)
                  {
                    cout << "too high p-index ...\n";
                    break;
                  }
		    
                if( (tcurrent <=800.7) && (pcurrent < vlh_liquid.Pressure()) )
                  {
                    storage_vector[t_dim*p_dim*temperature_index     + it*p_dim + ip] = tcurrent;
                    storage_vector[t_dim*p_dim*pressure_index        + it*p_dim + ip] = pcurrent;
                    storage_vector[t_dim*p_dim*composition_index     + it*p_dim + ip] = 0.0e0;
                    storage_vector[t_dim*p_dim*density_index         + it*p_dim + ip] = 0.0e0;
                    storage_vector[t_dim*p_dim*enthalpy_index        + it*p_dim + ip] = 0.0e0;
                    storage_vector[t_dim*p_dim*heatcapacity_index    + it*p_dim + ip] = 0.0e0;
                    storage_vector[t_dim*p_dim*compressibility_index + it*p_dim + ip] = 0.0e0;
                    storage_vector[t_dim*p_dim*viscosity_index       + it*p_dim + ip] = 0.0e0;
                    state_vector[it*p_dim + ip ] = V;
                  }
                else if(pcurrent > critcurve.Pressure())
                  {
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
                  {
                    xcurrent = twophase_l.MassFractionNaCl();
                    storage_vector[t_dim*p_dim*temperature_index     + it*p_dim + ip] = tcurrent;
                    storage_vector[t_dim*p_dim*pressure_index        + it*p_dim + ip] = pcurrent;
                    storage_vector[t_dim*p_dim*composition_index     + it*p_dim + ip] = xcurrent;
                    storage_vector[t_dim*p_dim*density_index         + it*p_dim + ip] = brine.Density();
                    storage_vector[t_dim*p_dim*enthalpy_index        + it*p_dim + ip] = brine.Enthalpy();
                    storage_vector[t_dim*p_dim*heatcapacity_index    + it*p_dim + ip] = brine.HeatCapacity();
                    storage_vector[t_dim*p_dim*compressibility_index + it*p_dim + ip] = brine.Compressibility();
                    storage_vector[t_dim*p_dim*viscosity_index       + it*p_dim + ip] = brine.Viscosity();
                    state_vector[it*p_dim + ip ] = VL;
                  }
              }
          }
	    
        FILE* outfile1;
        outfile1 = fopen( filename, "wb");
        if( outfile1 == NULL )
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
            fclose(outfile1);
            cout << "done!\n";
          }
	    
        FILE* outfile2;
        outfile2 = fopen( statefilename, "wb");
        if( outfile2 == NULL )
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
            fclose(outfile2);
            cout << "done!\n";
          }
	    
      }
	
    else
      {
        cout << "reading file " << filename << " ... ";
        skm_C_fread( infile1, storage_vector );
        fclose( infile1 );
        cout << "done!\n";
	    
        cout << "reading file " << statefilename << " ... ";
        skm_C_fread( infile2, state_vector );
        fclose( infile2 );
        cout << "done!\n";
	    
      }
    //d cout << "TwophaseLiquidLookup, leaving constructor ...\n\n";
  }
    
    
  TwophaseLiquidLookup::~TwophaseLiquidLookup()
  {
  }
    
    
    
  // The data interpolation routines
  double64 TwophaseLiquidLookup::Temperature(){     SetTemperatureAndPressure(); return Interpolate(temperature_index);    }
  double64 TwophaseLiquidLookup::Pressure(){        SetTemperatureAndPressure(); return Interpolate(pressure_index);       }
  double64 TwophaseLiquidLookup::MassFractionNaCl(){     SetTemperatureAndPressure(); return Interpolate(composition_index);    }
  double64 TwophaseLiquidLookup::Density(){         SetTemperatureAndPressure(); return Interpolate(density_index);        }
  double64 TwophaseLiquidLookup::Enthalpy(){        SetTemperatureAndPressure(); return Interpolate(enthalpy_index);       }
  double64 TwophaseLiquidLookup::HeatCapacity(){    SetTemperatureAndPressure(); return Interpolate(heatcapacity_index);   }
  double64 TwophaseLiquidLookup::Compressibility(){ SetTemperatureAndPressure(); return Interpolate(compressibility_index);}
  double64 TwophaseLiquidLookup::Viscosity(){       SetTemperatureAndPressure(); return Interpolate(viscosity_index);      }
    
  double64 TwophaseLiquidLookup::ReportMassFractionNaCl(){ return Interpolate(composition_index); }
  double64 TwophaseLiquidLookup::ReportEnthalpy(){    return Interpolate(enthalpy_index); }
    
  double64 TwophaseLiquidLookup::MinXResolution(){ return 1.0e-4; } //min x-res in H2ONaClLookup

  double64 TwophaseLiquidLookup::DEnthalpyDT()
  {
    SetTemperatureAndPressure(); 
    tdummy   = tcurrent;
	
    if(pcurrent >= vlh_pmax)
      {
        //d cout << "TwophaseLiquidLookup::DEnthalpyDT() - path 1\n";
        tcurrent += 0.1;
        x_high   = ReportEnthalpy();
        tcurrent = tdummy;
        return (x_high-ReportEnthalpy())/0.1;
      }
    else
      {
        if(tcurrent < vlh_tmax)
          {
            tvlh = vlh_liquid.TfromP(pcurrent,300.0);
            if(tcurrent < tvlh-0.101)
              {
                //d cout << "TwophaseLiquidLookup::DEnthalpyDT() - path 2\n";
                tcurrent += 0.1;
                x_high   = ReportEnthalpy();
                tcurrent = tdummy;
                return (x_high-ReportEnthalpy())/0.1;
              }
            else
              {
                //d cout << "TwophaseLiquidLookup::DEnthalpyDT() - path 3\n";
                tcurrent -= 0.1;
                x_low    = ReportEnthalpy();
                tcurrent = tdummy;
                return (ReportEnthalpy()-x_low)/0.1;
              }
          }
        else
          {
            tvlh = vlh_liquid.TfromP(pcurrent,700.0);
            if(tcurrent > tvlh+0.101)
              {
                //d cout << "TwophaseLiquidLookup::DEnthalpyDT() - path 4\n";
                tcurrent -= 0.1;
                x_low   = ReportEnthalpy();
                tcurrent = tdummy;
                return (ReportEnthalpy()-x_low)/0.1;
              }
            else
              {
                //d cout << "TwophaseLiquidLookup::DEnthalpyDT() - path 5\n";
                tcurrent += 0.1;
                x_high    = ReportEnthalpy();
                tcurrent = tdummy;
                return (x_high-ReportEnthalpy())/0.1;
              }
          }
      } 
  }
    
    
  double64 TwophaseLiquidLookup::DCompositionDT()
  {
    SetTemperatureAndPressure(); 
    tdummy   = tcurrent;
	
    if(pcurrent >= vlh_pmax)
      {
        //d cout << "TwophaseLiquidLookup::DCompositionDT() - path 1\n";
        tcurrent += 0.1;
        x_high   = ReportMassFractionNaCl();
        tcurrent = tdummy;
        return (x_high-ReportMassFractionNaCl())/0.1;
      }
    else
      {
        if(tcurrent < vlh_tmax)
          {
            tvlh = vlh_liquid.TfromP(pcurrent,300.0);
            if(tcurrent < tvlh-0.101)
              {
                //d cout << "TwophaseLiquidLookup::DCompositionDT() - path 2\n";
                tcurrent += 0.1;
                x_high   = ReportMassFractionNaCl();
                tcurrent = tdummy;
                return (x_high-ReportMassFractionNaCl())/0.1;
              }
            else
              {
                //d cout << "TwophaseLiquidLookup::DCompositionDT() - path 3\n";
                tcurrent -= 0.1;
                x_low    = ReportMassFractionNaCl();
                tcurrent = tdummy;
                return (ReportMassFractionNaCl()-x_low)/0.1;
              }
          }
        else
          {
            tvlh = vlh_liquid.TfromP(pcurrent,700.0);
            if(tcurrent > tvlh+0.101)
              {
                //d cout << "TwophaseLiquidLookup::DCompositionDT() - path 4\n";
                tcurrent -= 0.1;
                x_low   = ReportMassFractionNaCl();
                tcurrent = tdummy;
                return (ReportMassFractionNaCl()-x_low)/0.1;
              }
            else
              {
                //d cout << "TwophaseLiquidLookup::DCompositionDT() - path 5\n";
                tcurrent += 0.1;
                x_high    = ReportMassFractionNaCl();
                tcurrent = tdummy;
                return (x_high-ReportMassFractionNaCl())/0.1;
              }
          }
      } 
  }
    
    
  double64 TwophaseLiquidLookup::DSaltMassFractionDT()
  {
    SetTemperatureAndPressure(); 
    tdummy   = tcurrent;
	
    if(pcurrent >= vlh_pmax)
      {
        tcurrent += 0.1;
        x_high   = XNaCl2Massfraction(ReportMassFractionNaCl());
        tcurrent = tdummy;
        return (x_high-XNaCl2Massfraction(ReportMassFractionNaCl()))/0.1;
      }
    else
      {
        if(tcurrent < vlh_tmax)
          {
            tvlh = vlh_liquid.TfromP(pcurrent,300.0);
            if(tcurrent < tvlh-0.101)
              {
                tcurrent += 0.1;
                x_high   = XNaCl2Massfraction(ReportMassFractionNaCl());
                tcurrent = tdummy;
                return (x_high-XNaCl2Massfraction(ReportMassFractionNaCl()))/0.1;
              }
            else
              {
                tcurrent -= 0.1;
                x_low    = XNaCl2Massfraction(ReportMassFractionNaCl());
                tcurrent = tdummy;
                return (XNaCl2Massfraction(ReportMassFractionNaCl())-x_low)/0.1;
              }
          }
        else
          {
            tvlh = vlh_liquid.TfromP(pcurrent,700.0);
            if(tcurrent > tvlh+0.101)
              {
                tcurrent -= 0.1;
                x_low   = XNaCl2Massfraction(ReportMassFractionNaCl());
                tcurrent = tdummy;
                return (XNaCl2Massfraction(ReportMassFractionNaCl())-x_low)/0.1;
              }
            else
              {
                tcurrent += 0.1;
                x_high    = XNaCl2Massfraction(ReportMassFractionNaCl());
                tcurrent = tdummy;
                return (x_high-XNaCl2Massfraction(ReportMassFractionNaCl()))/0.1;
              }
          }
      } 
  }
    
    
    
  void TwophaseLiquidLookup::SetTemperatureAndPressure()
  {
    tcurrent = temperature;
    pcurrent = pressure;
    tdummy   = tcurrent;
    pdummy   = pcurrent;
    return;
  }
    
    
  void TwophaseLiquidLookup::GetIndex_iA(const int& property_index)
  {
    // CAUTION: tcurrent and pcurrent must be know before this function is called !!!
    GetTemperatureIndex(tcurrent);
    GetPressureIndex(pcurrent);
    iA  = t_dim*p_dim*property_index + it*p_dim + ip; //it*p_dim*(property_index+1) + ip;
    iB  = iA + p_dim;
    iC  = iB + 1;
    iD  = iA + 1;
	
    return;
  }
    

    
  double64 TwophaseLiquidLookup::Interpolate(const int& property_index)
  {
    // do NOT set tcurrent and pcurrent here, do it outside this function too keep it versatile !!!
	
    GetIndex_iA(property_index);
    //d cout << "\n\nTwophaseLiquidLookup::Interpolate for property_index " << property_index << endl;
    //d cout << "it = " << it << ", ip = " << ip << endl;

    //CheckForOutOfRange();
    state_iA = state_vector[ it*p_dim + ip ];
    state_iB = state_vector[ it*p_dim + ip + p_dim ];
    state_iC = state_vector[ it*p_dim + ip + p_dim +1 ];
    state_iD = state_vector[ it*p_dim + ip + 1 ];
	
    i_dummy  = t_dim*p_dim*temperature_index + it*p_dim + ip;
    //d cout << "i_dummy t = " << i_dummy << endl;
    t_iA     = storage_vector[i_dummy];
    t_iB     = storage_vector[i_dummy + p_dim];
    t_iC     = storage_vector[i_dummy +1 + p_dim];
    t_iD     = storage_vector[i_dummy + 1];
	
    i_dummy  = t_dim*p_dim*pressure_index + it*p_dim + ip;
    //d cout << "i_dummy p = " << i_dummy << endl;
    p_iA     = storage_vector[i_dummy];
    p_iB     = storage_vector[i_dummy + p_dim];
    p_iD     = storage_vector[i_dummy + 1];
    p_iC     = p_iD;
	
    //d cout << "Temperatures and Pressures iA-iD:\n";
    //d cout << t_iA << "\t" << t_iB << "\t" << t_iC << "\t" << t_iD << endl;
    //d cout << p_iA << "\t" << p_iB << "\t" << p_iC << "\t" << p_iD << endl;

    if( it == it_p_max && (ip == ip_p_max || ip == ip_p_max-1/*+1*/) )
      {
        eqtype = 99;
        return NearVLHMaxInterpolation(property_index);
      }

    // Now decisions based on topology
    //   P
    //   ^            critcurve
    //   |         | /
    //  -D---------C/-
    //   |         |
    //   |        /|
    //   |       / |
    //   |      /x |
    //  -A-----/---B--> T
    //        / 
    if( (state_iA == VL) && (state_iB == VL) && (state_iC == VL) && (state_iD == VL) )
      {
        //d cout << state_iA << "\t" << state_iB << "\t" << state_iC << "\t" << state_iD << endl;
        eqtype = 0;
        return NormalInterpolation(property_index); 
      }
	
    else
      {
        //d cout << state_iA << "\t" << state_iB << "\t" << state_iC << "\t" << state_iD << endl;
	
        // *** CAUTION: this applies to current lookup spacing only!!!
        if( tcurrent >= 373.9 && tcurrent <= 374.0 && pcurrent >= 220.5e5 && pcurrent <= 220.6e5 )
          {
            return NearCritpointInterpolation(property_index);
          }

        else if(tcurrent >= cp_h2o.Temperature())
          {
            eqtype = 1;
            if(      state_iD == F && state_iB == VL )                  return NearCritcurveInterpolation(property_index);
            else if( state_iB == V && state_iD == VL && it < it_p_max ) return NearVLHInterpolationLowT(property_index);
            else if( state_iA == V && state_iC == VL && it > it_p_max ) return NearVLHInterpolationHighT(property_index);
            else{
              eqtype = 2;
              //d cout << "TwophaseLiquidLookup::missed to find interpolation type in case 2, returning bogus value ...\n";
              return 9.9e99;
            }
          }
        else
          {
            // p_iD, p_iA, t_iA und t_iB muessen bekannt sein
            if(      state_iD == F && state_iB == VL ) return NearBoilingCurveInterpolation(property_index); // *
            else if( state_iB == V && state_iD == VL ) return NearVLHInterpolationLowT(property_index);
            else if( state_iD == F && state_iB == V 
                     && pcurrent <= water.SaturationPressureFromT(tcurrent)
                     && pcurrent >= vlh_liquid.Pressure() ) return InterpolateBetweenBoilingCurveAndVLH(property_index);
            else
              {
                eqtype = 3;
                //d cout << "TwophaseLiquidLookup::missed to find interpolation type in case 3, returning bogus value ...\n";
                return 9.9e99;
              }
          }
      }
  }
    
  double64 TwophaseLiquidLookup::NormalInterpolation( const int& property_index )
  {
    //d cout << "Using TwophaseLiquidLookup::NormalInterpolation( const int& property_index ) ...\n";
	
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
    //    //d cout << "values at iA-iD  = " << storage_vector[iA] << "\t" <<  storage_vector[iB] << "\t" << storage_vector[iC] << "\t" << storage_vector[iD] << endl;
    return value_interpolated;
  }
    
    
  double64 TwophaseLiquidLookup::NearVLHMaxInterpolation(const int& property_index)
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
    //   -D/A---/----\------C/B- ip_p_max //wrong: +1
    //     | _/       V\     |
    //     |/           L\   |
    //   _/|              H\ |  
    //  /  |       VH        \ 
    //     |                 | \
    //   --A-----------------B-- ip_p_max // needs:-1
    //     |                 |
    //   it_p_max
    //     

    // iA etc. have already been determined when this function is called
    if( ip == ip_p_max/*+1*/ ) // in this case, ABDC in the above figure are moved up by one cell
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
                csmp_error.notice( CSMP_FATAL_ERROR, 
                                   "TwophaseLiquidLookup::NearVLHMaxInterpolation(const int& property_index) - ",
                                   "p < vlh_pmax (in cell ip_p_max+1) but neither tcurrent<vlh_tmax nor tcurrent>vlh_tmax condition worked!\nReport issue to Thomas Driesner, thomas.driesner@erdw.ethz.ch"); 
                return std::numeric_limits<double>::quiet_NaN();
              }
          }
      }
    else if( ip == ip_p_max-1 )
      {
        if(tcurrent < vlh_tmax )
          {
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
            // first, value_before @vlh curve at p
            tdummy             =  vlh_liquid.TfromP(pcurrent,700.); // tdummy id the one that vlh_liquid uses, see constructor!
            tvlh               =  tdummy;
            value_before       =  vlh_liquid.ValueOf(property_index);
            // then value_behind 
            tdummy             =  t_iB;	    
            pnorm              = (pcurrent - vlh_liquid.Pressure()) / (p_iC - vlh_liquid.Pressure());
            value_behind       =  vlh_liquid.ValueOf(property_index) + pnorm*(storage_vector[iC]-vlh_liquid.ValueOf(property_index));
            // then interpolate
            tnorm              = (tcurrent - tvlh) / (t_iB - tvlh);
            value_interpolated =  value_before + tnorm*(value_behind-value_before);
            return                value_interpolated;		
          }
        else
          {
            csmp_error.notice( CSMP_FATAL_ERROR, 
                               "TwophaseLiquidLookup::NearVLHMaxInterpolation(const int& property_index) - ",
                               "p < vlh_pmax (in cell ip_p_max) but neither tcurrent<vlh_tmax nor tcurrent>vlh_tmax condition worked!\nReport issue to Thomas Driesner, thomas.driesner@erdw.ethz.ch"); 
            return std::numeric_limits<double>::quiet_NaN();
          }
      }
    else
      {
        csmp_error.notice( CSMP_FATAL_ERROR, 
                           "TwophaseLiquidLookup::NearVLHMaxInterpolation(const int& property_index) - ",
                           "Missed ALL if-statements!\nReport issue to Thomas Driesner, thomas.driesner@erdw.ethz.ch"); 
 
        return std::numeric_limits<double>::quiet_NaN();
     }
  }

  double64 TwophaseLiquidLookup::NearCritpointInterpolation( const int& property_index )
  {
    // This is  the topology:
    //                           
    //            critcurve         
    //   |        /                   
    //  -D-------/-C- 220.6 bar
    //   |      /  |                  
    //   |     /   |                  
    // --|----/----|------ pcrit_h2o   
    //   |   /     |                  
    //  -A--/------B- 220.5 bar        
    //   |         |                   
    // 373.9 C    374 C

    if( pcurrent == cp_h2o.Pressure() )
      {
        tdummy             =  cp_h2o.Temperature();	
        value_before       =  critcurve.ValueOf(property_index); // or create critprop_value vector?
        pnorm              = (cp_h2o.Pressure()-p_iB) / (p_iC-p_iB);
        value_behind       =  storage_vector[iB] + pnorm*(storage_vector[iC] - storage_vector[iB]);
        tnorm              = (tcurrent-cp_h2o.Temperature())/(t_iB-cp_h2o.Temperature());
        value_interpolated =  value_before+tnorm*(value_behind-value_before); 	
      }
    else if( pcurrent > cp_h2o.Pressure() ) // implies that also tcurrent > cp_h2o.Temperature(), must be made sure before getting here
      {
        tdummy             =  critcurve.TfromP(pcurrent);
        value_before       =  critcurve.ValueOf(property_index);
        pnorm              = (pcurrent-p_iB) / (p_iC-p_iB);
        value_behind       =  storage_vector[iB] + pnorm*(storage_vector[iC] - storage_vector[iB]);
        tnorm              = (tcurrent-tdummy) / (t_iB-tdummy);
        value_interpolated =  value_before+tnorm*(value_behind-value_before); 	
      } 
    else
      {
        tdummy             =  water.SaturationTemperatureFromP(pcurrent);
        value_boil         =  water.LiquidProperty(tdummy,property_index);
        pnorm              = (pcurrent-p_iB) / (p_iC-p_iB);
        value_behind       =  storage_vector[iB] + pnorm*( storage_vector[iC] - storage_vector[iB] );
        tnorm              = (tcurrent-tdummy) / (t_iB-tdummy);
        value_interpolated =  value_boil+tnorm*(value_behind-value_boil); 
      }
    return value_interpolated;

  }


  double64 TwophaseLiquidLookup::NearCritcurveInterpolation( const int& property_index )
  {
    //d cout << "Using TwophaseLiquidLookup::NearCritcurveInterpolation( const int& property_index ) ...\n";
    //d cout << "Critical Pressure at T of interest is " << critcurve.Pressure() << endl;
	
    // a few data needed to determine topology
    tcrit     = critcurve.TfromP(pcurrent);
    tdummy    = t_iB;
    pcrit_iB  = critcurve.Pressure();
    tdummy    = tcurrent;
	
    // There is a number of topologies possible (x is the value to be interpolated):
	
    if(tcrit >= t_iA)
      {

        if(pcrit_iB <= p_iC)
          { 
            //d cout << "TwophaseLiquidLookup::NearCritcurveInterpolation( const int& property_index ) ... eqtype 4\n";
            // case (1) and (2): Interpolation to be done between critcurve_at_pcurrent and B-critcurve_at_t_iB-segment
            //
            //  (1)                               (2)
            //
            //   P                                 P
            //   ^            critcurve            ^
            //   |         | /                     |         |   critcurve
            //  -D---------C/-                    -D---------C__/
            //   |         |                       |     ___/| 
            //   |        /|                     --|--__/--x-|------ pcurrent
            // --|-------/x|------ pcurrent      __|/        |
            //   |      /  |                    /  |         |
            //  -A-----/---B--> T                 -A---------B--> T
            //   |    /    |                       |         |
            //
            eqtype = 4;
            pnorm             = (pcurrent-p_iB) / (pcrit_iB-p_iB);
            value_iB              = storage_vector[iB];
		
            //this was missing(?)
            tdummy = t_iB;
            //
            value_crit_behind     = critcurve.ValueOf(property_index);
            tdummy            = tcrit;
            value_crit            = critcurve.ValueOf(property_index);

            value_behind          =  value_iB + pnorm*( value_crit_behind - value_iB );
            tnorm             = (tcurrent-tcrit) / (t_iB-tcrit);
            value_interpolated    =  value_crit + tnorm*(value_behind - value_crit);
          }
        else
          {
            //d cout << "TwophaseLiquidLookup::NearCritcurveInterpolation( const int& property_index ) ... eqtype 5\n";
		
            // case (3) and (4): Interpolation to be done between critcurve_at_pcurrent and B-C
            //
            //  (3)                               (4)
            //
            //   P                                 P
            //   ^         critcurve               ^    critcurve
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
            eqtype = 5;
            pnorm          = (pcurrent-p_iB) / (p_iC-p_iB);
            value_behind       =  storage_vector[iB] + pnorm*( storage_vector[iC] - storage_vector[iB] );
            tnorm          = (tcurrent-tcrit) / (t_iB-tcrit);
            tdummy         = tcrit;
            value_crit         = critcurve.ValueOf(property_index);
            value_interpolated =  value_crit + tnorm*(value_behind - value_crit);   
          }
      }
    else
      { // tcrit_pcurrent is < T-iA
	    
        if(pcrit_iB <= p_iC)
          {
            //d cout << "TwophaseLiquidLookup::NearCritcurveInterpolation( const int& property_index ) ... eqtype 6\n";
            // case (5): Interpolation to be done between A-critcurve_at_iA-segment and B-critcurve_at_t_iB-segment
            //
            //  (5)                      
            //
            //   P                        
            //   ^                        
            //   |                        
            //  -D---------C   critcurve  
            //   |       __|__/           
            //   | _____/  |              
            // __|/        |              
            //   |------x--|-pcurrent     
            //  -A---------B--> T         
            //   |         |              
            //
            //
            eqtype = 6;
            pnorm  = (pcurrent-p_iB) / (pcrit_iB-p_iB);
            value_iB   = storage_vector[iB];
            tdummy = t_iB;            
		
            value_crit_behind = critcurve.ValueOf(property_index);
            value_behind          =  value_iB + pnorm*( value_crit_behind - value_iB );
            // needs update


            tdummy    = t_iA;
            pcrit_iA  = critcurve.Pressure();
            pnorm     = (pcurrent-p_iA) / (pcrit_iA-p_iA);
            value_iA      = storage_vector[iA];
		
            value_crit    = critcurve.ValueOf(property_index);
            value_before           =  value_iA + pnorm*( value_crit - value_iA );
            tnorm             = (tcurrent-t_iA) / (t_iB-t_iA);
            value_interpolated    =  value_before + tnorm*(value_behind - value_before);
          }
        else
          {
            //d cout << "TwophaseLiquidLookup::NearCritcurveInterpolation( const int& property_index ) ... eqtype 7\n";
            // case (6): Interpolation to be done between A-critcurve_at_iA-segment and B-C 
            //
            //  (6)
            //
            //   P
            //   ^    critcurve
            //   |   /        
            //  -D--/------C-
            //   | /       | 
            //   |/        |
            //   /         |
            // -/|------x--|-pcurrent
            //  -A---------B--> T
            //   |         |        
            //
            eqtype = 7;
            pnorm          = (pcurrent-p_iB) / (p_iC-p_iB);
            value_behind       =  storage_vector[iB] + pnorm*( storage_vector[iC] - storage_vector[iB] );
            tnorm          = (tcurrent-t_iA) / (t_iB-t_iA);

            tdummy    = t_iA;
            pcrit_iA  = critcurve.Pressure();
            pnorm     = (pcurrent-p_iA) / (pcrit_iA-p_iA);
            value_iA      = storage_vector[iA];

            value_crit    = critcurve.ValueOf(property_index);
            value_before          =  value_iA + pnorm*( value_crit - value_iA );
            tnorm             = (tcurrent-t_iA) / (t_iB-t_iA);
            value_interpolated    =  value_before + tnorm*(value_behind - value_before);
          }
      }

    //d    cout << "value_behind = " << value_behind << ", value_crit = " << value_crit << ", value_interpolated = " << value_interpolated << endl;

    return value_interpolated;
  }



  double64 TwophaseLiquidLookup::NearBoilingCurveInterpolation( const int& property_index )
  {
    //d cout << "Using TwophaseLiquidLookup::NearBoilingCurveInterpolation( const int& property_index ) ...\n";


    // a few data needed to determine topology
    tboil     = water.SaturationTemperatureFromP(pcurrent);
    tdummy    = t_iB;
    pboil_iB  = water.SaturationPressureFromT(t_iB);
    tdummy    = tcurrent;

    // There is a number of topologies possible (x is the value to be interpolated):

    if(tboil >= t_iA)
      {
	
        if(pboil_iB <= p_iC)
          { 
            //d cout << "TwophaseLiquidLookup::NearBoilingCurveInterpolation( const int& property_index ) ... eqtype 8\n";
            // case (1) and (2): Interpolation to be done between boilcurve_at_pcurrent and B-boilcurve_at_t_iB-segment
            //
            //  (1)                               (2)
            //
            //   P                                 P
            //   ^            boilcurve            ^
            //   |         | /                     |         |   boilcurve
            //  -D---------C/-                    -D---------C__/
            //   |         |                       |     ___/| 
            //   |        /|                     --|--__/--x-|------ pcurrent
            // --|-------/x|------ pcurrent      __|/        |
            //   |      /  |                    /  |         |
            //  -A-----/---B--> T                 -A---------B--> T
            //   |    /    |                       |         |
            //
            eqtype = 8;
            pnorm             = (pcurrent-p_iB) / (pboil_iB-p_iB);
            value_iB              = storage_vector[iB];
	    
            // this was missing, I guess
            tdummy = t_iB;
            value_boil_behind     = water.LiquidProperty(t_iB,  property_index);
            value_boil            = water.LiquidProperty(tboil, property_index);

            value_behind          =  value_iB + pnorm*( value_boil_behind - value_iB );
            tnorm             = (tcurrent-tboil) / (t_iB-tboil);
            value_interpolated    =  value_boil + tnorm*(value_behind - value_boil);
          }
        else
          {
            //d cout << "TwophaseLiquidLookup::NearBoilingCurveInterpolation( const int& property_index ) ... eqtype 9\n";
            // case (3) and (4): Interpolation to be done between boiling_liquid_at_pcurrent and B-C
            //
            //  (3)                               (4)
            //
            //   P                                 P
            //   ^         boilcurve               ^    boilcurve
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
            eqtype = 9;
            pnorm          = (pcurrent-p_iB) / (p_iC-p_iB);
            value_behind       =  storage_vector[iB] + pnorm*( storage_vector[iC] - storage_vector[iB] );
            tnorm          = (tcurrent-tboil) / (t_iB-tboil);
            tdummy         = tboil;
            value_boil         = water.LiquidProperty(tboil,property_index);
	    
            value_interpolated =  value_boil + tnorm*(value_behind - value_boil);    
          }
      }
    else
      { // tboil_pcurrent is < T-iA
 
        if(pboil_iB <= p_iC)
          {
            //d cout << "TwophaseLiquidLookup::NearBoilingCurveInterpolation( const int& property_index ) ... eqtype 10\n";
            // case (5): Interpolation to be done between A-boiling_liquid_at_iA-segment and B-boiling_liquid_at_t_iB-segment
            //
            //  (5)                      
            //
            //   P                        
            //   ^                        
            //   |                        
            //  -D---------C   boilcurve  
            //   |       __|__/           
            //   | _____/  |              
            // __|/        |              
            //   |------x--|-pcurrent     
            //  -A---------B--> T         
            //   |         |              
            //
            //
            eqtype = 10;
            pnorm  = (pcurrent-p_iB) / (pboil_iB-p_iB);
            value_iB   = storage_vector[iB];
            tdummy = t_iB;            
	    
            value_boil_behind     = water.LiquidProperty(t_iB,property_index);
            value_behind          =  value_iB + pnorm*( value_boil_behind - value_iB );

            tdummy    = t_iA;
            pboil_iA  = water.SaturationPressureFromT(t_iA);
            pnorm     = (pcurrent-p_iA) / (pboil_iA-p_iA);
            value_iA      = storage_vector[iA];

            value_boil            = water.LiquidProperty(t_iA,property_index);
            value_before          =  value_iA + pnorm*( value_boil - value_iA );
            tnorm             = (tcurrent-t_iA) / (t_iB-t_iA);
            value_interpolated    =  value_before + tnorm*(value_behind - value_before);
          }
        else
          {
            //d cout << "TwophaseLiquidLookup::NearBoilingCurveInterpolation( const int& property_index ) ... eqtype 11\n";
            // case (6): Interpolation to be done between A-boiling_liquid_at_iA-segment and B-C
            //
            //  (6)
            //
            //   P
            //   ^    boilcurve
            //   |   /        
            //  -D--/------C-
            //   | /       | 
            //   |/        |
            //   /         |
            // -/|------x--|-pcurrent// needs update

            //  -A---------B--> T
            //   |         |        
            //
            eqtype = 11;
            pnorm          = (pcurrent-p_iB) / (p_iC-p_iB);
            value_behind       =  storage_vector[iB] + pnorm*( storage_vector[iC] - storage_vector[iB] );
            tnorm          = (tcurrent-t_iA) / (t_iB-t_iA);

            tdummy         = t_iA;
            pboil_iA       = water.SaturationPressureFromT(t_iA);
            pnorm          = (pcurrent-p_iA) / (pboil_iA-p_iA);
            value_iA           = storage_vector[iA];

            value_boil         = water.LiquidProperty(t_iA,property_index);
	    
            value_before           =  value_iA + pnorm*( value_boil - value_iA );
            tnorm             = (tcurrent-t_iA) / (t_iB-t_iA);
            value_interpolated    =  value_before + tnorm*(value_behind - value_before);
          }
      }

    //d    cout << "value_behind = " << value_behind << ", value_boil = " << value_boil << ", value_interpolated = " << value_interpolated << endl;

    return value_interpolated;
  }







  double64 TwophaseLiquidLookup::NearVLHInterpolationLowT( const int& property_index )
  {
    //d cout << "Using TwophaseLiquidLookup::NearVLHInterpolationLowT( const int& property_index ) ...\n";

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

            eqtype = 12;
            value_iA              = storage_vector[iA];
            value_iD              = storage_vector[iD];
	    
            pnorm             = (pcurrent-p_iD) / (p_iA-p_iD);
            value_before          = storage_vector[iD] + pnorm*( storage_vector[iA] - storage_vector[iD] );
            tdummy            = tvlh;
            value_vlh             = vlh_liquid.ValueOf(property_index);
            tnorm             = (tcurrent-t_iA) / (tvlh-t_iA);
            value_interpolated    =  value_before + tnorm*(value_vlh-value_before);
            //d cout << "computed via if-2 as value_before = " << value_before << ", " << "value_vlh = " 
            //d << value_vlh << ", " << value_interpolated << endl;
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
            eqtype = 13;
            pnorm             = (pcurrent-p_iD) / (vlh_liquid.Pressure()-p_iD);
            value_iD              = storage_vector[iD];
            value_vlh_before      = vlh_liquid.ValueOf(property_index);
            tdummy            = tvlh;
            value_vlh             = vlh_liquid.ValueOf(property_index);
            value_before          =  value_iD + pnorm*( value_vlh_before - value_iD );
            tnorm             = (tcurrent-t_iD) / (tvlh-t_iD);
            value_interpolated    =  value_before + tnorm*(value_vlh-value_before);
            //d cout << "value_iD = " << value_iD << ", tnorm = " << tnorm << endl;
            //d cout << "computed via if-1 as value_before = " << value_before << ", " 
            //d      << "value_vlh = " << value_vlh << ", " << value_interpolated << endl;
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
            eqtype = 14;
            pnorm             = (pcurrent-p_iD) / (vlh_liquid.Pressure()-p_iD);
            value_iC              = storage_vector[iC];
            value_iD              = storage_vector[iD];

            value_vlh_before  = vlh_liquid.ValueOf(property_index); // (between A and D)
            tdummy        = t_iB;
            value_vlh_behind  = vlh_liquid.ValueOf(property_index);
    
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
            eqtype = 15;
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



  // needs update
  double64 TwophaseLiquidLookup::NearVLHInterpolationHighT( const int& property_index )
  {
    //d cout << "Using TwophaseLiquidLookup::NearVLHInterpolationHighT( const int& property_index ) ...\n";

    tvlh      = vlh_liquid.TfromP(pcurrent,700.0e0);

    if(tvlh  >= t_iA)
      {
        //d cout << "cases 1/2\n";
        tdummy    = t_iB;
        if(vlh_liquid.Pressure() <= p_iB)
          {
	    
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
            eqtype = 16;
            pnorm          = (pcurrent-p_iB) / (p_iC-p_iB);
            value_behind       =  storage_vector[iB] + pnorm*( storage_vector[iC] - storage_vector[iB] );
            tnorm          = (tcurrent-tvlh) / (t_iB-tvlh);
            tdummy         =  tvlh;
            value_vlh          =  vlh_liquid.ValueOf(property_index);
	    
            value_interpolated =  value_vlh + tnorm*(value_behind - value_vlh);   
          }
        else
          {
            //d cout << "cases 3/4\n";
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
            pnorm              = (pcurrent-p_iC) / (vlh_liquid.Pressure()-p_iC);
            value_iC           = storage_vector[iC];
            tdummy             = t_iB; // anyhow the case but make sure
            value_vlh_behind   = vlh_liquid.ValueOf(property_index);
            tdummy             = tvlh;
            value_vlh          = vlh_liquid.ValueOf(property_index);
	    
            value_behind       =  value_iC + pnorm*( value_vlh_behind - value_iC );
            tnorm              = (tcurrent-tvlh) / (t_iB-tvlh);
            value_interpolated =  value_vlh + tnorm*(value_behind - value_vlh);
          }
      }
    else
      {
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
            tdummy            = t_iB; // as before: is the case anyhow but to make sure
            value_vlh_behind      = vlh_liquid.ValueOf(property_index);
            tdummy            = tvlh;
            value_vlh             = vlh_liquid.ValueOf(property_index);
            value_behind          =  value_iC + pnorm*( value_vlh_behind - value_iC );

            tdummy       = t_iD;
            pnorm        = (pcurrent-p_iD)/(vlh_liquid.Pressure()-p_iD);
            value_iD         = storage_vector[iD];
            value_vlh_before = vlh_liquid.ValueOf(property_index);

            value_before       = value_iD + pnorm*(value_vlh_before-value_iD);
            tnorm          = (tcurrent-t_iD) / (t_iB-t_iD);
            value_interpolated = value_before + tnorm*(value_behind-value_before);
          }
        else
          {
	    
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
            value_vlh_before = vlh_liquid.ValueOf(property_index);
            value_before       = value_iD + pnorm*(value_vlh_before-value_iD);
            tnorm          = (tcurrent-t_iD) / (t_iB-t_iD);
            value_interpolated = value_before + tnorm*(value_behind-value_before);
          }
      }

    return value_interpolated;
  }




  // needs update
  double64 TwophaseLiquidLookup::InterpolateBetweenBoilingCurveAndVLH( const int& property_index )
  {
    //d cout << "Using TwophaseLiquidLookup::InterpolateBetweenBoilingCurveAndVLH( const int& property_index ) ..\n";

    // ************************** this looks damn incomplete - also check in TwophaseVaporLookup!!!! ***************************

    eqtype = 14;
    tdummy = tcurrent;
    // new, complete(?) version, also simplified:
    pnorm = (pcurrent - vlh_liquid.Pressure())/(water.LiquidProperty(tcurrent,pressure_index) - vlh_liquid.Pressure());
    return vlh_liquid.ValueOf(property_index) 
      + pnorm * (water.LiquidProperty(tcurrent,property_index) - vlh_liquid.ValueOf(property_index));

    //     if(     property_index == temperature_index)  return water.LiquidProperty(tcurrent,temperature_index) - vlh_liquid.Temperature();
    //     else if(property_index == pressure_index)     return water.LiquidProperty(tcurrent,pressure_index) - vlh_liquid.Pressure();
    //     else if(property_index == composition_index){
    // 	pnorm = (pcurrent - vlh_liquid.Pressure())/(water.LiquidProperty(tcurrent,pressure_index) - vlh_liquid.Pressure());
    // 	return vlh_liquid.MassFractionNaCl() + pnorm*(water.LiquidProperty(tcurrent,composition_index) - vlh_liquid.MassFractionNaCl());
    //     }
    //     else if(property_index == density_index)         return water.LiquidProperty(tcurrent,density_index) - vlh_liquid.Density();
    //     else if(property_index == enthalpy_index)        return water.LiquidProperty(tcurrent,enthalpy_index) - vlh_liquid.Enthalpy();
    //     else if(property_index == heatcapacity_index)    return water.LiquidProperty(tcurrent,heatcapacity_index) - vlh_liquid.HeatCapacity();
    //     else if(property_index == compressibility_index) return water.LiquidProperty(tcurrent,compressibility_index) - vlh_liquid.Compressibility();
    //     else if(property_index == viscosity_index)       return water.LiquidProperty(tcurrent,viscosity_index) - vlh_liquid.Viscosity();
    //     else return 9.9e99;
  }

  void TwophaseLiquidLookup::GetTemperatureIndex(const double64& t)
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
    
  void TwophaseLiquidLookup::GetPressureIndex(const double64& p)
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
    else{ ip = 506; }//                   p_res = 100.0; ip = 476+static_cast<long>( (p-2000.0)/p_res );/* throw out of range ? */ }
    // p_dim is therefore (5000-2000)/100+476 + 1 = 507
  }


} // namespace csmp


