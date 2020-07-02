#include <limits>

#include "TwophaseVaporLookup.h"

#include "TwophaseVapor.h"
#include "Brine.h"
#include "LookupPropertyIndex.h"

#include "binaryReadWrite.h"

#include <cmath>

using namespace std;

namespace csmp
{
  TwophaseVaporLookup::TwophaseVaporLookup(const double64& externaltemperature, const double64& externalpressure)
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
    it_p_max(0),
    ip_p_max(0),
    state(none),
    state_iA(none),
    state_iB(none),
    state_iC(none),
    state_iD(none),
    critcurve(tdummy), 
    vlh_vapor(tdummy),
    water(),// initialization with tdummy, pdummy is ugly - could that be cleaned up?
    vlh_pmax(vlh_vapor.Pmax()),
    vlh_tmax(vlh_vapor.Tmax()),
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
    //d cerr << "TwophaseVaporLookup t_dim and p_dim = " << t_dim << "\t" << p_dim << endl;
	
    storage_vector.resize(t_dim*p_dim*max_index);
    state_vector.resize(t_dim*p_dim);

    char filename[60], statefilename[60];
    strcpy( filename, "TwophaseVaporPropertiesLookupTable.bin" );
    strcpy( statefilename, "TwophaseVaporStateLookupTable.bin" );
	fstream infile1(filename, ios::in | ios::binary);
	fstream infile2(statefilename, ios::in | ios::binary);

	if (!infile1.is_open() || !infile2.is_open())
      {
        cerr << "TwophaseVaporLookup : at least one lookup file missing, computing ...\n\n";
        it = 0;
	    
        state    = V; // for brine only!!!!!
	    
        TwophaseVapor                       twophase_v(tcurrent,pcurrent,ph2o);
        Brine                               brine(tcurrent,pcurrent,xcurrent);
	    
        //d  cerr << "TwophaseVaporLookup::TwophaseVaporLookup(...): computing subcritical data ...\n";
        for(tcurrent = 0.0e0; tcurrent < cp_h2o.Temperature(); tcurrent += t_res)
          {
            ph2o    = water.SaturationPressureFromT(tcurrent);
            GetTemperatureIndex(tcurrent+1.0e-3);
            cerr << "computing for t = " << tcurrent << ", data set " << it << endl << endl;
            if(it > t_dim-1)
              {
                cerr << "too high t-index ...\n";
                break;
              }
		
            tdummy = tcurrent;
		
            for(pcurrent = 1.0e5; pcurrent <= 5000.0e5; pcurrent += p_res)
              {
                GetPressureIndex(pcurrent+1.0e-3);
                if(ip > p_dim-1)
                  {
                    cerr << "too high p-index ...\n";
                    break;
                  }
		    
                if( (pcurrent < vlh_vapor.Pressure() ) )
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
                    xcurrent = twophase_v.MassFractionNaCl();
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
        //d  cerr << "TwophaseVaporLookup::TwophaseVaporLookup(...): computing supercritical data ...\n";
	    
        // ***** check if the usage of 374.0 causes problems
        for(tcurrent = 374.0e0; tcurrent < 1000.1e0; tcurrent += t_res)
          {
            GetTemperatureIndex(tcurrent+1.0e-3);
            cerr << "computing for t = " << tcurrent << ", data set " << it << endl << endl;
            if(it > t_dim-1)
              {
                cerr << "too high t-index ...\n";
                break;
              }
            tdummy = tcurrent; // for cirtcurve_lookup and vlh_vapor_lookup
		
            for(pcurrent = 0.5e5; pcurrent <= 5000.0e5; pcurrent += p_res)
              {
                GetPressureIndex(pcurrent+1.0e-3);
                if(ip > p_dim-1)
                  {
                    cerr << "too high p-index ...\n";
                    break;
                  }
		    
                if( (tcurrent <=800.7) && (pcurrent < vlh_vapor.Pressure()) )
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
                    xcurrent = twophase_v.MassFractionNaCl();
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
	    
		fstream outfile1(filename, ios::out | ios::binary);
		if (!outfile1.is_open())
          {
            cerr << "could not even it for writing, please stop program and debug !!!!\n";
            char yesno;
            cerr << "or enter any key to continue (simulation likely to crash or give wrong results!) :";
            cin  >> yesno;
          }
        else
          {
            cerr << "writing file " << filename << " ... ";
            binaryFileWrite( outfile1, storage_vector );
			      outfile1.close();
            cerr << "done!\n";
          }

		fstream outfile2(statefilename, ios::out | ios::binary);
		if (!outfile2.is_open())
          {
            cerr << "could not even it for writing, please stop program and debug !!!!\n";
            char yesno;
            cerr << "or enter any key to continue (simulation likely to crash or give wrong results!) :";
            cin  >> yesno;
          }
        else
          {
            cerr << "writing file " << statefilename << " ... ";
            binaryFileWrite( outfile2, state_vector );
			     outfile2.close();
            cerr << "done!\n";
          }

      }
	
    else
      {
        cerr << "reading file " << filename << " ... ";
        binaryFileRead( infile1, storage_vector );
		    infile1.close();		
        cerr << "done!\n";

        cerr << "reading file " << statefilename << " ... ";
        binaryFileRead( infile2, state_vector );
		    infile2.close();
        cerr << "done!\n";
      }
    cerr << "TwophaseVaporLookup, leaving constructor ...\n\n";
  }



  TwophaseVaporLookup::~TwophaseVaporLookup()
  {
  }



  // The data interpolation routines
  double64 TwophaseVaporLookup::Temperature(){     SetTemperatureAndPressure(); return Interpolate(temperature_index);    }
  double64 TwophaseVaporLookup::Pressure(){        SetTemperatureAndPressure(); return Interpolate(pressure_index);       }
  double64 TwophaseVaporLookup::MassFractionNaCl(){     SetTemperatureAndPressure(); return Interpolate(composition_index);    }
  double64 TwophaseVaporLookup::Density(){         SetTemperatureAndPressure(); return Interpolate(density_index);        }
  double64 TwophaseVaporLookup::Enthalpy(){        SetTemperatureAndPressure(); return Interpolate(enthalpy_index);       }
  double64 TwophaseVaporLookup::HeatCapacity(){    SetTemperatureAndPressure(); return Interpolate(heatcapacity_index);   }
  double64 TwophaseVaporLookup::Compressibility(){ SetTemperatureAndPressure(); return Interpolate(compressibility_index);}
  double64 TwophaseVaporLookup::Viscosity(){       SetTemperatureAndPressure(); return Interpolate(viscosity_index);      }
    
  // these are only needed for the DT function, consider cleaning up such that "SetTemperatureAndPressure" might
  // be moved into "Interpolate"(consistency with conventions in other Lookups !)
  double64 TwophaseVaporLookup::ReportMassFractionNaCl(){ return Interpolate(composition_index); }
  double64 TwophaseVaporLookup::ReportEnthalpy(){    return Interpolate(enthalpy_index); }
    


  double64 TwophaseVaporLookup::DEnthalpyDT()
  {
    SetTemperatureAndPressure(); 
    tdummy   = tcurrent;
	
    if(pcurrent >= vlh_vapor.Pmax())
      {
        //d  cout << "TwophaseVaporLookup::DEnthalpyDT() - path 1\n";
        tcurrent += 0.1;
        x_high   = ReportEnthalpy();
        tcurrent = tdummy;
        return (x_high-ReportEnthalpy())/0.1;
      }
    else
      {
        if(tcurrent < vlh_vapor.Tmax())
          {
            tvlh = vlh_vapor.TfromP(pcurrent,300.0);
            if(tcurrent < tvlh-0.101)
              {
                //d  cout << "TwophaseVaporLookup::DEnthalpyDT() - path 2\n";
                tcurrent += 0.1;
                x_high   = ReportEnthalpy();
                tcurrent = tdummy;
                return (x_high-ReportEnthalpy())/0.1;
              }
            else
              {
                //d  cout << "TwophaseVaporLookup::DEnthalpyDT() - path 3\n";
                tcurrent -= 0.1;
                x_low    = ReportEnthalpy();
                tcurrent = tdummy;
                return (ReportEnthalpy()-x_low)/0.1;
              }
          }
        else
          {
            tvlh = vlh_vapor.TfromP(pcurrent,700.0);
            if(tcurrent > tvlh+0.101)
              {
                //d  cout << "TwophaseVaporLookup::DEnthalpyDT() - path 4\n";
                tcurrent -= 0.1;
                x_low   = ReportEnthalpy();
                tcurrent = tdummy;
                return (ReportEnthalpy()-x_low)/0.1;
              }
            else
              {
                //d  cout << "TwophaseVaporLookup::DEnthalpyDT() - path 5\n";
                tcurrent += 0.1;
                x_high    = ReportEnthalpy();
                tcurrent = tdummy;
                return (x_high-ReportEnthalpy())/0.1;
              }
          }
      } 
  }
    
    
  double64 TwophaseVaporLookup::DCompositionDT()
  {
    SetTemperatureAndPressure(); 
    tdummy   = tcurrent;
	
    if(pcurrent >= vlh_vapor.Pmax())
      {
        //d  cout << "TwophaseVaporLookup::DCompositionDT() - path 1\n";
        tcurrent += 0.1;
        x_high   = ReportMassFractionNaCl();
        tcurrent = tdummy;
        return (x_high-ReportMassFractionNaCl())/0.1;
      }
    else
      {
        if(tcurrent < vlh_vapor.Tmax())
          {
            tvlh = vlh_vapor.TfromP(pcurrent,300.0);
            if(tcurrent < tvlh-0.101)
              {
                //d  cout << "TwophaseVaporLookup::DCompositionDT() - path 2\n";
                tcurrent += 0.1;
                x_high   = ReportMassFractionNaCl();
                tcurrent = tdummy;
                return (x_high-ReportMassFractionNaCl())/0.1;
              }
            else
              {
                //d  cout << "TwophaseVaporLookup::DCompositionDT() - path 3\n";
                tcurrent -= 0.1;
                x_low    = ReportMassFractionNaCl();
                tcurrent = tdummy;
                return (ReportMassFractionNaCl()-x_low)/0.1;
              }
          }
        else
          {
            tvlh = vlh_vapor.TfromP(pcurrent,700.0);
            if(tcurrent > tvlh+0.101)
              {
                //d  cout << "TwophaseVaporLookup::DCompositionDT() - path 4\n";
                tcurrent -= 0.1;
                x_low   = ReportMassFractionNaCl();
                tcurrent = tdummy;
                return (ReportMassFractionNaCl()-x_low)/0.1;
              }
            else
              {
                //d  cout << "TwophaseVaporLookup::DCompositionDT() - path 5\n";
                tcurrent += 0.1;
                x_high    = ReportMassFractionNaCl();
                tcurrent = tdummy;
                return (x_high-ReportMassFractionNaCl())/0.1;
              }
          }
      } 
  }
    

  double64 TwophaseVaporLookup::DSaltMassFractionDT()
  {
    SetTemperatureAndPressure(); 
    tdummy   = tcurrent;
	
    if(pcurrent >= vlh_vapor.Pmax())
      {
        tcurrent += 0.1;
        x_high   = XNaCl2Massfraction(ReportMassFractionNaCl());
        tcurrent = tdummy;
        return (x_high-XNaCl2Massfraction(ReportMassFractionNaCl()))/0.1;
      }
    else
      {
        if(tcurrent < vlh_vapor.Tmax())
          {
            tvlh = vlh_vapor.TfromP(pcurrent,300.0);
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
            tvlh = vlh_vapor.TfromP(pcurrent,700.0);
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
    
    
  void TwophaseVaporLookup::SetTemperatureAndPressure()
  {
    tcurrent = temperature;
    pcurrent = pressure;
    tdummy   = tcurrent;
    pdummy   = pcurrent;
    return;
  }
    
    
  void TwophaseVaporLookup::GetIndex_iA(const int& property_index)
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
    
    
  double64 TwophaseVaporLookup::Interpolate(const int& property_index)
  {
    // do NOT set tcurrent and pcurrent here, do it outside this function too keep it versatile !!!
	
    GetIndex_iA(property_index);
    //d  cout << "\n\nTwophaseVaporLookup::Interpolate for property_index " << property_index << endl;
    //d  cout << "it = " << it << ", ip = " << ip << endl;
	
    //CheckForOutOfRange();
    state_iA = state_vector[ it*p_dim + ip ];
    state_iB = state_vector[ it*p_dim + ip + p_dim ];
    state_iC = state_vector[ it*p_dim + ip + p_dim +1 ];
    state_iD = state_vector[ it*p_dim + ip + 1 ];
	
    i_dummy  = t_dim*p_dim*temperature_index + it*p_dim + ip;
    //d  cout << "i_dummy t = " << i_dummy << endl;
    t_iA     = storage_vector[i_dummy];
    t_iB     = storage_vector[i_dummy + p_dim];
    t_iC     = storage_vector[i_dummy +1 + p_dim];
    t_iD     = storage_vector[i_dummy + 1];
	
    i_dummy  = t_dim*p_dim*pressure_index + it*p_dim + ip;
    //d  cout << "i_dummy p = " << i_dummy << endl;
    p_iA     = storage_vector[i_dummy];
    p_iB     = storage_vector[i_dummy + p_dim];
    p_iD     = storage_vector[i_dummy + 1];
    p_iC     = p_iD;
	
    //d  cout << "Temperatures and Pressures iA-iD:\n";
    //d  cout << t_iA << "\t" << t_iB << "\t" << t_iC << "\t" << t_iD << endl;
    //d  cout << p_iA << "\t" << p_iB << "\t" << p_iC << "\t" << p_iD << endl;
	
    if( it == it_p_max && (ip == ip_p_max || ip == ip_p_max-1/*ip == ip_p_max+1*/) )
      {
        //d  cout << "TwophaseVaporLookup goes 99\n";
        return NearVLHMaxInterpolation(property_index);
      }

    // *** CAUTION: this applies to current lookup spacing only!!!
    if( tcurrent >= 373.9 && tcurrent <= 374.0 && pcurrent >= 220.5e5 && pcurrent <= 220.6e5 )
      {
        return NearCritpointInterpolation(property_index);
      }

    if( (state_iA == VL) && (state_iB == VL) && (state_iC == VL) && (state_iD == VL) )
      {
        //d  cout << state_iA << "\t" << state_iB << "\t" << state_iC << "\t" << state_iD << endl;
        return NormalInterpolation(property_index); 
      }
    //   P
    //   ^            phase boundary
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
        //d  cout << state_iA << "\t" << state_iB << "\t" << state_iC << "\t" << state_iD << endl;
	    
        // ***** careful: how to handle case very near tcrit ??? *****
        if(tcurrent >= cp_h2o.Temperature())
          {
            if(      state_iD == F && state_iB == VL ) return NearCritcurveInterpolation(property_index);
            else if( state_iB == V && state_iD == VL && it < it_p_max ) return NearVLHInterpolationLowT(property_index);
            else if( state_iA == V && state_iC == VL && it > it_p_max ) return NearVLHInterpolationHighT(property_index);
            else
              {
                //d  cout << "TwophaseVaporLookup::Interpolate: missed to find interpolation type in case 2, returning bogus value ...\n";
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
                     && pcurrent >= vlh_vapor.Pressure() ) return InterpolateBetweenBoilingCurveAndVLH(property_index);
            else
              {
                //d  cout << "TwophaseVaporLookup::Interpolate: missed to find interpolation type in case 3, returning bogus value ...\n";
                return 9.9e99;
              }
          }
      }
  }
    
  double64 TwophaseVaporLookup::NormalInterpolation( const int& property_index )
  {
    //d  cout << "Using TwophaseVaporLookup::NormalInterpolation( const int& property_index ) ...\n";
	
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
    

  double64 TwophaseVaporLookup::NearVLHMaxInterpolation(const int& property_index)
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
    //   -D/A---/----\------C/B- ip_p_max+1 // isn't it rather ip_p_max ?
    //     | _/       V\     |
    //     |/           L\   |
    //   _/|              H\ |  
    //  /  |       VH          \ 
    //     |                 | \
    //   --A-----------------B-- ip_p_max // and this ip_p_max-1 ?
    //     |                 |
    //   it_p_max
    //     

    // iA etc. have already been determined when this function is called
    if( ip == ip_p_max/*+1*/ ) // in this case, ABDC in the above figure are moved up by one cell
      {
        if( pcurrent >= vlh_pmax ) 
          {
            //d  cout << "p is >= pmax \n";
            // this may seem overly complicated but is required to maintain consistency with
            // VLH_VaporLookup, which provides exact values tmax, pmax
            if( tcurrent <= vlh_tmax )
              {
                //d cout << "TwophaseVaporLookup::NearVLHMaxInterpolation - case1\n";
                // first, value_before
                pnorm              = (pcurrent - p_iA) / (p_iD - p_iA);
                value_before       =  storage_vector[iA] + pnorm*(storage_vector[iD]-storage_vector[iA]);
                // then value_behind = @vlh_tmax
                tnorm              = (vlh_tmax - t_iD) / (t_iC - t_iD);
                value_top          =  storage_vector[iD] + tnorm*(storage_vector[iC] - storage_vector[iD]);
                pnorm              = (pcurrent - vlh_pmax) / (p_iD - vlh_pmax);
                tdummy             =  vlh_tmax; // needed for vlh_vapor in following line
                value_behind       =  vlh_vapor.ValueOf(property_index) + pnorm*(value_top - vlh_vapor.ValueOf(property_index));
                // then interpolate
                tnorm              = (tcurrent - t_iD) / (vlh_tmax-t_iD);
                value_interpolated =  value_before + tnorm*(value_behind-value_before);
                return                value_interpolated;
              }
            else
              {
                //d cout << "TwophaseVaporLookup::NearVLHMaxInterpolation - case2\n";
                // first, value_before
                tnorm              = (vlh_tmax - t_iD) / (t_iC - t_iD);
                value_top          =  storage_vector[iD] + tnorm*(storage_vector[iC] - storage_vector[iD]);
                pnorm              = (pcurrent - vlh_pmax) / (p_iD - vlh_pmax);
                tdummy             =  vlh_tmax; // needed for vlh_vapor in following line
                value_before       =  vlh_vapor.ValueOf(property_index) + pnorm*(value_top - vlh_vapor.ValueOf(property_index));
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
                //d cout << "TwophaseVaporLookup::NearVLHMaxInterpolation - case3\n";
                // first, value_before
                pnorm              = (pcurrent - p_iA) / (p_iD - p_iA);
                value_before       =  storage_vector[iA] + pnorm*(storage_vector[iD]-storage_vector[iA]);
                // then value_behind @vlh curve at p
                tdummy             =  vlh_vapor.TfromP(pcurrent,300.); // needed for vlh_vapor in following line
                value_behind       =  vlh_vapor.ValueOf(property_index);
                // then interpolate
                tnorm              = (tcurrent - t_iA) / (tdummy - t_iA);
                value_interpolated =  value_before + tnorm*(value_behind-value_before);
                return                value_interpolated;		
              }
            else if( tcurrent > vlh_tmax )
              {
                //d cout << "TwophaseVaporLookup::NearVLHMaxInterpolation - case4\n";
                // first, value_before @vlh curve at p
                tdummy             =  vlh_vapor.TfromP(pcurrent,700.); // needed for vlh_vapor in following line
                value_before       =  vlh_vapor.ValueOf(property_index);
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
                                   "TwophaseVaporLookup::NearVLHMaxInterpolation(const int& property_index) - ",
                                   "p < vlh_pmax (in cell ip_p_max+1) but neither tcurrent<vlh_tmax nor tcurrent>vlh_tmax condition worked!\nReport issue to Thomas Driesner, thomas.driesner@erdw.ethz.ch"); 
                return std::numeric_limits<double>::signaling_NaN();
             }
          }
      }
    else if( ip == ip_p_max-1 )
      {
        if(tcurrent < vlh_tmax )
          {
            //d cout << "TwophaseVaporLookup::NearVLHMaxInterpolation - case5\n";
            // first, value_before
            tdummy             =  t_iA;	    
            pnorm              = (pcurrent - vlh_vapor.Pressure()) / (p_iD -  vlh_vapor.Pressure());
            value_before       =  vlh_vapor.ValueOf(property_index) + pnorm*(storage_vector[iD]-vlh_vapor.ValueOf(property_index));
            // then value_behind @vlh curve at p
            tdummy             =  vlh_vapor.TfromP(pcurrent,300.); // tdummy id the one that vlh_vapor uses, see constructor!
            value_behind       =  vlh_vapor.ValueOf(property_index);
            // then interpolate
            tnorm              = (tcurrent - t_iA) / (tdummy - t_iA);
            value_interpolated =  value_before + tnorm*(value_behind-value_before);
            return                value_interpolated;		
          }
        else if( tcurrent > vlh_tmax )
          {
            //d cout << "TwophaseVaporLookup::NearVLHMaxInterpolation - case6\n";
            // first, value_before @vlh curve at p
            tdummy             =  vlh_vapor.TfromP(pcurrent,700.); // tdummy id the one that vlh_vapor uses, see constructor!
            value_before       =  vlh_vapor.ValueOf(property_index);
            //d cout << "tdummy       = " << tdummy << endl;
            //d cout << "value_before = " << value_before << endl;
            // then value_behind 
            tdummy             =  t_iB;	    
            pnorm              = (pcurrent - vlh_vapor.Pressure()) / (p_iC - vlh_vapor.Pressure());
            value_behind       =  vlh_vapor.ValueOf(property_index) + pnorm*(storage_vector[iC]-vlh_vapor.ValueOf(property_index));
            //d cout << "tdummy       = " << tdummy << endl;
            //d cout << "pnorm        = " << pnorm << endl;
            //d cout << "value_behind = " << value_behind << endl;
            // then interpolate
            tnorm              = (tcurrent - tvlh) / (t_iB - tvlh);
            //d cout << "tnorm        = " << tnorm << endl;
            value_interpolated =  value_before + tnorm*(value_behind-value_before);
            return                value_interpolated;		
          }
        else
          {
            csmp_error.notice( FATAL_ERROR, 
                               "TwophaseVaporLookup::NearVLHMaxInterpolation(const int& property_index) - ",
                               "p < vlh_pmax (in cell ip_p_max) but neither tcurrent<vlh_tmax nor tcurrent>vlh_tmax condition worked!\nReport issue to Thomas Driesner, thomas.driesner@erdw.ethz.ch"); 
            return std::numeric_limits<double>::signaling_NaN();
         }
      }
    else
      {
        csmp_error.notice( FATAL_ERROR, 
                           "TwophaseVaporLookup::NearVLHMaxInterpolation(const int& property_index) - ",
                           "Missed ALL if-statements!\nReport issue to Thomas Driesner, thomas.driesner@erdw.ethz.ch"); 

        return std::numeric_limits<double>::signaling_NaN();
      }
  }
    

  double64 TwophaseVaporLookup::NearCritpointInterpolation( const int& property_index )
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
        value_boil         =  water.VaporProperty(tdummy,property_index);
        pnorm              = (pcurrent-p_iB) / (p_iC-p_iB);
        value_behind       =  storage_vector[iB] + pnorm*( storage_vector[iC] - storage_vector[iB] );
        tnorm              = (tcurrent-tdummy) / (t_iB-tdummy);
        value_interpolated =  value_boil+tnorm*(value_behind-value_boil); 
      }
    return value_interpolated;

  }



    
  double64 TwophaseVaporLookup::NearCritcurveInterpolation( const int& property_index )
  {
    //d  cout << "Using TwophaseVaporLookup::NearCritcurveInterpolation( const int& property_index ) ...\n";
    //d  cout << "Critical Pressure at T of interest is " << critcurve.Pressure() << endl;
	
	
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

            // I think this was forgotten before
            tdummy = t_iB;
            //

            pnorm             = (pcurrent-p_iB) / (pcrit_iB-p_iB);
            value_iB              = storage_vector[iB];
		
            value_crit_behind     = critcurve.ValueOf(property_index);
            tdummy            = tcrit;
            value_crit            = critcurve.ValueOf(property_index);

            value_behind          =  value_iB + pnorm*( value_crit_behind - value_iB );
            tnorm             = (tcurrent-tcrit) / (t_iB-tcrit);
            value_interpolated    =  value_crit + tnorm*(value_behind - value_crit);
          }
        else
          {
		
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
            pnorm          = (pcurrent-p_iB) / (p_iC-p_iB);
            value_behind       =  storage_vector[iB] + pnorm*( storage_vector[iC] - storage_vector[iB] );
            tnorm          = (tcurrent-t_iA) / (t_iB-t_iA);
		
            tdummy         = t_iA;
            pcrit_iA       = critcurve.Pressure();
            pnorm          = (pcurrent-p_iA) / (pcrit_iA-p_iA);
            value_iA           = storage_vector[iA];
		
            value_crit         = critcurve.ValueOf(property_index);
            value_before       =  value_iA + pnorm*( value_crit - value_iA );
            tnorm          = (tcurrent-t_iA) / (t_iB-t_iA);
            value_interpolated =  value_before + tnorm*(value_behind - value_before);
          }
      }
	
    //	cout << pcurrent << "\t" << "value_before = " << value_before << "\t" << "value_behind = " << value_behind << ", value_crit = " << value_crit << ", value_interpolated = " << value_interpolated << "tnorm = " << tnorm << endl;
	
    return value_interpolated;
  }
    
    
    
    
  double64 TwophaseVaporLookup::NearBoilingCurveInterpolation( const int& property_index )
  {
    //d  cout << "Using TwophaseVaporLookup::NearBoilingCurveInterpolation( const int& property_index ) ...\n";
    //d  cout << "Boiling Pressure at T of interest is " << water.SaturationPressureFromT(tcurrent) << endl;
	
    // a few data needed to determine topology
    tboil     = water.SaturationTemperatureFromP(pcurrent);
    tdummy    = t_iB;
    pboil_iB  = water.SaturationPressureFromT(tdummy);
    tdummy    = tcurrent;
	
    // There is a number of topologies possible (x is the value to be interpolated):
	
    if(tboil >= t_iA)
      {
	    
        if(pboil_iB <= p_iC)
          { 
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
            pnorm             = (pcurrent-p_iB) / (pboil_iB-p_iB);
            value_iB              = storage_vector[iB];
		
            value_boil_behind = water.VaporProperty(t_iB, property_index);
            value_boil        = water.VaporProperty(tboil,property_index);
            value_behind          =  value_iB + pnorm*( value_boil_behind - value_iB );
            tnorm             = (tcurrent-tboil) / (t_iB-tboil);
            value_interpolated    =  value_boil + tnorm*(value_behind - value_boil);
          }
        else
          {
		
            // case (3) and (4): Interpolation to be done between boiling_vapor_at_pcurrent and B-C
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
            pnorm          = (pcurrent-p_iB) / (p_iC-p_iB);
            value_behind       =  storage_vector[iB] + pnorm*( storage_vector[iC] - storage_vector[iB] );
            tnorm          = (tcurrent-tboil) / (t_iB-tboil);
            tdummy         = tboil;
            value_boil         = water.VaporProperty(tboil,property_index);
            value_interpolated =  value_boil + tnorm*(value_behind - value_boil);    
          }
      }
    else
      { // tboil_pcurrent is < T-iA
	    
        if(pboil_iB <= p_iC)
          {
            // case (5): Interpolation to be done between A-boiling_vapor_at_iA-segment and B-boiling_vapor_at_t_iB-segment
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
            pnorm  = (pcurrent-p_iB) / (pboil_iB-p_iB);
            value_iB   = storage_vector[iB];
            tdummy = t_iB;            
		
            value_boil_behind = water.VaporProperty(t_iB, property_index);
            value_behind          =  value_iB + pnorm*( value_boil_behind - value_iB );
		
		
            tdummy    = t_iA;
            pboil_iA  = water.VaporProperty(t_iA, pressure_index);
            pnorm     = (pcurrent-p_iA) / (pboil_iA-p_iA);
            value_iA      = storage_vector[iA];
		
            value_boil = water.VaporProperty(t_iA,property_index);
            value_before          =  value_iA + pnorm*( value_boil - value_iA );
            tnorm             = (tcurrent-t_iA) / (t_iB-t_iA);
            value_interpolated    =  value_before + tnorm*(value_behind - value_before);
          }
        else
          {
            // case (6): Interpolation to be done between A-boiling_vapor_at_iA-segment and B-C
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
            pnorm          = (pcurrent-p_iB) / (p_iC-p_iB);
            value_behind       =  storage_vector[iB] + pnorm*( storage_vector[iC] - storage_vector[iB] );
            tnorm          = (tcurrent-t_iA) / (t_iB-t_iA);
		
            tdummy    = t_iA;
            pboil_iA  = water.VaporProperty(t_iA,pressure_index);
            pnorm     = (pcurrent-p_iA) / (pboil_iA-p_iA);
            value_iA      = storage_vector[iA];
		
            value_boil = water.VaporProperty(t_iA,property_index);
            value_before           =  value_iA + pnorm*( value_boil - value_iA );
            tnorm             = (tcurrent-t_iA) / (t_iB-t_iA);
            value_interpolated    =  value_before + tnorm*(value_behind - value_before);
          }
      }
	
    //d     cout << "value_behind = " << value_behind << ", value_boil = " << value_boil << ", value_interpolated = " << value_interpolated << endl;
	
    return value_interpolated;
	
  }
    
    
    
  double64 TwophaseVaporLookup::NearVLHInterpolationLowT( const int& property_index )
  {
    //d  cout << "Using TwophaseVaporLookup::NearVLHInterpolationLowT( const int& property_index ) ...\n";
	
    tvlh      = vlh_vapor.TfromP(pcurrent,300.0e0);
    if(tvlh  <= t_iC)
      {
	    
        tdummy    = t_iD;
        if(vlh_vapor.Pressure() <= p_iA)
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
		
            value_vlh             = vlh_vapor.ValueOf(property_index);
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
            pnorm             = (pcurrent-p_iD) / (vlh_vapor.Pressure()-p_iD);
            value_iD              =  storage_vector[iD];
            value_vlh_before      =  vlh_vapor.ValueOf(property_index);
            value_before          =  value_iD + pnorm*( value_vlh_before - value_iD );
            tdummy            =  tvlh;
            value_vlh             =  vlh_vapor.ValueOf(property_index);
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
        if(vlh_vapor.Pressure() >= p_iA)
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
            pnorm             = (pcurrent-p_iD) / (vlh_vapor.Pressure()-p_iD);
            value_iC              = storage_vector[iC];
            value_iD              = storage_vector[iD];
            value_vlh_before      = vlh_vapor.ValueOf(property_index); // (between A and D)
            tdummy            = t_iB;
            value_vlh_behind      = vlh_vapor.ValueOf(property_index);

            value_before          =  value_iD + pnorm*( value_vlh_before - value_iD );
            tdummy            =  t_iC;
            pnorm             = (pcurrent-p_iC) / (vlh_vapor.Pressure()-p_iC);
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
            value_vlh_behind      = vlh_vapor.ValueOf(property_index);
            value_before          =  value_iD + pnorm*( value_iA - value_iD );
            tdummy            =  t_iC;
            pnorm             = (pcurrent-p_iC) / (vlh_vapor.Pressure()-p_iC);
            value_behind          =  value_iC + pnorm*( value_vlh_behind - value_iC); 
            tnorm             = (tcurrent-t_iD) / (t_iB-t_iD);
            value_interpolated    =  value_before + tnorm*(value_behind-value_before);
          }
      }
    return value_interpolated;
	
  }




  double64 TwophaseVaporLookup::NearVLHInterpolationHighT( const int& property_index )
  {
    //d  cout << "Using TwophaseVaporLookup::NearVLHInterpolationHighT( const int& property_index ) ...\n";
    tdummy = tcurrent;
	
    tvlh      = vlh_vapor.TfromP(pcurrent,700.0e0);
	
    if(tvlh  >= t_iA)
      {
	    
        tdummy    = t_iB;
        if(vlh_vapor.Pressure() <= p_iB)
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
            pnorm          = (pcurrent-p_iB) / (p_iC-p_iB);
            value_behind       =  storage_vector[iB] + pnorm*( storage_vector[iC] - storage_vector[iB] );
            tnorm          = (tcurrent-tvlh) / (t_iB-tvlh);
            tdummy         =  tvlh;
            value_vlh          = vlh_vapor.ValueOf(property_index);
            value_interpolated =  value_vlh + tnorm*(value_behind - value_vlh);   
          }
        else
          {
            // vlh_vapor.Pressure() > p_iB
		
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
            pnorm             = (pcurrent-p_iC) / (vlh_vapor.Pressure()-p_iC);
            value_iC              = storage_vector[iC];
            // tdummy    = t_iB !
            value_vlh_behind      = vlh_vapor.ValueOf(property_index);
            tdummy            = tvlh;
            value_vlh             = vlh_vapor.ValueOf(property_index);
            value_behind          =  value_iC + pnorm*( value_vlh_behind - value_iC );
            tnorm             = (tcurrent-tvlh) / (t_iB-tvlh);
            value_interpolated    =  value_vlh + tnorm*(value_behind - value_vlh);
          }
      }
    else
      {
        // tvlh < t_iA
        tdummy = t_iB;
        if(vlh_vapor.Pressure() >= p_iB)
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
		
		
            pnorm             = (pcurrent-p_iC) / (vlh_vapor.Pressure()-p_iC);
            value_iC              = storage_vector[iC];
            // tdummy    = t_iB !
            value_vlh_behind      = vlh_vapor.ValueOf(property_index);
            tdummy            = tvlh;
            value_vlh             = vlh_vapor.ValueOf(property_index);
            value_behind          =  value_iC + pnorm*( value_vlh_behind - value_iC );
		
            tdummy       = t_iD;
            pnorm        = (pcurrent-p_iD)/(vlh_vapor.Pressure()-p_iD);
            value_iD         = storage_vector[iD];

            value_vlh_before = vlh_vapor.ValueOf(property_index);
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
            pnorm        = (pcurrent-p_iD)/(vlh_vapor.Pressure()-p_iD);
            value_iD         = storage_vector[iD];
            value_vlh_before = vlh_vapor.ValueOf(property_index);
            value_before       = value_iD + pnorm*(value_vlh_before-value_iD);
            tnorm          = (tcurrent-t_iD) / (t_iB-t_iD);
            value_interpolated = value_before + tnorm*(value_behind-value_before);
          }
      }
	
    return value_interpolated;
  }
    
  double64 TwophaseVaporLookup::InterpolateBetweenBoilingCurveAndVLH( const int& property_index )
  {
    //d  cout << "Using TwophaseVaporLookup::InterpolateBetweenBoilingCurveAndVLH( const int& property_index ) ..\n";
    tdummy = tcurrent;

    pnorm = (pcurrent - vlh_vapor.Pressure())/(water.VaporProperty(tcurrent,pressure_index) - vlh_vapor.Pressure());
    return vlh_vapor.ValueOf(property_index) 
      + pnorm * (water.VaporProperty(tcurrent,property_index) - vlh_vapor.ValueOf(property_index));
  }
    
    
    
  void TwophaseVaporLookup::GetTemperatureIndex(const double64& t)
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
    
  void TwophaseVaporLookup::GetPressureIndex(const double64& p)
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
