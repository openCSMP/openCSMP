#include "VLH_VaporLookup.h"

#include "Brine.h"
#include "ThreephaseHLV.h"
#include "NaClSaturatedVapor.h"
#include "H2OLookup.h"
#include "LookupPropertyIndex.h"
#include "compareFloats.h"
#include "binaryReadWrite.h"

#include <cmath>


using namespace std;

namespace csmp
{
  extern csmp::ErrorHandler           csmp_error;

  VLH_VaporLookup::VLH_VaporLookup(const double& externaltemperature)
    : temperature(externaltemperature),
      tcurrent(0.0),
      pcurrent(0.0),
      xcurrent(0.0),
      tmax(0.0),
      pmax(0.0),
      t_res(0.0),
      tnorm(0.0),
      pnorm(0.0),
      dp(0.0),
      it(0),
      t_dim(0),
      i_guess(0),
      i_max(0),
      i_min(0),
      it_p_max(0),
      tdim_times_pindex(0),
      state(none),
      csmp_error( ErrorHandler::Instance() )
  {
    // ---------------------------------------------------------------
    // 0. Set up external objects needed to compute lookup tables etc.
    ThreephaseHLV                vlh(tcurrent);
    NaClSaturatedVapor           naclsatvap(tcurrent,pcurrent);
    Brine                        brine(tcurrent,pcurrent,xcurrent);
    // ---------------------------------------------------------------
    
    // ---------------------------------------------------------------
    // 1. Preparing vector for properties at Tmax, Pmax;
    tmax                       = vlh.Tmax();
    pmax                       = vlh.Pmax();
    GetTemperatureIndex(tmax);
    it_p_max                   = it;
    tcurrent                   = tmax; // needs to be intialized before calling naclsatvap!
    pcurrent                   = pmax; // dito;
    xcurrent                   = naclsatvap.MassFractionNaCl();
    // what follows is supposed to be a solid check for NaN (http://bytes.com/topic/c/answers/588254-how-check-double-inf-nan)
    // reason to include it is that I recently experienced NaN for naclsatvap.MassFractionNaCl() at exactly these
    // conditions
    bool my_isnan;
    my_isnan = ( xcurrent != xcurrent); 
    if( my_isnan )
      {
        cout << "xcurrent = " << xcurrent << endl;
        csmp_error.notice( FATAL_ERROR, 
                           "Constructor VLH_VaporLookup::VLH_VaporLookup(const double& externaltemperature) -",
                           "FATAL_ERROR: when constructing the properties_at _tmax vector, x at Tmax, Pmax was calaculated as NaN!\ncontact developer.\nResponsible developer: Thomas Driesner, thomas.driesner@erdw.ethz.ch");
        return;
      }
    if( max_index-1 != 7 )
      {
        cerr << "max_index is " << max_index << endl;
        csmp_error.notice( FATAL_ERROR, 
                           "Constructor VLH_VaporLookup::VLH_VaporLookup(const double& externaltemperature) -",
                           "FATAL_ERROR: wrong size of properties_at_tmax vector compared to LookupPropertyIndex.h!\ncontact developer.\nResponsible developer: Thomas Driesner, thomas.driesner@erdw.ethz.ch");
        return;
      }	
    properties_at_tmax.resize(   max_index);
    properties_at_tmax[0]      = tmax;
    properties_at_tmax[1]      = pmax;
    properties_at_tmax[2]      = xcurrent;
    properties_at_tmax[3]      = brine.Density();
    properties_at_tmax[4]      = brine.Enthalpy();
    properties_at_tmax[5]      = brine.HeatCapacity();
    properties_at_tmax[6]      = brine.Compressibility();
    properties_at_tmax[7]      = brine.Viscosity();
    // ---------------------------------------------------------------

    // ---------------------------------------------------------------
    // 2. Check for already exisiting lookup table:
    //    - read in if existing
    //    - compute, write and establish if not yet existing
    GetTemperatureIndex(         1000.0001);
    t_dim                      = it+1;
    storage_vector.resize(       t_dim*max_index);
    tdim_times_pindex          = t_dim*pressure_index;

    char filename[60];
    strcpy(filename,"VLH_VaporLookupTable.bin");
	fstream infile(filename, ios::in | ios::binary);

	if (!infile.is_open())
      {
        cout << "VLH_VaporLookup : Lookup file missing, computing ...\n\n";
        it                     = 0;
        state                  = V;

        H2OLookup              water;

        for(tcurrent = 0.0e0; tcurrent < 800.7e0; tcurrent += t_res)
          {
            if(essentiallyEqual( tcurrent, tmax, numeric_limits<double>::epsilon() ) )
              {
                csmp_error.notice( FATAL_ERROR, 
                                   "Constructor VLH_VaporLookup::VLH_VaporLookup(const double& externaltemperature) -",
                                   "while building lookup table, tcurrent was == tmax of vlh curve.\nThis can potentially mess up computations of fluid\nproperties during simulations.\nProbably you or somebody else changed the t-resolution of lookup tables - re-think those.\nElse: if you have source code access you might set the error level associated\nwith this message to WARNING and pray ;-) but better contact developer.\nResponsible developer: Thomas Driesner, thomas.driesner@erdw.ethz.ch");
                return;
              }

            GetTemperatureIndex(tcurrent+1.0e-3);
            cout << "VLH_VaporLookup computing for t = " << tcurrent << ", data set " << it << endl;
            if(it > t_dim-1)
              {
                cout << "too high index ...\n";
                break;
              }
            pcurrent  = vlh.Pressure();
            xcurrent  = naclsatvap.MassFractionNaCl();
            storage_vector[t_dim*temperature_index+it]     = tcurrent;
            storage_vector[tdim_times_pindex+it]           = pcurrent;
            storage_vector[t_dim*composition_index+it]     = xcurrent;
            if(xcurrent < 1.0e-6)
              {
                // this section is a work-around since I encountered problems with SoWat brine
                // that weren't there before ...
                storage_vector[t_dim*density_index+it]         = water.Density(tcurrent,pcurrent);
                storage_vector[t_dim*enthalpy_index+it]        = water.Enthalpy(tcurrent,pcurrent);
                storage_vector[t_dim*heatcapacity_index+it]    = water.HeatCapacity(tcurrent,pcurrent);
                storage_vector[t_dim*compressibility_index+it] = water.Compressibility(tcurrent,pcurrent);
                storage_vector[t_dim*viscosity_index+it]       = water.Viscosity(tcurrent,pcurrent);
                if(pcurrent < 1.1)
                  {
                    storage_vector[t_dim*density_index+it]         = water.VaporDensityFromT(tcurrent);
                    storage_vector[t_dim*enthalpy_index+it]        = water.VaporEnthalpyFromT(tcurrent);
                    storage_vector[t_dim*heatcapacity_index+it]    = water.VaporHeatCapacityFromT(tcurrent);
                    storage_vector[t_dim*compressibility_index+it] = water.VaporCompressibilityFromT(tcurrent);
                    storage_vector[t_dim*viscosity_index+it]       = water.VaporViscosityFromT(tcurrent);
                  }
              }
            else
              {
                storage_vector[t_dim*density_index+it]         = brine.Density();
                storage_vector[t_dim*enthalpy_index+it]        = brine.Enthalpy();
                storage_vector[t_dim*heatcapacity_index+it]    = brine.HeatCapacity();
                storage_vector[t_dim*compressibility_index+it] = brine.Compressibility();
                storage_vector[t_dim*viscosity_index+it]       = brine.Viscosity();
              }
          }
	    
        tcurrent = tp_nacl.Temperature();
        pcurrent = tp_nacl.Pressure();
        xcurrent = tp_nacl.CompositionVapor();
        GetTemperatureIndex(tcurrent+1.0e-3);
        cout << "computing for t, p  = " << tcurrent << ", " << pcurrent << endl;
        storage_vector[t_dim*temperature_index+it]     = tcurrent;
        storage_vector[tdim_times_pindex+it]           = pcurrent;
        storage_vector[t_dim*composition_index+it]     = xcurrent;
        storage_vector[t_dim*density_index+it]         = brine.Density();
        storage_vector[t_dim*enthalpy_index+it]        = brine.Enthalpy();
        storage_vector[t_dim*heatcapacity_index+it]    = brine.HeatCapacity();
        storage_vector[t_dim*compressibility_index+it] = brine.Compressibility();
        storage_vector[t_dim*viscosity_index+it]       = brine.Viscosity();
	    
        // write file
		fstream outfile(filename, ios::out | ios::binary);
		if (!outfile.is_open())
          {
            cout << "could not even open it for writing, please stop program and debug !!!!\n";
            char yesno;
            cout << "or enter any key to continue (simulation likely to crash or give wrong results!) :";
            cin  >> yesno;
          }
        else
          {
            cout << "writing file " << filename << " ... ";
            binaryFileWrite( outfile, storage_vector );
			      outfile.close();
            cout << "done!\n";
          }
      }
    
    else
      {
        cout << "reading file " << filename << " ... ";
        binaryFileRead( infile, storage_vector );
		    infile.close();
        cout << "done!\n";
	
        // finding value for Tmax, Pmax
        ThreephaseHLV  vlh(tcurrent);
        GetTemperatureIndex(vlh.Tmax());
        it_p_max = it;
        if(storage_vector[tdim_times_pindex+it_p_max+1] > storage_vector[tdim_times_pindex+it_p_max]) it_p_max += 1;
        if(storage_vector[tdim_times_pindex+it_p_max-1] > storage_vector[tdim_times_pindex+it_p_max]) it_p_max -= 1;
        if(storage_vector[tdim_times_pindex+it_p_max-2] > storage_vector[tdim_times_pindex+it_p_max]) it_p_max -= 2;
        cout << "it_p_max found as " << it_p_max << ", with Tmax = " << storage_vector[t_dim*temperature_index+it_p_max] << endl;	
      }
    // ---------------------------------------------------------------

    cout << "VLH_VaporLookup, leaving constructor ...\n\n";
  }
  

  VLH_VaporLookup::~VLH_VaporLookup()
  {
  }
  

  double VLH_VaporLookup::DPressureDT()
  {
    // May be in error near Tmax, Pmax!!!
    GetTemperatureIndex(tcurrent);
    if(it == t_dim-1)
      {
        dp  = (storage_vector[tdim_times_pindex+it]-storage_vector[tdim_times_pindex+it-1]);
        dp /= (storage_vector[t_dim*temperature_index+it]-storage_vector[t_dim*temperature_index+it-1]);
      }
    else if(it == it_p_max)
      {
        if(tcurrent < tmax)
          {
            dp  =  pmax - storage_vector[tdim_times_pindex+it];
            dp /= (tmax - storage_vector[t_dim*temperature_index+it]);
          }
        else if(tcurrent > tmax)
          {
            dp  =  storage_vector[tdim_times_pindex+it+1]    - pmax;
            dp /= (storage_vector[t_dim*temperature_index+it+1] - tmax);
          }
        else
          {
            dp = 1.0e-10; // arbitrary - in H2ONaClFluidproperties (threephase compressibility) there is a division
            // by  this value, hence, zero would be fatal
          }
      }
    else
      {
        dp  = (storage_vector[tdim_times_pindex+it+1]-storage_vector[tdim_times_pindex+it]);
        dp /= (storage_vector[t_dim*temperature_index+it+1]-storage_vector[t_dim*temperature_index+it]);
      }
    return dp;
  }
  
  double VLH_VaporLookup::DEnthalpyDT()
  {
    // May be in error near Tmax, Pmax!!!
    GetTemperatureIndex(tcurrent);
    if(it == t_dim-1)
      {
        dp  = (storage_vector[t_dim*enthalpy_index+it]-storage_vector[t_dim*enthalpy_index+it-1]);
        dp /= (storage_vector[t_dim*temperature_index+it]-storage_vector[t_dim*temperature_index+it-1]);
      }
    else if(it == it_p_max)
      {
        if(tcurrent < tmax)
          {
            dp  =  properties_at_tmax[4] - storage_vector[t_dim*enthalpy_index+it];
            dp /= (tmax - storage_vector[t_dim*temperature_index+it]);
          }
        else if(tcurrent > tmax)
          {
            dp  =  storage_vector[t_dim*enthalpy_index+it+1]    - properties_at_tmax[4];
            dp /= (storage_vector[t_dim*temperature_index+it+1] - tmax);
          }
        else
          {
            dp = 1.0e-10; // arbitrary 
          }
      }
    else
      {
        dp  = (storage_vector[t_dim*enthalpy_index+it+1]-storage_vector[t_dim*enthalpy_index+it]);
        dp /= (storage_vector[t_dim*temperature_index+it+1]-storage_vector[t_dim*temperature_index+it]);
      }
    return dp;
  }
  
  double VLH_VaporLookup::TfromP(const double& press, const double& t_estimate)
  {
    // to catch all possibilities, first check for press == pmax
    // this should have been caught by code using VLH_VaporLookup but since that is not fool-proof, here's another
    // check; this, however, doens't prevent the user from running into problems if he/she wanted to avoid
    // ambiguities regarding on which side of the pressure maximum this pressure was hit 
    // if( definitelyGreaterThan( press, pmax, numeric_limits<double>::epsilon() ) )
    //   {
    // 	cerr << "pmax = " << pmax << ", pcurrent was " << press << endl;
    // 	csmp_error.notice( FATAL_ERROR, 
    // 			"VLH_VaporLookup::TfromP(const double& press, const double& t_estimate) -",
    // 			"FATAL_ERROR: you tried to invoke this function at p>pmax, this makes no sense, terminating!\nThomas Driesner, thomas.driesner@erdw.ethz.ch");
    // 	return 9.9e99;
    //   }

    if( essentiallyEqual( press, pmax, numeric_limits<double>::epsilon() ) )
      return tmax;

    else if( definitelyGreaterThan( t_estimate, tmax, numeric_limits<double>::epsilon() )
             &&
             !definitelyLessThan( press, storage_vector[it_p_max+1+tdim_times_pindex], numeric_limits<double>::epsilon() ) )
      {
        // we are in the lookup cell that contains pmax at t > tmax. If so, compute directly and return
        //
        //       |                 |              ( /--  represents VLH curve in T-P space )
        //       |        o (Pmax) |              ( *    is tcurrent on lookup curve       )
        //       |       / \       |
        //       |   /--/   \+\    |
        //      /|--/          \-\ |
        //   --/ |                \|
        //    it_p_max         it_p_max+1
        //
	
        //	cout << "VLH_VaporLookup::TfromP(const double& press, const double& t_estimate) uses direct computation at t > tmax!\n";
        pnorm  = press - pmax;
        pnorm /= storage_vector[it_p_max+1+tdim_times_pindex] - pmax;
        return  tmax + pnorm * ( storage_vector[it_p_max+1]-tmax );
      }

    else if( definitelyLessThan( t_estimate, tmax, numeric_limits<double>::epsilon() )
             &&
             !definitelyLessThan( press, storage_vector[it_p_max+tdim_times_pindex], numeric_limits<double>::epsilon() ) )
      // we are in the lookup cell that contains pmax at t<tmax. If so, compute directly and return
      //
      //       |                 |              ( /--  represents VLH curve in T-P space )
      //       |        - (Pmax) |              ( +    is tcurrent on lookup curve       )
      //       |       / \       |
      //       |   /+-/   \-\    |
      //      /|--/          \-\ |
      //   --/ |                \|
      //    it_p_max         it_p_max+1
      //
      {
        //	cout << "VLH_VaporLookup::TfromP(const double& press, const double& t_estimate) uses direct computation at t < tmax!\n";
        pnorm  = press - storage_vector[it_p_max+tdim_times_pindex];
        pnorm /= pmax  - storage_vector[it_p_max+tdim_times_pindex];
        return storage_vector[it_p_max] + pnorm * ( tmax-storage_vector[it_p_max] );
      }

    else
      {
        // use clear t_estimates, e.g. 200 and 700
        // here clarify which of two t's should be found for given p
        // since t(pmax) in the table near 590C with it=it_p_max, that's the boundary    
        GetTemperatureIndex(t_estimate+numeric_limits<double>::epsilon());

        if(it < it_p_max)
          {
            i_max = it_p_max; i_min = 0; i_guess = it_p_max/2;
            int        myi = 0;
            cout.setf(ios::scientific);
            while( i_max-i_guess > 1)
              {	
                if(press == storage_vector[i_guess+tdim_times_pindex])
                  {
                    return storage_vector[i_guess]; // if p is exactly correct at i_guess
                  }
                else if(press < storage_vector[i_guess+tdim_times_pindex])
                  { 
                    i_max   = i_guess; 
                    i_guess = (i_guess+i_min)/2; 
                  }
                else
                  { 
                    i_min = i_guess; 
                    i_guess = (i_max+i_guess)/2; 
                  }
                myi ++;
                //cout << myi << "\t" << storage_vector[i_min] << "\t" << storage_vector[i_max] << endl;
                if(myi > 50)
                  {
                    cerr << "VLH_VaporLookup::TfromP(const double& press) not converged ..." << endl;
                    return 0.0e0;
                  }
              }
            if(press < storage_vector[i_guess+tdim_times_pindex]) i_guess -= 1;
            pnorm = ( press-storage_vector[i_guess+tdim_times_pindex] ) / (storage_vector[i_guess+1+tdim_times_pindex]-storage_vector[i_guess+tdim_times_pindex]);
            return storage_vector[i_guess] + pnorm*(storage_vector[i_guess+1]-storage_vector[i_guess]);
          }
      	else
          {
            i_max    = t_dim-1;  i_guess  = (it_p_max+i_max)/2;  i_min    = it_p_max+1;
            int        myi = 0;
            cout.setf(ios::scientific);
            while( i_max-i_guess > 1)
              {	
                if(press == storage_vector[i_guess+tdim_times_pindex])
                  {
                    return storage_vector[i_guess]; // if p is exactly correct at i_guess
                  }
                else if(press > storage_vector[i_guess+tdim_times_pindex])
                  { 
                    i_max   = i_guess; 
                    i_guess = (i_guess+i_min)/2; 
                  }
                else
                  { 
                    i_min = i_guess; 
                    i_guess = (i_max+i_guess)/2; 
                  }
                myi ++;
                if(myi > 50)
                  {
                    cerr << "VLH_VaporLookup::TfromP(const double& press) not converged ..." << endl;
                    return 0.0e0;
                  }
              }
            if(press > storage_vector[i_guess+tdim_times_pindex]) i_guess -= 1;
            pnorm = ( press-storage_vector[i_guess+tdim_times_pindex] ) / (storage_vector[i_guess+1+tdim_times_pindex]-storage_vector[i_guess+tdim_times_pindex]);
            return storage_vector[i_guess] + pnorm*(storage_vector[i_guess+1]-storage_vector[i_guess]);
          }
      }
  } 

  // The data interpolation routines
  double VLH_VaporLookup::Temperature(){     return ValueOf(temperature_index); }
  double VLH_VaporLookup::Pressure(){        return ValueOf(pressure_index); }
  double VLH_VaporLookup::MassFractionNaCl(){     return ValueOf(composition_index); }
  double VLH_VaporLookup::Density(){         return ValueOf(density_index); }
  double VLH_VaporLookup::Enthalpy(){        return ValueOf(enthalpy_index); }
  double VLH_VaporLookup::HeatCapacity(){    return ValueOf(heatcapacity_index); }
  double VLH_VaporLookup::Compressibility(){ return ValueOf(compressibility_index); }
  double VLH_VaporLookup::Viscosity(){       return ValueOf(viscosity_index); }
    
    
  double VLH_VaporLookup::ValueOf(const int& property_index)
  {
    tcurrent = temperature;
    GetTemperatureIndex(tcurrent);
    if( it == it_p_max ) // we are in the lookup cell that contains the pressure-maximum
      {
        if( definitelyLessThan( tcurrent, tmax, numeric_limits<double>::epsilon() ) )
          // we are at a temperature LOWER than that of the maximum, i.e., interpolate between the
          // value at the LOWER temperature-index and that at the maximum
          {
            //	    cout << "VLH_VaporLookup::ValueOf(const int& property_index) at t < tmax\n";
            tnorm    = ( tcurrent - storage_vector[it] )/( tmax - storage_vector[it] );
            it      += t_dim*property_index;
            return storage_vector[it]+tnorm*( properties_at_tmax[property_index]-storage_vector[it] );
          }
        else if( definitelyGreaterThan( tcurrent, tmax, numeric_limits<double>::epsilon() ) )
          {
            //	    cout << "VLH_VaporLookup::ValueOf(const int& property_index) at t > tmax\n";
            // we are at a temperature HIGHER than that of the maximum, i.e., interpolate between the
            // value at the next HIGHER temperature-index and that at the maximum
            tnorm    = ( tcurrent - tmax ) / ( storage_vector[it+1] - tmax );
            it      += t_dim*property_index;
            return properties_at_tmax[property_index] + tnorm*( storage_vector[it+1] - properties_at_tmax[property_index] );
          }
        else
          // if at exactly the maximum coordinates, return pre-calculated values at that point
          {
            //	    cout << "VLH_VaporLookup::ValueOf(const int& property_index) at t == tmax\n";
            return properties_at_tmax[property_index];
          }
      }
    else
      // no danger of being near the maximum, proceed normally
      {
        //	cout << "VLH_VaporLookup::ValueOf(const int& property_index) at normal condition\n";
        tnorm    = (tcurrent - storage_vector[it])/(storage_vector[it+1]-storage_vector[it]);
        it      += t_dim*property_index; // convert to enthalpy section of storage vector
        return storage_vector[it]+tnorm*(storage_vector[it+1]-storage_vector[it]);
      }
  }
    
    
  // data indexing
  void VLH_VaporLookup::GetTemperatureIndex(const double& t)
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
    else{ it = 330; }//                   t_res = 10.0; it = 285+static_cast<long>( (t-550.0)/t_res ); /* throw out of range ? */ }
  }
  
  
} // namespace csmp
