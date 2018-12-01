#include <cmath>

#include "NaClMeltingCurveHaliteLookup.h"
#include "binaryReadWrite.h"
#include "NaClMeltingCurve.h"
#include "Halite.h"
#include "LookupPropertyIndex.h"
#include "States.h"

using namespace std;

namespace csmp
{
  NaClMeltingCurveHaliteLookup::NaClMeltingCurveHaliteLookup(const double64& externaltemperature,
                                                             const double64& externalpressure)
    : temperature(externaltemperature),
      pressure(externalpressure),
      tcurrent(0.0),
      pcurrent(0.0),
      xcurrent(0.0),
      t_res(0.0),
      tnorm(0.0),
      pnorm(0.0),
      tmin(0.0),
      it(0),
      t_dim(0),
      i_guess(0),
      i_max(0),
      i_min(0),
      it_max(0),
      it_min(0),
      state(0)
  {
    // this sequence MUST be done as the very first thing!
    // also ***** NEVER ***** change it_min's value anywhere if you do changes to this file
    it_min   = 0;
    GetTemperatureIndex(tp_nacl.Temperature());
    it_min   = it;
	
    GetTemperatureIndex(1000.0001);
    it_max   = it;
    t_dim    = it_max+1;

    storage_vector.resize(t_dim*max_index);
	
    char filename[60];
    strcpy(filename,"NaClMeltingCurveHaliteLookupTable.bin");
    
	fstream infile(filename, ios::in | ios::binary);
	if (!infile.is_open())
      {
        cerr << "NaClMeltingCurveHaliteLookup: file " << filename;
        cerr << " could not be openend, computing ...\n";
        NaClMeltingCurve      naclmelt(tcurrent, pcurrent);
        Halite                halite(tcurrent,pcurrent);                            
        tcurrent            = tp_nacl.Temperature();
        pcurrent            = tp_nacl.Pressure();
        xcurrent            = tp_nacl.CompositionHalite();
        state               = H;
        cerr << "computing for t = " << tcurrent << endl;
        it = 0;
        storage_vector[it]                             = tcurrent;
        storage_vector[it+t_dim*pressure_index]        = pcurrent;
        storage_vector[it+t_dim*composition_index]     = xcurrent;
        storage_vector[it+t_dim*density_index]         = halite.Density();
        storage_vector[it+t_dim*enthalpy_index]        = halite.Enthalpy();
        storage_vector[it+t_dim*heatcapacity_index   ] = halite.HeatCapacity();
        storage_vector[it+t_dim*compressibility_index] = halite.Compressibility();
        storage_vector[it+t_dim*viscosity_index      ] = 1.0e20;

        // Get the minimum temperature for looping, ASSUMING there is not another point 
        // between 800.0 C and and tp_nacl.Temperature()
        // If you change t_res to much smaller values, you must adjust this
        GetTemperatureIndex(tp_nacl.Temperature());
        tmin = 800.0e0+t_res; 

        for(tcurrent = tmin; tcurrent < 1000.0e0; tcurrent += t_res)
          {
            GetTemperatureIndex(tcurrent+1.0e-3);
            cerr << "NaClMeltingCurveHaliteLookup computing table entry for t = " << tcurrent << ", data set " << it << endl;
            if(it > it_max){
              cerr << "too high index ...\n";
              break;
            }
            pcurrent             = naclmelt.PmeltFromT();
            xcurrent             = 1.0e0;
            storage_vector[t_dim*temperature_index+it]     = tcurrent;
            storage_vector[t_dim*pressure_index+it]        = pcurrent;
            storage_vector[t_dim*composition_index+it]     = xcurrent;
            if(pcurrent > 5000.0e5)
              {
                //d cerr << "NaClMeltingCurveHaliteLookup constructor, highP entry for t = " << tcurrent << " : survived tpx\n";
                storage_vector[t_dim*density_index         + it] = 9.e99;
                //d cerr << "NaClMeltingCurveHaliteLookup constructor, highP entry for t = " << tcurrent << " : survived rho\n";
                storage_vector[t_dim*enthalpy_index        + it] = 9.e99;
                //d  cerr << "NaClMeltingCurveHaliteLookup constructor, highP entry for t = " << tcurrent << " : survived h\n";
                storage_vector[t_dim*heatcapacity_index    + it] = 9.e99;
                //d  cerr << "NaClMeltingCurveHaliteLookup constructor, highP entry for t = " << tcurrent << " : survived cp\n";
                storage_vector[t_dim*compressibility_index + it] = 9.e99;
                //d  cerr << "NaClMeltingCurveHaliteLookup constructor, highP entry for t = " << tcurrent << " : survived beta\n";
                storage_vector[t_dim*viscosity_index       + it] = 9.e99;
                //d  cerr << "NaClMeltingCurveHaliteLookup constructor, highP entry for t = " << tcurrent << " : survived all\n";
              }
            else
              {
                //d cerr << "NaClMeltingCurveHaliteLookup constructor, entry for t = " << tcurrent << " : survived tpx\n";
                storage_vector[t_dim*density_index         + it] = halite.Density();
                //d cerr << "Tstar_V and Tstar_H were " << halite.Tstar_V() << "\t" << halite.Tstar_H() << endl;
                //d cerr << "NaClMeltingCurveHaliteLookup constructor, entry for t = " << tcurrent << " : survived rho\n";
                storage_vector[t_dim*enthalpy_index        + it] = halite.Enthalpy();
                //d  cerr << "Tstar_V and Tstar_H were " << halite.Tstar_V() << "\t" << halite.Tstar_H() << endl;
                //d  cerr << "NaClMeltingCurveHaliteLookup constructor, entry for t = " << tcurrent << " : survived h\n";
                storage_vector[t_dim*heatcapacity_index    + it] = halite.HeatCapacity();
                //d cerr << "Tstar_V and Tstar_H were " << halite.Tstar_V() << "\t" << halite.Tstar_H() << endl;
                //d cerr << "NaClMeltingCurveHaliteLookup constructor, entry for t = " << tcurrent << " : survived cp\n";
                storage_vector[t_dim*compressibility_index + it] = halite.Compressibility();
                //d cerr << "Tstar_V and Tstar_H were " << halite.Tstar_V() << "\t" << halite.Tstar_H() << endl;
                //d cerr << "NaClMeltingCurveHaliteLookup constructor, entry for t = " << tcurrent << " : survived beta\n";
                storage_vector[t_dim*viscosity_index       + it] = 1.0e20;
                //d cerr << "Tstar_V and Tstar_H were " << halite.Tstar_V() << "\t" << halite.Tstar_H() << endl;
                //d cerr << "NaClMeltingCurveHaliteLookup constructor, entry for t = " << tcurrent << " : survived all\n";
              }
          }
        cerr << "vector computed, now writing file NaClMeltingCurveHaliteLookupTable.bin...\n";
        // write file
        fstream outfile(filename, ios::out | ios::binary);
        if( !outfile.is_open() )
          {
            cerr << "could not even it for writing, please stop program and debug !!!!\n";
            char yesno;
            cerr << "or enter any key to continue (simulation likely to crash or give wrong results!) :";
            cin  >> yesno;
          }
        else
          {
            cerr << "writing file " << filename << " ... ";
            skm_C_fwrite( outfile, storage_vector );
            outfile.close();
            cerr << "done!\n";
          }
      }
    else
      {
        cerr << "reading file " << filename << " ... ";
        skm_C_fread( infile, storage_vector );
		infile.close();
        cerr << "done!\n";
      }
    cerr << "NaClMeltingCurveHaliteLookup, leaving constructor ...\n\n";
  }

  NaClMeltingCurveHaliteLookup::~NaClMeltingCurveHaliteLookup()
  {
  }
    


  // The data interpolation routines
  double64 NaClMeltingCurveHaliteLookup::TmeltFromP(){ pcurrent = pressure; return TfromP(pcurrent); }
  double64 NaClMeltingCurveHaliteLookup::PmeltFromT(){      return ValueOf(pressure_index); }
  double64 NaClMeltingCurveHaliteLookup::MassFractionNaCl(){return ValueOf(composition_index); }
  double64 NaClMeltingCurveHaliteLookup::Density(){         return ValueOf(density_index); }
  double64 NaClMeltingCurveHaliteLookup::Enthalpy(){        return ValueOf(enthalpy_index); }
  double64 NaClMeltingCurveHaliteLookup::HeatCapacity(){    return ValueOf(heatcapacity_index); }
  double64 NaClMeltingCurveHaliteLookup::Compressibility(){ return ValueOf(compressibility_index); }
  double64 NaClMeltingCurveHaliteLookup::Viscosity(){       return ValueOf(viscosity_index); }
    
    
  double64 NaClMeltingCurveHaliteLookup::ValueOf(const int& property_index)
  {
    tcurrent = temperature;
    GetTemperatureIndex(tcurrent);
    tnorm    = (tcurrent - storage_vector[it])/(storage_vector[it+1]-storage_vector[it]);
    it      += t_dim*property_index; // convert to enthalpy section of storage vector
    return storage_vector[it]+tnorm*(storage_vector[it+1]-storage_vector[it]);
  }
    
    
  // data indexing
  void NaClMeltingCurveHaliteLookup::GetTemperatureIndex(const double64& t)
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
    it -= it_min;
    // t_dim is therefore (1000-550)/10+285 + 1 = 331
  }


  double64 NaClMeltingCurveHaliteLookup::TfromP(const double64& press)
  {
    i_max    = it_max;
    i_guess  = it_max/2;
    i_min    = 0;
    int        myi = 0;
    cerr.setf(ios::scientific);
    while( i_max-i_guess > 1)
      {	
        if(press == storage_vector[i_guess+t_dim*pressure_index])
          {
            return storage_vector[i_guess]; // if p is exactly correct at i_guess
          }
        else if(press < storage_vector[i_guess+t_dim*pressure_index])
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
            cerr << "NaClMeltingCurveLiquidLookup::TfromP - not converged ..." << endl;
            return 0.0e0;
          }
      }
    if(press < storage_vector[i_guess+t_dim*pressure_index]) i_guess -= 1;
    pnorm = ( press-storage_vector[i_guess+t_dim*pressure_index] ) / (storage_vector[i_guess+1+t_dim*pressure_index]-storage_vector[i_guess+t_dim*pressure_index]);
    return storage_vector[i_guess] + pnorm*(storage_vector[i_guess+1]-storage_vector[i_guess]); 
  }
    
    
} // namespace csmp
