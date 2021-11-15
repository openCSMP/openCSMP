#include <limits>

#include "ConvertConcentrationUnitsNaCl.h"
#include "CriticalCurveLookup.h"
#include "CriticalCurve.h"
#include "compareFloats.h"
#include "binaryReadWrite.h"
#include "Brine.h"
#include "LookupPropertyIndex.h"
#include "States.h"

using namespace std;

namespace csmp
{

  CriticalCurveLookup::CriticalCurveLookup(const double& externaltemperature)
    : temperature(externaltemperature),
      tcurrent(0.0),
      pcurrent(0.0),
      t_res(0.0),
      tnorm(0.0),
      pnorm(0.0),
      it(0),
      it_min(0),
      it_max(0),
      t_dim(0),
      i_guess(0),
      i_max(0),
      i_min(0),
      csmp_error( ErrorHandler::Instance() )
  {
    // the following sequence MUST be done as the very first thing!
    // also NEVER change it_min somewhere if you do changes to this file
    it_min   = 0;
    GetTemperatureIndex( cp_h2o.Temperature() );
    it_min   = it;
    cout << "it_min = " << it_min << endl;
    
    GetTemperatureIndex(1000.0001);
    it_max   = it;
    cout << "it_max = " << it_max << endl;
    t_dim    = it_max+1-it_min;
    cout << "t_dim  = " << t_dim << endl;
    //    char yesno; cout << "continue [y/n]:"; cin >> yesno;

    storage_vector.resize(t_dim*max_index);

    char filename[60];
    strcpy(filename,"CriticalCurveLookupTable.bin");
    fstream infile (filename, ios::in | ios::binary);
    cout << "point1\n"; 
    if ( !infile.is_open() )
      {
        cout <<"\nCriticalCurveLookup: file "<< filename;
        cout <<" could not be opened, computing ..."<< endl;
	
        double        xcurrent(0.0);

        CriticalCurve   critcurve(tcurrent);
        Brine           brine(tcurrent,pcurrent,xcurrent); 
                           
        // keep, needed for values from Brine object
        tcurrent = cp_h2o.Temperature();
        pcurrent = cp_h2o.Pressure();
        xcurrent = cp_h2o.MassFractionNaCl();

        cout << "computing for t = " << tcurrent << endl;
        it = 0;
        storage_vector[it]                             = tcurrent;
        storage_vector[it+t_dim*pressure_index]        = pcurrent;
        storage_vector[it+t_dim*composition_index]     = xcurrent;
        storage_vector[it+t_dim*density_index]         = cp_h2o.Density();
        storage_vector[it+t_dim*enthalpy_index]        = cp_h2o.Enthalpy();
        storage_vector[it+t_dim*heatcapacity_index]    = brine.HeatCapacity();    // value?
        storage_vector[it+t_dim*compressibility_index] = brine.Compressibility(); // value?
        storage_vector[it+t_dim*viscosity_index]       = brine.Viscosity();       // value?


        int myit = 1;
        GetTemperatureIndex(374.0001);
        for(tcurrent = 374.0e0; tcurrent < 1000.1e0; tcurrent += t_res)
          {
            GetTemperatureIndexCriticalCurve(tcurrent+2.0*numeric_limits<double>::epsilon());
            cout << "computing for t = " << tcurrent << ", t_res = " << t_res << ", data set " << it << "\t" << myit << endl;
            if(it > t_dim+it_min-1) break;

            pcurrent                                        = critcurve.Pressure();
            xcurrent                                        = critcurve.MassFractionNaCl();
            storage_vector[t_dim*temperature_index     +it] = tcurrent;
            storage_vector[t_dim*pressure_index        +it] = pcurrent;
            storage_vector[t_dim*composition_index     +it] = xcurrent;
            storage_vector[t_dim*density_index         +it] = brine.Density();
            storage_vector[t_dim*enthalpy_index        +it] = brine.Enthalpy();
            storage_vector[t_dim*heatcapacity_index    +it] = brine.HeatCapacity();
            storage_vector[t_dim*compressibility_index +it] = brine.Compressibility();
            storage_vector[t_dim*viscosity_index       +it] = brine.Viscosity();
            myit++;

          }
        cout << "vector computed, now writing file CriticalCurveLookupTable.bin...\n";
        // write file

        fstream outfile ( filename, ios::out | ios::binary );
        if( !outfile.is_open() )
          {
            cout << "could not even it for writing, please stop program and debug !!!!\n";
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
      }
    cout << "CriticalCurveLookup, leaving constructor ...\n\n";
  }


  CriticalCurveLookup::~CriticalCurveLookup()
  {
  }


  // The public interface
  double CriticalCurveLookup::Temperature(){      return ValueOf(temperature_index);                     }
  double CriticalCurveLookup::Pressure(){         return ValueOf(pressure_index);                        }
  double CriticalCurveLookup::MassFractionNaCl(){ return ValueOf(composition_index);                     }
  //  double CriticalCurveLookup::MoleFractionNaCl(){ return Massfraction2XNaCl(ValueOf(composition_index)); }
  double CriticalCurveLookup::Density(){          return ValueOf(density_index);                         }
  double CriticalCurveLookup::Enthalpy(){         return ValueOf(enthalpy_index);                        }
  double CriticalCurveLookup::HeatCapacity(){     return ValueOf(heatcapacity_index);                    }
  double CriticalCurveLookup::Compressibility(){  return ValueOf(compressibility_index);                 }
  double CriticalCurveLookup::Viscosity(){        return ValueOf(viscosity_index);                       }


  double CriticalCurveLookup::ValueOf( const int& property_index)
  {
    tcurrent = temperature;
    GetTemperatureIndexCriticalCurve(tcurrent);
    tnorm    = (tcurrent - storage_vector[it])/(storage_vector[it+1]-storage_vector[it]);
    it      +=  t_dim*property_index; 
    return storage_vector[it]+tnorm*(storage_vector[it+1]-storage_vector[it]);
  }    

    
  void CriticalCurveLookup::GetTemperatureIndex(const double& t)
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

  void CriticalCurveLookup::GetTemperatureIndexCriticalCurve(const double& t)
  {
    // the if-statements could be enclosed with preprocessor-idefs (e.g. referring to "verbose")
    if(!definitelyLessThan(    t, cp_h2o.Temperature(), numeric_limits<double>::epsilon() ) && 
       !definitelyGreaterThan( t, 1000.0e0, 10.0*numeric_limits<double>::epsilon() ) )
      {
        GetTemperatureIndex(t);
        it -= it_min;
        return;
      }
    else
      {
        if( definitelyLessThan (t, cp_h2o.Temperature(),numeric_limits<double>::epsilon() ) )
          {
            csmp_error.notice( FATAL_ERROR, 
                             "CriticalCurveLookup::GetTemperatureIndexCriticalCurve(const double& t) -",
                             "temperature is less than critical temperature for H2O and therefore out of range, terminating ...\n");
          }
        if( definitelyGreaterThan (t, 1000.0,10.0*numeric_limits<double>::epsilon() ) )
          {
            cerr << "delta is " << t-1000.0e0 << " for t = " << t << endl;
            csmp_error.notice( FATAL_ERROR, 
                             "CriticalCurveLookup::GetTemperatureIndexCriticalCurve(const double& t) -",
                             "temperature is higher than 1000 C and therefore out of range, terminating ...\n");
          }
      }
    return;
  }

  double CriticalCurveLookup::TfromP(const double& press)
  {
    double mypress = press;
    tcurrent = temperature; // debugging only
    i_max    = it_max;
    i_guess  = it_max/2;
    i_min    = 0;
    int        myi = 0;
    cout.setf(ios::scientific);
    while( i_max-i_guess > 1)
      {
        if(mypress == storage_vector[i_guess+t_dim*pressure_index])
          {
            return storage_vector[i_guess]; // if p is exactly correct at i_guess
          }
        else if(mypress < storage_vector[i_guess+t_dim*pressure_index])
          { i_max   = i_guess; 
            i_guess = (i_guess+i_min)/2; 
          }
        else{ i_min = i_guess; i_guess = (i_max+i_guess)/2; }
        myi ++;
        if(myi > 50)
          {
            cout << "CriticalCurveLookup::TfromP - not converged ..." << endl;
            return 0.0e0;
          }
      }
    if(mypress < storage_vector[i_guess+t_dim*pressure_index]) i_guess -= 1;
    pnorm = ( mypress-storage_vector[i_guess+t_dim*pressure_index] ) / 
      (storage_vector[i_guess+1+t_dim*pressure_index]-storage_vector[i_guess+t_dim*pressure_index]);
    return storage_vector[i_guess] + pnorm*(storage_vector[i_guess+1]-storage_vector[i_guess]);
  }
    
} // namespace csmp
