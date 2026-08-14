#include "HaliteLookup.h"
#include "binaryReadWrite.h"
#include "TriplePointNaCl.h"
#include "Halite.h"
#include "LookupPropertyIndex.h"
#include "States.h"

using namespace std;

namespace csmp
{
  HaliteLookup::HaliteLookup(const double& externaltemperature, const double& externalpressure)
    : 
    temperature_(externaltemperature), 
    pressure_(externalpressure), 
    tcurrent_(0.0), 
    pcurrent_(0.0), 
    tdummy_(0.0), 
    pdummy_(0.0),
    t_res_(0.0), 
    p_res_(0.0), 
    tnorm_(0.0), 
    pnorm_(0.0), 
    t_iA_(0.0), 
    t_iB_(0.0), 
    t_iC_(0.0), 
    t_iD_(0.0), 
    p_iA_(0.0), 
    p_iB_(0.0), 
    p_iC_(0.0), 
    p_iD_(0.0),
    v_bottom_(0.0), 
    v_top_(0.0), 
    v_interpolated_(0.0), 
    v_iA_(0.0), 
    v_iB_(0.0), 
    v_iC_(0.0), 
    v_iD_(0.0), 
    v_before_(0.0), 
    v_behind_(0.0), 
    v_vlh_(0.0), 
    v_vlh_behind_(0.0), 
    v_vlh_before_(0.0),
    tvlh_(0.0),
    it_(0), 
    ip_(0), 
    t_dim_(0), 
    p_dim_(0), 
    iA_(0), 
    iB_(0), 
    iC_(0), 
    iD_(0), 
    i_dummy_(0),
    state_(0), 
    state_iA_(0), 
    state_iB_(0), 
    state_iC_(0), 
    state_iD_(0),
    naclmelt_h_lookup(tdummy_,pdummy_)
  {
    GetTemperatureIndex(1000.0001);
    t_dim_ = it_+1;
    GetPressureIndex(5000.0e5+0.001);
    p_dim_ = ip_+1;
    //d cout << "HaliteLookup t_dim_ and p_dim_ = " << t_dim_ << "\t" << p_dim_ << endl;
	
    storage_vector.resize(t_dim_*p_dim_*max_index);
    state_vector.resize(t_dim_*p_dim_);

    char filename[60], statefilename[60];
    strcpy( filename, "HalitePropertiesLookupTable.bin" );
    strcpy( statefilename, "HaliteStateLookupTable.bin" );
    
	fstream infile1(filename, ios::in | ios::binary);
	fstream infile2(statefilename, ios::in | ios::binary);
	
    if( !infile1.is_open() || !infile2.is_open() )
      {
        cout << "HaliteLookup : at least one lookup file missing, computing ...\n\n";
        it_ = 0;
	    
        state_    = H;
	    
        Halite                       halite(tcurrent_,pcurrent_);
	    
        cout << "HaliteLookup::HaliteLookup(...): computing data ...\n";
        for(tcurrent_ = 0.0e0; tcurrent_ < 1000.1e0; tcurrent_ += t_res_)
          {
            GetTemperatureIndex(tcurrent_+1.0e-3);
            cout << "HaliteLookup computing for t = " << tcurrent_ << ", data set " << it_ << endl << endl;
            if(it_ > t_dim_-1)
              {
                cout << "too high t-index ...\n";
                break;
              }
		
            tdummy_ = tcurrent_;
		
            // ********************* if-statement for melting curve *************************
		
            for(pcurrent_ = 1.0e5; pcurrent_ <= 5000.0e5; pcurrent_ += p_res_)
              {
                pdummy_ = pcurrent_;
                GetPressureIndex(pcurrent_+1.0e-3);
                if(ip_ > p_dim_-1)
                  {
                    cout << "too high p-index ...\n";
                    break;
                  }
		    
                if( tcurrent_ > naclmelt_h_lookup.TmeltFromP() )
                  {
                    storage_vector[t_dim_*p_dim_*temperature_index     + it_*p_dim_ + ip_] = tcurrent_;
                    storage_vector[t_dim_*p_dim_*pressure_index        + it_*p_dim_ + ip_] = pcurrent_;
                    storage_vector[t_dim_*p_dim_*composition_index     + it_*p_dim_ + ip_] = 1.0; // XNaCl
                    storage_vector[t_dim_*p_dim_*density_index         + it_*p_dim_ + ip_] = halite.Density();
                    storage_vector[t_dim_*p_dim_*enthalpy_index        + it_*p_dim_ + ip_] = halite.Enthalpy();
                    storage_vector[t_dim_*p_dim_*heatcapacity_index    + it_*p_dim_ + ip_] = halite.HeatCapacity();
                    storage_vector[t_dim_*p_dim_*compressibility_index + it_*p_dim_ + ip_] = halite.Compressibility();
                    storage_vector[t_dim_*p_dim_*viscosity_index       + it_*p_dim_ + ip_] = 1.0e20;
                    state_vector[it_*p_dim_ + ip_ ] = L;
                  }
                else
                  {
                    storage_vector[t_dim_*p_dim_*temperature_index     + it_*p_dim_ + ip_] = tcurrent_;
                    storage_vector[t_dim_*p_dim_*pressure_index        + it_*p_dim_ + ip_] = pcurrent_;
                    storage_vector[t_dim_*p_dim_*composition_index     + it_*p_dim_ + ip_] = 1.0; // XNaCl
                    storage_vector[t_dim_*p_dim_*density_index         + it_*p_dim_ + ip_] = halite.Density();
                    storage_vector[t_dim_*p_dim_*enthalpy_index        + it_*p_dim_ + ip_] = halite.Enthalpy();
                    storage_vector[t_dim_*p_dim_*heatcapacity_index    + it_*p_dim_ + ip_] = halite.HeatCapacity();
                    storage_vector[t_dim_*p_dim_*compressibility_index + it_*p_dim_ + ip_] = halite.Compressibility();
                    storage_vector[t_dim_*p_dim_*viscosity_index       + it_*p_dim_ + ip_] = 1.0e20;
                    state_vector[it_*p_dim_ + ip_ ] = H;
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
    cout << "HaliteLookup, leaving constructor ...\n\n";
  }
    
    
  HaliteLookup::~HaliteLookup()
  {
  }

  // The data interpolation routines
  double HaliteLookup::MassFractionNaCl(){     SetTemperatureAndPressure(); return ValueOf(composition_index);    }
  double HaliteLookup::Density(){         SetTemperatureAndPressure(); return ValueOf(density_index);        }
  double HaliteLookup::Enthalpy(){        SetTemperatureAndPressure(); return ValueOf(enthalpy_index);       }
  double HaliteLookup::HeatCapacity(){    SetTemperatureAndPressure(); return ValueOf(heatcapacity_index);   }
  double HaliteLookup::Compressibility(){ SetTemperatureAndPressure(); return ValueOf(compressibility_index);}
  double HaliteLookup::Viscosity(){       SetTemperatureAndPressure(); return ValueOf(viscosity_index);      }



  void HaliteLookup::GetIndex_iA(const int& property_index)
  {
    // CAUTION: tcurrent_ and pcurrent_ must be know before this function is called !!!
    GetTemperatureIndex(tcurrent_);
    GetPressureIndex(pcurrent_);
    iA_  = t_dim_*p_dim_*property_index + it_*p_dim_ + ip_; 
    iB_  = iA_ + p_dim_;
    iC_  = iB_ + 1;
    iD_  = iA_ + 1;
	
    return;
  }
    
  void HaliteLookup::SetTemperatureAndPressure()
  {
    tcurrent_ = temperature_;
    pcurrent_ = pressure_;
    tdummy_   = tcurrent_;
    pdummy_   = pcurrent_;
    return;
  }
    
  double HaliteLookup::ValueOf(const int& property_index)
  {
    GetIndex_iA(property_index);
    //   cout << "it = " << it << ", ip_ = " << ip_ << endl;
	
    //CheckForOutOfRange();
    state_iA_ = state_vector[ it_*p_dim_ + ip_ ];
    state_iB_ = state_vector[ it_*p_dim_ + ip_ + p_dim_ ];
    state_iC_ = state_vector[ it_*p_dim_ + ip_ + p_dim_ +1 ];
    state_iD_ = state_vector[ it_*p_dim_ + ip_ + 1 ];
	
    i_dummy_  = t_dim_*p_dim_*temperature_index + it_*p_dim_ + ip_;
    //    cout << "i_dummy_ t = " << i_dummy_ << endl;
    t_iA_     = storage_vector[i_dummy_];
    t_iB_     = storage_vector[i_dummy_ + p_dim_];
    t_iC_     = storage_vector[i_dummy_ +1 + p_dim_];
    t_iD_     = storage_vector[i_dummy_ + 1];
	
    i_dummy_  = t_dim_*p_dim_*pressure_index + it_*p_dim_ + ip_;
    //    cout << "i_dummy_ p = " << i_dummy_ << endl;
    p_iA_     = storage_vector[i_dummy_];
    p_iB_     = storage_vector[i_dummy_ + p_dim_];
    p_iD_     = storage_vector[i_dummy_ + 1];
    p_iC_     = p_iD_;
	
    //    cout << "Temperatures and Pressures iA-iD:\n";
    //     cout << t_iA_ << "\t" << t_iB_ << "\t" << t_iC_ << "\t" << t_iD_ << endl;
    //     cout << p_iA_ << "\t" << p_iB_ << "\t" << p_iC_ << "\t" << p_iD_ << endl;
	
    // *****************************************************************************
    // in principle, NormalInterpolation( const int& property_index ) is sufficient
    // since halite properties have also been calculated beyond the melting curve
    // *****************************************************************************
	
    //    if( (state_iA_ == H) && (state_iB_ == H) && (state_iC_ == H) && (state_iD_ == H) ){
    //	cout << state_iA_ << "\t" << state_iB_ << "\t" << state_iC_ << "\t" << state_iD_ << endl;
    return NormalInterpolation(property_index); 
    //    }
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
	
    //     else if( (state_iD_ == H) && (state_iB_ == M) ){
    // 	return NearNaClMeltInterpolation(property_index);
    //     }
    //     else{
    // 	cout << "HaliteLookup::ValueOf missed to find proper interpolation type, returning bogus value ...\n";
    // 	return 9.9e99;
    //     }
  }
    
    
  double HaliteLookup::NormalInterpolation( const int& property_index )
  {
    //    cout << "Using HaliteLookup::NormalInterpolation( const int& property_index ) ...\n";
	
    //   P
    //   ^
    //   |         |
    //  -D---------C-
    //   |  x      |
    //   |         |
    //  -A---------B--> T
    //
	
    // ****** now compute iA  ******
	
    tnorm_          = (tcurrent_ - t_iA_ ) / ( t_iB_ - t_iA_ );
    pnorm_          = (pcurrent_ - p_iA_ ) / ( p_iD_ - p_iA_ );
    v_bottom_       =  storage_vector[iA_] + tnorm_*(storage_vector[iB_] - storage_vector[iA_]);
    v_top_          =  storage_vector[iD_] + tnorm_*(storage_vector[iC_] - storage_vector[iD_]);
    v_interpolated_ =  v_bottom_ + pnorm_*(v_top_-v_bottom_);
    return v_interpolated_;
  }
    
    
    
    
    
  double HaliteLookup::NearNaClMeltInterpolation( const int& property_index )
  {
    //    cout << "Using HaliteLookup::NearNaClMeltInterpolation( const int& property_index ) ...\n";
    pdummy_ = pcurrent_;
    tvlh_      = naclmelt_h_lookup.TmeltFromP(); 
    if(tvlh_  <= t_iC_)
      {
	    
        tdummy_    = t_iD_;
        if(naclmelt_h_lookup.PmeltFromT() <= p_iA_)
          {
            // case (1) and (2): Interpolation to be done between A-D and vlh_curve_at_pcurrent_
            //
            //  (1)                                (2)                          
            //                                   
            //   P                                  P                             
            //   ^            meltcurve             ^         meltcurve          
            //   |         | /                      |        /                    
            //  -D---------C/-                     -D-------/-C-              
            //   |         |                        |      /  |           
            //   |        /|                        |     /   |        
            // --|--x----/-|------ pcurrent_      --|-x--/----|------ pcurrent_  
            //   |      /  |                        |   /     |         
            //  -A-----/---B--> T                  -A--/------B--> T      
            //   |    /    |                        |         |   
            //
		
            v_iA_              = storage_vector[iA_];
            v_iD_              = storage_vector[iD_];
		
            pnorm_             = (pcurrent_-p_iD_) / (p_iA_-p_iD_);
            v_before_          = storage_vector[iD_] + pnorm_*( storage_vector[iA_] - storage_vector[iD_] );
            tdummy_            = tvlh_;
            v_vlh_             = naclmelt_h_lookup.ValueOf(property_index);
            tnorm_             = (tcurrent_-t_iA_) / (tvlh_-t_iA_);
            v_interpolated_    =  v_before_ + tnorm_*(v_vlh_-v_before_);
            //d cout << "computed via if-2 as v_before_ = " << v_before_ << ", " << "v_vlh = " 
            //d << v_vlh << ", " << v_interpolated_ << endl;
          }
        else
          {
            // case (3) and (4): Interpolation to be done between vlh_at_t_iD_-D-segment 
            // and vlh_curve_at_pcurrent_
            //
            //      (4)                              (2)
            //                                                               
            //   P                                P
            //   ^    meltcurve                   ^
            //   |   /                            |         |   meltcurve
            //  -D--/------C-                    -D---------C__/
            // --|x/-------|------ pcurrent_      |     ___/| 
            //   |/        |                    --|-x__/----|------ pcurrent_
            //   /         |                    __|/        |
            //  /|         |                   /  |         |
            //  -A---------B--> T                -A---------B--> T
            //   |         |                      |         |  
            //
	    
            // notice: it is still tdummy_ = t_iD_;
            pnorm_             = (pcurrent_-p_iD_) / (naclmelt_h_lookup.PmeltFromT()*1.0e-5-p_iD_);
            v_iD_              = storage_vector[iD_];
            v_vlh_before_      = naclmelt_h_lookup.ValueOf(property_index); // (between A and D)
            tdummy_            = tvlh_;
            v_vlh_             = naclmelt_h_lookup.ValueOf(property_index);
            v_before_          =  v_iD_ + pnorm_*( v_vlh_before_ - v_iD_ );
            tnorm_             = (tcurrent_-t_iD_) / (tvlh_-t_iD_);
            v_interpolated_    =  v_before_ + tnorm_*(v_vlh_-v_before_);
            //d cout << "v_iD_ = " << v_iD_ << ", tnorm_ = " << tnorm_ << endl;
            //d cout << "computed via if-1 as v_before_ = " << v_before_ << ", " 
            //d      << "v_vlh = " << v_vlh << ", " << v_interpolated_ << endl;
          }
      }
    else
      {
        // tvlh_ is > t_iB_ 
        tdummy_    = t_iD_;
        if(naclmelt_h_lookup.PmeltFromT()*1.0e-5 >= p_iA_)
          {
            // case (5): Interpolation to be done between D-vlh_at_t_iD_-segment and C-vlh_at_t_iC_-segment
            //
            //  (5)                      
            //
            //   P                        
            //   ^                        
            //   |                        
            //  -D---------C ___  meltcurve  
            // --|-----x--_|/-----pcurrent_           
            //   | _____/  |              
            // __|/        |              
            //   |         |     
            //  -A---------B--> T         
            //   |         |              
            //
            //
            pnorm_             = (pcurrent_-p_iD_) / (naclmelt_h_lookup.PmeltFromT()*1.0e-5-p_iD_);
            v_iC_              = storage_vector[iC_];
            v_iD_              = storage_vector[iD_];
            v_vlh_before_      = naclmelt_h_lookup.ValueOf(property_index); // (between A and D)
            tdummy_            = t_iB_;
            v_vlh_behind_      = naclmelt_h_lookup.ValueOf(property_index);

            v_before_          =  v_iD_ + pnorm_*( v_vlh_before_ - v_iD_ );
            tdummy_            =  t_iC_;
            pnorm_             = (pcurrent_-p_iC_) / (naclmelt_h_lookup.PmeltFromT()*1.0e-5-p_iC_);
            v_behind_          =  v_iC_ + pnorm_*( v_vlh_behind_ - v_iC_); 
            tnorm_             = (tcurrent_-t_iD_) / (t_iB_-t_iD_);
            v_interpolated_    =  v_before_ + tnorm_*(v_behind_-v_before_);
          }
        else
          {
	    
            //  (6)
            //
            //   P
            //   ^            meltcurve
            //   |           /
            //  -D--------C-/
            //  -|--x-----|/--- pcurrent_
            //   |        /
            //   |       /|
            //   |      / |
            //  -A-----/--B--> T
            //   |    /   |        
            //
            pnorm_             = (pcurrent_-p_iD_) / (p_iA_-p_iD_);
            v_iD_              = storage_vector[iD_];
            v_iC_              = storage_vector[iC_];
            v_iA_              = storage_vector[iA_];
	    
            tdummy_            = t_iC_;
            v_vlh_behind_      = naclmelt_h_lookup.ValueOf(property_index);
	    
            v_before_          =  v_iD_ + pnorm_*( v_iA_ - v_iD_ );
            tdummy_            =  t_iC_;
            pnorm_             = (pcurrent_-p_iC_) / (naclmelt_h_lookup.PmeltFromT()*1.0e-5-p_iC_);
            v_behind_          =  v_iC_ + pnorm_*( v_vlh_behind_ - v_iC_); 
            tnorm_             = (tcurrent_-t_iD_) / (t_iB_-t_iD_);
            v_interpolated_    =  v_before_ + tnorm_*(v_behind_-v_before_);
          }
      }
    return v_interpolated_;
  }
    
  void HaliteLookup::GetTemperatureIndex(const double& t)
  {
    // new version
    if(     t <    0.0e0){  t_res_ =  5.0; it_ = 0; }
    else if(t <= 250.0e0){  t_res_ =  5.0; it_ =     static_cast<long>( (t-  0.0)/t_res_ ); }
    else if(t <= 350.0e0){  t_res_ =  2.0; it_ =  50+static_cast<long>( (t-250.0)/t_res_ ); }
    else if(t <= 360.0e0){  t_res_ =  1.0; it_ = 100+static_cast<long>( (t-350.0)/t_res_ ); }
    else if(t <= 370.0e0){  t_res_ =  0.5; it_ = 110+static_cast<long>( (t-360.0)/t_res_ ); }
    else if(t <= 372.0e0){  t_res_ =  0.2; it_ = 130+static_cast<long>( (t-370.0)/t_res_ ); }
    else if(t <= 378.0e0){  t_res_ =  0.1; it_ = 140+static_cast<long>( (t-372.0)/t_res_ ); }
    else if(t <= 380.0e0){  t_res_ =  0.2; it_ = 200+static_cast<long>( (t-378.0)/t_res_ ); }
    else if(t <= 390.0e0){  t_res_ =  0.5; it_ = 210+static_cast<long>( (t-380.0)/t_res_ ); }
    else if(t <= 400.0e0){  t_res_ =  1.0; it_ = 230+static_cast<long>( (t-390.0)/t_res_ ); }
    else if(t <= 450.0e0){  t_res_ =  2.0; it_ = 240+static_cast<long>( (t-400.0)/t_res_ ); }
    else if(t <= 550.0e0){  t_res_ =  5.0; it_ = 265+static_cast<long>( (t-450.0)/t_res_ ); }
    else if(t <= 1000.0e0){ t_res_ = 10.0; it_ = 285+static_cast<long>( (t-550.0)/t_res_ ); }
    else{ it_ = 330; }//                  t_res_ = 10.0; it_ = 285+static_cast<long>( (t-550.0)/t_res_ ); /* throw out of range ? */ }
    // t_dim_ is therefore (1000-550)/10+285 + 1 = 331
  }
    
  void HaliteLookup::GetPressureIndex(const double& p)
  {
    // new version
    if(     p <=   20.0e5){ p_res_ =   0.5e5; ip_ =     static_cast<long>( (p-   0.5e5)/p_res_ ); }    
    else if(p <=  210.0e5){ p_res_ =   1.0e5; ip_ =  39+static_cast<long>( (p-  20.0e5)/p_res_ ); }
    else if(p <=  215.0e5){ p_res_ =   0.5e5; ip_ = 229+static_cast<long>( (p- 210.0e5)/p_res_ ); }
    else if(p <=  225.0e5){ p_res_ =   0.1e5; ip_ = 239+static_cast<long>( (p- 215.0e5)/p_res_ ); }
    else if(p <=  230.0e5){ p_res_ =   0.5e5; ip_ = 339+static_cast<long>( (p- 225.0e5)/p_res_ ); }
    else if(p <=  250.0e5){ p_res_ =   1.0e5; ip_ = 349+static_cast<long>( (p- 230.0e5)/p_res_ ); }
    else if(p <=  300.0e5){ p_res_ =   2.0e5; ip_ = 369+static_cast<long>( (p- 250.0e5)/p_res_ ); }
    else if(p <=  400.0e5){ p_res_ =   5.0e5; ip_ = 394+static_cast<long>( (p- 300.0e5)/p_res_ ); }
    else if(p <=  700.0e5){ p_res_ =  10.0e5; ip_ = 414+static_cast<long>( (p- 400.0e5)/p_res_ ); }
    else if(p <= 1000.0e5){ p_res_ =  25.0e5; ip_ = 444+static_cast<long>( (p- 700.0e5)/p_res_ ); }
    else if(p <= 2000.0e5){ p_res_ =  50.0e5; ip_ = 456+static_cast<long>( (p-1000.0e5)/p_res_ ); }
    else if(p <= 5000.0e5){ p_res_ = 100.0e5; ip_ = 476+static_cast<long>( (p-2000.0e5)/p_res_ ); }
    else{ ip_ = 506; }//                   p_res_ = 100.0; ip_ = 476+static_cast<long>( (p-2000.0)/p_res_ );/* throw out of range ? */ }
    // p_dim_ is therefore (5000-2000)/100+476 + 1 = 507
  }


} // namespace csmp
