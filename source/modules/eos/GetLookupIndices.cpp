#include "GetLookupIndices.h"

using namespace std;

namespace csmp
{
  long GetTemperatureIndex(const double& t)
  {
    // new version
    double t_res;
    long      it;
    if(     t <    0.0e0)
      {
	cerr << "H2OLookup::GetTemperatureIndex(const double& t) : t < 0 (t = " << t << "), better terminate ...\n";
	cerr << "or to continue, enter any key : ";
	char yesno;
	cin >> yesno;
      }
    else if(t <= 250.0e0){  t_res =  5.0; it =     static_cast<long>( (t-   0.0)/t_res ); }
    else if(t <= 350.0e0){  t_res =  2.0; it =  50+static_cast<long>( (t- 250.0)/t_res ); }
    else if(t <= 360.0e0){  t_res =  1.0; it = 100+static_cast<long>( (t- 350.0)/t_res ); }
    else if(t <= 370.0e0){  t_res =  0.5; it = 110+static_cast<long>( (t- 360.0)/t_res ); }
    else if(t <= 372.0e0){  t_res =  0.2; it = 130+static_cast<long>( (t- 370.0)/t_res ); }
    else if(t <= 378.0e0){  t_res =  0.1; it = 140+static_cast<long>( (t- 372.0)/t_res ); }
    else if(t <= 380.0e0){  t_res =  0.2; it = 200+static_cast<long>( (t- 378.0)/t_res ); }
    else if(t <= 390.0e0){  t_res =  0.5; it = 210+static_cast<long>( (t- 380.0)/t_res ); }
    else if(t <= 400.0e0){  t_res =  1.0; it = 230+static_cast<long>( (t- 390.0)/t_res ); }
    else if(t <= 450.0e0){  t_res =  2.0; it = 240+static_cast<long>( (t- 400.0)/t_res ); }
    else if(t <= 550.0e0){  t_res =  5.0; it = 265+static_cast<long>( (t- 450.0)/t_res ); }
    else if(t <= 1000.0e0){ t_res = 10.0; it = 285+static_cast<long>( (t- 550.0)/t_res ); }
    else if(t <= 2000.0e0){ t_res = 50.0; it = 330+static_cast<long>( (t-1000.0)/t_res ); }
    else{ it = 350; }; //                  t_res = 10.0; it = 285+static_cast<long>( (t-550.0)/t_res ); /* throw out of range ? */ }
    // t_dim is therefore (2000-1000)/50+330 + 1 = 351
    return it;
  }
    
  long GetPressureIndex(const double& p)
  {
    double p_res;
    long      ip;
    /*!
      New version, August 24, 2009, Thomas Driesner

      In order to resolve low-Temperature vapor reasonably well, the resolution
      in that range could either be increased or a fake "0 bar" value be invented.
      I decided to try the latter, i.e., i will use - in Table2 - a value at 0.005 bar, 
      which is always below psat. This increases p_dim by 1 compared to earlier versions.
      This should probably also be invented in the H2O-NaCl Tables ...
    */ 
    if(     p <=   20.0e0){  p_res =   0.5; ip =     static_cast<long>( (p       )/p_res ); }    
    else if(p <=  210.0e0){  p_res =   1.0; ip =  40+static_cast<long>( (p-  20.0)/p_res ); }
    else if(p <=  215.0e0){  p_res =   0.5; ip = 230+static_cast<long>( (p- 210.0)/p_res ); }
    else if(p <=  225.0e0){  p_res =   0.1; ip = 240+static_cast<long>( (p- 215.0)/p_res ); }
    else if(p <=  230.0e0){  p_res =   0.5; ip = 340+static_cast<long>( (p- 225.0)/p_res ); }
    else if(p <=  250.0e0){  p_res =   1.0; ip = 350+static_cast<long>( (p- 230.0)/p_res ); }
    else if(p <=  300.0e0){  p_res =   2.0; ip = 370+static_cast<long>( (p- 250.0)/p_res ); }
    else if(p <=  400.0e0){  p_res =   5.0; ip = 395+static_cast<long>( (p- 300.0)/p_res ); }
    else if(p <=  700.0e0){  p_res =  10.0; ip = 415+static_cast<long>( (p- 400.0)/p_res ); }
    else if(p <= 1000.0e0){  p_res =  25.0; ip = 445+static_cast<long>( (p- 700.0)/p_res ); }
    else if(p <= 2000.0e0){  p_res =  50.0; ip = 457+static_cast<long>( (p-1000.0)/p_res ); }
    else if(p <= 5000.0e0){  p_res = 100.0; ip = 477+static_cast<long>( (p-2000.0)/p_res ); }
    else if(p <= 10000.0e0){ p_res = 500.0; ip = 507+static_cast<long>( (p-5000.0)/p_res ); }
    else{ ip = 517; }//                  p_res = 100.0; ip = 476+static_cast<long>( (p-2000.0)/p_res );/* throw out of range ? */ }
    // p_dim is therefore (10000-5000)/500+506 + 1 = 518
    return ip;
  }


  double GetTemperatureResolution(const double& t)
  {
    // new version
    double t_res;
    if(     t <= 250.0e0){  t_res =  5.0; }
    else if(t <= 350.0e0){  t_res =  2.0; }
    else if(t <= 360.0e0){  t_res =  1.0; }
    else if(t <= 370.0e0){  t_res =  0.5; }
    else if(t <= 372.0e0){  t_res =  0.2; }
    else if(t <= 378.0e0){  t_res =  0.1; }
    else if(t <= 380.0e0){  t_res =  0.2; }
    else if(t <= 390.0e0){  t_res =  0.5; }
    else if(t <= 400.0e0){  t_res =  1.0; }
    else if(t <= 450.0e0){  t_res =  2.0; }
    else if(t <= 550.0e0){  t_res =  5.0; }
    else if(t <= 1000.0e0){ t_res = 10.0; }
    else{                   t_res = 50.0; }
    // t_dim is therefore (1000-550)/10+285 + 1 = 331
    return t_res;
  }
    
  double GetPressureResolution(const double& p)
  {
    // new version
    double p_res;
    if(     p <=   20.0e0){ p_res =   0.5; }    
    else if(p <=  210.0e0){ p_res =   1.0; }
    else if(p <=  215.0e0){ p_res =   0.5; }
    else if(p <=  225.0e0){ p_res =   0.1; }
    else if(p <=  230.0e0){ p_res =   0.5; }
    else if(p <=  250.0e0){ p_res =   1.0; }
    else if(p <=  300.0e0){ p_res =   2.0; }
    else if(p <=  400.0e0){ p_res =   5.0; }
    else if(p <=  700.0e0){ p_res =  10.0; }
    else if(p <= 1000.0e0){ p_res =  25.0; }
    else if(p <= 2000.0e0){ p_res =  50.0; }
    else if(p <= 5000.0e0){ p_res = 100.0; }
    else{                   p_res = 500.0; }
    return p_res;
  }

}
