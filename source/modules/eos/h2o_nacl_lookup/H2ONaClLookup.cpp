#include <iostream>
#include <cmath>
#include <string>

#include "H2ONaClLookup.h"
#include "compareFloats.h"
#include "binaryReadWrite.h"
#include "LookupPropertyIndex.h"
#include "ConvertConcentrationUnitsNaCl.h"
#include "H2OLookup.h"
#include "CriticalCurveLookup.h"
#include "CriticalPointH2O.h"
#include "VLH_VaporLookup.h"
#include "VLH_LiquidLookup.h"
#include "Brine.h"
#include "HaliteLiquidusLookup.h"
#include "HaliteLiquidus.h"
#include "LookupPropertyIndex.h"
#include "TwophaseLiquidLookup.h"
#include "TwophaseVaporLookup.h"
#include "NaClSaturatedVaporLookup.h"


using namespace std;

namespace csmp
{
    
  H2ONaClLookup::H2ONaClLookup(const double64& externaltemperature, 
                               const double64& externalpressure, 
                               const double64& externalcomposition)
    : t_dim(0), // to be re-computed in constructor body
      p_dim(0), // to be re-computed in constructor body
      x_dim(83),
      x_dim200(37),
      x_dim400(44),
      x_dim600(57),
      x_dim700(69),
      it_200(39),
      it_400(251),
      it_600(301),
      it_700(311),
      n_tables(6),
      table_id(0),
      i(0),
      it(0),
      ip(0),
      ix(0),
      iA(0),
      iB(0),
      iC(0),
    iD(0),
    iE(0),
    iF(0),
    iG(0),
    iH(0),
    block_id(0),
    iA_t(0),
    iB_t(0),
    iA_p(0),
    iD_p(0),
    iA_x(0),
    iE_x(0),
    storage_index(0),
    ip_boiling(0),
    storageindex_lowp1(0),
    storageindex_lowp2(0),
    storageindex_highp1(0),
    storageindex_highp2(0),
    ix_twophase(0),
    storage_index_left(0),
    ix_old(0),
    ip_low(0),
    ip_high(0),
    state(L),
    temperature(externaltemperature),
    pressure(externalpressure),
    composition(externalcomposition),
    tcurrent(temperature-1.0),
    pcurrent(pressure-1.0),
    xcurrent(composition+1.0),
    xdummy(0.0),
    pdummy(0.0),
    tdummy(0.0),
    t_res(0.0),
    p_res(0.0),
    x_res(0.0),
    t_iA(0.0),
    t_iB(0.0),
    p_iA(0.0),
    p_iD(0.0),
    x_iA(0.0),
    x_iE(0.0),
    tnorm(0.0),
    pnorm(0.0),
    xnorm(0.0),
    v_bottom(0.0),
    v_top(0.0),
    v_right(0.0),
    v_left(0.0),
    v_point(0.0),
    dh(0.0),
    DP_extrapol(0.0),
    DP_lowp(0.0),
    DDensityDP(0.0),
    DEnthalpyDP(0.0),
    DHeatCapacityDP(0.0),
    DCompressibilityDP(0.0),
    DViscosityDP(0.0),
    DP_highp(0.0),
  // density_left(0.0),
  // density_right(0.0),
  // enthalpy_left(0.0),
  // enthalpy_right(0.0),
  // heatcapacity_left(0.0),
  // heatcapacity_right(0.0),
  // compressibility_left(0.0),
  // compressibility_right(0.0),
  // viscosity_left(0.0),
  // viscosity_right(0.0),
  // dx_left(0.0),
  // dx_right(0.0),
  // x_left(0.0),
  // x_right(0.0),
    p_boil(0.0)      
  { 

    // cout << "H2ONaClLookupDevel.cpp::Constructor - after initialization list\n";
    t_dim = ComputeTemperatureIndex(1000.0001)+1;
    p_dim = ComputePressureIndex(5000.e5+0.01)+1;
    ip_vl.resize(x_dim);

    strcpy(filename[0],"H2ONaClPropertiesLookupTableSubcriticalVapor.bin");
    strcpy(filename[1],"H2ONaClPropertiesLookupTableLiquidTo200C.bin");
    strcpy(filename[2],"H2ONaClPropertiesLookupTableLiquidTo400C.bin");
    strcpy(filename[3],"H2ONaClPropertiesLookupTableTo600.bin");
    strcpy(filename[4],"H2ONaClPropertiesLookupTableTo700.bin");
    strcpy(filename[5],"H2ONaClPropertiesLookupTableTo1000.bin");

    t_min_table.resize(n_tables);
    t_min_table[0] =    0.0e0;
    t_min_table[1] =    0.0e0;
    t_min_table[2] =  200.0e0;
    t_min_table[3] =  400.0e0;
    t_min_table[4] =  600.0e0;
    t_min_table[5] =  700.0e0;

    t_max_table.resize(n_tables);
    t_max_table[0] =  376.0e0;
    t_max_table[1] =  200.0e0;
    t_max_table[2] =  400.0e0;
    t_max_table[3] =  600.0e0;
    t_max_table[4] =  700.0e0;
    t_max_table[5] = 1000.0e0;

    t_dim_table.resize(n_tables);
    it_offset.resize(n_tables);
    it_max_table.resize(n_tables);

    for(i = 0; i < n_tables; i++) t_dim_table[i]  = ComputeTemperatureIndex(t_max_table[i])-ComputeTemperatureIndex(t_min_table[i])+1;
    for(i = 0; i < n_tables; i++) it_offset[i]    = ComputeTemperatureIndex(t_min_table[i]);
    for(i = 0; i < n_tables; i++) it_max_table[i] = ComputeTemperatureIndex(t_max_table[i]);

    x_dim_table.resize(n_tables);
    x_dim_table[0]  = 4;  // H2ONaClPropertiesLookupTableSubcriticalVapor.bin
    x_dim_table[1]  = 37; // H2ONaClPropertiesLookupTableLiquidTo200C.bin
    x_dim_table[2]  = 44; // H2ONaClPropertiesLookupTableLiquidTo400C.bin
    x_dim_table[3]  = 57; // H2ONaClPropertiesLookupTableTo600.bin
    x_dim_table[4]  = 69; // H2ONaClPropertiesLookupTableTo700.bin
    x_dim_table[5]  = 83; // H2ONaClPropertiesLookupTableTo1000.bin

    tablesize.resize(n_tables);
    tablesize[0] = t_dim_table[0] * p_dim * x_dim_table[0]; 
    tablesize[1] = t_dim_table[1] * p_dim * x_dim_table[1]; 
    tablesize[2] = t_dim_table[2] * p_dim * x_dim_table[2]; 
    tablesize[3] = t_dim_table[3] * p_dim * x_dim_table[3]; 
    tablesize[4] = t_dim_table[4] * p_dim * x_dim_table[4]; 
    tablesize[5] = t_dim_table[5] * p_dim * x_dim_table[5]; 

    table0.resize(tablesize[0]*max_index);
    table1.resize(tablesize[1]*max_index);
    table2.resize(tablesize[2]*max_index);
    table3.resize(tablesize[3]*max_index);
    table4.resize(tablesize[4]*max_index);
    table5.resize(tablesize[5]*max_index);
      	
    BuildTable0();
    BuildTable1();
    BuildTable2();
    BuildSimpleTable( 3, table3 );
    BuildSimpleTable( 4, table4 );
    BuildSimpleTable( 5, table5 );
    cout << "H2ONaClLookup, leaving constructor ...\n\n";

    return;
  }


    
    
  H2ONaClLookup::~H2ONaClLookup()
  {
  }
    
  // need if statements to determine table_id
  double64 H2ONaClLookup::Temperature(){                      return Interpolate( temperature, temperature_index     );}
  double64 H2ONaClLookup::Pressure(){                         return Interpolate( temperature, pressure_index        );}
  double64 H2ONaClLookup::MassFractionNaCl(){                      return Interpolate( temperature, composition_index     );}
  double64 H2ONaClLookup::Density(){                          return Interpolate( temperature, density_index         );}
  double64 H2ONaClLookup::Enthalpy(){                         return Interpolate( temperature, enthalpy_index        );}
  double64 H2ONaClLookup::HeatCapacity(){                     return Interpolate( temperature, heatcapacity_index    );}
  double64 H2ONaClLookup::Compressibility(){                  return Interpolate( temperature, compressibility_index );}
  double64 H2ONaClLookup::Viscosity(){                        return Interpolate( temperature, viscosity_index       );}

  double64 H2ONaClLookup::SubcriticalVaporDensity(){          return TrilinearInterpolation( 0, table0, density_index         );}
  double64 H2ONaClLookup::SubcriticalVaporEnthalpy(){         return TrilinearInterpolation( 0, table0, enthalpy_index        );}
  double64 H2ONaClLookup::SubcriticalVaporHeatCapacity(){     return TrilinearInterpolation( 0, table0, heatcapacity_index    );}
  double64 H2ONaClLookup::SubcriticalVaporCompressibility(){  return TrilinearInterpolation( 0, table0, compressibility_index );}
  double64 H2ONaClLookup::SubcriticalVaporViscosity(){        return TrilinearInterpolation( 0, table0, viscosity_index       );}

  double64 H2ONaClLookup::SubcriticalLiquidDensity()
  {         
    if(tcurrent < 200.0e0) return TrilinearInterpolation( 1, table1, density_index         );
    else                   return TrilinearInterpolation( 2, table2, density_index         );
  }
  double64 H2ONaClLookup::SubcriticalLiquidEnthalpy()
  {
    if(tcurrent < 200.0e0) return TrilinearInterpolation( 1, table1, enthalpy_index        );
    else                   return TrilinearInterpolation( 2, table2, enthalpy_index        );
  }
  double64 H2ONaClLookup::SubcriticalLiquidHeatCapacity()
  {
    if(tcurrent < 200.0e0) return TrilinearInterpolation( 1, table1, heatcapacity_index    );
    else                   return TrilinearInterpolation( 2, table2, heatcapacity_index    );
  }
  double64 H2ONaClLookup::SubcriticalLiquidCompressibility()
  {
    if(tcurrent < 200.0e0) return TrilinearInterpolation( 1, table1, compressibility_index );
    else                   return TrilinearInterpolation( 2, table2, compressibility_index );
  }
  double64 H2ONaClLookup::SubcriticalLiquidViscosity()
  {
    if(tcurrent < 200.0e0) return TrilinearInterpolation( 1, table1, viscosity_index       );
    else                   return TrilinearInterpolation( 2, table2, viscosity_index       );
  }
    
  double64 H2ONaClLookup::Interpolate( const double64& t,  const int& property_index )
  {
    if(     t < 200.0e0){ table_id = 1; return TrilinearInterpolation( table_id, table1, property_index ); }
    else if(t < 400.0e0){ table_id = 2; return TrilinearInterpolation( table_id, table2, property_index ); }
    else if(t < 600.0e0){ table_id = 3; return TrilinearInterpolation( table_id, table3, property_index ); }
    else if(t < 700.0e0){ table_id = 4; return TrilinearInterpolation( table_id, table4, property_index ); }
    else                { table_id = 5; return TrilinearInterpolation( table_id, table5, property_index ); }
  }
  // ***** the "throw out of range ?" comment below needs really the "?" because at least the t_res info is queried when constructing the table


  void H2ONaClLookup::GetTemperatureIndex(const double64& t)
  {
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
    return;
  }
    
  void H2ONaClLookup::GetPressureIndex(const double64& p)
  {
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
    return;
  }


  long H2ONaClLookup::ComputeTemperatureIndex(const double64& t)
  {
    // new version
    double64 t_res;
    long      it;
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
    return it;
  }
    
  long H2ONaClLookup::ComputePressureIndex(const double64& p)
  {
    // new version
    double64 p_res;
    long      ip;
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
    return ip;
  }
    
  //CAUTION, THE NEXT TO ARE STILL MOLE-FRACTION-BASED!!!
  void H2ONaClLookup::GetCompositionIndex(const double64& x)
  {
    if(     x <= 1.0e-3){ x_res = 1.0e-4; ix =    static_cast<long>((x-0.0e0 )/x_res); }
    else if(x <= 1.0e-2){ x_res = 1.0e-3; ix = 10+static_cast<long>((x-1.0e-3)/x_res); }
    else if(x <= 5.0e-2){ x_res = 5.0e-3; ix = 19+static_cast<long>((x-1.0e-2)/x_res); }
    else if(x <= 2.0e-1){ x_res = 1.0e-2; ix = 27+static_cast<long>((x-5.0e-2)/x_res); }
    else if(x <= 1.0e0 ){ x_res = 0.02e0; ix = 42+static_cast<long>((x-2.0e-1)/x_res); } 
    else{ ix = 82; }//                x_res = 0.02e0; ix = 42+static_cast<long>((x-2.0e-1)/x_res); /* throw out of range ? */ }
    return;
  }

  long H2ONaClLookup::ComputeCompositionIndex(const double64& x)
  {
    long myix;
    if(     x <= 1.0e-3){ x_res = 1.0e-4; myix =    static_cast<long>((x-0.0e0 )/x_res); }
    else if(x <= 1.0e-2){ x_res = 1.0e-3; myix = 10+static_cast<long>((x-1.0e-3)/x_res); }
    else if(x <= 5.0e-2){ x_res = 5.0e-3; myix = 19+static_cast<long>((x-1.0e-2)/x_res); }
    else if(x <= 2.0e-1){ x_res = 1.0e-2; myix = 27+static_cast<long>((x-5.0e-2)/x_res); }
    else if(x <= 1.0e0 ){ x_res = 0.02e0; myix = 42+static_cast<long>((x-2.0e-1)/x_res); } 
    else{ myix =  82; }//                x_res = 0.02e0; myix = 42+static_cast<long>((x-2.0e-1)/x_res); /* throw out of range ? */ }
    return myix;
  }
  
  
  
  
  double64 H2ONaClLookup::TrilinearInterpolation( const long& table_id, 
                                                  const std::vector<double64>& storage_vector, 
                                                  const int& property_index )
  {
    // maybe get rid of if-statement anyhow
    // if( !(tcurrent == temperature && pcurrent == pressure && xcurrent == Massfraction2XNaCl(composition)) )
    //   {
    tcurrent = temperature;
    pcurrent = pressure;
    xcurrent = Massfraction2XNaCl(composition);
    GetTemperatureIndex(tcurrent);
    GetPressureIndex(pcurrent);
    GetCompositionIndex(xcurrent);
    if(pcurrent <= 1.01325e5)
      {
        pcurrent = 1.01325e5;
      }
	
    storage_index = (it-it_offset[table_id]) * x_dim_table[table_id] * p_dim + x_dim_table[table_id]*ip + ix;
    // storage_index here refers always to the it-coordinate of A in the sketch below
    // axes are t,p,x
    //
    //     p
    //     ^
    //     |/    |/
    //    -C-----G-
    //    /|    /|  t
    //  |/ |  |/ | ^
    // -D-----H- |/
    // /| -B-/|--F-
    //  | /|  | /|
    //  |/    |/
    // -A-----E-> x
    // /|    /|
    //
    // coordinates tA, tB, ..., pA, pB, as well as values vA, vB, ... at A,B,... have been determined 
    // prior to this step; t, p, x  and vpoint are the coordinates and value at the point to be interpolated
    // t = temperature coordinate, p = pressure coordinate
	
    // *** begin checks for values on outer table surfaces ***
    // on backside
    if(it == it_max_table[table_id])
      {
        //d	    cout << tcurrent << ": data sit on backside of table " << table_id << "...\n";
        if(ip == p_dim-1)
          {
            if(ix == x_dim_table[table_id]-1)
              {
                // conditions are tmax_table, 5000bar, xmax_table
                return storage_vector[tablesize[table_id]*property_index+tablesize[table_id]];
              }
            else
              { // tmax_table, 5000 bar, X<xmax_table
                iA    = tablesize[table_id]*property_index + storage_index;
                iE    = iA+1;
                iA_x  = tablesize[table_id]*composition_index + storage_index;
                x_iA  = storage_vector[iA_x];
                x_iE  = storage_vector[iA_x+1];
                xnorm = (xcurrent-x_iA)/(x_iE-x_iA);
                return storage_vector[iA]+xnorm*(storage_vector[iE]-storage_vector[iA]);
              }
          } // tmax_table, <5000bar, X<xmax_table
        iA      = tablesize[table_id]*property_index + storage_index;
        iD      = iA + x_dim;
        iE      = iA+1;
        iH      = iD +1;
        iA_x    = tablesize[table_id]*composition_index + storage_index;
        iE_x    = iA_x+1;
        iA_p    = tablesize[table_id]*pressure_index + storage_index;
        iD_p    = tablesize[table_id]*pressure_index + storage_index + x_dim;
        p_iA    = storage_vector[iA_p];
        p_iD    = storage_vector[iD_p];
        x_iA    = storage_vector[iA_x];
        x_iE    = storage_vector[iE_x];
        pnorm   = (pcurrent-p_iA)/(p_iD-p_iA);
        xnorm   = (xcurrent-x_iA)/(x_iE-x_iA);
        v_left  = storage_vector[iA]+pnorm*(storage_vector[iD]-storage_vector[iA]);
        v_right = storage_vector[iE]+pnorm*(storage_vector[iH]-storage_vector[iE]);
        v_point = v_left + xnorm*(v_right-v_left);
        return    v_point;
      }
	
    // on top side
    if(ip == p_dim-1)
      { // the case with both p_dim-1 and tmax_table has already been treated ...
        //d	    cout << "data sit on top side of table ...\n";
        if(ix == x_dim_table[table_id]-1)
          {
            // conditions are 5000bar, X=xmax_table at some t
            iA      = tablesize[table_id]*property_index + storage_index;
            iB      = iA + p_dim*x_dim_table[table_id];
            iB_t    = iA_t  + p_dim*x_dim_table[table_id];
            iA_t    = tablesize[table_id]*temperature_index + storage_index;
            tnorm   = (tcurrent-t_iA)/(t_iB-t_iA);
            return     storage_vector[iA]+tnorm*(storage_vector[iB]-storage_vector[iA]);
          }
        else
          { // 5000bar, X<xmax_table
            iA      = tablesize[table_id]*property_index + storage_index;
            iE      = iA + x_dim_table[table_id];
            iB      = iA + p_dim*x_dim_table[table_id];
            iB_t    = iA_t  + p_dim*x_dim_table[table_id];
            iF      = iB + 1;
            iA_x    = tablesize[table_id]*composition_index + storage_index;
            x_iA    = storage_vector[iA_x];
            x_iE    = storage_vector[iA_x+1];
            xnorm   = (xcurrent-x_iA)/(x_iE-x_iA);
            iA_t    = tablesize[table_id]*temperature_index + storage_index;
            tnorm   = (tcurrent-t_iA)/(t_iB-t_iA);
            v_left  = storage_vector[iA]+tnorm*(storage_vector[iB]-storage_vector[iA]);
            v_right = storage_vector[iE]+tnorm*(storage_vector[iF]-storage_vector[iE]);
            v_point =  v_left + xnorm*(v_right-v_left);
            return     v_point;
          }
      }
	
    // on right side
    if(ix == x_dim_table[table_id]-1)
      {
        //d	    cout << "data sit on right side of table ...\n";
        //d     cout << "tcurrent, pcurrent, xcurrent : " << tcurrent << ", " << pcurrent << ", " << xcurrent << endl;
        iA    = tablesize[table_id]*property_index + storage_index;
        iA_t  = tablesize[table_id]*temperature_index + storage_index;
        iB    = iA    + p_dim*x_dim_table[table_id];
        iB_t  = iA_t  + p_dim*x_dim_table[table_id];
        t_iA  = storage_vector[iA_t];
        t_iB  = storage_vector[iB_t];
        tnorm = (tcurrent-t_iA)/(t_iB-t_iA);
        return storage_vector[iA]+tnorm*(storage_vector[iB]-storage_vector[iA]); 
      }
    // *** end checks for values on outer table surfaces ***
	
    // Now determine all parameters that are independet of property_index
    iA_t = tablesize[table_id]*temperature_index + storage_index;
    iB_t = iA_t  + p_dim*x_dim_table[table_id];
    iA_p = tablesize[table_id]*pressure_index    + storage_index;
    iD_p = tablesize[table_id]*pressure_index    + storage_index + x_dim_table[table_id];
    iA_x = tablesize[table_id]*composition_index + storage_index;
    iE_x = iA_x+1;
	
    // iA = tablesize[table_id]*property_index + storage_index;
    // iE = iA + 1;
    // iD = iA + x_dim_table[table_id];
    // iH = iD + 1;
	
    // iB = iA + p_dim*x_dim_table[table_id];
    // iC = iB + x_dim_table[table_id];
    // iF = iB + 1;
    // iG = iC + 1;
	
    t_iA = storage_vector[iA_t];
    t_iB = storage_vector[iB_t];
    p_iA = storage_vector[iA_p];
    p_iD = storage_vector[iD_p];
    x_iA = storage_vector[iA_x];
    x_iE = storage_vector[iE_x];
	
    //d     cout << "ip is found as " << ip << endl;
    //d     cout << "it is found as " << it << endl;
    //d     cout << "ix is found as " << ix << endl;
    //d     cout << "iA...iH     as " << iA << ", "  << iB << ", "  << iC << ", "  << iD << ", "  << iE << ", "  << iF << ", "  << iG << ", "  << iH << endl;
    //d     cout << "t_iA = " << t_iA << endl; 
    //d     cout << "t_iB = " << t_iB << endl; 
    //d     cout << "p_iA = " << p_iA << endl; 
    //d     cout << "p_iD = " << p_iD << endl; 
    //d     cout << "x_iA = " << x_iA << endl; 
    //d     cout << "x_iE = " << x_iE << endl; 
    //d     cout << "storage_vector[iA] = " << storage_vector[iA] << endl;
    //d     cout << "storage_vector[iB] = " << storage_vector[iB] << endl;
    //d     cout << "storage_vector[iC] = " << storage_vector[iC] << endl;
    //d     cout << "storage_vector[iD] = " << storage_vector[iD] << endl;
    //d     cout << "storage_vector[iE] = " << storage_vector[iE] << endl;
    //d     cout << "storage_vector[iF] = " << storage_vector[iF] << endl;
    //d     cout << "storage_vector[iG] = " << storage_vector[iG] << endl;
    //d     cout << "storage_vector[iH] = " << storage_vector[iH] << endl << endl;
	
    tnorm    = (tcurrent - t_iA ) / ( t_iB - t_iA );	    
    pnorm    = (pcurrent - p_iA ) / ( p_iD - p_iA );
    xnorm    = (xcurrent - x_iA ) / ( x_iE - x_iA );
	
    //d     if(property_index == viscosity_index)
    //d  {
    //d 	cout << "Viscosity interpolation iA-iH: \n";
    //d 	cout << storage_vector[iA] << endl;
    //d 	cout << storage_vector[iB] << endl;
    //d 	cout << storage_vector[iC] << endl;
    //d 	cout << storage_vector[iD] << endl;
    //d 	cout << storage_vector[iE] << endl;
    //d 	cout << storage_vector[iF] << endl;
    //d 	cout << storage_vector[iG] << endl;
    //d 	cout << storage_vector[iH] << endl;
    //d     
    //d  }
    //      }

    // evtl noch einen check property_index vs "last_property_index" einfuehren?

    iA = tablesize[table_id]*property_index + storage_index;
    iE = iA + 1;
    iD = iA + x_dim_table[table_id];
    iH = iD + 1;
    
    iB = iA + p_dim*x_dim_table[table_id];
    iC = iB + x_dim_table[table_id];
    iF = iB + 1;
    iG = iC + 1;
    
    // Interpolation left side in t,p
    v_bottom =  storage_vector[iA] + tnorm*(storage_vector[iB] - storage_vector[iA]);
    v_top    =  storage_vector[iD] + tnorm*(storage_vector[iC] - storage_vector[iD]);
    v_left   =  v_bottom + pnorm*(v_top-v_bottom);
    
    // Interpolation right side in t,p
    v_bottom =  storage_vector[iE] + tnorm*(storage_vector[iF] - storage_vector[iE]);
    v_top    =  storage_vector[iH] + tnorm*(storage_vector[iG] - storage_vector[iH] );
    v_right  =  v_bottom + pnorm*(v_top-v_bottom);
    v_point  =  v_left + xnorm*(v_right-v_left);

    return  v_point;
          
  }




  //   //d cout << "H2ONaClLookup::TrilinearInterpolation is being carried out for table_id = "
  //     //d		       << table_id << " and property_index = " << property_index << endl;
  //   tcurrent = temperature;
  //   pcurrent = pressure;
  //   xcurrent = composition;
  //   GetTemperatureIndex(tcurrent);
  //   GetPressureIndex(pcurrent);
  //   GetCompositionIndex(xcurrent);
  //   if(pcurrent <= 1.01325e5)
  //     {
  // 	pcurrent = 1.01325e5;
  //     }
    
	
  //   // ********************* Implement general range check? ***********************************
	
  //   // storage_index here refers always to the it-coordinate of A in the sketch below
  //   storage_index = (it-it_offset[table_id]) * x_dim_table[table_id] * p_dim + x_dim_table[table_id]*ip + ix;
	
  //   // *** begin checks for values on outer table surfaces ***
  //   // on backside
  //   if(it == it_max_table[table_id])
  //     {
  // 	//d	    cout << tcurrent << ": data sit on backside of table " << table_id << "...\n";
  // 	if(ip == p_dim-1)
  // 	  {
  // 	  if(ix == x_dim_table[table_id]-1)
  // 	    {
  // 	      // conditions are tmax_table, 5000bar, xmax_table
  // 	      return storage_vector[tablesize[table_id]*property_index+tablesize[table_id]];
  // 	    }
  // 	  else
  // 	    { // tmax_table, 5000 bar, X<xmax_table
  // 	      iA    = tablesize[table_id]*property_index + storage_index;
  // 	      iE    = iA+1;
  // 	      iA_x  = tablesize[table_id]*composition_index + storage_index;
  // 	      x_iA  = storage_vector[iA_x];
  // 	      x_iE  = storage_vector[iA_x+1];
  // 	      xnorm = (xcurrent-x_iA)/(x_iE-x_iA);
  // 	      return storage_vector[iA]+xnorm*(storage_vector[iE]-storage_vector[iA]);
  // 	    }
  // 	  } // tmax_table, <5000bar, X<xmax_table
  // 	iA      = tablesize[table_id]*property_index + storage_index;
  // 	iD      = iA + x_dim;
  // 	iE      = iA+1;
  // 	iH      = iD +1;
  // 	iA_x    = tablesize[table_id]*composition_index + storage_index;
  // 	iE_x    = iA_x+1;
  // 	iA_p    = tablesize[table_id]*pressure_index + storage_index;
  // 	iD_p    = tablesize[table_id]*pressure_index + storage_index + x_dim;
  // 	p_iA    = storage_vector[iA_p];
  // 	p_iD    = storage_vector[iD_p];
  // 	x_iA    = storage_vector[iA_x];
  // 	x_iE    = storage_vector[iE_x];
  // 	pnorm   = (pcurrent-p_iA)/(p_iD-p_iA);
  // 	xnorm   = (xcurrent-x_iA)/(x_iE-x_iA);
  // 	v_left  = storage_vector[iA]+pnorm*(storage_vector[iD]-storage_vector[iA]);
  // 	v_right = storage_vector[iE]+pnorm*(storage_vector[iH]-storage_vector[iE]);
  // 	v_point = v_left + xnorm*(v_right-v_left);
  // 	return    v_point;
  //     }
    
  //   // on top side
  //   if(ip == p_dim-1)
  //     { // the case with both p_dim-1 and tmax_table has already been treated ...
  // 	//d	    cout << "data sit on top side of table ...\n";
  // 	if(ix == x_dim_table[table_id]-1)
  // 	  {
  // 	    // conditions are 5000bar, X=xmax_table at some t
  // 	    iA      = tablesize[table_id]*property_index + storage_index;
  // 	    iB      = iA + p_dim*x_dim_table[table_id];
  // 	    iB_t    = iA_t  + p_dim*x_dim_table[table_id];
  // 	    iA_t    = tablesize[table_id]*temperature_index + storage_index;
  // 	    tnorm   = (tcurrent-t_iA)/(t_iB-t_iA);
  // 	    return     storage_vector[iA]+tnorm*(storage_vector[iB]-storage_vector[iA]);
  // 	  }
  // 	else
  // 	  { // 5000bar, X<xmax_table
  // 	    iA      = tablesize[table_id]*property_index + storage_index;
  // 	    iE      = iA + x_dim_table[table_id];
  // 	    iB      = iA + p_dim*x_dim_table[table_id];
  // 	    iB_t    = iA_t  + p_dim*x_dim_table[table_id];
  // 	    iF      = iB + 1;
  // 	    iA_x    = tablesize[table_id]*composition_index + storage_index;
  // 	    x_iA    = storage_vector[iA_x];
  // 	    x_iE    = storage_vector[iA_x+1];
  // 	    xnorm   = (xcurrent-x_iA)/(x_iE-x_iA);
  // 	    iA_t    = tablesize[table_id]*temperature_index + storage_index;
  // 	    tnorm   = (tcurrent-t_iA)/(t_iB-t_iA);
  // 	    v_left  = storage_vector[iA]+tnorm*(storage_vector[iB]-storage_vector[iA]);
  // 	    v_right = storage_vector[iE]+tnorm*(storage_vector[iF]-storage_vector[iE]);
  // 	    v_point =  v_left + xnorm*(v_right-v_left);
  // 	    return     v_point;
  // 	  }
  //     }

  //   // on right side
  //   if(ix == x_dim_table[table_id]-1)
  //     {
  // 	//d	    cout << "data sit on right side of table ...\n";
  // 	cout << "tcurrent, pcurrent, xcurrent : " << tcurrent << ", " << pcurrent << ", " << xcurrent << endl;
  // 	iA    = tablesize[table_id]*property_index + storage_index;
  // 	iA_t  = tablesize[table_id]*temperature_index + storage_index;
  // 	iB    = iA    + p_dim*x_dim_table[table_id];
  // 	iB_t  = iA_t  + p_dim*x_dim_table[table_id];
  // 	t_iA  = storage_vector[iA_t];
  // 	t_iB  = storage_vector[iB_t];
  // 	tnorm = (tcurrent-t_iA)/(t_iB-t_iA);
  // 	return storage_vector[iA]+tnorm*(storage_vector[iB]-storage_vector[iA]); 
  //     }
    
  //   // *** end checks for values on outer table surfaces ***
    
    
    
  //   // axes are t,p,x
  //   //
  //   //     p
  //   //     ^
  //   //     |/    |/
  //   //    -C-----G-
  //   //    /|    /|  t
  //   //  |/ |  |/ | ^
  //   // -D-----H- |/
  //   // /| -B-/|--F-
  //   //  | /|  | /|
  //   //  |/    |/
  //   // -A-----E-> x
  //   // /|    /|
  //   //
  //   // coordinates tA, tB, ..., pA, pB, as well as values vA, vB, ... at A,B,... have been determined 
  //   // prior to this step; t, p, x  and vpoint are the coordinates and value at the point to be interpolated
  //   // t = temperature coordinate, p = pressure coordinate
    
    
  //   iA_t = tablesize[table_id]*temperature_index + storage_index;
  //   iA_p = tablesize[table_id]*pressure_index    + storage_index;
  //   iD_p = tablesize[table_id]*pressure_index    + storage_index + x_dim_table[table_id];
  //   iA_x = tablesize[table_id]*composition_index + storage_index;
  //   iE_x = iA_x+1;
    
  //   iA = tablesize[table_id]*property_index + storage_index;
  //   iE = iA + 1;
  //   iD = iA + x_dim_table[table_id];
  //   iH = iD + 1;
    
  //   iB = iA + p_dim*x_dim_table[table_id];
  //   iC = iB + x_dim_table[table_id];
  //   iF = iB + 1;
  //   iG = iC + 1;
  //   iB_t = iA_t  + p_dim*x_dim_table[table_id];

  //   t_iA = storage_vector[iA_t];
  //   t_iB = storage_vector[iB_t];
  //   p_iA = storage_vector[iA_p];
  //   p_iD = storage_vector[iD_p];
  //   x_iA = storage_vector[iA_x];
  //   x_iE = storage_vector[iE_x];

  //   //d     cout << "ip is found as " << ip << endl;
  //   //d     cout << "it is found as " << it << endl;
  //   //d     cout << "ix is found as " << ix << endl;
  //   //d     cout << "iA...iH     as " << iA << ", "  << iB << ", "  << iC << ", "  << iD << ", "  << iE << ", "  << iF << ", "  << iG << ", "  << iH << endl;
  //   //d     cout << "t_iA = " << t_iA << endl; 
  //   //d     cout << "t_iB = " << t_iB << endl; 
  //   //d     cout << "p_iA = " << p_iA << endl; 
  //   //d     cout << "p_iD = " << p_iD << endl; 
  //   //d     cout << "x_iA = " << x_iA << endl; 
  //   //d     cout << "x_iE = " << x_iE << endl; 
  //   //d     cout << "storage_vector[iA] = " << storage_vector[iA] << endl;
  //   //d     cout << "storage_vector[iB] = " << storage_vector[iB] << endl;
  //   //d     cout << "storage_vector[iC] = " << storage_vector[iC] << endl;
  //   //d     cout << "storage_vector[iD] = " << storage_vector[iD] << endl;
  //   //d     cout << "storage_vector[iE] = " << storage_vector[iE] << endl;
  //   //d     cout << "storage_vector[iF] = " << storage_vector[iF] << endl;
  //   //d     cout << "storage_vector[iG] = " << storage_vector[iG] << endl;
  //   //d     cout << "storage_vector[iH] = " << storage_vector[iH] << endl << endl;

  //   // Interpolation left side in t,p
  //   tnorm    = (tcurrent - t_iA ) / ( t_iB - t_iA );	    

  //   pnorm    = (pcurrent - p_iA ) / ( p_iD - p_iA );
  //   v_bottom =  storage_vector[iA] + tnorm*(storage_vector[iB] - storage_vector[iA]);
  //   v_top    =  storage_vector[iD] + tnorm*(storage_vector[iC] - storage_vector[iD]);
  //   v_left   =  v_bottom + pnorm*(v_top-v_bottom);

  //   // Interpolation right side in t,p
  //   v_bottom =  storage_vector[iE] + tnorm*(storage_vector[iF] - storage_vector[iE]);
  //   v_top    =  storage_vector[iH] + tnorm*(storage_vector[iG] - storage_vector[iH] );
  //   v_right  =  v_bottom + pnorm*(v_top-v_bottom);
    
  //   // Interpolation between the two sides along x
  //   xnorm   = (xcurrent - x_iA ) / ( x_iE - x_iA );
    
  //   v_point =  v_left + xnorm*(v_right-v_left);
    
  //   //d     if(property_index == viscosity_index)
  //   //d  {
  // 	//d 	cout << "Viscosity interpolation iA-iH: \n";
  // 	//d 	cout << storage_vector[iA] << endl;
  // 	//d 	cout << storage_vector[iB] << endl;
  // 	//d 	cout << storage_vector[iC] << endl;
  // 	//d 	cout << storage_vector[iD] << endl;
  // 	//d 	cout << storage_vector[iE] << endl;
  // 	//d 	cout << storage_vector[iF] << endl;
  // 	//d 	cout << storage_vector[iG] << endl;
  // 	//d 	cout << storage_vector[iH] << endl;
  // 	//d     
  //   //d  }
    
  //   return  v_point;
  // }



  void H2ONaClLookup::BuildSimpleTable( const long& i, std::vector<double64>& table )
  {
    table_id   = i;
    
    FILE* infile;
    infile = fopen(filename[table_id],"rb");
	
    if(infile == NULL)
      {
        cerr << "H2ONaClLookup : Lookup file \"" << filename[table_id] << "\" missing, computing ...\n\n";
        it = 0;

        double64 xbrine(XNaCl2Massfraction(xcurrent));
	
        H2OLookup                              water;
        Brine                                  brine(tcurrent,pcurrent,xbrine);
        CriticalCurveLookup                    critcurve(tcurrent);
	
        for(tcurrent = t_min_table[table_id]; tcurrent < t_max_table[table_id]+1.0e-3; tcurrent += t_res)
          {
            cerr << "H2ONaClLookup::BuildSimpleTable computing for t = " << tcurrent << endl;
            GetTemperatureIndex(tcurrent+1.0e-3);
            if(it > it_max_table[table_id])
              {
                break;
              }
		
            for(pcurrent = 0.5e5; pcurrent <= 5000.1e5; pcurrent += p_res)
              {
                GetPressureIndex(pcurrent+1.0e-3);
                //		    cout << "H2ONaClLookup::BuildSimpleTable" << filename[table_id] << " tcurrent, pcurrent = " << tcurrent << "," << pcurrent << endl;
                if(pcurrent > 5000.e5+0.001)
                  {
                    break;
                  }
                if(pcurrent > 5000.0e5) pcurrent = 5000.0e5;
                for(xcurrent = 0.0e0; xcurrent <= 1.00001e0; xcurrent += x_res)
                  {
                    xbrine = XNaCl2Massfraction(xcurrent);
                    GetCompositionIndex(xcurrent + 1.0e-6);
                    if(ix > x_dim_table[table_id]-1)
                      {
                        break;
                      }
			
                    if(pcurrent < critcurve.Pressure() && xbrine < critcurve.MassFractionNaCl()) state = V;
                    else state = F; // sufficient?
                    storage_index = (it-it_offset[table_id]) * x_dim_table[table_id]*p_dim + x_dim_table[table_id]*ip + ix;
                    table[tablesize[table_id]*temperature_index     + storage_index ] = tcurrent;
                    table[tablesize[table_id]*pressure_index        + storage_index ] = pcurrent;
                    table[tablesize[table_id]*composition_index     + storage_index ] = xcurrent; //xbrine?
                    if(ix != 0)
                      {
                        table[tablesize[table_id]*density_index         + storage_index ] = brine.Density();
                        table[tablesize[table_id]*enthalpy_index        + storage_index ] = brine.Enthalpy();
                        table[tablesize[table_id]*heatcapacity_index    + storage_index ] = brine.HeatCapacity();
                        table[tablesize[table_id]*compressibility_index + storage_index ] = brine.Compressibility();
                        table[tablesize[table_id]*viscosity_index       + storage_index ] = brine.Viscosity();
                      }
                    else
                      {
                        table[tablesize[table_id]*density_index         + storage_index ] = water.Density(tcurrent, pcurrent);
                        table[tablesize[table_id]*enthalpy_index        + storage_index ] = water.Enthalpy(tcurrent, pcurrent);
                        table[tablesize[table_id]*heatcapacity_index    + storage_index ] = water.HeatCapacity(tcurrent, pcurrent);
                        table[tablesize[table_id]*compressibility_index + storage_index ] = water.Compressibility(tcurrent, pcurrent);
                        table[tablesize[table_id]*viscosity_index       + storage_index ] = water.Viscosity(tcurrent, pcurrent);
                      }
                  } // x-loop
              } // p-loop
          } // t_loop
        FILE* outfile;
        outfile = fopen( filename[table_id], "wb");
        if( outfile == NULL )
          {
            cerr << "could not even find it for writing, please stop program and debug !!!!\n";
            char yesno;
            cerr << "or enter any key to continue (simulation likely to crash or give wrong results!) :";
            cin  >> yesno;
          }
        else
          {
            cerr << "writing file " << filename[table_id] << " ... ";
            skm_C_fwrite( outfile, table );
            fclose(outfile);
            cerr << "done!\n";
          }
      }
    else
      {
        cerr << "reading file " << filename[table_id] << " ... ";
        skm_C_fread( infile, table );
        fclose( infile );
        cerr << "done!\n";
      }
    cerr << "H2ONaClLookup, leaving BuildSimpleTable" << table_id << " ...\n\n";
  }



  void H2ONaClLookup::BuildTable1()
  {
    
    table_id   = 1;
    state      = L;
    
    FILE* infile1;
    infile1=fopen( filename[1], "rb");
    
    if ( infile1 == NULL ) 
      {
        cerr << "H2ONaClLookup : Lookup file " << filename[1] << " missing, computing ...\n\n";
        it = 0;

        double64                               xbrine(XNaCl2Massfraction(xcurrent));

        H2OLookup                              water;
        CriticalPointH2O                       cp_h2o;
        //	    WaterBoilingCurveLiquidLookup          boiling_l(tcurrent);
        Brine                                  brine(tcurrent,pcurrent,xbrine);
        HaliteLiquidusLookup                   liquidus(tcurrent,pcurrent);
        VLH_LiquidLookup                       vlh_l(tcurrent);
        TwophaseLiquidLookup                   twophase_l(tcurrent,pcurrent);
	
        for(tcurrent = 0.0e0; tcurrent < 200.01e0; tcurrent += t_res)
          {
            //		cout << "H2ONaClLookup::BuildTable1() : before GetTemperatureIndex 1 at tcurrent = " << tcurrent << endl;
            GetTemperatureIndex(tcurrent+1.0e-3);
            if(it > it_max_table[1])
              {
                //		    cout << "reached temperature limit at t = " << tcurrent << " and it = " << it << endl;
                break;
              }
            cerr << "H2ONaClLookup::BuildTable1() : t = " << tcurrent << endl;
		
            if(tcurrent < cp_h2o.Temperature())
              {
                //		    cout << "H2ONaClLookup::BuildTable1() : Point 1\n";
		
                // ********** search routine for VL boundary **********
                // ip_vl[ix] is the ip of the last point below the VL surface before going into the L field
                // for x > x_liquidus this applies to the last point in VH before entering LH
		
                for(xcurrent = 0.0e0; xcurrent <= 1.00001e0; xcurrent += x_res)
                  {
                    xbrine = XNaCl2Massfraction(xcurrent);
                    //		    cout << "H2ONaClLookup::BuildTable1() : Point 2\n";
                    GetCompositionIndex(xcurrent + 1.0e-6);
                    if(ix > x_dim_table[1]-1)
                      {
                        break;
                      }
                    else if(xbrine > liquidus.MassFractionNaCl())
                      { 
                        GetPressureIndex(vlh_l.Pressure());
                        ip_vl[ix] = ip;
                        GetPressureIndex(pcurrent);
                      }
                    else
                      {
                        ip_vl[ix] = 0;
                        for(pcurrent = 0.5e5; pcurrent <= 5000.e5+0.1; pcurrent += p_res)
                          {
                            GetPressureIndex(pcurrent+1.0e-3);
                            if(ip > p_dim-1)
                              {
                                break;
                              }
                            if(pcurrent < water.SaturationPressureFromT(tcurrent))
                              {
                                if(xbrine < twophase_l.MassFractionNaCl())
                                  {
                                    ip_vl[ix] = ip;
                                  }
                              } 
                            else break;
                          }
                      }
                  }
		
                // here's the problem - determination of p_res upon downward looping
                for(pcurrent = 5000.0e5; pcurrent >= 0.5e5; pcurrent -= p_res)
                  {
                    //			cout << "H2ONaClLookup::BuildTable1() : Point 3\n";
                    if(pcurrent < 1.0e5) break;
		    
                    GetPressureIndex(pcurrent-1.0e-3);
                    ip += 1; // because determined it at too low P, but should now have correct p_res;
                    if(ip > p_dim-1)
                      {
                        //			    cout << "H2ONaClLookup::BuildTable1() : Point 3a\n";
                        break;
                      }
                    if(pcurrent < 1.0e5) break;
		    
                    for(xcurrent = 0.0e0; xcurrent <= 1.00001e0; xcurrent += x_res)
                      {
                        xbrine = XNaCl2Massfraction(xcurrent);
                        //			    cout << "H2ONaClLookup::BuildTable1() : Point 3b\n";
                        GetCompositionIndex(xcurrent + 1.0e-6);
                        if(ix > x_dim_table[1]-1)
                          {
                            break;
                          }
			
                        // ********** p-storage is wrong !!! **********
                        // much better now but still not perfect
                        //			    cout << "H2ONaClLookup::BuildTable1() : Point 3c\n";
                        storage_index = (it-it_offset[1])  *x_dim_table[1]*p_dim + x_dim_table[1]*ip + ix;
                        table1[tablesize[1]*temperature_index     + storage_index ] = tcurrent;
                        table1[tablesize[1]*pressure_index        + storage_index ] = pcurrent;
                        table1[tablesize[1]*composition_index     + storage_index ] = xcurrent;
			
                        if(ip > ip_vl[ix])
                          {
                            if( ix != 0)
                              {
                                table1[tablesize[1]*density_index         + storage_index ] = brine.Density();
                                table1[tablesize[1]*enthalpy_index        + storage_index ] = brine.Enthalpy();
                                table1[tablesize[1]*heatcapacity_index    + storage_index ] = brine.HeatCapacity();
                                table1[tablesize[1]*compressibility_index + storage_index ] = brine.Compressibility();
                                table1[tablesize[1]*viscosity_index       + storage_index ] = brine.Viscosity();
                              }
                            else
                              {
                                table1[tablesize[1]*density_index         + storage_index ] = water.Density(tcurrent, pcurrent);
                                table1[tablesize[1]*enthalpy_index        + storage_index ] = water.Enthalpy(tcurrent, pcurrent);
                                table1[tablesize[1]*heatcapacity_index    + storage_index ] = water.HeatCapacity(tcurrent, pcurrent);
                                table1[tablesize[1]*compressibility_index + storage_index ] = water.Compressibility(tcurrent, pcurrent);
                                table1[tablesize[1]*viscosity_index       + storage_index ] = water.Viscosity(tcurrent, pcurrent);
                              }
                          }			    
                        else
                          {
                            storageindex_highp1 = (it-it_offset[1]) * x_dim_table[1]*p_dim + x_dim_table[1]*(ip_vl[ix]+1) + ix;
                            storageindex_highp2 = (it-it_offset[1]) * x_dim_table[1]*p_dim + x_dim_table[1]*(ip_vl[ix]+2) + ix;
                            DP_extrapol         =  table1[tablesize[1]*pressure_index + storageindex_highp1];
                            DP_extrapol        -=  table1[tablesize[1]*pressure_index + storage_index];
                            DP_highp            =  table1[tablesize[1]*pressure_index + storageindex_highp2];
                            DP_highp           -=  table1[tablesize[1]*pressure_index + storageindex_highp1];
                            DDensityDP          =  table1[tablesize[1]*density_index         + storageindex_highp2];
                            DDensityDP         -=  table1[tablesize[1]*density_index         + storageindex_highp1];
                            DDensityDP         /=  DP_highp;
                            DEnthalpyDP         =  table1[tablesize[1]*enthalpy_index        + storageindex_highp2];
                            DEnthalpyDP        -=  table1[tablesize[1]*enthalpy_index        + storageindex_highp1];
                            DEnthalpyDP        /=  DP_highp;
                            DHeatCapacityDP     =  table1[tablesize[1]*heatcapacity_index    + storageindex_highp2];
                            DHeatCapacityDP    -=  table1[tablesize[1]*heatcapacity_index    + storageindex_highp1];
                            DHeatCapacityDP    /=  DP_highp;
                            DCompressibilityDP  =  table1[tablesize[1]*compressibility_index + storageindex_highp2];
                            DCompressibilityDP -=  table1[tablesize[1]*compressibility_index + storageindex_highp1];
                            DCompressibilityDP /=  DP_highp;
                            DViscosityDP        =  table1[tablesize[1]*viscosity_index       + storageindex_highp2];
                            DViscosityDP       -=  table1[tablesize[1]*viscosity_index       + storageindex_highp1];
                            DViscosityDP       /=  DP_highp;
			    
                            table1[tablesize[1]*density_index             + storage_index ] = 
                              table1[tablesize[1]*density_index         + storageindex_highp1] - DDensityDP * DP_extrapol;
                            table1[tablesize[1]*enthalpy_index            + storage_index ] = 				
                              table1[tablesize[1]*enthalpy_index        + storageindex_highp1] - DEnthalpyDP * DP_extrapol;
                            table1[tablesize[1]*heatcapacity_index        + storage_index ] = 
                              table1[tablesize[1]*heatcapacity_index    + storageindex_highp1] - DHeatCapacityDP * DP_extrapol;
                            table1[tablesize[1]*compressibility_index + storage_index ] = 
                              table1[tablesize[1]*compressibility_index + storageindex_highp1] - DCompressibilityDP * DP_extrapol;
                            table1[tablesize[1]*viscosity_index           + storage_index ] = 
                              table1[tablesize[1]*viscosity_index       + storageindex_highp1] - DViscosityDP * DP_extrapol;
                          }
                      }
                  }
              }
            else
              { //tcurrent >= cp_h2o.Temperature()){
                //p and x-loop!!!!
                pcurrent = 5000.0e5;
                GetPressureIndex(pcurrent+1.0e-3);
                xcurrent = 0.0e0;
                GetCompositionIndex(xcurrent + 1.0e-6);
		
                for(pcurrent = 0.5e5; pcurrent <= 5000.e5+0.1; pcurrent += p_res)
                  {
                    GetPressureIndex(pcurrent+1.0e-3);
                    if(ip > p_dim-1)
                      {
                        break;
                      }
		    
                    for(xcurrent = 0.0e0; xcurrent <= 1.00001e0; xcurrent += x_res)
                      {
                        xbrine = XNaCl2Massfraction(xcurrent);
                        GetCompositionIndex(xcurrent + 1.0e-6);
                        if(ix > x_dim_table[1]-1)
                          {
                            break;
                          }
			
                        storage_index = (it-it_offset[1]) * x_dim_table[1]*p_dim + x_dim_table[1]*ip + ix;
                        table1[tablesize[1]*temperature_index     + storage_index ] = tcurrent;
                        table1[tablesize[1]*pressure_index        + storage_index ] = pcurrent;
                        table1[tablesize[1]*composition_index     + storage_index ] = xcurrent;
                        if(ix != 0)
                          {
                            table1[tablesize[1]*density_index         + storage_index ] = brine.Density();
                            table1[tablesize[1]*enthalpy_index        + storage_index ] = brine.Enthalpy();
                            table1[tablesize[1]*heatcapacity_index    + storage_index ] = brine.HeatCapacity();
                            table1[tablesize[1]*compressibility_index + storage_index ] = brine.Compressibility();
                            table1[tablesize[1]*viscosity_index       + storage_index ] = brine.Viscosity();
                          }
                        else
                          {
                            table1[tablesize[1]*density_index         + storage_index ] = water.Density(tcurrent, pcurrent);
                            table1[tablesize[1]*enthalpy_index        + storage_index ] = water.Enthalpy(tcurrent, pcurrent);
                            table1[tablesize[1]*heatcapacity_index    + storage_index ] = water.HeatCapacity(tcurrent, pcurrent);
                            table1[tablesize[1]*compressibility_index + storage_index ] = water.Compressibility(tcurrent, pcurrent);
                            table1[tablesize[1]*viscosity_index       + storage_index ] = water.Viscosity(tcurrent, pcurrent);
                          }
                      } // x-loop
                  } // p-loop
              }
          }
        FILE* outfile1;
        outfile1 = fopen( filename[1], "wb");
        if( outfile1 == NULL )
          {
            cerr << "could not even it for writing, please stop program and debug !!!!\n";
            char yesno;
            cerr << "or enter any key to continue (simulation likely to crash or give wrong results!) :";
            cin  >> yesno;
          }
        else
          {
            cerr << "writing file " << filename[1] << " ... ";
            skm_C_fwrite( outfile1, table1 );
            fclose(outfile1);
            cerr << "done!\n";
          }
      }
    else
      {
        cerr << "reading file " << filename[1] << " ... ";
        skm_C_fread( infile1, table1 );
        fclose( infile1 );
        cerr << "done!\n";
      }
    cerr << "H2ONaClLookup, leaving BuildTable1 ...\n\n";
  }
    



  void H2ONaClLookup::BuildTable2()
  {
    table_id   = 2;
    state      = L;
    
    FILE* infile2;
    infile2=fopen( filename[2], "rb");
    
    if ( infile2 == NULL ) 
      {
        cerr << "H2ONaClLookup : Lookup file " << filename[2] << " missing, computing ...\n\n";
	
        it = 0;

        double64                               xbrine(0.0);
	
        H2OLookup                              water;
        CriticalPointH2O                       cp_h2o;
        Brine                                  brine(tcurrent,pcurrent,xbrine  );
        HaliteLiquidusLookup                   liquidus(tcurrent,pcurrent);
        HaliteLiquidus                         liquidus_orig(tcurrent,pcurrent);
        VLH_LiquidLookup                       vlh_l(tcurrent);
        TwophaseLiquidLookup                   twophase_l(tcurrent,pcurrent);
	
        for(tcurrent = 200.0e0; tcurrent < 400.01e0; tcurrent += t_res)
          {
            GetTemperatureIndex(tcurrent+1.0e-3);
            if(it > it_max_table[2])
              {
                cerr << "reached temperature limit at t = " << tcurrent << " and it = " << it << endl;
                break;
              }
            cerr << "BuildTable2 : t = " << tcurrent << endl;
	    
            if(tcurrent < cp_h2o.Temperature())
              {
		
                // ********** search routine for VL boundary **********
                // ip_vl[ix] is the ip of the last point below the VL surface before going into the L field
                // for x > x_liquidus this applies to the last point in VH before entering LH
		
                for(xcurrent = 0.0e0; xcurrent < 1.00001e0; xcurrent += x_res)
                  {
                    xbrine = XNaCl2Massfraction(xcurrent);
                    GetCompositionIndex(xcurrent + 1.0e-6);
                    //d			cerr << "doing xcurrent = " << xcurrent << endl;
                    if(ix > x_dim_table[2]-1)
                      {
                        //d			    cerr << "H2ONaClLookup::BuildTable2() : point 1 - pcurrent is " << pcurrent << "\n";
                        break;
                      }
                    else if(xbrine > liquidus.MassFractionNaCl())
                      { 
                        //d			    cerr << "H2ONaClLookup::BuildTable2() : point 2 - pcurrent is " << pcurrent << "\n";
                        GetPressureIndex(vlh_l.Pressure());
                        ip_vl[ix] = ip;
                        GetPressureIndex(pcurrent); // ********** Why this??
                      }
                    else
                      {
                        //d			    cerr << "H2ONaClLookup::BuildTable2() : point 3\n";
                        ip_vl[ix] = 0;
                        for(pcurrent = 0.5e5; pcurrent < 5000.e5+0.1e0; pcurrent += p_res)
                          {
                            //d				cerr << "H2ONaClLookup::BuildTable2() : point 4 - pcurrent is " << pcurrent << "\n";
                            GetPressureIndex(pcurrent+1.0e-3);
                            if(ip > p_dim-1)
                              {
                                break;
                              }
                            if(pcurrent < water.SaturationPressureFromT(tcurrent))
                              {
                                if(xbrine < twophase_l.MassFractionNaCl())
                                  {
                                    ip_vl[ix] = ip;
                                  }
                              } 
                            else break;
                          }
                      }
                  }
                // debug only
                // 		    if(tcurrent > 373.89e0 && tcurrent < 373.91e0){
                // 			int myerror;
                // 			cerr << "for 373 C: " << ip_vl[0] << endl;
                // 			cin >> myerror;
                // 		    }
                // here's the problem - determination of p_res upon downward looping
                for(pcurrent = 5000.0e5; pcurrent >= 0.5e5; pcurrent -= p_res)
                  {
                    //d			cerr << "H2ONaClLookup::BuildTable2() : point 5 - pcurrent is " << pcurrent << "\n";
                    if(pcurrent < 1.0e5) break;
		    
                    GetPressureIndex(pcurrent-1.0e-3);
                    ip += 1; // because determined it at too low P, but should now have correct p_res;
                    if(ip > p_dim-1)
                      {
                        break;
                      }
                    if(pcurrent < 1.0e5) break;
		    
                    for(xcurrent = 0.0e0; xcurrent <= 1.00001e0; xcurrent += x_res)
                      {
                        xbrine = XNaCl2Massfraction(xcurrent);
                        GetCompositionIndex(xcurrent + 1.0e-6);
                        if(ix > x_dim_table[2]-1)
                          {
                            break;
                          }
			
                        // ********** p-storage is wrong !!! **********
                        // much better now but still not perfect
                        storage_index = (it-it_offset[2])  *x_dim_table[2]*p_dim + x_dim_table[2]*ip + ix;
                        table2[tablesize[2]*temperature_index     + storage_index ] = tcurrent;
                        table2[tablesize[2]*pressure_index        + storage_index ] = pcurrent;
                        table2[tablesize[2]*composition_index     + storage_index ] = xcurrent;
			
                        if(ip > ip_vl[ix])
                          {
                            if( ix != 0)
                              {
                                table2[tablesize[2]*density_index         + storage_index ] = brine.Density();
                                table2[tablesize[2]*enthalpy_index        + storage_index ] = brine.Enthalpy();
                                table2[tablesize[2]*heatcapacity_index    + storage_index ] = brine.HeatCapacity();
                                if(brine.HeatCapacity() < 0.0)
                                  {
                                    cerr << "H2ONaClLookup::BuildTable2() : encountered negative liquid heat capacity for\n";
                                    cerr << "tcurrent = " << tcurrent << endl;
                                    cerr << "pcurrent = " << pcurrent << endl;
                                    cerr << "xcurrent = " << xcurrent << endl;
                                    // char yesno;
                                    // cerr << "enter any char to continue : "; cin >> yesno;
                                  }
                                table2[tablesize[2]*compressibility_index + storage_index ] = brine.Compressibility();
                                //d if( (tcurrent > 300.0 && 
                                //d      xcurrent < liquidus_orig.MassFractionNaCl() &&
                                //d	 table2[tablesize[2]*compressibility_index + storage_index ] < 0.0))
                                //d  {
                                //d int myerror;
                                //d    cerr << "H2ONaClLookup::BuildTable2(): Error is in brine at point 1a, at conditions"
                                //d	 << "\ntcurrent = " << tcurrent
                                //d	 << "\npcurrent = " << pcurrent
                                //d	 << "\nxcurrent = " << xcurrent
                                //d	 << "\ncompress = " << brine.Compressibility() << endl;
                                //d cin >> myerror;
                                //d }
                                table2[tablesize[2]*viscosity_index       + storage_index ] = brine.Viscosity();
                              }
                            else
                              {
                                table2[tablesize[2]*density_index         + storage_index ] = water.Density(tcurrent, pcurrent);
                                table2[tablesize[2]*enthalpy_index        + storage_index ] = water.Enthalpy(tcurrent, pcurrent);
                                table2[tablesize[2]*heatcapacity_index    + storage_index ] = water.HeatCapacity(tcurrent, pcurrent);
                                if(water.HeatCapacity(tcurrent, pcurrent) < 0.0)
                                  {
                                    cerr << "H2ONaClLookup::BuildTable2() : encountered negative liquid heat capacity water for\n";
                                    cerr << "tcurrent = " << tcurrent << endl;
                                    cerr << "pcurrent = " << pcurrent << endl;
                                    cerr << "xcurrent = " << xcurrent << endl;
                                    // char yesno;
                                    // cerr << "enter any char to continue : "; cin >> yesno;
                                  }
                                table2[tablesize[2]*compressibility_index + storage_index ] = water.Compressibility(tcurrent, pcurrent);
                                //d if( (tcurrent > 300.0 && 
                                //d	       xcurrent < liquidus_orig.MassFractionNaCl() &&
                                //d	       table2[tablesize[2]*compressibility_index + storage_index ] < 0.0))
                                //d  {
                                //d int myerror;
                                //d    cerr << "H2ONaClLookup::BuildTable2(): Error is in brine at point 1b, at conditions"
                                //d	 << "\ntcurrent = " << tcurrent
                                //d	 << "\npcurrent = " << pcurrent
                                //d	 << "\nxcurrent = " << xcurrent
                                //d	 << "\ncompress = " << brine.Compressibility()
                                //d  << "\nbetah2o  = " << water.Compressibility(tcurrent,pcurrent) << endl;
                                //d cin >> myerror;
                                //d }
                                table2[tablesize[2]*viscosity_index       + storage_index ] = water.Viscosity(tcurrent, pcurrent);
                              }
                          }			    
                        else
                          {
                            storageindex_highp1 = (it-it_offset[2]) * x_dim_table[2]*p_dim + x_dim_table[2]*(ip_vl[ix]+1) + ix;
                            storageindex_highp2 = (it-it_offset[2]) * x_dim_table[2]*p_dim + x_dim_table[2]*(ip_vl[ix]+2) + ix;
                            DP_extrapol         =  table2[tablesize[2]*pressure_index + storageindex_highp1];
                            DP_extrapol        -=  table2[tablesize[2]*pressure_index + storage_index];
                            DP_highp            =  table2[tablesize[2]*pressure_index + storageindex_highp2];
                            DP_highp           -=  table2[tablesize[2]*pressure_index + storageindex_highp1];
                            DDensityDP          =  table2[tablesize[2]*density_index         + storageindex_highp2];
                            DDensityDP         -=  table2[tablesize[2]*density_index         + storageindex_highp1];
                            DDensityDP         /=  DP_highp;
                            DEnthalpyDP         =  table2[tablesize[2]*enthalpy_index        + storageindex_highp2];
                            DEnthalpyDP        -=  table2[tablesize[2]*enthalpy_index        + storageindex_highp1];
                            DEnthalpyDP        /=  DP_highp;
                            DHeatCapacityDP     =  table2[tablesize[2]*heatcapacity_index    + storageindex_highp2];
                            DHeatCapacityDP    -=  table2[tablesize[2]*heatcapacity_index    + storageindex_highp1];
                            DHeatCapacityDP    /=  DP_highp;
                            DCompressibilityDP  =  table2[tablesize[2]*compressibility_index + storageindex_highp2];
                            DCompressibilityDP -=  table2[tablesize[2]*compressibility_index + storageindex_highp1];
                            DCompressibilityDP /=  DP_highp;
                            DViscosityDP        =  table2[tablesize[2]*viscosity_index       + storageindex_highp2];
                            DViscosityDP       -=  table2[tablesize[2]*viscosity_index       + storageindex_highp1];
                            DViscosityDP       /=  DP_highp;
			    
                            table2[tablesize[2]*density_index             + storage_index ] = 
                              table2[tablesize[2]*density_index         + storageindex_highp1] - DDensityDP * DP_extrapol;
                            if(table2[tablesize[2]*density_index+ storage_index ] < 322.0)
                              {
                                cerr << "H2ONaClLookup::BuildTable2() : encountered too low liquid density for\n";
                                cerr << "tcurrent = " << tcurrent << endl;
                                cerr << "pcurrent = " << pcurrent << endl;
                                cerr << "xcurrent = " << xcurrent << endl;
                                cerr << table2[tablesize[2]*density_index         + storageindex_highp1] - DDensityDP * DP_extrapol << " changed to value just above boiling surface! : ";
                                // char yesno;
                                // cerr << "enter any char to continue : "; cin >> yesno;
				
                                table2[tablesize[2]*density_index+ storage_index] =
                                  table2[tablesize[2]*density_index    + storageindex_highp1];
                                cerr << table2[tablesize[2]*density_index+ storage_index] << endl;
                              }
                            table2[tablesize[2]*enthalpy_index            + storage_index ] = 				
                              table2[tablesize[2]*enthalpy_index        + storageindex_highp1] - DEnthalpyDP * DP_extrapol;
                            if(table2[tablesize[2]*enthalpy_index+ storage_index ] < 0.0)
                              {
                                cerr << "H2ONaClLookup::BuildTable2() : encountered negative liquid enthalpy for\n";
                                cerr << "tcurrent = " << tcurrent << endl;
                                cerr << "pcurrent = " << pcurrent << endl;
                                cerr << "xcurrent = " << xcurrent << endl;
                                cerr << "changed to value just above boiling surface! : ";
                                table2[tablesize[2]*enthalpy_index+ storage_index] =
                                  table2[tablesize[2]*enthalpy_index    + storageindex_highp1];
                                cerr << table2[tablesize[2]*enthalpy_index+ storage_index] << endl;
                              }
                            table2[tablesize[2]*heatcapacity_index        + storage_index ] = 
                              table2[tablesize[2]*heatcapacity_index    + storageindex_highp1] - DHeatCapacityDP * DP_extrapol;
                            if(table2[tablesize[2]*heatcapacity_index+ storage_index ] < 0.0)
                              {
                                cerr << "H2ONaClLookup::BuildTable2() : encountered negative liquid heat capacity for\n";
                                cerr << "tcurrent = " << tcurrent << endl;
                                cerr << "pcurrent = " << pcurrent << endl;
                                cerr << "xcurrent = " << xcurrent << endl;
                                cerr << "changed to value just above boiling surface! : ";
                                table2[tablesize[2]*heatcapacity_index+ storage_index] =
                                  table2[tablesize[2]*heatcapacity_index    + storageindex_highp1];
                                cerr << table2[tablesize[2]*heatcapacity_index+ storage_index] << endl;
                                // 				    char yesno;
                                // 				    cerr << "enter any char to continue : "; cin >> yesno;
                              }
			    
			    
                            // the following is a work-around to avoid extrapolation artifacts of DCompressibilityDP being very negative
                            // near the critical point; for highly accurate simulations in this region, you would have to set up a 
                            // different table
                            if(tcurrent >= 373.79 && tcurrent <= 374.0 && ix <= 1)
                              table2[tablesize[2]*compressibility_index + storage_index ] = 
                                water.Compressibility(tcurrent,pcurrent);
                            else table2[tablesize[2]*compressibility_index + storage_index ] = 
                                   table2[tablesize[2]*compressibility_index + storageindex_highp1] - DCompressibilityDP * DP_extrapol;
			    
			    
                            //d if( tcurrent > 300.0 &&  
                            //d	      xcurrent < liquidus_orig.MassFractionNaCl() &&
                            //d	      table2[tablesize[2]*compressibility_index + storage_index ] < 0.0)
                            //d {
                            //d int myerror;
                            //d    cerr << "H2ONaClLookup::BuildTable2(): Error is in extrapolation at point 2, at conditions"
                            //d	 << "\ntcurrent = " << tcurrent
                            //d	 << "\npcurrent = " << pcurrent
                            //d	 << "\nxcurrent = " << xcurrent
                            //d        << "\nxsat     = " << liquidus.MassFractionNaCl() 
                            //d	 << "\nbrine = " << brine.Compressibility()
                            //d        << "\nextrapol = " <<  table2[tablesize[2]*compressibility_index + storageindex_highp1] - DCompressibilityDP * DP_extrapol << endl;
                            //d cerr << "table2[tablesize[2]*compressibility_index + storageindex_highp2] = " << table2[tablesize[2]*compressibility_index + storageindex_highp2] << endl;
                            //d cerr << "table2[tablesize[2]*compressibility_index + storageindex_highp1] = " << table2[tablesize[2]*compressibility_index + storageindex_highp1] << endl;
                            //d cerr << "DP_highp = " << DP_highp << endl;
                            //d cerr << "DCompressibilityDP = " << DCompressibilityDP << endl;
                            //dcin >> myerror;
                            //d }
                            table2[tablesize[2]*viscosity_index           + storage_index ] = 
                              table2[tablesize[2]*viscosity_index       + storageindex_highp1] - DViscosityDP * DP_extrapol;
                          }
                      }
                  }
              }
            else
              { //tcurrent >= cp_h2o.Temperature()){
                //p and x-loop!!!!
                state = L; // sufficient?
                pcurrent = 5000.0e5;
                GetPressureIndex(pcurrent+1.0e-3);
                xcurrent = 0.0e0;
                GetCompositionIndex(xcurrent + 1.0e-6);
		
                for(pcurrent = 0.5e5; pcurrent <= 5000.e5+0.1e0; pcurrent += p_res)
                  {
                    GetPressureIndex(pcurrent+1.0e-3);
                    if(ip > p_dim-1)
                      {
                        break;
                      }
		    
                    for(xcurrent = 0.0e0; xcurrent <= 1.00001e0; xcurrent += x_res)
                      {
                        xbrine = XNaCl2Massfraction(xcurrent);
                        GetCompositionIndex(xcurrent + 1.0e-6);
                        if(ix > x_dim_table[2]-1)
                          {
                            break;
                          }
			
                        storage_index = (it-it_offset[2]) * x_dim_table[2]*p_dim + x_dim_table[2]*ip + ix;
                        table2[tablesize[2]*temperature_index     + storage_index ] = tcurrent;
                        table2[tablesize[2]*pressure_index        + storage_index ] = pcurrent;
                        table2[tablesize[2]*composition_index     + storage_index ] = xcurrent;
                        if( ix != 0 )
                          {
                            table2[tablesize[2]*density_index         + storage_index ] = brine.Density();
                            table2[tablesize[2]*enthalpy_index        + storage_index ] = brine.Enthalpy();
                            table2[tablesize[2]*heatcapacity_index    + storage_index ] = brine.HeatCapacity();
                            table2[tablesize[2]*compressibility_index + storage_index ] = brine.Compressibility();
                            // 			    //d if( tcurrent > 300.0 && 
                            // 					  xcurrent < liquidus_orig.MassFractionNaCl() &&
                            // 					  table2[tablesize[2]*compressibility_index + storage_index ] < 0.0){
                            // 				//d int myerror;
                            // 				//d    cerr << "H2ONaClLookup::BuildTable2(): Error is in brine at point 3, at conditions"
                            // 				    //d	 << "\ntcurrent = " << tcurrent
                            // 				    //d	 << "\npcurrent = " << pcurrent
                            // 				    //d	 << "\nxcurrent = " << xcurrent
                            // 				    //d	 << "\ncompress = " << brine.Compressibility() << endl;
                            // 				//d cin >> myerror;}
                            table2[tablesize[2]*viscosity_index       + storage_index ] = brine.Viscosity();
                          }
                        else
                          {
                            table2[tablesize[2]*density_index         + storage_index ] = water.Density(tcurrent, pcurrent);
                            table2[tablesize[2]*enthalpy_index        + storage_index ] = water.Enthalpy(tcurrent, pcurrent);
                            table2[tablesize[2]*heatcapacity_index    + storage_index ] = water.HeatCapacity(tcurrent, pcurrent);
                            table2[tablesize[2]*compressibility_index + storage_index ] = water.Compressibility(tcurrent, pcurrent);
                            // 			    //d if( tcurrent > 300.0 && 
                            // 					  xcurrent < liquidus_orig.MassFractionNaCl() &&
                            // 					  table2[tablesize[2]*compressibility_index + storage_index ] < 0.0){
                            // 				//d int myerror;
                            // 				//d    cerr << "H2ONaClLookup::BuildTable2(): Error is in brine at point 3, at conditions"
                            // 				    //d	 << "\ntcurrent = " << tcurrent
                            // 				    //d	 << "\npcurrent = " << pcurrent
                            // 				    //d	 << "\nxcurrent = " << xcurrent
                            // 				    //d	 << "\ncompress = " << water.Compressibility() << endl;
                            // 				//d cin >> myerror;}
                            table2[tablesize[2]*viscosity_index       + storage_index ] = water.Viscosity(tcurrent, pcurrent);
                          }
                      } // x-loop
                  } // p-loop    
              }
          }
        FILE* outfile2;
        outfile2 = fopen( filename[2], "wb");
        if( outfile2 == NULL )
          {
            cerr << "could not even open file for writing, please stop program and debug !!!!\n";
            char yesno;
            cerr << "or enter any key to continue (simulation likely to crash or give wrong results!) :";
            cin  >> yesno;
          }
        else
          {
            cerr << "writing file " << filename[2] << " ... ";
            skm_C_fwrite( outfile2, table2 );
            fclose(outfile2);
            cerr << "done!\n";
          }
      }
    else
      {
        cerr << "reading file " << filename[2] << " ... ";
        skm_C_fread( infile2, table2 );
        fclose( infile2 );
        cerr << "done!\n";
      }
    cerr << "H2ONaClLookup, leaving BuildTable2 ...\n\n";
  }
  




  void H2ONaClLookup::BuildTable0()
  {
    /* 
       Table0 provides data for the subcritical vapor, i.e., a fluid for the following conditions combined:
       
       - at temperatures below the critical temperature of water
       - at pressures below the boiling pressure of water for that temperature
       - and has a composition smaller than or equal to Twophase Vapor or halite-saturated vapor, respectively
       
       In order to avoid mismatches with adjoining regions - and in particular to avoid that when using the TrilinearInterpolation 
       member one or more of the eight corner points is NOT in the vapor field - data are 
       
       a) extended to 376C by using direct values from Brine 
       b) extrapolated into the L or F fields when encountering the VL boundaries
       
       Maximum XNaCl is 0.0003 and was determined by checking the vapor composition on the 376C V+L isotherm.
       From this, maximum ix is 3.
       
       Implementation: Thomas Driesner, April 25 and 28, 2008
    */
    
    table_id   = 0;
    state      = V;
    
    FILE* infile0;
    infile0=fopen( filename[0], "rb");
    
    if ( infile0 == NULL ) 
      {
        cerr << "H2ONaClLookup : Lookup file " << filename[0] << " missing, computing ...\n\n";
        it = 0;
	
        double64                               xbrine(0.0);

        H2OLookup                              water;
        CriticalPointH2O                       cp_h2o;
        Brine                                  brine(tcurrent,pcurrent,xbrine);
        VLH_VaporLookup                        vlh_v(tcurrent);
        NaClSaturatedVaporLookup               naclsatvap(tcurrent,pcurrent);
        TwophaseVaporLookup                    twophase_v(tcurrent,pcurrent);
        TwophaseLiquidLookup                   twophase_l(tcurrent,pcurrent);
	
        for(tcurrent = 0.0e0; tcurrent < 376.001e0; tcurrent += t_res)
          {
            cerr << "H2ONaClLookup::BuildTable0() : t = " << tcurrent << endl;

            GetTemperatureIndex(tcurrent+1.0e-3);
            if(it > it_max_table[0])
              {
                break;
              }
            // get the data for extrapolation coordinates
            // as the L part of the VL surface is bending downwards in T-X, we take
            // pressures low enough to not get into any L-like field
            ip_boiling  = ComputePressureIndex(water.SaturationPressureFromT(tcurrent));
            ip_high     = ip_boiling-1; // upper p
            ip_low      = ip_high-1;    // lower p
	    
            for(xcurrent = 0.0e0; xcurrent <= 1.00001e0; xcurrent += x_res)
              {
                xbrine = XNaCl2Massfraction(xcurrent);
                GetCompositionIndex(xcurrent + 1.0e-6);
                if(ix > x_dim_table[0]-1)
                  {
                    break;
                  }

                for(pcurrent = 0.5e5; pcurrent <= 5000.e5+0.1e0; pcurrent += p_res)
                  {
                    GetPressureIndex(pcurrent+1.0e-3);
                    if(ip > p_dim-1)
                      {
                        break;
                      }

                    storage_index = (it-it_offset[0]) * x_dim_table[0]*p_dim + x_dim_table[0]*ip + ix;
                    table0[tablesize[0]*temperature_index     + storage_index ] = tcurrent;
                    table0[tablesize[0]*pressure_index        + storage_index ] = pcurrent;
                    table0[tablesize[0]*composition_index     + storage_index ] = xcurrent;
	      
	      
                    if(tcurrent < cp_h2o.Temperature())
                      {
		  
                        if(tcurrent <= 150.0)
                          {
                            // use boiling_v properties
                            table0[tablesize[0]*density_index         + storage_index ] = water.VaporDensityFromT(tcurrent);
                            table0[tablesize[0]*enthalpy_index        + storage_index ] = water.VaporEnthalpyFromT(tcurrent);
                            table0[tablesize[0]*heatcapacity_index    + storage_index ] = water.VaporHeatCapacityFromT(tcurrent);
                            table0[tablesize[0]*compressibility_index + storage_index ] = water.VaporCompressibilityFromT(tcurrent);
                            table0[tablesize[0]*viscosity_index       + storage_index ] = water.VaporViscosityFromT(tcurrent);
                          }
                        // crude fix for near-critical region
                        else if( tcurrent > 373.5 )
                          {
                            table0[tablesize[0]*density_index         + storage_index ] = water.Density(tcurrent,pcurrent) * (1.0-xcurrent);
                            table0[tablesize[0]*enthalpy_index        + storage_index ] = water.Enthalpy(tcurrent,pcurrent) * (1.0-xcurrent);
                            table0[tablesize[0]*heatcapacity_index    + storage_index ] = water.HeatCapacity(tcurrent,pcurrent) * (1.0-xcurrent);
                            table0[tablesize[0]*compressibility_index + storage_index ] = water.Compressibility(tcurrent,pcurrent) * (1.0-xcurrent);
                            table0[tablesize[0]*viscosity_index       + storage_index ] = water.Viscosity(tcurrent,pcurrent) * (1.0-xcurrent);
                          }
                        else if( ip < ip_boiling)
                          {
                            if(ix != 0)
                              {
                                table0[tablesize[0]*density_index         + storage_index ] = brine.Density();
                                table0[tablesize[0]*enthalpy_index        + storage_index ] = brine.Enthalpy();
                                table0[tablesize[0]*heatcapacity_index    + storage_index ] = brine.HeatCapacity();
                                table0[tablesize[0]*compressibility_index + storage_index ] = brine.Compressibility();
                                table0[tablesize[0]*viscosity_index       + storage_index ] = brine.Viscosity();
                              }
                            else
                              {
                                table0[tablesize[0]*density_index         + storage_index ] = water.Density(tcurrent,pcurrent);
                                table0[tablesize[0]*enthalpy_index        + storage_index ] = water.Enthalpy(tcurrent,pcurrent);
                                table0[tablesize[0]*heatcapacity_index    + storage_index ] = water.HeatCapacity(tcurrent,pcurrent);
                                table0[tablesize[0]*compressibility_index + storage_index ] = water.Compressibility(tcurrent,pcurrent);
                                table0[tablesize[0]*viscosity_index       + storage_index ] = water.Viscosity(tcurrent,pcurrent);
                              }
                            // ***************** Baustelle ***************
                            // 				if(brine.Density() > 322.0){
                            // 				    char mychar;
                            // 				    cerr << "H2ONaClLookup::H2ONaClLookup() - when computing Table 0:\n"
                            // 					 << "potentially found liquid-like vapor density!!!\n";
                            // 				    cerr << "At " << tcurrent << "C, " << pcurrent << " bar, and " << xcurrent << " x,\n";
                            // 				    cerr << "with boiling pressure being " << water.SaturationPressureFromT(tcurrent) << endl;
                            // 				    cerr << "density is " << brine.Density() << endl;
                            // 				    cin  >> block_id;
                            // 				}
                          }
                        else
                          {
                            storageindex_lowp2 = (it-it_offset[0]) * x_dim_table[0]*p_dim + x_dim_table[0]*(ip_high)     + ix;
                            storageindex_lowp1 = (it-it_offset[0]) * x_dim_table[0]*p_dim + x_dim_table[0]*(ip_low)      + ix;
                            DP_lowp      = table0[tablesize[0]*pressure_index + storageindex_lowp2 ];
                            DP_lowp     -= table0[tablesize[0]*pressure_index + storageindex_lowp1 ];
                            DP_extrapol  = pcurrent;
                            DP_extrapol -= table0[tablesize[0]*pressure_index + storageindex_lowp2 ];
                            DDensityDP          =  table0[tablesize[0]*density_index         + storageindex_lowp2];
                            DDensityDP         -=  table0[tablesize[0]*density_index         + storageindex_lowp1];
                            DDensityDP         /=  DP_lowp;
                            DEnthalpyDP         =  table0[tablesize[0]*enthalpy_index        + storageindex_lowp2];
                            DEnthalpyDP        -=  table0[tablesize[0]*enthalpy_index        + storageindex_lowp1];
                            DEnthalpyDP        /=  DP_lowp;
                            DHeatCapacityDP     =  table0[tablesize[0]*heatcapacity_index    + storageindex_lowp2];
                            DHeatCapacityDP    -=  table0[tablesize[0]*heatcapacity_index    + storageindex_lowp1];
                            DHeatCapacityDP    /=  DP_lowp; 
                            DCompressibilityDP  =  table0[tablesize[0]*compressibility_index + storageindex_lowp2];
                            DCompressibilityDP -=  table0[tablesize[0]*compressibility_index + storageindex_lowp1];
                            DCompressibilityDP /=  DP_lowp;
                            DViscosityDP        =  table0[tablesize[0]*viscosity_index       + storageindex_lowp2];
                            DViscosityDP       -=  table0[tablesize[0]*viscosity_index       + storageindex_lowp1];
                            DViscosityDP       /=  DP_lowp;
		      
                            table0[tablesize[0]*density_index         + storage_index ] = 
                              table0[tablesize[0]*density_index         + storageindex_lowp2] + DDensityDP         * DP_extrapol;
                            table0[tablesize[0]*enthalpy_index        + storage_index ] = 
                              table0[tablesize[0]*enthalpy_index        + storageindex_lowp2] + DEnthalpyDP        * DP_extrapol;
                            table0[tablesize[0]*heatcapacity_index    + storage_index ] = 
                              table0[tablesize[0]*heatcapacity_index    + storageindex_lowp2] + DHeatCapacityDP    * DP_extrapol;
                            table0[tablesize[0]*compressibility_index + storage_index ] = 
                              table0[tablesize[0]*compressibility_index + storageindex_lowp2] + DCompressibilityDP * DP_extrapol;
                            table0[tablesize[0]*viscosity_index       + storage_index ] = 
                              table0[tablesize[0]*viscosity_index       + storageindex_lowp2] + DViscosityDP       * DP_extrapol;
                          }
                      }
                    else
                      {
                        // supercritical
                        storage_index = (it-it_offset[0]) * x_dim_table[0]*p_dim + x_dim_table[0]*ip + ix;
                        table0[tablesize[0]*temperature_index     + storage_index ] = tcurrent;
                        table0[tablesize[0]*pressure_index        + storage_index ] = pcurrent;
                        table0[tablesize[0]*composition_index     + storage_index ] = xcurrent;
                        if(ix != 0)
                          {
                            table0[tablesize[0]*density_index         + storage_index ] = brine.Density();
                            table0[tablesize[0]*enthalpy_index        + storage_index ] = brine.Enthalpy();
                            table0[tablesize[0]*heatcapacity_index    + storage_index ] = brine.HeatCapacity();
                            table0[tablesize[0]*compressibility_index + storage_index ] = brine.Compressibility();
                            table0[tablesize[0]*viscosity_index       + storage_index ] = brine.Viscosity();
                          }
                        else
                          {
                            table0[tablesize[0]*density_index         + storage_index ] = water.Density(tcurrent, pcurrent);
                            table0[tablesize[0]*enthalpy_index        + storage_index ] = water.Enthalpy(tcurrent, pcurrent);
                            table0[tablesize[0]*heatcapacity_index    + storage_index ] = water.HeatCapacity(tcurrent, pcurrent);
                            table0[tablesize[0]*compressibility_index + storage_index ] = water.Compressibility(tcurrent, pcurrent);
                            table0[tablesize[0]*viscosity_index       + storage_index ] = water.Viscosity(tcurrent, pcurrent);
                          }   
                      }
                  }
              }
          }
        FILE* outfile0;
        outfile0 = fopen( filename[0], "wb");
        if( outfile0 == NULL )
          {
            cerr << "could not even it for writing, please stop program and debug !!!!\n";
            char yesno;
            cerr << "or enter any key to continue (simulation likely to crash or give wrong results!) :";
            cin  >> yesno;
          }
        else
          {
            cerr << "writing file " << filename[0] << " ... ";
            skm_C_fwrite( outfile0, table0 );
            fclose(outfile0);
            cerr << "done!\n";
          }
      }
    else
      {
        cerr << "reading file " << filename[0] << " ... ";
        skm_C_fread( infile0, table0 );
        fclose( infile0 );
        cerr << "done!\n";
      }
    cerr << "H2ONaClLookup, leaving BuildTable0 ...\n\n";
  }


  void H2ONaClLookup::PrintAll()
  {
    cout <<    Temperature() << endl;
    cout <<    Pressure() << endl;
    cout <<    MassFractionNaCl() << endl;
    cout <<    Density() << endl;
    cout <<    Enthalpy() << endl;
    cout <<    HeatCapacity() << endl;
    cout <<    Compressibility() << endl;
    cout <<    Viscosity() << endl;
    return;
  }


} // end namespace csmp
