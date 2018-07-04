#ifndef H2ONACLLOOKUP_H
#define H2ONACLLOOKUP_H

#include "CSMP_definitions.h"
#include "States.h"

namespace csmp
{

  // 2D grid

  class H2ONaClLookup
  {

  public:
    H2ONaClLookup(const double64& externaltemperature, 
                  const double64& externalpressure, 
                  const double64& externalcomposition);
    ~H2ONaClLookup();

    double64   Temperature();
    double64   Pressure();
    double64   MassFractionNaCl();
    double64   Density();
    double64   Enthalpy();
    double64   HeatCapacity();
    double64   Compressibility();
    double64   Viscosity();
    //    double64   DEnthalpyDX();
    void        PrintAll();

    double64 SubcriticalVaporDensity();
    double64 SubcriticalVaporEnthalpy();
    double64 SubcriticalVaporHeatCapacity();
    double64 SubcriticalVaporCompressibility();
    double64 SubcriticalVaporViscosity();

    double64 SubcriticalLiquidDensity();
    double64 SubcriticalLiquidEnthalpy();
    double64 SubcriticalLiquidHeatCapacity();
    double64 SubcriticalLiquidCompressibility();
    double64 SubcriticalLiquidViscosity();
   
    long GetTableID(const double64& t);

  private:
    long t_dim;
    long p_dim;
    const long x_dim;
    const long x_dim200;
    const long x_dim400;
    const long x_dim600;
    const long x_dim700;
    const long it_200;
    const long it_400;
    const long it_600;
    const long it_700;

    long  n_tables, table_id;
    long  i,it,ip,ix,iA,iB,iC,iD,iE,iF,iG,iH,block_id;
    long  iA_t, iB_t, iA_p, iD_p, iA_x, iE_x;
    long  storage_index;
    long  ip_boiling,storageindex_lowp1,storageindex_lowp2,storageindex_highp1,storageindex_highp2;
    long  ix_twophase, storage_index_left, ix_old,ip_low,ip_high;

    States   state;
    char  filename[6][60];

    bool  BinaryOut( const char* bin_name, long table_id ) const;
    bool  BinaryIn( const char* fname, long dimension, long table_id ); 
    bool file_found;

    const double64& temperature;
    const double64& pressure;
    const double64& composition;

    double64    tcurrent,pcurrent,xcurrent,xdummy,pdummy,tdummy,t_res,p_res,x_res;
    double64    t_iA, t_iB, p_iA, p_iD, x_iA, x_iE;
    double64    tnorm,pnorm,xnorm;
    double64    v_bottom,v_top,v_right,v_left,v_point;
    double64    dh;
    double64    DP_extrapol,DP_lowp,DDensityDP,DEnthalpyDP,DHeatCapacityDP,DCompressibilityDP,DViscosityDP,DP_highp;
    double64    p_boil;
    std::vector<double64>   table0,table1,table2,table3,table4,table5;
    std::vector<double64>   t_min_table, t_max_table;
    std::vector<long>        ip_vl,it_max_table;
    std::vector<long>        it_offset,x_dim_block,t_dim_table,x_dim_table,tablesize;

    void BuildTable0();
    //    void Table0LeftRightInterpolation();
    void BuildTable1();
    void BuildTable2();
    void BuildSimpleTable( const long& i, 
                           std::vector<double64>& table );

    void GetTemperatureIndex(const double64& t);
    void GetPressureIndex(const double64& p);    
    void GetCompositionIndex(const double64& p);    

    long ComputeTemperatureIndex(const double64& t);    
    long ComputePressureIndex(const double64& p);    
    long ComputeCompositionIndex(const double64& p);    
    double64 Interpolate( const double64& t,  
                          const int& property_index );
    double64 TrilinearInterpolation( const long& table_id, 
                                     const std::vector<double64>& storage_vector, 
                                     const int& property_index );
    double64 PropertyAtPointH( const int& property_index );


  };

} // csmp

#endif

/**


   @todo
   - mv table building into external class!
 */
