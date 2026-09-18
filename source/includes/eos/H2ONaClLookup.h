// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef H2ONACLLOOKUP_H
#define H2ONACLLOOKUP_H

#include "CSMP_definitions.h"
#include "States.h"
#include "ErrorHandler.h"

namespace csmp
{

  // 2D grid

  class H2ONaClLookup
  {

  public:
    H2ONaClLookup(const double& externaltemperature,
                  const double& externalpressure,
                  const double& externalcomposition);
    ~H2ONaClLookup();

    double   Temperature();
    double   Pressure();
    double   MassFractionNaCl();
    double   Density();
    double   Enthalpy();
    double   HeatCapacity();
    double   Compressibility();
    double   Viscosity();
    //    double   DEnthalpyDX();
    void        PrintAll();

    double SubcriticalVaporDensity();
    double SubcriticalVaporEnthalpy();
    double SubcriticalVaporHeatCapacity();
    double SubcriticalVaporCompressibility();
    double SubcriticalVaporViscosity();

    double SubcriticalLiquidDensity();
    double SubcriticalLiquidEnthalpy();
    double SubcriticalLiquidHeatCapacity();
    double SubcriticalLiquidCompressibility();
    double SubcriticalLiquidViscosity();
   
    long GetTableID(const double& t);

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

    const double& temperature;
    const double& pressure;
    const double& composition;

    double    tcurrent,pcurrent,xcurrent,xdummy,pdummy,tdummy,t_res,p_res,x_res;
    double    t_iA, t_iB, p_iA, p_iD, x_iA, x_iE;
    double    tnorm,pnorm,xnorm;
    double    v_bottom,v_top,v_right,v_left,v_point;
    double    dh;
    double    DP_extrapol,DP_lowp,DDensityDP,DEnthalpyDP,DHeatCapacityDP,DCompressibilityDP,DViscosityDP,DP_highp;
    double    p_boil;
    std::vector<double>   table0,table1,table2,table3,table4,table5;
    std::vector<double>   t_min_table, t_max_table;
    std::vector<long>        ip_vl,it_max_table;
    std::vector<long>        it_offset,x_dim_block,t_dim_table,x_dim_table,tablesize;

    void BuildTable0();
    //    void Table0LeftRightInterpolation();
    void BuildTable1();
    void BuildTable2();
    void BuildSimpleTable( const long& i, 
                           std::vector<double>& table );

    void GetTemperatureIndex(const double& t);
    void GetPressureIndex(const double& p);
    void GetCompositionIndex(const double& p);

    long ComputeTemperatureIndex(const double& t);
    long ComputePressureIndex(const double& p);
    long ComputeCompositionIndex(const double& p);
    double Interpolate( const double& t,
                          const int& property_index );
    double TrilinearInterpolation( const long& table_id,
                                     const std::vector<double>& storage_vector,
                                     const int& property_index );
    double PropertyAtPointH( const int& property_index );

    ErrorHandler&           csmp_error;
  };

} // csmp

#endif

/**


   @todo
   - mv table building into external class!
 */
