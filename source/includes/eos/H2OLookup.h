#ifndef  H2O_LOOKUP_H
#define  H2O_LOOKUP_H

#include "CSMP_definitions.h"
#include "ErrorHandler.h"

// TODO: move includes to cpp file
//#include "ConvertConcentrationUnitsNaCl.h" //! included in CriticalPointH2O.h
#include "CriticalPointH2O.h"
#include "LookupPropertyIndex.h"
#include "GetLookupIndices.h"
//#include "steam4.h" //! included in iaps.h
//#include "decl.h" //! included in iaps.h
#include "iaps.h"

namespace csmp {

  // 2D grid
  class H2OLookup 
  {
  public:
    H2OLookup();
    ~H2OLookup();

    double    VaporProperty(const long& it, const int& property_index); 
    double    LiquidProperty(const long& it, const int& property_index); 
    double    VaporProperty(const double& t, const int& property_index); 
    double    LiquidProperty(const double& t, const int& property_index); 
    double    SinglePhaseProperty(const long& it, const long& ip, const int& property_index);
    double    SinglePhaseProperty(const double& t, const double& p, const int& property_index);
    double    TfromP(const double& press);

    // for given p
    double    SaturationTemperatureFromP( double p );
    double    LiquidDensityFromP( double p );
    double    VaporDensityFromP( double p );
    double    LiquidEnthalpyFromP( double p );
    double    VaporEnthalpyFromP( double p );
    double    LiquidHeatCapacityFromP( double p );
    double    VaporHeatCapacityFromP( double p );
    double    LiquidViscosityFromP( double p );
    double    VaporViscosityFromP( double p );
    double    LiquidCompressibilityFromP( double p );
    double    VaporCompressibilityFromP( double p );
    
    // for given t
    double    SaturationPressureFromT( double t );
    double    DSaturationPressureFromTDT( double t );
    double    LiquidDensityFromT( double t );
    double    VaporDensityFromT( double t );
    double    LiquidEnthalpyFromT( double t );
    double    VaporEnthalpyFromT( double t );
    double    LiquidHeatCapacityFromT( double t );
    double    VaporHeatCapacityFromT( double t );
    double    LiquidViscosityFromT( double t );
    double    VaporViscosityFromT( double t );
    double    LiquidCompressibilityFromT( double t );
    double    VaporCompressibilityFromT( double t );
    
    // for given p, t
    double    Density( double t, double p );
    double    HeatCapacity( double t, double p );
    double    Enthalpy( double t, double p );
    double    Viscosity( double t, double p );
    double    Compressibility( double t, double p );
    
    // TD's additions from Feb 2007
    double    MolarVolume(double t, double p);
    double    LiquidMolarVolumeFromP(double p);
    double    VaporMolarVolumeFromP(double p);
    double    Dvdt(double t, double p);
    double    Dvdpbar(double t, double p);

  private:

    CriticalPointH2O  cp_h2o;

    void        BuildTable0And1();
    void        BuildTable2();

    double    DynamicViscosity( const double& T, const double& rho )  const;

    long        GetTemperatureIndex(const double& t);
    long        GetPressureIndex(const double& p);

    bool        BoilingCurveInInterpolationCell();
    void        GetIndex_iA(const double& tcurrent, const double& pcurrent, const int& property_index);
    double    NormalInterpolation(const double& tcurrent, const double& pcurrent, const int& property_index );
    double    LiquidInterpolationNearBoilingCurve(const double& t, const double& p, const int& property_index);
    double    VaporInterpolationNearBoilingCurve(const double& t, const double& p, const int& property_index);
    double    NearCriticalInterpolation(const double& tcurrent, const double& pcurrent , const int& property_index);
    double    AccidentalBoilingCurveEncounter(const double& t, const double& p, const int& property_index);

    long        t_dim, p_dim;
    long        n_tables, table_id;
    long        i,it,ip,iA,iB,iC,iD,i_dummy,i_min,i_max,i_guess;
    long        it_A, it_B, ip_A, ip_D;
    char        filename[6][60];
    double    tcurrent,pcurrent,pdummy;
    double    t_res,p_res,T,tboil,myt;
    double    v_top,v_bottom,v_before,v_behind,v_boil,v_boil_before,v_boil_behind;
    double    v_iA,v_iB,v_iC,v_iD,v_interpolated;
    double    pboil,p_boil_low,p_boil_high,t_iA,t_iB,t_iC,t_iD,p_iA,p_iB,p_iC,p_iD,pboil_iB,pboil_iA;
    double    dvdt;
    double    ak[4], bij[6][5];


    std::vector<double>   satvap; 
    std::vector<double>   satliq;
    std::vector<double>   singlephase;
    std::vector<long>       t_dim_table;

  };
  /**
     author: Thomas Driesner, ETH Zuerich
     contact: thomas.driesner@erdw.ethz.ch
     latest modifications: 
     - Oct  4, 2012: moved pressure units from bar to Pa and ported to CSMP++

     Licensing issues:

     Data contained in these tables are based on the Haar-Gallagher-Kell equation of state (also known as "IAPS84", check www.iapws.org for more information). In order to compute the correct values when the tables are first created, H2OLookup makes use of the "PROST4" library (as seen from the #include "steam4.h"). PROST is distributed under GPL (which version, please re-check) and, hence, PROST4 is distributed as stand-alone external library as we do currently NOT intend to distribute H2OLookup under GPL. This needs to be sorted out!

     Technically:
   
     H2OLookup generates and queries three lookup tables for the properties of pure water in binary format:
     - H2OPropertiesLookupTableSaturatedVapor.bin  ("Table0" in H2OLookup.cpp)
     - H2OPropertiesLookupTableSaturatedLiquid.bin ("Table1" in H2OLookup.cpp)
     - H2OPropertiesLookupTableSinglePhase.bin     ("Table2" in H2OLookup.cpp)

     Properties currently stored are:
     temperature [C], pressure [Pa], composition [mole fraction NaCl], density [kg m^-3], specific enthalpy [J kg^-1], specific heat capacity [J kg^-1 C^-1], compressibility [Pa^-1] and dynamic viscosity [Pa s] as defined in file "LookupPropertyIndex.h"

     Saturation properties:
     "H2OPropertiesLookupTableSaturatedVapor.bin" and "H2OPropertiesLookupTableSaturatedLiquid.bin" store the properties  of the respective phase along the boiling curve of pure water. The range of validity for these two lookup tables is 0 to 373.976 Celsius (critical temperature of H2O, see CriticalPointH2O.cpp). Temperature spacing is in accordance with other lookup classes to 1000C and 500MPa and is currently as follows: for t <= 250.0 C: 5.0 C; for t <= 350.0 C: 2.0 C; for t <= 360.0 C: 1.0 C; for t <= 370.0 C: 0.5 C; for t <= 372.0 C: 0.2 C; for t <= 373.9: 0.1 C; and 0.076 C between 373.9 C and the critical point. For the extended range that is not covered by the H2O-NaCl lookup tables (but used by them in property calculations!!!) between 1000 and 2000 C, the spacing is 50 C, and between 500 to 1000 MPa, it is 50 MPa.

     The tables are written as binary files, and each stores the properties listed above as a single, 1-dimensional std::vector<double>. The vector comprises subsequent blocks in each of which the values for the one property are stored in sequence of ascending temperature. Each block is of size t_dim (the number of entries, the first thing that is computed in the constructor), the sequence of properties is defined in the file "LookupPropertyIndex.h".
   
     Single phase properties:
     "H2OPropertiesLookupTableSinglePhase.bin" stores the properties in the single phase regions of the pure water phase diagram in a temperature-pressure grid from 0 to 2000 C and "0" to 1000 MPa ("0" means that the respective pressure is assigned 0 bar in the lookup table but values are actually computed at 5.0e2 Pa, i.e. just below the vapor pressure at the triple point. While p < 0.1 MPa is usually irrelevant in geologic applications, it may nevertheless be relevant to compute vapor properties at temperatures below 100 C where water vapor pressures ARE below 0.1 MPa. Also, the near-surface parts of hydrothermal systems at high elevation may be operating under conditions below 0.1 MPa). The grid has variable resolution (***insert details***). The table is written as a binary file and stores the properties in a single, 1-dimensional std::vector<double>. The vector comprises subsequent blocks in each of which the values for the one property are stored in sequence of ascending temperature (with all pressure entries per temperature being filled before the next temperature is accessed, i.e., storage sequence is: property0 for t = 0 C and p = 0...1000 MPa, property0 for t = 5C and p = 0...1000 MPa,...,property0 for t = 2000 C and p = 0...1000 MPa, then the same for property1 etc.). Each block is of size t_dim*p_dim (the number of entries, the first thing that is computed in the constructor), the sequence of properties is defined in the file "LookupPropertyIndex.h".

     Usage:
     In this version, only default construction of the object is possible. Public member names are self-explanatory from the header file. If temperature and/or pressure are needed as function arguments, they should be passed with units of Celsius and Pascal, respectively. 

     Requirements:
     The PROST4 library is needed. The header files "CSMP_definitions.h", "ConvertConcentrationUnitsNaCl.h", "CriticalPointH2O.h", "States.h", "LookupPropertyIndex.h", "steam4.h" are required.
   
     Known issues: 
     - Interpolations between grid points are linear, hence, if very high accuracy in the critical region is required, one needs to check the applicability
     - Possible bugs due to the move from bar-based to pasacl based have not yet been tested for as of the date stated above 
  */

} // csmp


#endif
