#include "H2OPropertiesLookUpTable.h"

using namespace std;

namespace csmp {

template<typename fT>
H2OPropertiesLookUpTable<fT>::H2OPropertiesLookUpTable( bool use_lookup_tables )
  : pcrit(22.054915e+06), // PROST P CRIT 
    tcrit(373.976), // PROST T CRIT
    dp(1.0e-08),
    kelvin(273.15), // PROST KELVIN
    dP(2.5e+05),
    dT(1.0),
    dT_fine(0.1),
    dP_fine(1.0e+04),
    d(1.0),
    patm(101325.0),
    tatm(5.0),
    tsat_atm(100.001),
    tmax_c(400.0),
    tmin_c(350.0),
    pmax_c(2.7e+07),
    pmin_c(1.7e+07),
    p_sat_min(1.0e+05),
    p_cols(1601.0), // p(j)
    t_rows(996.0), // t(i)
    pc_cols(1001.0),
    tc_rows(501.0),
    cols(static_cast<int32_t>(p_cols)),
    rows(static_cast<int32_t>(t_rows)),
    c_cols(static_cast<int32_t>(pc_cols)),
    c_rows(static_cast<int32_t>(tc_rows)),
    t_sat_rows(3691.0), // 5.0 oC to 373.9 oC in 0.1 oC steps plus 2 entries for t_atm and tcrit
    p_sat_rows(2197.0), // 1.05e+05 to 22.05e+06 Pa in 5.0e+03 Pa stest plus 2 entries for p_atm and pcrit
    sat_cols(17.0)  // p or t, rl, rv, hl, hv, cpl, cpv, ml, muv, bl, bv, al, av, dp_d_CTl, dp_d_CTv, dp_T_Cdl, dp_T_Cdv
  { 
    
    // variables for viscosity
    ak[0] =  0.0181583;
    ak[1] =  0.0177624;
    ak[2] =  0.0105287;
    ak[3] = -0.0036744;

    bij[0][0] =  0.501938;
    bij[1][0] =  0.162888;
    bij[2][0] = -0.130356;
    bij[3][0] =  0.907919;
    bij[4][0] = -0.551119;
    bij[5][0] =  0.146543;
    bij[0][1] =  0.235622;
    bij[1][1] =  0.789393;
    bij[2][1] =  0.673665;
    bij[3][1] =  1.207552;
    bij[4][1] =  0.0670665;
    bij[5][1] = -0.084337;
    bij[0][2] = -0.274637;
    bij[1][2] = -0.743539;
    bij[2][2] = -0.959456;
    bij[3][2] = -0.687343;
    bij[4][2] = -0.497089;
    bij[5][2] =  0.195286;
    bij[0][3] =  0.145831;
    bij[1][3] =  0.263129;
    bij[2][3] =  0.347247;
    bij[3][3] =  0.213486;
    bij[4][3] =  0.100754;
    bij[5][3] = -0.032932;
    bij[0][4] = -0.0270448;
    bij[1][4] = -0.0253093;
    bij[2][4] = -0.0267758;
    bij[3][4] = -0.0822904;
    bij[4][4] =  0.0602253;
    bij[5][4] = -0.0202595;
    
    // max T and P for lookup table
    pmax = dP * static_cast<fT>(cols-1);
    tmax = dT * static_cast<fT>(rows-1) + tatm;
        
    // Trying to read lookup table from file
    if ( use_lookup_tables ) {
        if ( LoadLookupTables() ) {
            cout << "\nH2OPropertiesLookUpTable(constructor): At least one lookup table is missing ";
            cout << "\nMake sure that lookup tables are in same directory as exectuable";
            cout << "\nStarting to recompute lookup tables... " << endl;
            ComputeLookupTables();
          }    
    
        if ( LoadLookupTablesCriticalPoint() ) {
            cout << "\nH2OPropertiesLookUpTable(constructor): At least one lookup table for properties around critical point is missing ";
            cout << "\nMake sure that lookup tables are in same directory as exectuable";
            cout << "\nStarting to recompute lookup tables... " << endl;
            ComputeLookupTablesCriticalPoint();
          }    

        if ( LoadLookupTablesTwoPhase() ) {
            cout << "\nH2OPropertiesLookUpTable(constructor): At least one lookup table for properties along 2-phase curve is missing ";
            cout << "\nMake sure that lookup tables are in same directory as exectuable";
            cout << "\nStarting to recompute lookup tables... " << endl;
            ComputeLookupTablesTwoPhase();
          }
      }  

   
  }

template<typename fT>
H2OPropertiesLookUpTable<fT>::H2OPropertiesLookUpTable( const H2OPropertiesLookUpTable<fT>& lookup ) 
  : pcrit(22.054915e+06), // PROST P CRIT 
    tcrit(373.976), // PROST T CRIT
    dp(1.0e-08),
    kelvin(273.15), // PROST KELVIN
    dP(2.5e+05),
    dT(1.0),
    dT_fine(0.1),
    dP_fine(1.0e+04),
    d(1.0),
    patm(101325.0),
    tatm(5.0),
    tsat_atm(100.001),
    tmax_c(400.0),
    tmin_c(350.0),
    pmax_c(2.7e+07),
    pmin_c(1.7e+07),
    p_sat_min(1.0e+05),
    p_cols(1601.0), // p(j)
    t_rows(996.0), // t(i)
    pc_cols(1001.0),
    tc_rows(501.0),
    cols(static_cast<int32_t>(p_cols)),
    rows(static_cast<int32_t>(t_rows)),
    c_cols(static_cast<int32_t>(pc_cols)),
    c_rows(static_cast<int32_t>(tc_rows)),
    t_sat_rows(3691.0), // 5.0 oC to 373.9 oC in 0.1 oC steps plus 2 entries for t_atm and tcrit
    p_sat_rows(2197.0), // 1.05e+05 to 22.05e+06 Pa in 5.0e+03 Pa stest plus 2 entries for p_atm and pcrit
    sat_cols(17.0)  // p or t, rl, rv, hl, hv, cpl, cpv, ml, muv, bl, bv, al, av, dp_d_CTl, dp_d_CTv, dp_T_Cdl, dp_T_Cdv
{

    // variables for viscosity
    ak[0] =  0.0181583;
    ak[1] =  0.0177624;
    ak[2] =  0.0105287;
    ak[3] = -0.0036744;

    bij[0][0] =  0.501938;
    bij[1][0] =  0.162888;
    bij[2][0] = -0.130356;
    bij[3][0] =  0.907919;
    bij[4][0] = -0.551119;
    bij[5][0] =  0.146543;
    bij[0][1] =  0.235622;
    bij[1][1] =  0.789393;
    bij[2][1] =  0.673665;
    bij[3][1] =  1.207552;
    bij[4][1] =  0.0670665;
    bij[5][1] = -0.084337;
    bij[0][2] = -0.274637;
    bij[1][2] = -0.743539;
    bij[2][2] = -0.959456;
    bij[3][2] = -0.687343;
    bij[4][2] = -0.497089;
    bij[5][2] =  0.195286;
    bij[0][3] =  0.145831;
    bij[1][3] =  0.263129;
    bij[2][3] =  0.347247;
    bij[3][3] =  0.213486;
    bij[4][3] =  0.100754;
    bij[5][3] = -0.032932;
    bij[0][4] = -0.0270448;
    bij[1][4] = -0.0253093;
    bij[2][4] = -0.0267758;
    bij[3][4] = -0.0822904;
    bij[4][4] =  0.0602253;
    bij[5][4] = -0.0202595;

    *this = lookup;
}

template<typename fT>
H2OPropertiesLookUpTable<fT>& H2OPropertiesLookUpTable<fT>::operator=( const H2OPropertiesLookUpTable<fT>& lookup )
{
    if ( &lookup == this ) return *this;
    
    // properties
    d          = lookup.d;
    cols       = lookup.cols;
    rows       = lookup.rows;
    c_cols     = lookup.c_cols;
    c_rows     = lookup.c_rows;
    pmax       = lookup.pmax;
    tmax       = lookup.tmax;
    pcurrent   = -1.0e+300;
    tcurrent   = -1.0e+300;
    // binary lookup tables
    density              = lookup.density;
    enthalpy             = lookup.enthalpy;
    heat_capacity        = lookup.heat_capacity;
    viscosity            = lookup.viscosity; 
    expansivity          = lookup.expansivity;
    compressibility      = lookup.compressibility;
    dp_dd_CT             = lookup.dp_dd_CT;
    dp_dT_Cd             = lookup.dp_dT_Cd;
    t_sat                = lookup.t_sat;
    p_sat                = lookup.p_sat;
    crit_density         = lookup.crit_density;
    crit_enthalpy        = lookup.crit_enthalpy;
    crit_heat_capacity   = lookup.crit_heat_capacity;
    crit_viscosity       = lookup.crit_viscosity;
    crit_expansivity     = lookup.crit_expansivity;
    crit_compressibility = lookup.crit_compressibility;
    crit_dp_dd_CT        = lookup.crit_dp_dd_CT;
    crit_dp_dT_Cd        = lookup.crit_dp_dT_Cd;
    
    return *this;
}


template<typename fT>
H2OPropertiesLookUpTable<fT>::~H2OPropertiesLookUpTable() 
 {}


template<typename fT>
bool H2OPropertiesLookUpTable<fT>::LoadLookupTables() 
 {
    bool lookup_missing(false);
    
    cout <<"\nH2OPropertiesLookUpTable<fT, dim>::LoadLookUpTables: Trying to load lookup tables into memory... " << endl;
    
    // density
    if ( density.BinaryIn("H2O_Density_LookupTable") )
        cout <<"H2OPropertiesLookUpTable<fT, dim>::LoadLookUpTables: 'H2O_Density_LookupTable' table was found... " << endl;
    else {
        lookup_missing = true;
        cout <<"H2OPropertiesLookUpTable<fT, dim>::LoadLookUpTables: 'H2O_Density_LookupTable' was not found... " << endl;
        density.Initialize( p_cols, t_rows, 1.0, 1.0 );
      }
    // enthalpy
    if ( enthalpy.BinaryIn("H2O_Enthalpy_LookupTable") ) 
        cout <<"H2OPropertiesLookUpTable<fT, dim>::LoadLookUpTables: 'H2O_Enthalpy_LookupTable' table was found... " << endl;
    else {
        lookup_missing = true;
        cout <<"H2OPropertiesLookUpTable<fT, dim>::LoadLookUpTables: 'H2O_Enthalpy_LookupTable' was not found... " << endl;
        enthalpy.Initialize( p_cols, t_rows, 1.0, 1.0 );
      }
    // heat capacity
    if ( heat_capacity.BinaryIn("H2O_HeatCapacity_LookupTable") ) 
        cout <<"H2OPropertiesLookUpTable<fT, dim>::LoadLookUpTables: 'H2O_HeatCapacity_LookupTable' was found... " << endl;
    else {
        lookup_missing = true;
        cout <<"H2OPropertiesLookUpTable<fT, dim>::LoadLookUpTables: 'H2O_HeatCapacity_LookupTable' was not found... " << endl;
        heat_capacity.Initialize( p_cols, t_rows, 1.0, 1.0 );
      }
    // viscosity
    if ( viscosity.BinaryIn("H2O_Viscosity_LookupTable") ) 
        cout <<"H2OPropertiesLookUpTable<fT, dim>::LoadLookUpTables: 'H2O_Viscosity_LookupTable' was found... " << endl;
    else {
        lookup_missing = true;
        cout <<"H2OPropertiesLookUpTable<fT, dim>::LoadLookUpTables: 'H2O_Viscosity_LookupTable' was not found... " << endl;
        viscosity.Initialize( p_cols, t_rows, 1.0, 1.0 );
      }
    // thermal expansivity
    if ( expansivity.BinaryIn("H2O_Expansivity_LookupTable") ) 
        cout <<"H2OPropertiesLookUpTable<fT, dim>::LoadLookUpTables: 'H2O_Expansivity_LookupTable' was found... " << endl;
    else {
        lookup_missing = true;
        cout <<"H2OPropertiesLookUpTable<fT, dim>::LoadLookUpTables: 'H2O_Expansivity_LookupTable' was not found... " << endl;
        expansivity.Initialize( p_cols, t_rows, 1.0, 1.0 );
      }
    // compressibility
    if ( compressibility.BinaryIn("H2O_Compressibility_LookupTable") ) 
        cout <<"H2OPropertiesLookUpTable<fT, dim>::LoadLookUpTables: 'H2O_Compressibility_LookupTable' was found... " << endl;
    else {
        lookup_missing = true;
        cout <<"H2OPropertiesLookUpTable<fT, dim>::LoadLookUpTables: 'H2O_Compressibility_LookupTable' was not found... " << endl;
        compressibility.Initialize( p_cols, t_rows, 1.0, 1.0 );
      }
    // first derivatives
    if ( dp_dd_CT.BinaryIn("H2O_dp_dd_CT_LookupTable") ) 
        cout <<"H2OPropertiesLookUpTable<fT, dim>::LoadLookUpTables: 'H2O_dp_dd_CT_LookupTable' was found... " << endl;
    else {
        lookup_missing = true;
        cout <<"H2OPropertiesLookUpTable<fT, dim>::LoadLookUpTables: 'H2O_dp_d_CT_LookupTable' was not found... " << endl;
        dp_dd_CT.Initialize( p_cols, t_rows, 1.0, 1.0 );
      }
    
    if ( dp_dT_Cd.BinaryIn("H2O_dp_dT_Cd_LookupTable") ) 
        cout <<"H2OPropertiesLookUpTable<fT, dim>::LoadLookUpTables: 'H2O_dp_dT_Cd_LookupTable' was found... " << endl;
    else {
        lookup_missing = true;
        cout <<"H2OPropertiesLookUpTable<fT, dim>::LoadLookUpTables: 'H2O_dp_dT_Cd_LookupTable' was not found... " << endl;
        dp_dT_Cd.Initialize( p_cols, t_rows, 1.0, 1.0 );
      }
      
    return lookup_missing;  
 
 }


template<typename fT>
bool H2OPropertiesLookUpTable<fT>::LoadLookupTablesCriticalPoint() 
 {
    bool lookup_missing(false);
    
    cout <<"\nH2OPropertiesLookUpTable<fT, dim>::LoadLookUpTablesCriticalPoint: Trying to load lookup tables for critical point into memory... " << endl;
    
    // density
    if ( crit_density.BinaryIn("H2O_Density_LookupTableCriticalPoint") )
        cout <<"H2OPropertiesLookUpTable<fT, dim>::LoadLookUpTablesCriticalPoint: 'H2O_Density_LookupTableCriticalPoint' table was found... " << endl;
    else {
        lookup_missing = true;
        cout <<"H2OPropertiesLookUpTable<fT, dim>::LoadLookUpTablesCriticalPoint: 'H2O_Density_LookupTableCriticalPoint' was not found... " << endl;
        crit_density.Initialize( pc_cols, tc_rows, 1.0, 1.0 );
      }
    // enthalpy
    if ( crit_enthalpy.BinaryIn("H2O_Enthalpy_LookupTableCriticalPoint") ) 
        cout <<"H2OPropertiesLookUpTable<fT, dim>::LoadLookUpTablesCriticalPoint: 'H2O_Enthalpy_LookupTableCriticalPoint' table was found... " << endl;
    else {
        lookup_missing = true;
        cout <<"H2OPropertiesLookUpTable<fT, dim>::LoadLookUpTablesCriticalPoint: 'H2O_Enthalpy_LookupTableCriticalPoint' was not found... " << endl;
        crit_enthalpy.Initialize( pc_cols, tc_rows, 1.0, 1.0 );
      }
    // heat capacity
    if ( crit_heat_capacity.BinaryIn("H2O_HeatCapacity_LookupTableCriticalPoint") ) 
        cout <<"H2OPropertiesLookUpTable<fT, dim>::LoadLookUpTablesCriticalPoint: 'H2O_HeatCapacity_LookupTableCriticalPoint' was found... " << endl;
    else {
        lookup_missing = true;
        cout <<"H2OPropertiesLookUpTable<fT, dim>::LoadLookUpTablesCriticalPoint: 'H2O_HeatCapacity_LookupTableCriticalPoint' was not found... " << endl;
        crit_heat_capacity.Initialize( pc_cols, tc_rows, 1.0, 1.0 );
      }
    // viscosity
    if ( crit_viscosity.BinaryIn("H2O_Viscosity_LookupTableCriticalPoint") ) 
        cout <<"H2OPropertiesLookUpTable<fT, dim>::LoadLookUpTablesCriticalPoint: 'H2O_Viscosity_LookupTableCriticalPoint' was found... " << endl;
    else {
        lookup_missing = true;
        cout <<"H2OPropertiesLookUpTable<fT, dim>::LoadLookUpTablesCriticalPoint: 'H2O_Viscosity_LookupTableCriticalPoint' was not found... " << endl;
        crit_viscosity.Initialize( pc_cols, tc_rows, 1.0, 1.0 );
      }
    // thermal expansivity
    if ( crit_expansivity.BinaryIn("H2O_Expansivity_LookupTableCriticalPoint") ) 
        cout <<"H2OPropertiesLookUpTable<fT, dim>::LoadLookUpTablesCriticalPoint: 'H2O_Expansivity_LookupTableCriticalPoint' was found... " << endl;
    else {
        lookup_missing = true;
        cout <<"H2OPropertiesLookUpTable<fT, dim>::LoadLookUpTablesCriticalPoint: 'H2O_Expansivity_LookupTableCriticalPoint' was not found... " << endl;
        crit_expansivity.Initialize( pc_cols, tc_rows, 1.0, 1.0 );
      }
    // compressibility
    if ( crit_compressibility.BinaryIn("H2O_Compressibility_LookupTableCriticalPoint") ) 
        cout <<"H2OPropertiesLookUpTable<fT, dim>::LoadLookUpTablesCriticalPoint: 'H2O_Compressibility_LookupTableCriticalPoint' was found... " << endl;
    else {
        lookup_missing = true;
        cout <<"H2OPropertiesLookUpTable<fT, dim>::LoadLookUpTablesCriticalPoint: 'H2O_Compressibility_LookupTableCriticalPoint' was not found... " << endl;
        crit_compressibility.Initialize( pc_cols, tc_rows, 1.0, 1.0 );
      }
    // first derivatives
    if ( crit_dp_dd_CT.BinaryIn("H2O_dp_dd_CT_LookupTableCriticalPoint") ) 
        cout <<"H2OPropertiesLookUpTable<fT, dim>::LoadLookUpTables: 'H2O_dp_dd_CT_LookupTableCriticalPoint' was found... " << endl;
    else {
        lookup_missing = true;
        cout <<"H2OPropertiesLookUpTable<fT, dim>::LoadLookUpTables: 'H2O_dp_d_CT_LookupTableCriticalPoint' was not found... " << endl;
        crit_dp_dd_CT.Initialize( pc_cols, tc_rows, 1.0, 1.0 );
      }
    
    if ( crit_dp_dT_Cd.BinaryIn("H2O_dp_dT_Cd_LookupTableCriticalPoint") ) 
        cout <<"H2OPropertiesLookUpTable<fT, dim>::LoadLookUpTables: 'H2O_dp_dT_Cd_LookupTableCriticalPoint' was found... " << endl;
    else {
        lookup_missing = true;
        cout <<"H2OPropertiesLookUpTable<fT, dim>::LoadLookUpTables: 'H2O_dp_dT_Cd_LookupTableCriticalPoint' was not found... " << endl;
        crit_dp_dT_Cd.Initialize( pc_cols, tc_rows, 1.0, 1.0 );
      }
      
    return lookup_missing;  
 
 }

 
template<typename fT>
bool H2OPropertiesLookUpTable<fT>::LoadLookupTablesTwoPhase() 
 {
    bool lookup_missing(false);
    
    cout <<"\nH2OPropertiesLookUpTable<fT, dim>::LoadLookUpTablesTwoPhaseCurve: Trying to load lookup tables for two-phase curve into memory... " << endl;
    
    // saturation properties for given p
    if ( p_sat.BinaryIn("H2O_P_Saturation_LookupTable") )
        cout <<"H2OPropertiesLookUpTable<fT, dim>::LoadLookUpTablesTwoPhaseCurve: 'H2O_P_Saturation_LookupTable' table was found... " << endl;
    else {
        lookup_missing = true;
        cout <<"H2OPropertiesLookUpTable<fT, dim>::LoadLookUpTablesTwoPhaseCurve: 'H2O_P_Saturation_LookupTable' was not found... " << endl;
        p_sat.Initialize( sat_cols, p_sat_rows, 1.0, 1.0 );
      }
    // saturation properties for given t
    if ( t_sat.BinaryIn("H2O_T_Saturation_LookupTable") ) 
        cout <<"H2OPropertiesLookUpTable<fT, dim>::LoadLookUpTablesTwoPhaseCurve: 'H2O_T_Saturation_LookupTable' table was found... " << endl;
    else {
        lookup_missing = true;
        cout <<"H2OPropertiesLookUpTable<fT, dim>::LoadLookUpTablesTwoPhaseCurve: 'H2O_T_Saturation_LookupTable' was not found... " << endl;
        t_sat.Initialize( sat_cols, t_sat_rows, 1.0, 1.0 );
      }  
    return lookup_missing;  
 
 }

 
template<typename fT>
void H2OPropertiesLookUpTable<fT>::ComputeLookupTables()   
 { 
    Prop *properties, *liqprops, *vapprops;;
    properties = newProp('T', 'p', 1);
    liqprops   = newProp('T', 'p', 1);
    vapprops   = newProp('T', 'p', 1);
    
    // p increment 1bar, t increment 1 oC
    fT p0, t, t0, t1, r0, r1, alpha, mu1, dtemp;

    cout << "\nH2OPropertiesLookUpTable<fT, dim>::ComputeLookupTables: ";
    cout << "\nComputing H2O properties in the range from " << tatm << " to " << tmax << " oC ";
    cout << "and from " << patm << " to " << pmax << " Pa " << endl;

    // loop over rows (temperature)
    for ( int32_t i=0; i<rows; i++ ) {
        t   = tatm + (dT * static_cast<fT>(i));
        t0  = t + kelvin;
        // loop over columns (pressure)
        for ( int32_t j=0; j<cols; j++ ) {
            if ( j == 0 ) p0 = patm;
            else          p0 = dP * static_cast<fT>(j);
            // properties at given t and p
            water_tp( t0, p0, d, dp, properties );
            CheckError( properties );
            
            // write properties to lookup tables
            density(i,j)         = properties->d;
            enthalpy(i,j)        = properties->h;
            heat_capacity(i,j)   = properties->cp;
            compressibility(i,j) = 1.0/(properties->dp->d_CT*properties->d);
            dp_dd_CT(i,j)        = properties->dp->d_CT;
            dp_dT_Cd(i,j)        = properties->dp->T_Cd;
            
            // checl if viscosity is computed correct (Prost has smaller range of validity) if necessary use other equation
            //mu = viscos( properties );
            //if ( mu <= 0.0 ) mu = DynamicViscosity( t, properties->d );
            mu1 = DynamicViscosity( t, properties->d );
            viscosity(i,j) = mu1;
            
            //if ( t >= 370.0 ) 
            //cout << "\nP: " << p0 << ", T: " << t << ", r: " << density(i,j) << ", h: " << enthalpy(i,j) << ", cp: " << heat_capacity(i,j) << ", mu: " << viscosity(i,j) << ", b: " << compressibility(i,j);

            // compute thermal expansivity
            r0    = properties->d;
            dtemp = 0.1;
            t1    = t0 + dtemp;
                
            // check if two-phase curve is crossed
            if ( p0 < pcrit ) {
                sat_p( p0, liqprops, vapprops );
                CheckError( liqprops );
                CheckError( vapprops );

                if ( t0 < liqprops->T && t1 > liqprops->T ) {
                    r1 = liqprops->d;
                    dtemp = liqprops->T - t0;
                  }
                else {
                    water_tp( t1, p0, d, dp, properties );
                    CheckError( properties );
                    r1 = properties->d;
                  }
              }
            else {
                water_tp( t1, p0, d, dp, properties ); 
                CheckError( properties );
                r1 = properties->d;
              }
                  
            alpha  = r1 - r0;
            alpha /= -r0;
            alpha /= dtemp;
            
            if ( alpha > 1.0 || alpha < 0.0 )
            cout << "\np: " << p0 << ", t: " << t << ", r0: " << r0 << ", r1: " << r1 << ", a: " << alpha << ", dt: " << dtemp;
                                           
            // store thermal expansivity      
            expansivity(i,j) = alpha;
         
         } // end p
       } // end t  
     
     freeProp( properties );
     freeProp( liqprops );
     freeProp( vapprops );
     
     cout << "\nH2OPropertiesLookUpTable<fT, dim>::ComputeLookupTable: Lookup Tables successfully computed, saving as binary files... "; 
     
     // saving FD Grids to Binary Files
     density.BinaryOut("H2O_Density_LookupTable");
     enthalpy.BinaryOut("H2O_Enthalpy_LookupTable"); 
     heat_capacity.BinaryOut("H2O_HeatCapacity_LookupTable"); 
     viscosity.BinaryOut("H2O_Viscosity_LookupTable");   
     expansivity.BinaryOut("H2O_Expansivity_LookupTable"); 
     compressibility.BinaryOut("H2O_Compressibility_LookupTable"); 
     dp_dd_CT.BinaryOut("H2O_dp_dd_CT_LookupTable"); 
     dp_dT_Cd.BinaryOut("H2O_dp_dT_Cd_LookupTable"); 
  }  



template<typename fT>
void H2OPropertiesLookUpTable<fT>::ComputeLookupTablesCriticalPoint()   
 { 

    Prop *properties, *liqprops, *vapprops;
    properties = newProp('T', 'p', 1);
    liqprops   = newProp('T', 'p', 1);
    vapprops   = newProp('T', 'p', 1);

    // p increment 1bar, t increment 1 oC
    fT p0, t, t0, t1, r0, r1, alpha, mu1, dtemp, dT_correc(0.015);

    cout << "\nH2OPropertiesLookUpTable<fT, dim>::ComputeLookupTablesCriticalPoint: ";
    cout << "\nComputing refined H2O properties around critical point in the range from " << tmin_c << " to " << tmax_c << " oC ";
    cout << "and from " << pmin_c << " to " << pmax_c << " Pa " << endl;

    // loop over rows (temperature)
    for ( int32_t i=0; i<c_rows; i++ ) {
        t   = tmin_c + (dT_fine * static_cast<fT>(i));
        t0  = t + kelvin;
        // loop over columns (pressure)
        for ( int32_t j=0; j<c_cols; j++ ) {
            p0 = pmin_c + (dP_fine * static_cast<fT>(j));
            // properties at given t and p
            if ( t == 374.0 || t == 374.1 ) water_tp( t0, p0, 500.0, dp, properties );
            else                            water_tp( t0, p0, d, dp, properties );
            //CheckError( properties );
            if ( properties->error == 1 ) {
                  //cout << "\nP: " << p0 << ", T: " << t0-kelvin << endl;
                  // check that T is moved away from saturation curve
                  sat_p( p0, liqprops, vapprops );                  
                  // decreae T if T is below saturation T
                  if ( t0 < liqprops->T ) t0 -= dT_correc;
                  // increase T if T is above saturation T
                  else                    t0 += dT_correc;
                  //cout << "\nCorrected P: " << p0 << ", T: " << t0-kelvin << endl;
                  if ( t == 374.0 || t == 374.1 ) water_tp( t0, p0, 500.0, dp, properties );
                  else                            water_tp( t0, p0, d, dp, properties );
                  CheckError( properties );
              }
                  
            
            // write properties to lookup tables
            crit_density(i,j)         = properties->d;
            crit_enthalpy(i,j)        = properties->h;
            crit_heat_capacity(i,j)   = properties->cp;
            crit_compressibility(i,j) = 1.0/(properties->dp->d_CT*properties->d);
            crit_dp_dd_CT(i,j)        = properties->dp->d_CT;
            crit_dp_dT_Cd(i,j)        = properties->dp->T_Cd;
            
            // checl if viscosity is computed correct (Prost has smaller range of validity) if necessary use other equation
            //mu = viscos( properties );
            //if ( mu <= 0.0 ) mu = DynamicViscosity( t, properties->d );
            mu1 = DynamicViscosity( t, properties->d );
            crit_viscosity(i,j) = mu1;
            
            //if ( crit_compressibility(i,j) <= 0.0 || crit_heat_capacity(i,j) <= 0.0 | crit_density(i,j) <= 0.0 | crit_enthalpy(i,j ) <= 0.0 ) 
            //cout << "\nP: " << p0 << ", T: " << t << ", r: " << crit_density(i,j) << ", h: " << crit_enthalpy(i,j) << ", cp: " << crit_heat_capacity(i,j) << ", mu: " << crit_viscosity(i,j) << ", b: " << crit_compressibility(i,j);

            // compute thermal expansivity
            r0    = properties->d;
            dtemp = 0.1;
            t1    = t0 + dtemp;
                
            // check if two-phase curve is crossed
            if ( p0 < pcrit ) {
                sat_p( p0, liqprops, vapprops );
                //CheckError( liqprops );
                //CheckError( vapprops );

                if ( t0 < liqprops->T && t1 > liqprops->T ) {
                    r1 = liqprops->d;
                    dtemp = liqprops->T - t0;
                  }
                else {
                    if ( t1 == (374.0+kelvin) || t1 == (374.1+kelvin) ) water_tp( t1, p0, 500.0, dp, properties );
                    else                                                water_tp( t1, p0, d, dp, properties );
                    //CheckError( properties );
                    r1 = properties->d;                  
                    if ( properties->error == 1 ) {
                        //cout << "\nP: " << p0 << ", T: " << t1-kelvin << endl;
                        t1 -= 0.01;
                        dtemp -= 0.01;
                        //cout << "\nCorrected P: " << p0 << ", T: " << t0-kelvin << endl;
                        if ( t == 374.0 || t == 374.1 ) water_tp( t0, p0, 500.0, dp, properties );
                        else                            water_tp( t0, p0, d, dp, properties );
                        CheckError( properties );
                    }
                 }
              }
            else {
                if ( t1 == (374.0+kelvin) || t1 == (374.1+kelvin) ) water_tp( t1, p0, 500.0, dp, properties );
                else                                                water_tp( t1, p0, d, dp, properties );
                //CheckError( properties );
                r1 = properties->d;
              }
                  
            alpha  = r1 - r0;
            alpha /= -r0;
            alpha /= dtemp;
            
            if ( alpha < 0.0 )
            cout << "\np: " << p0 << ", t: " << t << ", r0: " << r0 << ", r1: " << r1 << ", a: " << alpha << ", dt: " << dtemp;
                                           
            // store thermal expansivity      
            crit_expansivity(i,j) = alpha;
            if ( crit_compressibility(i,j) <= 0.0 )
            cout << "\nP: " << p0 << ", T: " << t << ", a: " << crit_expansivity(i,j);
            t0  = t + kelvin; // reset temperature if corrected previously
         
         } // end p
       } // end t  
     
     freeProp( properties );
     freeProp( liqprops );
     freeProp( vapprops );
     
     cout << "\nH2OPropertiesLookUpTable<fT, dim>::ComputeTableCriticalPoint: Lookup Tables fro critical point successfully computed, saving as binary files... "; 
     
     // saving FD Grids to Binary Files
     crit_density.BinaryOut("H2O_Density_LookupTableCriticalPoint");
     crit_enthalpy.BinaryOut("H2O_Enthalpy_LookupTableCriticalPoint"); 
     crit_heat_capacity.BinaryOut("H2O_HeatCapacity_LookupTableCriticalPoint"); 
     crit_viscosity.BinaryOut("H2O_Viscosity_LookupTableCriticalPoint");   
     crit_expansivity.BinaryOut("H2O_Expansivity_LookupTableCriticalPoint"); 
     crit_compressibility.BinaryOut("H2O_Compressibility_LookupTableCriticalPoint"); 
     crit_dp_dd_CT.BinaryOut("H2O_dp_dd_CT_LookupTableCriticalPoint"); 
     crit_dp_dT_Cd.BinaryOut("H2O_dp_dT_Cd_LookupTableCriticalPoint"); 
      
  }  


template<typename fT>
void H2OPropertiesLookUpTable<fT>::ComputeLookupTablesTwoPhase()   
 {
    
    Prop *properties, *liqprops, *vapprops;
    properties = newProp('T', 'p', 1);
    liqprops   = newProp('T', 'p', 1);
    vapprops   = newProp('T', 'p', 1); 
            
    fT mu1, dtemp, r, alpha, t, p, dp2(5.0e-07), plow(1.0e+04); 
       
    cout << "\nH2OPropertiesLookUpTable<fT>::ComputeLookupTablesTwoPhase(): Computing lookup tables along 2-phase curve. " << endl;
    
    for ( int32_t i=0; i<t_sat.Rows(); i++ ) {
        if ( i < t_sat.Rows()-1) t = tatm + (dT_fine * static_cast<fT>(i)) + kelvin;
        else                     t = tcrit + kelvin;
        
        // get saturation properties for given and check for Prost errors
        sat_t( t, liqprops, vapprops );
        CheckError( liqprops );
        CheckError( vapprops );
                
        // save 2-phase properties in the following order p, rl, rv, hl, hv, cpl, cpv, ml, mv, bl, bv, al, av
        t_sat(i,0)  = liqprops->p;
        t_sat(i,1)  = liqprops->d;
        t_sat(i,2)  = vapprops->d;
        t_sat(i,3)  = liqprops->h;
        t_sat(i,4)  = vapprops->h;
        t_sat(i,5)  = liqprops->cp;
        t_sat(i,6)  = vapprops->cp;
        //mu = viscos( liqprops );
        mu1 = DynamicViscosity( t-kelvin, liqprops->d );
        t_sat(i,7) = mu1;
        //if ( mu <= 0.0 ) t_sat(i,7) = DynamicViscosity( t, liqprops->d );
        //mu = viscos( vapprops );
        mu1 = DynamicViscosity( t-kelvin, vapprops->d );
        t_sat(i,8) = mu1;
        //if ( mu <= 0.0 ) t_sat(i,8) = DynamicViscosity( t, vapprops->d );
        t_sat(i,9)  = 1.0/(liqprops->dp->d_CT*liqprops->d);
        t_sat(i,10) = 1.0/(vapprops->dp->d_CT*vapprops->d);
        t_sat(i,13) = liqprops->dp->d_CT;
        t_sat(i,14) = vapprops->dp->d_CT;
        t_sat(i,15) = liqprops->dp->T_Cd;
        t_sat(i,16) = vapprops->dp->T_Cd;
        
        
        // check for negative heat capacities or compressibilities (can happen in Prost)
        if ( t_sat(i,5) <= 0.0 )  t_sat(i,5)  = NumericalHeatCapacity( liqprops->h, t, liqprops->p, -0.01 );  
        if ( t_sat(i,6) <= 0.0 )  t_sat(i,6)  = NumericalHeatCapacity( vapprops->h, t, vapprops->p,  0.01 );
        if ( t_sat(i,9) <= 0.0 )  t_sat(i,9)  = NumericalCompressibility( liqprops->d, t, liqprops->p,  1000.0 );  
        if ( t_sat(i,10) <= 0.0 ) t_sat(i,10) = NumericalCompressibility( vapprops->d, t, vapprops->p, -1000.0 );
        if ( t_sat(i,13) <= 0.0 ) t_sat(i,13) = NumericalPressureDensityDerivative( liqprops->d, t, liqprops->p,  1000.0 );  
        if ( t_sat(i,14) <= 0.0 ) t_sat(i,14) = NumericalPressureDensityDerivative( vapprops->d, t, vapprops->p, -1000.0 );

        //cout <<  "\nT: " << t-kelvin << ", i: " << i << ", dp->d_CTl: " << t_sat(i,13) << ", dp->d_CTv: " << t_sat(i,14) << ", dp->T_Cdl: " << t_sat(i,15) << ", dp->T_Cdv: " << t_sat(i,16);
        
        //if ( t-kelvin > 370.0 ) {
        //   cout << "\nT: " << t-kelvin << ", p: " << t_sat(i,0) << ", rl: " << t_sat(i,1) << ", rv: " << t_sat(i,2) << ", hl: " << t_sat(i,3) << ", hv: " << t_sat(i,4);
        //   cout << ", cpl: " << t_sat(i,5) << ", cpv: " << t_sat(i,6)  << ", ml: " << t_sat(i,7) << ", mv: " << t_sat(i,8) << ", bl: " << t_sat(i,9);
        //   cout << ", bv: " << t_sat(i,10);
        //  }
        
        // thermal expansivity liquid phase
        if ( liqprops->p < plow ) dtemp = 0.5; // prost needs larger delta T and less precision for very low pressures close at saturation curve
        else                      dtemp = 0.1;
        if ( liqprops->p < plow ) water_tp( (t-dtemp), liqprops->p, d, dp2, properties );
        else                      water_tp( (t-dtemp), liqprops->p, d, dp,  properties );
        CheckError( properties ); 
        r = properties->d;
        alpha  = liqprops->d - r;
        alpha /= -r;
        alpha /= dtemp;
        t_sat(i,11) = alpha;
        //cout << ", t: " << t-dtemp-kelvin << " , r: " << r << ", al: " << t_sat(i,11);

        // thermal expansivity vapor phase
        if ( i < t_sat.Rows()-1) {
            if ( liqprops->p < plow ) water_tp( (t+dtemp), liqprops->p, d, dp2, properties );
            else                      water_tp( (t+dtemp), liqprops->p, d, dp,  properties );
            CheckError( properties ); 
            r = properties->d;
            alpha  =  r - vapprops->d;
            alpha /= -r;
            alpha /= dtemp;
            t_sat(i,12) = alpha;
          }
        else t_sat(i,12) = t_sat(i,11); // at critical point only one thermal expansivity (ste to the one of 'liquid' phase)
        //cout << ", t: " << t+dtemp-kelvin << " , r: " << r << ", av: " << t_sat(i,12);
        
        
      }  // end t loop  
       

    // get properties for pcrit from above 
    for ( int32_t i=0; i<p_sat.Rows()-1; i++ ) {
        if ( i == 0 )  p = patm;
        else           p = p_sat_min + (dP_fine * static_cast<fT>(i));
        
        // get saturation properties for given and check for Prost errors
        sat_p( p, liqprops, vapprops );
        CheckError( liqprops );
        CheckError( vapprops );
        
        // save 2-phase properties in the following order p, rl, rv, hl, hv, cpl, cpv, ml, mv, bl, bv, al, av
        p_sat(i,0)  = liqprops->T-kelvin;
        p_sat(i,1)  = liqprops->d;
        p_sat(i,2)  = vapprops->d;
        p_sat(i,3)  = liqprops->h;
        p_sat(i,4)  = vapprops->h;
        p_sat(i,5)  = liqprops->cp;
        p_sat(i,6)  = vapprops->cp;
        //mu = viscos( liqprops );
        mu1 = DynamicViscosity( t-kelvin, liqprops->d );
        p_sat(i,7) = mu1;
        //if ( mu <= 0.0 ) p_sat(i,7) = DynamicViscosity( t, liqprops->d );
        //mu = viscos( vapprops );
        mu1 = DynamicViscosity( t-kelvin, vapprops->d );
        p_sat(i,8) = mu1;
        //if ( mu <= 0.0 ) p_sat(i,8) = DynamicViscosity( t, vapprops->d );
        p_sat(i,9)  = 1.0/(liqprops->dp->d_CT*liqprops->d);
        p_sat(i,10) = 1.0/(vapprops->dp->d_CT*vapprops->d);
        p_sat(i,13) = liqprops->dp->d_CT;
        p_sat(i,14) = vapprops->dp->d_CT;
        p_sat(i,15) = liqprops->dp->T_Cd;
        p_sat(i,16) = vapprops->dp->T_Cd;
        

        // check for negative heat capacities or compressibilities (can happen in Prost)
        if ( p_sat(i,5) <= 0.0 )  p_sat(i,5)  = NumericalHeatCapacity( liqprops->h, liqprops->T, p, -0.01 );  
        if ( p_sat(i,6) <= 0.0 )  p_sat(i,6)  = NumericalHeatCapacity( vapprops->h, vapprops->T, p,  0.01 );
        if ( p_sat(i,9) <= 0.0 )  p_sat(i,9)  = NumericalCompressibility( liqprops->d, liqprops->T, p,  2500.0 );  
        if ( p_sat(i,10) <= 0.0 ) p_sat(i,10) = NumericalCompressibility( vapprops->d, vapprops->T, p, -2500.0 );
        if ( p_sat(i,13) <= 0.0 ) p_sat(i,13) = NumericalPressureDensityDerivative( liqprops->d, liqprops->T, p,  2500.0 );  
        if ( p_sat(i,14) <= 0.0 ) p_sat(i,14) = NumericalPressureDensityDerivative( vapprops->d, vapprops->T, p, -2500.0 );
        
        //if ( p > 22.0e+06 ) {
        //    cout << "\nP: " << p << ", T: " << p_sat(i,0) << ", rl: " << p_sat(i,1) << ", rv: " << p_sat(i,2) << ", hl: " << p_sat(i,3) << ", hv: " << p_sat(i,4);
        //    cout << ", cpl: " << p_sat(i,5) << ", cpv: " << p_sat(i,6)  << ", ml: " << p_sat(i,7) << ", mv: " << p_sat(i,8) << ", bl: " << p_sat(i,9);
        //    cout << ", bv: " << p_sat(i,10);
        //  }
        
        // thermal expansivity liquid phase
        water_tp( (liqprops->T-dtemp), p, d, dp,  properties );
        CheckError( properties ); 
        r = properties->d;
        alpha  = liqprops->d - r;
        alpha /= -r;
        alpha /= dtemp;
        p_sat(i,11) = alpha;
        //cout << ", t: " << t-dtemp-kelvin << " , r: " << r << ", al: " << t_sat(i,11);

        // thermal expansivity vapor phase
        water_tp( (liqprops->T+dtemp), p, d, dp,  properties );
        CheckError( properties ); 
        r = properties->d;
        alpha  =  r - vapprops->d;
        alpha /= -r;
        alpha /= dtemp;
        p_sat(i,12) = alpha;
        //cout << ", t: " << t+dtemp-kelvin << " , r: " << r << ", av: " << t_sat(i,12);
        
        
      }  // end p loop 
      
      // properties at pcrit from t_sat
      p_sat(p_sat.Rows()-1,0) = tcrit;
      for ( int32_t i=1; i<sat_cols; i++ ) p_sat(p_sat.Rows()-1,i) = t_sat(t_sat.Rows()-1,i);
      //for ( int32_t i=0; i<sat_cols; i++ ) cout << "\nValue: " << p_sat(p_sat.Rows()-1,i);
      
      freeProp( properties );
      freeProp( liqprops );
      freeProp( vapprops );

      cout << "\nH2OPropertiesLookUpTable<fT, dim>::ComputeLookupTableTwoPhase: Lookup Tables at two-phase curve successfully computed, saving as binary files... ";

      // save lookup tables as binaries
      t_sat.BinaryOut("H2O_T_Saturation_LookupTable");
      p_sat.BinaryOut("H2O_P_Saturation_LookupTable"); 
      

 }


template<typename fT> 
void H2OPropertiesLookUpTable<fT>::GetLookupTableEntries( fT t, fT p )
 {
   // 0. only if p or t has changed find new entries
   if ( p == pcurrent && t == tcurrent ) return;
   
   // 1. check that p and t are witin lookup table range
   // using a ccw numbering of points, point 1 is at t_min, p_min
   if ( t < tmax && t > tatm && p < pmax && p > patm ) {
       // 1.a get row and column entries for bilinear interpolation
       // first column entry -> t at first node (i1)
       ij[0][0] = static_cast<int32_t>((t-tatm)/dT);
       // second column entry -> p at first node (j1)
       if ( p < dP ) ij[0][1] = 0;
       else          ij[0][1] = static_cast<int32_t>(p/dP);
       // now increment ij in a counter-clockwise fashion (see NumRec in C p. 123)
       ij[1][0] = ij[0][0]+1; // i2
       ij[1][1] = ij[0][1];   // j2
       ij[2][0] = ij[1][0];   // i3
       ij[2][1] = ij[0][1]+1; // j3
       ij[3][0] = ij[0][0];   // i4
       ij[3][1] = ij[2][1];   // j4
       // 1.b check if entries lie around 2-phase curve and properties from steam and liquid field are read in at the same time
       // note that 2-phase curve could also be crossed if t < tcrit and p < pcrit but only by a little. in this case, however,
       // use above bilinear interpolation as properties *very* close to critical point are similar on liquid and vapor side,
       // other interpolation schemes are time consuming and prost does not give the most exact values anyway.
       if ( t < tcrit && p < pcrit ) {
           fT p1(static_cast<fT>(ij[0][1])*dP),
              p2(p1),
              p3(static_cast<fT>(ij[2][1])*dP),
              p4(p3),
              t1(static_cast<fT>(ij[0][0])*dT+tatm),
              t2(static_cast<fT>(ij[1][0])*dT+tatm),
              psat1(SaturationPressureFromT(t1)),
              psat2(SaturationPressureFromT(t2));
           
           // for t at ij[0][0] and ij[1][0] test if the corresponding p is either always above
           // or always below the saturation pressure (i.e., all in liquid or all in vapor field)
           // in this case, the p-t positions identified above are o.k.
           if ( ( p1 < psat1 && p4 < psat1 && p2 < psat2 && p3 < psat2 ) || // all vapor
                ( p1 > psat1 && p4 > psat1 && p2 > psat2 && p3 > psat2 ) )  // all liquid
                {
                  interpolation_type = regular;
                  return;
                }
           // 4 interpolation points sourround ciritcal point, treat as regular interpolation as not much
           // is won if a more sophisticated interpolation scheme is used     
           if ( (static_cast<fT>(ij[2][0])*dT+tatm) > tcrit && (static_cast<fT>(ij[2][1])*dP) > pcrit ) {
               interpolation_type = regular;
               return;
             }
                      
           // otherwise if 2-phase curve is crossed identify which points are on vapor or liquid side
           bool vapor;
           if ( p < SaturationPressureFromT(t) ) vapor = true; // point in vapor field
           else                                  vapor = false; // point in liquid field
           
           // point 4 in liquid field, all other in vapor field
           if ( ( p1 < psat1 && p4 > psat1 && p2 < psat2 && p3 < psat2 ) ) {
               if ( !vapor ) {
                   // move entire grid cell up by one p increment such that point 2 is in vapor field
                   ij[0][1] = ij[2][1];
                   ij[1][0] = ij[1][1] = -1;
                   ij[3][1] = ij[2][1] = ij[0][1] + 1;
                   interpolation_type = three_pt_liquid; 
                   return;
                 }
               else {
                   // point 4 must be replaced by properties at 2-phase curve
                   ij[3][0] = ij[3][1] = -1;
                   interpolation_type = three_pt_vapor;
                   return;
                 }
             }
           // point 2 in vapor field, all other in liquid field
           if ( ( p1 > psat1 && p4 > psat1 && p2 < psat2 && p3 > psat2 ) ) {
               if ( !vapor ) {
                   // point 2 must be replaced by properties at 2-phase curve
                   ij[1][0] = ij[1][1] = -1;
                   interpolation_type = three_pt_liquid;
                   return;
                 }
               else {
                   // move entire grid cell down by one p increment such thatpoint 4 is in liquid vield
                   ij[2][1] = ij[0][1]; 
                   ij[3][0] = ij[3][1] = -1;
                   ij[0][1] = ij[1][1] = ij[2][1] -1;
                   interpolation_type = three_pt_vapor;
                   return;
                 }
             }
           // 2 points in vapor and 2 in liquid field  
           if ( p1 < psat1 && p4 > psat1 && p2 < psat2 && p3 > psat2 ) {
               if ( !vapor ) {
                   ij[0][0] = ij[0][1] = ij[1][0] = ij[1][1] = -1;
                   interpolation_type = two_pt_liquid;
                   return;
                 }
               else {
                   ij[2][0] = ij[2][1] = ij[3][0] = ij[3][1] = -1;
                   interpolation_type = two_pt_vapor;
                   return;
                 }
             }
           
         
         } // end if ( t < tcrit && p < pcrit )
       else {
           interpolation_type = regular;
           return;
         }
       
     }
   // 2. p-t at boundary or outside of lookup table     
   else {
       interpolation_type = boundary;
       if ( t < tmax && t > tatm && p >= pmax ) {
           ij[0][0] = ij[3][0] = static_cast<int32_t>((t-tatm)/dT);  
           ij[1][0] = ij[2][0] = ij[0][0]+1;  
           ij[0][1] = ij[1][1] = ij[2][1] = ij[3][1] = cols-1;
           return;
         }
       if ( t < tmax && t > tatm && p <= patm ) {
           ij[0][0] = ij[3][0] = static_cast<int32_t>((t-tatm)/dT);  
           ij[1][0] = ij[2][0] = ij[0][0]+1;  
           ij[0][1] = ij[1][1] = ij[2][1] = ij[3][1] = 0;
           // needs a check to make sure that 2-phase curve is not crossed at 1bar
           fT temp1(static_cast<fT>(ij[0][0])*dT+tatm),
              temp2(static_cast<fT>(ij[1][0])*dT+tatm);
           if ( ( temp1 > tsat_atm && temp2 > tsat_atm ) || ( temp1 < tsat_atm && temp2 < tsat_atm ) ) return;
           else {
               if ( t > tsat_atm ) ij[0][0] = -1; // vapor
               else                ij[1][0] = -1; // liquid
               return;
             }

         }
       if ( t >= tmax && p > patm && p < pmax ) {
           ij[0][1] = ij[1][1] = static_cast<int32_t>(p/dP);  
           ij[2][1] = ij[3][1] = ij[0][1]+1;  
           ij[0][0] = ij[1][0] = ij[2][0] = ij[3][0] = rows-1;
           return;
         }
       if ( t <= tatm && p > patm && p < pmax ) {
           ij[0][1] = ij[1][1] = static_cast<int32_t>(p/dP);  
           ij[2][1] = ij[3][1] = ij[0][1]+1;  
           ij[0][0] = ij[1][0] = ij[2][0] = ij[3][0] = 0;
           return;
         }
       if ( t >= tmax && p >= pmax ) {
           ij[0][1] = ij[1][1] = ij[2][1] = ij[3][1] = cols-1;
           ij[0][0] = ij[1][0] = ij[2][0] = ij[3][0] = rows-1;
           return;
         }
       if ( t <= tatm && p <= patm ) {
           ij[0][1] = ij[1][1] = ij[2][1] = ij[3][1] = 0;
           ij[0][0] = ij[1][0] = ij[2][0] = ij[3][0] = 0;
           return;
         }
       if ( t <= tatm && p >= pmax ) {
           ij[0][1] = ij[1][1] = ij[2][1] = ij[3][1] = cols-1;
           ij[0][0] = ij[1][0] = ij[2][0] = ij[3][0] = 0;
           return;
         }
       if ( t >= tmax && p <= patm ) {
           ij[0][1] = ij[1][1] = ij[2][1] = ij[3][1] = 0;
           ij[0][0] = ij[1][0] = ij[2][0] = ij[3][0] = rows-1;
           return;
         }
         
     }
 }


template<typename fT>
void H2OPropertiesLookUpTable<fT>::GetLookupTableEntriesCriticalPoint( fT t, fT p )
 {
   // 0. only if p or t has changed find new entries
   if ( p == pcurrent && t == tcurrent ) return;
   
   // 1. check that p and t are witin lookup table range
   // using a ccw numbering of points, point 1 is at t_min, p_min
   // 1.a get row and column entries for bilinear interpolation
   // first column entry -> t at first node (i1)
   ij[0][0] = static_cast<int32_t>((t-tmin_c)/dT_fine);
   // second column entry -> p at first node (j1)
   ij[0][1] = static_cast<int32_t>((p-pmin_c)/dP_fine);
   // now increment ij in a counter-clockwise fashion (see NumRec in C p. 123)
   ij[1][0] = ij[0][0]+1; // i2
   ij[1][1] = ij[0][1];   // j2
   ij[2][0] = ij[1][0];   // i3
   ij[2][1] = ij[0][1]+1; // j3
   ij[3][0] = ij[0][0];   // i4
   ij[3][1] = ij[2][1];   // j4
   // 1.b check if entries lie around 2-phase curve and properties from steam and liquid field are read in at the same time
   // note that 2-phase curve could also be crossed if t < tcrit and p < pcrit but only by a little. in this case, however,
   // use above bilinear interpolation as properties *very* close to critical point are similar on liquid and vapor side,
   // other interpolation schemes are time consuming and prost does not give the most exact values anyway.
   if ( t < tcrit && p < pcrit ) {
       fT p1(static_cast<fT>(ij[0][1])*dP_fine+pmin_c),
          p2(p1),
          p3(static_cast<fT>(ij[2][1])*dP_fine+pmin_c),
          p4(p3),
          t1(static_cast<fT>(ij[0][0])*dT_fine+tmin_c),
          t2(static_cast<fT>(ij[1][0])*dT_fine+tmin_c),
          psat1(SaturationPressureFromT(t1)),
          psat2(SaturationPressureFromT(t2));
           
       // for t at ij[0][0] and ij[1][0] test if the corresponding p is either always above
       // or always below the saturation pressure (i.e., all in liquid or all in vapor field)
       // in this case, the p-t positions identified above are o.k.
       if ( ( p1 < psat1 && p4 < psat1 && p2 < psat2 && p3 < psat2 ) || // all vapor
            ( p1 > psat1 && p4 > psat1 && p2 > psat2 && p3 > psat2 ) )  // all liquid
            {
              interpolation_type = regular;
              return;
            }
       // 4 interpolation points sourround ciritcal point, treat as regular interpolation as not much
       // is won if a more sophisticated interpolation scheme is used     
       if ( (static_cast<fT>(ij[2][0])*dT_fine+tmin_c) > tcrit && (static_cast<fT>(ij[2][1])*dP_fine+pmin_c) > pcrit ) {
           interpolation_type = regular;
           return;
         }
                      
       // otherwise if 2-phase curve is crossed identify which points are on vapor or liquid side
       bool vapor;
       if ( p < SaturationPressureFromT(t) ) vapor = true; // point in vapor field
       else                                  vapor = false; // point in liquid field
           
       // point 4 in liquid field, all other in vapor field
       if ( ( p1 < psat1 && p4 > psat1 && p2 < psat2 && p3 < psat2 ) ) {
           if ( !vapor ) {
               // move entire grid cell up by one p increment such that point 2 is in vapor field
               ij[0][1] = ij[2][1];
               ij[1][0] = ij[1][1] = -1;
               ij[3][1] = ij[2][1] = ij[0][1] + 1;
               interpolation_type = three_pt_liquid; 
               return;
             }
           else {
               // point 4 must be replaced by properties at 2-phase curve
               ij[3][0] = ij[3][1] = -1;
               interpolation_type = three_pt_vapor;
               return;
             }
         }
       // point 2 in vapor field, all other in liquid field
       if ( ( p1 > psat1 && p4 > psat1 && p2 < psat2 && p3 > psat2 ) ) {
           if ( !vapor ) {
               // point 2 must be replaced by properties at 2-phase curve
               ij[1][0] = ij[1][1] = -1;
               interpolation_type = three_pt_liquid;
               return;
             }
           else {
               // move entire grid cell down by one p increment such thatpoint 4 is in liquid vield
               ij[2][1] = ij[0][1]; 
               ij[3][0] = ij[3][1] = -1;
               ij[0][1] = ij[1][1] = ij[2][1] -1;
               interpolation_type = three_pt_vapor;
               return;
             }
         }
       // 2 points in vapor and 2 in liquid field, 2-phase curve crosses grid horizontally 
       if ( p1 < psat1 && p4 > psat1 && p2 < psat2 && p3 > psat2 ) {
           if ( !vapor ) {
               ij[0][0] = ij[0][1] = ij[1][0] = ij[1][1] = -1;
               interpolation_type = two_pt_liquid;
               return;
             }
           else {
               ij[2][0] = ij[2][1] = ij[3][0] = ij[3][1] = -1;
               interpolation_type = two_pt_vapor;
               return;
             }
         }
       // 2 points in vapor and 2 in liquid field, 2-phase curve crosses grid vertically 
       if ( p1 > psat1 && p4 > psat1 && p2 < psat2 && p3 < psat2 ) {
           if ( !vapor ) {
               ij[1][0] = ij[1][1] = ij[2][0] = ij[2][1] = -1;
               interpolation_type = two_pt_liquid;
               return;
             }
           else {
               ij[0][0] = ij[0][1] = ij[3][0] = ij[3][1] = -1;
               interpolation_type = two_pt_vapor;
               return;
             }
         }
       
         
     } // end if ( t < tcrit && p < pcrit )
   else {
     interpolation_type = regular;
     return;
   }
       
 }



template<typename fT>
void H2OPropertiesLookUpTable<fT>::GetDistances( fT t, fT p )  
 {
   fT press;
   if ( p == pcurrent && t == tcurrent ) return;
   
   if ( interpolation_type == regular || interpolation_type == three_pt_vapor || interpolation_type == three_pt_liquid ) {
       xy[0] = ( t - static_cast<fT>(ij[0][0])*dT - tatm ) / dT;
       if ( ij[2][1] != 1 ) xy[1] = ( p - static_cast<fT>(ij[0][1])*dP ) / dP;
       else                 xy[1] = ( p - patm ) / ( dP - patm );
       pcurrent = p;
       tcurrent = t;
       return;
     }
   if ( interpolation_type == boundary ) {
       if      ( t <= tatm ) xy[0] = 0.0;
       else if ( t >= tmax ) xy[0] = 1.0;
       else {
           if ( p > patm ) xy[0] = ( t - static_cast<fT>(ij[0][0])*dT - tatm ) / dT;
           else {
               if ( ij[0][0] != -1 && ij[1][0] != -1 ) xy[0] = ( t - static_cast<fT>(ij[0][0])*dT - tatm ) / dT;
               else {
                   if ( ij[0][0] == -1 ) xy[0] = ( t - tsat_atm ) / ( static_cast<fT>(ij[1][0])*dT + tatm - tsat_atm );  //vapor
                   else                  xy[0] = ( t - static_cast<fT>(ij[0][0])*dT - tatm ) / ( tsat_atm - static_cast<fT>(ij[0][0])*dT - tatm ); // liquid
                 }
             }
         }
       
       if      ( p <= patm ) xy[1] = 0.0;  
       else if ( p >= pmax ) xy[1] = 1.0;  
       else {
           if ( ij[2][1] != 1 ) xy[1] = ( p - static_cast<fT>(ij[0][1])*dP ) / dP;
           else                 xy[1] = ( p - patm ) / ( dP - patm );
         }
       pcurrent = p;
       tcurrent = t;
       return;
     }
   if ( interpolation_type == two_pt_vapor ) {
       if ( p < dP ) press = patm;
       else          press = static_cast<fT>(ij[0][1])*dP;
       xy[0] = ( t - static_cast<fT>(ij[0][0])*dT - tatm ) / dT;
       xy[1] = ( p - press ) / ( SaturationPressureFromT( static_cast<fT>(ij[1][0])*dT + tatm ) - press );
       pcurrent = p;
       tcurrent = t;
       return;
     }
   if ( interpolation_type == two_pt_liquid ) {
       press = SaturationPressureFromT( static_cast<fT>(ij[3][0])*dT + tatm );
       xy[0] = ( t - static_cast<fT>(ij[3][0])*dT - tatm ) / dT;
       xy[1] = ( p - press ) / ( static_cast<fT>(ij[3][1])*dP - press );
       pcurrent = p;
       tcurrent = t;
       return;
     }
   
   
 }


template<typename fT>
void H2OPropertiesLookUpTable<fT>::GetDistancesCriticalPoint( fT t, fT p )  
 {
   fT press, temp;
   if ( p == pcurrent && t == tcurrent ) return;
   
   if ( interpolation_type == regular || interpolation_type == three_pt_vapor || interpolation_type == three_pt_liquid ) {
       xy[0] = ( t - static_cast<fT>(ij[0][0])*dT_fine - tmin_c ) / dT_fine;
       xy[1] = ( p - static_cast<fT>(ij[0][1])*dP_fine - pmin_c ) / dP_fine;
       pcurrent = p;
       tcurrent = t;
       return;
     }
   if ( interpolation_type == two_pt_vapor ) {
       // 2-phase curve crosses horizontally
       if ( ij[0][1] != -1 ) {
           press = static_cast<fT>(ij[0][1])*dP_fine + pmin_c;
           xy[0] = ( t - static_cast<fT>(ij[0][0])*dT_fine - tmin_c ) / dT_fine;
           xy[1] = ( p - press ) / ( SaturationPressureFromT( static_cast<fT>(ij[1][0])*dT_fine + tmin_c ) - press );
         }
       // 2-phase curve crosses vertically
       else {
           press = static_cast<fT>(ij[1][1])*dP_fine + pmin_c;
           temp  = SaturationTemperatureFromP( press );
           xy[0] = ( t - temp ) / ( static_cast<fT>(ij[1][1])*dT_fine + tmin_c - temp );
           xy[1] = ( p - press ) / dP_fine; 
         }
       pcurrent = p;
       tcurrent = t;
       return;
     }
   if ( interpolation_type == two_pt_liquid ) {
       // 2-phase curve crosses horizontally
       if ( ij[2][0] != -1 ) {
           temp  = static_cast<fT>(ij[3][0])*dT_fine + tmin_c;
           press = SaturationPressureFromT( temp );
           xy[0] = ( t - temp ) / dT_fine;
           xy[1] = ( p - press ) / ( static_cast<fT>(ij[3][1])*dP_fine + pmin_c - press );
         }
       // 2-phase curve crosses vertically
       else {
           temp = static_cast<fT>(ij[0][0])*dT_fine + tmin_c;
           xy[0] = ( t - temp  ) / ( SaturationTemperatureFromP( static_cast<fT>(ij[3][1])*dP_fine + pmin_c ) - temp );
           xy[1] = ( p - static_cast<fT>(ij[0][1])*dP_fine - pmin_c ) / dP_fine;
         }
       pcurrent = p;
       tcurrent = t;
       return;
     }
   
   
 }



template<typename fT>
void H2OPropertiesLookUpTable<fT>::GetPropertyValues( property p )  
 {
   fT t, press, psat1, psat2, dist;
   for ( int32_t i=0; i<4; i++ ) {
       if ( ij[i][0] != -1 && ij[i][1] != -1 ) {
           if      ( p == cp )     prop[i] = heat_capacity( ij[i][0], ij[i][1] ); 
           else if ( p == h )      prop[i] = enthalpy( ij[i][0], ij[i][1] ); 
           else if ( p == rho )    prop[i] = density( ij[i][0], ij[i][1] ); 
           else if ( p == mu )     prop[i] = viscosity( ij[i][0], ij[i][1] ); 
           else if ( p == a )      prop[i] = expansivity( ij[i][0], ij[i][1] ); 
           else if ( p == b )      prop[i] = compressibility( ij[i][0], ij[i][1] ); 
           else if ( p == dpddCT ) prop[i] = dp_dd_CT( ij[i][0], ij[i][1] ); 
           else if ( p == dpdTCd ) prop[i] = dp_dT_Cd( ij[i][0], ij[i][1] ); 
         }
       else {
           if ( interpolation_type == two_pt_vapor ) {
               if ( i == 2 ) {
                   t = static_cast<fT>(ij[1][0])*dT+tatm;
                   if      ( p == cp )     prop[i] = VaporHeatCapacityFromT( t ); 
                   else if ( p == h )      prop[i] = VaporEnthalpyFromT( t ); 
                   else if ( p == rho )    prop[i] = VaporDensityFromT( t ); 
                   else if ( p == mu )     prop[i] = VaporViscosityFromT( t ); 
                   else if ( p == a )      prop[i] = VaporThermalExpansivityFromT( t ); 
                   else if ( p == b )      prop[i] = VaporCompressibilityFromT( t ); 
                   else if ( p == dpddCT ) prop[i] = VaporPressureDensityDerivativeFromT( t );
                   else if ( p == dpdTCd ) prop[i] = VaporPressureTemperatureDerivativeFromT( t );
                 }
               if ( i == 3 ) {
                   t = static_cast<fT>(ij[0][0])*dT+tatm;
                   if ( ij[0][1] != 0 ) press = static_cast<fT>(ij[0][1])*dP;
                   else                 press = patm;
                   psat1 = SaturationPressureFromT( t );
                   psat2 = SaturationPressureFromT( static_cast<fT>(ij[1][0])*dT+tatm );
                   
                   if      ( p == cp )       prop[i] = prop[0] + ( psat2 - press ) * ( VaporHeatCapacityFromT( t ) - prop[0] ) / ( psat1 - press ); 
                   else if ( p == h )        prop[i] = prop[0] + ( psat2 - press ) * ( VaporEnthalpyFromT( t ) - prop[0] ) / ( psat1 - press ); 
                   else if ( p == rho )      prop[i] = prop[0] + ( psat2 - press ) * ( VaporDensityFromT( t ) - prop[0] ) / ( psat1 - press ); 
                   else if ( p == mu )       prop[i] = prop[0] + ( psat2 - press ) * ( VaporViscosityFromT( t ) - prop[0] ) / ( psat1 - press ); 
                   else if ( p == a )        prop[i] = prop[0] + ( psat2 - press ) * ( VaporThermalExpansivityFromT( t ) - prop[0] ) / ( psat1 - press ); 
                   else if ( p == b )        prop[i] = prop[0] + ( psat2 - press ) * ( VaporCompressibilityFromT( t ) - prop[0] ) / ( psat1 - press ); 
                   else if ( p == dpddCT )   prop[i] = prop[0] + ( psat2 - press ) * ( VaporPressureDensityDerivativeFromT( t ) - prop[0] ) / ( psat1 - press ); 
                   else if ( p == dpdTCd )   prop[i] = prop[0] + ( psat2 - press ) * ( VaporPressureTemperatureDerivativeFromT( t ) - prop[0] ) / ( psat1 - press ); 
                 }
             }
           else if ( interpolation_type == two_pt_liquid ) {
               if ( i == 0 ) {
                   t = static_cast<fT>(ij[3][0])*dT+tatm;
                   if      ( p == cp )     prop[i] = LiquidHeatCapacityFromT( t ); 
                   else if ( p == h )      prop[i] = LiquidEnthalpyFromT( t ); 
                   else if ( p == rho )    prop[i] = LiquidDensityFromT( t ); 
                   else if ( p == mu )     prop[i] = LiquidViscosityFromT( t ); 
                   else if ( p == a )      prop[i] = LiquidThermalExpansivityFromT( t ); 
                   else if ( p == b )      prop[i] = LiquidCompressibilityFromT( t ); 
                   else if ( p == dpddCT ) prop[i] = LiquidPressureDensityDerivativeFromT( t );
                   else if ( p == dpdTCd ) prop[i] = LiquidPressureTemperatureDerivativeFromT( t );
                }
               if ( i == 1 ) {
                   t = static_cast<fT>(ij[2][0])*dT+tatm;
                   press = static_cast<fT>(ij[3][1])*dP;
                   psat1 = SaturationPressureFromT( t );
                   psat2 = SaturationPressureFromT( static_cast<fT>(ij[3][0])*dT+tatm );
                   
                   if      ( p == cp )       prop[i] = heat_capacity( ij[2][0], ij[2][1] ) - ( press - psat2 ) * ( heat_capacity( ij[2][0], ij[2][1] ) - LiquidHeatCapacityFromT( t ) ) / ( press - psat1 ); 
                   else if ( p == h )        prop[i] = enthalpy( ij[2][0], ij[2][1] ) - ( press - psat2 ) * ( enthalpy( ij[2][0], ij[2][1] ) - LiquidEnthalpyFromT( t ) ) / ( press - psat1 ); 
                   else if ( p == rho )      prop[i] = density( ij[2][0], ij[2][1] ) - ( press - psat2 ) * ( density( ij[2][0], ij[2][1] ) - LiquidDensityFromT( t ) ) / ( press - psat1 ); 
                   else if ( p == mu )       prop[i] = viscosity( ij[2][0], ij[2][1] ) - ( press - psat2 ) * ( viscosity( ij[2][0], ij[2][1] ) - LiquidViscosityFromT( t ) ) / ( press - psat1 ); 
                   else if ( p == a )        prop[i] = expansivity( ij[2][0], ij[2][1] ) - ( press - psat2 ) * ( expansivity( ij[2][0], ij[2][1] ) - LiquidThermalExpansivityFromT( t ) ) / ( press - psat1 ); 
                   else if ( p == b )        prop[i] = compressibility( ij[2][0], ij[2][1] ) - ( press - psat2 ) * ( compressibility( ij[2][0], ij[2][1] ) - LiquidCompressibilityFromT( t ) ) / ( press - psat1 ); 
                   else if ( p == dpddCT )   prop[i] = dp_dd_CT( ij[2][0], ij[2][1] ) - ( press - psat2 ) * ( dp_dd_CT( ij[2][0], ij[2][1] ) - LiquidPressureDensityDerivativeFromT( t ) ) / ( press - psat1 ); 
                   else if ( p == dpdTCd )   prop[i] = dp_dT_Cd( ij[2][0], ij[2][1] ) - ( press - psat2 ) * ( dp_dT_Cd( ij[2][0], ij[2][1] ) - LiquidPressureTemperatureDerivativeFromT( t ) ) / ( press - psat1 ); 
                 }
             }
           else if ( interpolation_type == three_pt_vapor ) {
               if ( i == 3 ) {
                   t = static_cast<fT>(ij[0][0])*dT+tatm;
                   if ( ij[2][1] != 1 ) press = static_cast<fT>(ij[0][1])*dP;
                   else                 press = patm;
                   psat1 = SaturationPressureFromT( t );
                   if ( ij[2][1] != 1 ) dist = dP;
                   else                 dist = dP - patm;
                   
                   if      ( p == cp )       prop[i] = prop[0] + dist * ( VaporHeatCapacityFromT( t ) - prop[0] ) / ( psat1 - press ); 
                   else if ( p == h )        prop[i] = prop[0] + dist * ( VaporEnthalpyFromT( t ) - prop[0] ) / ( psat1 - press ); 
                   else if ( p == rho )      prop[i] = prop[0] + dist * ( VaporDensityFromT( t ) - prop[0] ) / ( psat1 - press ); 
                   else if ( p == mu )       prop[i] = prop[0] + dist * ( VaporViscosityFromT( t ) - prop[0] ) / ( psat1 - press ); 
                   else if ( p == a )        prop[i] = prop[0] + dist * ( VaporThermalExpansivityFromT( t ) - prop[0] ) / ( psat1 - press ); 
                   else if ( p == b )        prop[i] = prop[0] + dist * ( VaporCompressibilityFromT( t ) - prop[0] ) / ( psat1 - press ); 
                   else if ( p == dpddCT )   prop[i] = prop[0] + dist * ( VaporPressureDensityDerivativeFromT( t ) - prop[0] ) / ( psat1 - press ); 
                   else if ( p == dpdTCd )   prop[i] = prop[0] + dist * ( VaporPressureTemperatureDerivativeFromT( t ) - prop[0] ) / ( psat1 - press ); 
                }
             }
           else if ( interpolation_type == three_pt_liquid ) {
               if ( i == 1 ) {
                   t = static_cast<fT>(ij[2][0])*dT+tatm;
                   press = static_cast<fT>(ij[2][1])*dP;
                   psat1 = SaturationPressureFromT( t );
                   if ( ij[2][1] != 1 ) dist = dP;
                   else                 dist = dP - patm;
                   
                   if      ( p == cp )       prop[i] = heat_capacity( ij[2][0], ij[2][1] ) - dist * ( heat_capacity( ij[2][0], ij[2][1] ) - LiquidHeatCapacityFromT( t ) ) / ( press - psat1 ); 
                   else if ( p == h )        prop[i] = enthalpy( ij[2][0], ij[2][1] ) - dist * ( enthalpy( ij[2][0], ij[2][1] ) - LiquidEnthalpyFromT( t ) ) / ( press - psat1 ); 
                   else if ( p == rho )      prop[i] = density( ij[2][0], ij[2][1] ) - dist * ( density( ij[2][0], ij[2][1] ) - LiquidDensityFromT( t ) ) / ( press - psat1 ); 
                   else if ( p == mu )       prop[i] = viscosity( ij[2][0], ij[2][1] ) - dist * ( viscosity( ij[2][0], ij[2][1] ) - LiquidViscosityFromT( t ) ) / ( press - psat1 ); 
                   else if ( p == a )        prop[i] = expansivity( ij[2][0], ij[2][1] ) - dist * ( expansivity( ij[2][0], ij[2][1] ) - LiquidThermalExpansivityFromT( t ) ) / ( press - psat1 ); 
                   else if ( p == b )        prop[i] = compressibility( ij[2][0], ij[2][1] ) - dist * ( compressibility( ij[2][0], ij[2][1] ) - LiquidCompressibilityFromT( t ) ) / ( press - psat1 ); 
                   else if ( p == dpddCT )   prop[i] = dp_dd_CT( ij[2][0], ij[2][1] ) - dist * ( dp_dd_CT( ij[2][0], ij[2][1] ) - LiquidPressureDensityDerivativeFromT( t ) ) / ( press - psat1 ); 
                   else if ( p == dpdTCd )   prop[i] = dp_dT_Cd( ij[2][0], ij[2][1] ) - dist * ( dp_dT_Cd( ij[2][0], ij[2][1] ) - LiquidPressureTemperatureDerivativeFromT( t ) ) / ( press - psat1 ); 
                 }
             }
           else if ( interpolation_type == boundary ) {
               if ( i == 0 ) {
                   if      ( p == cp )     prop[i] = VaporHeatCapacityFromP( patm ); 
                   else if ( p == h )      prop[i] = VaporEnthalpyFromP( patm ); 
                   else if ( p == rho )    prop[i] = VaporDensityFromP( patm ); 
                   else if ( p == mu )     prop[i] = VaporViscosityFromP( patm ); 
                   else if ( p == a )      prop[i] = VaporThermalExpansivityFromP( patm ); 
                   else if ( p == b )      prop[i] = VaporCompressibilityFromP( patm ); 
                   else if ( p == dpddCT ) prop[i] = VaporPressureDensityDerivativeFromP( patm );
                   else if ( p == dpdTCd ) prop[i] = VaporPressureTemperatureDerivativeFromP( patm );
                 }
               if ( i == 1 ) {
                   if      ( p == cp )     prop[i] = LiquidHeatCapacityFromP( patm ); 
                   else if ( p == h )      prop[i] = LiquidEnthalpyFromP( patm ); 
                   else if ( p == rho )    prop[i] = LiquidDensityFromP( patm ); 
                   else if ( p == mu )     prop[i] = LiquidViscosityFromP( patm ); 
                   else if ( p == a )      prop[i] = LiquidThermalExpansivityFromP( patm ); 
                   else if ( p == b )      prop[i] = LiquidCompressibilityFromP( patm ); 
                   else if ( p == dpddCT ) prop[i] = LiquidPressureDensityDerivativeFromP( patm );
                   else if ( p == dpdTCd ) prop[i] = LiquidPressureTemperatureDerivativeFromP( patm );
                 }
             }

         }
     }
 }



template<typename fT>
void H2OPropertiesLookUpTable<fT>::GetPropertyValuesCriticalPoint( property pr )  
 {
   fT t, p, press, temp, psat1, psat2, tsat1, tsat2;
   for ( int32_t i=0; i<4; i++ ) {
       if ( ij[i][0] != -1 && ij[i][1] != -1 ) {
           if      ( pr == cp )       prop[i] = crit_heat_capacity( ij[i][0], ij[i][1] ); 
           else if ( pr == h )        prop[i] = crit_enthalpy( ij[i][0], ij[i][1] ); 
           else if ( pr == rho )      prop[i] = crit_density( ij[i][0], ij[i][1] ); 
           else if ( pr == mu )       prop[i] = crit_viscosity( ij[i][0], ij[i][1] ); 
           else if ( pr == a )        prop[i] = crit_expansivity( ij[i][0], ij[i][1] ); 
           else if ( pr == b )        prop[i] = crit_compressibility( ij[i][0], ij[i][1] ); 
           else if ( pr == dpddCT )   prop[i] = crit_dp_dd_CT( ij[i][0], ij[i][1] ); 
           else if ( pr == dpdTCd )   prop[i] = crit_dp_dT_Cd( ij[i][0], ij[i][1] ); 
         }
       else {
           if ( interpolation_type == two_pt_vapor ) {
               // vertical
               if ( i == 0 ) {
                   p = static_cast<fT>(ij[1][1])*dP_fine+pmin_c;
                   if      ( pr == cp )       prop[i] = VaporHeatCapacityFromP( p ); 
                   else if ( pr == h )        prop[i] = VaporEnthalpyFromP( p ); 
                   else if ( pr == rho )      prop[i] = VaporDensityFromP( p ); 
                   else if ( pr == mu )       prop[i] = VaporViscosityFromP( p ); 
                   else if ( pr == a )        prop[i] = VaporThermalExpansivityFromP( p ); 
                   else if ( pr == b )        prop[i] = VaporCompressibilityFromP( p ); 
                   else if ( pr == dpddCT )   prop[i] = VaporPressureDensityDerivativeFromP( p ); 
                   else if ( pr == dpdTCd )   prop[i] = VaporPressureTemperatureDerivativeFromP( p ); 
                 }
               // horizontal
               if ( i == 2 ) {
                   t = static_cast<fT>(ij[1][0])*dT_fine+tmin_c;
                   if      ( pr == cp )       prop[i] = VaporHeatCapacityFromT( t ); 
                   else if ( pr == h )        prop[i] = VaporEnthalpyFromT( t ); 
                   else if ( pr == rho )      prop[i] = VaporDensityFromT( t ); 
                   else if ( pr == mu )       prop[i] = VaporViscosityFromT( t ); 
                   else if ( pr == a )        prop[i] = VaporThermalExpansivityFromT( t ); 
                   else if ( pr == b )        prop[i] = VaporCompressibilityFromT( t ); 
                   else if ( pr == dpddCT )   prop[i] = VaporPressureDensityDerivativeFromT( t ); 
                   else if ( pr == dpdTCd )   prop[i] = VaporPressureTemperatureDerivativeFromT( t ); 
                 }
               // horizontal
               if ( i == 3 && ij[2][0] == -1 ) {
                   t = static_cast<fT>(ij[0][0])*dT_fine+tmin_c;
                   press = static_cast<fT>(ij[0][1])*dP_fine+pmin_c;
                   psat1 = SaturationPressureFromT( t );
                   psat2 = SaturationPressureFromT( static_cast<fT>(ij[1][0])*dT_fine+tmin_c );
                   
                   if      ( pr == cp )       prop[i] = prop[0] + ( psat2 - press ) * ( VaporHeatCapacityFromT( t ) - prop[0] ) / ( psat1 - press ); 
                   else if ( pr == h )        prop[i] = prop[0] + ( psat2 - press ) * ( VaporEnthalpyFromT( t ) - prop[0] ) / ( psat1 - press ); 
                   else if ( pr == rho )      prop[i] = prop[0] + ( psat2 - press ) * ( VaporDensityFromT( t ) - prop[0] ) / ( psat1 - press ); 
                   else if ( pr == mu )       prop[i] = prop[0] + ( psat2 - press ) * ( VaporViscosityFromT( t ) - prop[0] ) / ( psat1 - press ); 
                   else if ( pr == a )        prop[i] = prop[0] + ( psat2 - press ) * ( VaporThermalExpansivityFromT( t ) - prop[0] ) / ( psat1 - press ); 
                   else if ( pr == b )        prop[i] = prop[0] + ( psat2 - press ) * ( VaporCompressibilityFromT( t ) - prop[0] ) / ( psat1 - press ); 
                   else if ( pr == dpddCT )   prop[i] = prop[0] + ( psat2 - press ) * ( VaporPressureDensityDerivativeFromT( t ) - prop[0] ) / ( psat1 - press ); 
                   else if ( pr == dpdTCd )   prop[i] = prop[0] + ( psat2 - press ) * ( VaporPressureTemperatureDerivativeFromT( t ) - prop[0] ) / ( psat1 - press ); 
                 }
               // vertical
               if ( i == 3 && ij[0][0] == -1 ) {
                   p = static_cast<fT>(ij[2][1])*dP_fine+pmin_c;
                   temp  = static_cast<fT>(ij[2][0])*dT_fine+tmin_c;
                   tsat1 = SaturationTemperatureFromP( p );
                   tsat2 = SaturationTemperatureFromP( static_cast<fT>(ij[1][1])*dP_fine+pmin_c );
                   
                   if      ( pr == cp )       prop[i] = prop[2] + ( temp - tsat2 ) * ( VaporHeatCapacityFromP( p ) - prop[2] ) / ( temp - tsat1 ); 
                   else if ( pr == h )        prop[i] = prop[2] + ( temp - tsat2 ) * ( VaporEnthalpyFromP( p ) - prop[2] ) / ( temp - tsat1 ); 
                   else if ( pr == rho )      prop[i] = prop[2] + ( temp - tsat2 ) * ( VaporDensityFromP( p ) - prop[2] ) / ( temp - tsat1 ); 
                   else if ( pr == mu )       prop[i] = prop[2] + ( temp - tsat2 ) * ( VaporViscosityFromP( p ) - prop[2] ) / ( temp - tsat1 ); 
                   else if ( pr == a )        prop[i] = prop[2] + ( temp - tsat2 ) * ( VaporThermalExpansivityFromP( p ) - prop[2] ) / ( temp - tsat1 ); 
                   else if ( pr == b )        prop[i] = prop[2] + ( temp - tsat2 ) * ( VaporCompressibilityFromP( p ) - prop[2] ) / ( temp - tsat1 ); 
                   else if ( pr == dpddCT )   prop[i] = prop[2] + ( temp - tsat2 ) * ( VaporPressureDensityDerivativeFromP( p ) - prop[2] ) / ( temp - tsat1 ); 
                   else if ( pr == dpdTCd )   prop[i] = prop[2] + ( temp - tsat2 ) * ( VaporPressureTemperatureDerivativeFromP( p ) - prop[2] ) / ( temp - tsat1 ); 
                 }
             }
           else if ( interpolation_type == two_pt_liquid ) {
               // horizontal
               if ( i == 0 ) {
                   t = static_cast<fT>(ij[3][0])*dT_fine+tmin_c;
                   if      ( pr == cp )       prop[i] = LiquidHeatCapacityFromT( t ); 
                   else if ( pr == h )        prop[i] = LiquidEnthalpyFromT( t ); 
                   else if ( pr == rho )      prop[i] = LiquidDensityFromT( t ); 
                   else if ( pr == mu )       prop[i] = LiquidViscosityFromT( t ); 
                   else if ( pr == a )        prop[i] = LiquidThermalExpansivityFromT( t ); 
                   else if ( pr == b )        prop[i] = LiquidCompressibilityFromT( t ); 
                   else if ( pr == dpddCT )   prop[i] = LiquidPressureDensityDerivativeFromT( t ); 
                   else if ( pr == dpdTCd )   prop[i] = LiquidPressureTemperatureDerivativeFromT( t ); 
                 }
               // horizontal
               if ( i == 1 && ij[0][0] == -1 ) {
                   t = static_cast<fT>(ij[2][0])*dT_fine+tmin_c;
                   press = static_cast<fT>(ij[3][1])*dP_fine+pmin_c;
                   psat1 = SaturationPressureFromT( t );
                   psat2 = SaturationPressureFromT( static_cast<fT>(ij[3][0])*dT_fine+tmin_c );
                   
                   if      ( pr == cp )       prop[i] = crit_heat_capacity( ij[2][0], ij[2][1] ) - ( press - psat2 ) * ( crit_heat_capacity( ij[2][0], ij[2][1] ) - LiquidHeatCapacityFromT( t ) ) / ( press - psat1 ); 
                   else if ( pr == h )        prop[i] = crit_enthalpy( ij[2][0], ij[2][1] ) - ( press - psat2 ) * ( crit_enthalpy( ij[2][0], ij[2][1] ) - LiquidEnthalpyFromT( t ) ) / ( press - psat1 ); 
                   else if ( pr == rho )      prop[i] = crit_density( ij[2][0], ij[2][1] ) - ( press - psat2 ) * ( crit_density( ij[2][0], ij[2][1] ) - LiquidDensityFromT( t ) ) / ( press - psat1 ); 
                   else if ( pr == mu )       prop[i] = crit_viscosity( ij[2][0], ij[2][1] ) - ( press - psat2 ) * ( crit_viscosity( ij[2][0], ij[2][1] ) - LiquidViscosityFromT( t ) ) / ( press - psat1 ); 
                   else if ( pr == a )        prop[i] = crit_expansivity( ij[2][0], ij[2][1] ) - ( press - psat2 ) * ( crit_expansivity( ij[2][0], ij[2][1] ) - LiquidThermalExpansivityFromT( t ) ) / ( press - psat1 ); 
                   else if ( pr == b )        prop[i] = crit_compressibility( ij[2][0], ij[2][1] ) - ( press - psat2 ) * ( crit_compressibility( ij[2][0], ij[2][1] ) - LiquidCompressibilityFromT( t ) ) / ( press - psat1 ); 
                   else if ( pr == dpddCT )   prop[i] = crit_dp_dd_CT( ij[2][0], ij[2][1] ) - ( press - psat2 ) * ( crit_dp_dd_CT( ij[2][0], ij[2][1] ) - LiquidPressureDensityDerivativeFromT( t ) ) / ( press - psat1 ); 
                   else if ( pr == dpdTCd )   prop[i] = crit_dp_dT_Cd( ij[2][0], ij[2][1] ) - ( press - psat2 ) * ( crit_dp_dT_Cd( ij[2][0], ij[2][1] ) - LiquidPressureTemperatureDerivativeFromT( t ) ) / ( press - psat1 ); 
                 }
               // vertical
               if ( i == 1 && ij[2][0] == -1 ) {
                   p = static_cast<fT>(ij[0][1])*dP_fine+pmin_c;
                   temp = static_cast<fT>(ij[0][0])*dT_fine+tmin_c;
                   tsat1 = SaturationTemperatureFromP( p );
                   tsat2 = SaturationTemperatureFromP( static_cast<fT>(ij[3][1])*dP_fine+pmin_c );
                   if      ( pr == cp )       prop[i] = prop[0] - ( tsat2 - temp ) * ( prop[0] - LiquidHeatCapacityFromP( p ) ) / ( tsat1 - temp  ); 
                   else if ( pr == h )        prop[i] = prop[0] - ( tsat2 - temp ) * ( prop[0] - LiquidEnthalpyFromP( p ) ) / ( tsat1 - temp ); 
                   else if ( pr == rho )      prop[i] = prop[0] - ( tsat2 - temp ) * ( prop[0] - LiquidDensityFromP( p ) ) / ( tsat1 - temp  ); 
                   else if ( pr == mu )       prop[i] = prop[0] - ( tsat2 - temp ) * ( prop[0] - LiquidViscosityFromP( p ) ) / ( tsat1 - temp ); 
                   else if ( pr == a )        prop[i] = prop[0] - ( tsat2 - temp ) * ( prop[0] - LiquidThermalExpansivityFromP( p ) ) / ( tsat1 - temp ); 
                   else if ( pr == b )        prop[i] = prop[0] - ( tsat2 - temp ) * ( prop[0] - LiquidCompressibilityFromP( p ) ) / ( tsat1 - temp ); 
                   else if ( pr == dpddCT )   prop[i] = prop[0] - ( tsat2 - temp ) * ( prop[0] - LiquidPressureDensityDerivativeFromP( p ) ) / ( tsat1 - temp ); 
                   else if ( pr == dpdTCd )   prop[i] = prop[0] - ( tsat2 - temp ) * ( prop[0] - LiquidPressureTemperatureDerivativeFromP( p ) ) / ( tsat1 - temp ); 
                 }
               // vertical
               if ( i == 2 ) {
                   p = static_cast<fT>(ij[3][1])*dP_fine+pmin_c;
                   if      ( pr == cp )       prop[i] = LiquidHeatCapacityFromP( p ); 
                   else if ( pr == h )        prop[i] = LiquidEnthalpyFromP( p ); 
                   else if ( pr == rho )      prop[i] = LiquidDensityFromP( p ); 
                   else if ( pr == mu )       prop[i] = LiquidViscosityFromP( p ); 
                   else if ( pr == a )        prop[i] = LiquidThermalExpansivityFromP( p ); 
                   else if ( pr == b )        prop[i] = LiquidCompressibilityFromP( p ); 
                   else if ( pr == dpddCT )   prop[i] = LiquidPressureDensityDerivativeFromP( p ); 
                   else if ( pr == dpdTCd )   prop[i] = LiquidPressureTemperatureDerivativeFromP( p ); 
                 }
             }
           else if ( interpolation_type == three_pt_vapor ) {
               if ( i == 3 ) {
                   t = static_cast<fT>(ij[0][0])*dT_fine+tmin_c;
                   press = static_cast<fT>(ij[0][1])*dP_fine+pmin_c;
                   psat1 = SaturationPressureFromT( t );
                   
                   if      ( pr == cp )       prop[i] = prop[0] + dP_fine * ( VaporHeatCapacityFromT( t ) - prop[0] ) / ( psat1 - press ); 
                   else if ( pr == h )        prop[i] = prop[0] + dP_fine * ( VaporEnthalpyFromT( t ) - prop[0] ) / ( psat1 - press ); 
                   else if ( pr == rho )      prop[i] = prop[0] + dP_fine * ( VaporDensityFromT( t ) - prop[0] ) / ( psat1 - press ); 
                   else if ( pr == mu )       prop[i] = prop[0] + dP_fine * ( VaporViscosityFromT( t ) - prop[0] ) / ( psat1 - press ); 
                   else if ( pr == a )        prop[i] = prop[0] + dP_fine * ( VaporThermalExpansivityFromT( t ) - prop[0] ) / ( psat1 - press ); 
                   else if ( pr == b )        prop[i] = prop[0] + dP_fine * ( VaporCompressibilityFromT( t ) - prop[0] ) / ( psat1 - press ); 
                   else if ( pr == dpddCT )   prop[i] = prop[0] + dP_fine * ( VaporPressureDensityDerivativeFromT( t ) - prop[0] ) / ( psat1 - press ); 
                   else if ( pr == dpdTCd )   prop[i] = prop[0] + dP_fine * ( VaporPressureTemperatureDerivativeFromT( t ) - prop[0] ) / ( psat1 - press ); 
                }
             }
           else if ( interpolation_type == three_pt_liquid ) {
               if ( i == 1 ) {
                   t = static_cast<fT>(ij[2][0])*dT_fine+tmin_c;
                   press = static_cast<fT>(ij[2][1])*dP_fine+pmin_c;
                   psat1 = SaturationPressureFromT( t );

                   if      ( pr == cp )       prop[i] = crit_heat_capacity( ij[2][0], ij[2][1] ) - dP_fine * ( crit_heat_capacity( ij[2][0], ij[2][1] ) - LiquidHeatCapacityFromT( t ) ) / ( press - psat1 ); 
                   else if ( pr == h )        prop[i] = crit_enthalpy( ij[2][0], ij[2][1] ) - dP_fine * ( crit_enthalpy( ij[2][0], ij[2][1] ) - LiquidEnthalpyFromT( t ) ) / ( press - psat1 ); 
                   else if ( pr == rho )      prop[i] = crit_density( ij[2][0], ij[2][1] ) - dP_fine * ( crit_density( ij[2][0], ij[2][1] ) - LiquidDensityFromT( t ) ) / ( press - psat1 ); 
                   else if ( pr == mu )       prop[i] = crit_viscosity( ij[2][0], ij[2][1] ) - dP_fine * ( crit_viscosity( ij[2][0], ij[2][1] ) - LiquidViscosityFromT( t ) ) / ( press - psat1 ); 
                   else if ( pr == a )        prop[i] = crit_expansivity( ij[2][0], ij[2][1] ) - dP_fine * ( crit_expansivity( ij[2][0], ij[2][1] ) - LiquidThermalExpansivityFromT( t ) ) / ( press - psat1 ); 
                   else if ( pr == b )        prop[i] = crit_compressibility( ij[2][0], ij[2][1] ) - dP_fine * ( crit_compressibility( ij[2][0], ij[2][1] ) - LiquidCompressibilityFromT( t ) ) / ( press - psat1 ); 
                   else if ( pr == dpddCT )   prop[i] = crit_dp_dd_CT( ij[2][0], ij[2][1] ) - dP_fine * ( crit_dp_dd_CT( ij[2][0], ij[2][1] ) - LiquidPressureDensityDerivativeFromT( t ) ) / ( press - psat1 ); 
                   else if ( pr == dpdTCd )   prop[i] = crit_dp_dT_Cd( ij[2][0], ij[2][1] ) - dP_fine * ( crit_dp_dT_Cd( ij[2][0], ij[2][1] ) - LiquidPressureTemperatureDerivativeFromT( t ) ) / ( press - psat1 ); 
                 }
             }

         }
     }
 }




template<typename fT>
fT H2OPropertiesLookUpTable<fT>::NumericalPressureDensityDerivative( fT r1, fT t, fT p, fT dpress )  
 {
    Prop *properties;
    properties = newProp('T', 'p', 1);

    fT r2;
    water_tp( t, p+dpress, d, dp, properties );
    CheckError( properties );
    r2 = properties->d;
    freeProp( properties );
    return ( dpress / ( r2 - r1 ) ); 


 }   

template<typename fT>
fT H2OPropertiesLookUpTable<fT>::NumericalCompressibility( fT r1, fT t, fT p, fT dpress )  
 {
    Prop *properties;
    properties = newProp('T', 'p', 1);
    
    fT r2, invr;
    water_tp( t, p+dpress, d, dp, properties );
    CheckError( properties );
    r2 = properties->d;
    if ( r1 > r2 ) invr = 1.0/r1; // use larger density
    else           invr = 1.0/r2;
    freeProp( properties );
    return ( invr * ( r2 - r1 ) / dpress ); 
 }   


template<typename fT>
fT H2OPropertiesLookUpTable<fT>::NumericalHeatCapacity( fT h1, fT t, fT p, fT dtemp )  
 {
    Prop *properties;
    properties = newProp('T', 'p', 1);
    
    fT h2;
    water_tp( t+dtemp, p, d, dp, properties );
    CheckError( properties );
    h2 = properties->h;
    freeProp( properties );
    return ( ( h2 - h1 ) / dtemp ); 
 }   



template<typename fT>
void H2OPropertiesLookUpTable<fT>::CheckError( Prop* test_property )  const   
 {
     if ( test_property->error == 1 ) {
         cout << "\nH2OPropertiesLookUpTable<fT>::CheckError: Error in PROST, duming output " << endl; 
         dumpProp( stdout, test_property );
         cout << endl;
       }
 }


template<typename fT>
fT H2OPropertiesLookUpTable<fT>::DynamicViscosity( fT t, fT input_density )  const  
 {
   const double tstar(647.27), rhostar(317.763);
   double trat,trat1,rhorat,rhorat1,n0,n;
   int i,j,k;

   trat    = (t+kelvin)/tstar;
   trat1   = 1.0/trat-1.0;
   rhorat  = input_density/rhostar;
   rhorat1 = rhorat-1.0;

   n0=0;
   for(k=0;k<4;++k){
     n0 += (ak[k]*pow((1.0/trat),k));
   }  
   n0 = 1.0/n0;
   n0 *= (sqrt(trat));

   n=0;
   for(i=0;i<6;++i){
     for(j=0;j<5;++j){
       n += (bij[i][j]*pow(trat1,i)*pow(rhorat1,j));
     }
   }
   n = n0*exp(rhorat*n);

   return n*1.0e-06; // converts to kg/s/m
}


template class H2OPropertiesLookUpTable<double>;

} // end name space csmp
