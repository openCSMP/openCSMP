#include "IAPWS_H2OPropertiesVisitor.h"
#include "Model.h"
#include "Element.h"
#include "Exception.h"
#include <sstream>

using namespace std;

namespace csmp {

template<uint32_t dim>
IAPWS_H2OPropertiesVisitor<dim>::IAPWS_H2OPropertiesVisitor( Model<dim>& sg,
                                                             const char* fluidPressure,
                                                             const char* fluidDensity,
                                                             const char* fluidViscosity,
                                                             bool with_steam, bool verbose )
    : Visitor<dim>( MODEL, NODE ),
      verbose(verbose),
      with_steam(with_steam),
      kelvin(273.15),
      bar(1.0e+05),
      dt(1.0e-02),
      dp(1.0e-07), // convergence criterion for PROST
      // getting csmp::Index values for the involved Node and integration point variables
      P_key(sg.Database().StorageKey(fluidPressure)),
      T_key(sg.Database().StorageKey("temperature")),
      rho_key(sg.Database().StorageKey(fluidDensity)),
      cp_key(sg.Database().StorageKey("fluid heat capacity")),
      h_key(sg.Database().StorageKey("fluid enthalpy")),
      mu_key(sg.Database().StorageKey(fluidViscosity)),
      alpha_key(sg.Database().StorageKey("fluid expansivity")),
      beta_key(sg.Database().StorageKey("fluid compressibility")),
      s_key((with_steam==true)? sg.Database().StorageKey("saturation steam") : csmp::Index())
{
    if (rho_key.place == ELEMENT_INTEGRATION_POINT ) Visitor<dim>::ApplicationTarget(ELEMENT);

    this->PlacementChecks();
    // allocating the memory for the PROST output structures
    props = newProp('p', 't', 2);

    if ( with_steam ) {
        sprops = newProp('p', 't', 2);
        lprops = newProp('p', 't', 2);
    }
}




template<uint32_t dim>
IAPWS_H2OPropertiesVisitor<dim>::IAPWS_H2OPropertiesVisitor( Model<dim>& sg,
                                                             Index& fluidPressure,
                                                             Index& fluidTemperature,
                                                             Index& fluidDensity,
                                                             Index& fluidViscosity,
                                                             Index& fluidEnthalpy,
                                                             Index& fluidHeatCapacity,
                                                             Index& fluidExpansivity,
                                                             Index& fluidCompressibility,
                                                             bool with_steam,
                                                             bool verbose )
    : Visitor<dim>( MODEL, NODE ),
      verbose(verbose),
      with_steam(with_steam),
      kelvin(273.15),
      bar(1.0e+05),
      dt(1.0e-02),
      dp(1.0e-07), // convergence criterion for PROST
      // getting csmp::Index values for the involved Node and integration point variables
      P_key(fluidPressure),
      T_key(fluidTemperature),
      rho_key(fluidDensity),
      cp_key(fluidHeatCapacity),
      h_key(fluidEnthalpy),
      mu_key(fluidViscosity),
      alpha_key(fluidExpansivity),
      beta_key(fluidCompressibility),
      s_key((with_steam==true)? sg.Database().StorageKey("saturation steam") : csmp::Index())
{
    if (rho_key.place == ELEMENT_INTEGRATION_POINT ) Visitor<dim>::ApplicationTarget(ELEMENT);
    this->PlacementChecks();
    // allocating the memory for the PROST output structures
    props = newProp('p', 't', 2);

    if ( with_steam ) {
        sprops = newProp('p', 't', 2);
        lprops = newProp('p', 't', 2);
    }
}

template<uint32_t dim>
IAPWS_H2OPropertiesVisitor<dim>::IAPWS_H2OPropertiesVisitor( Model<dim>& sg,
                                                             const char* fluidPressure,
                                                             const char* fluidTemperature,
                                                             const char* fluidDensity,
                                                             const char* fluidViscosity,
                                                             const char* fluidEnthalpy,
                                                             const char* fluidHeatCapacity,
                                                             const char* fluidExpansivity,
                                                             const char* fluidCompressibility,
                                                             bool with_steam,
                                                             bool verbose )
    : Visitor<dim>( MODEL, NODE ),
      verbose(verbose),
      with_steam(with_steam),
      kelvin(273.15),
      bar(1.0e+05),
      dt(1.0e-02),
      dp(1.0e-07), // convergence criterion for PROST
      // getting csmp::Index values for the involved Node and integration point variables
      P_key(sg.Database().StorageKey(fluidPressure)),
      T_key(sg.Database().StorageKey(fluidTemperature)),
      rho_key(sg.Database().StorageKey(fluidDensity)),
      cp_key(sg.Database().StorageKey(fluidHeatCapacity)),
      h_key(sg.Database().StorageKey(fluidEnthalpy)),
      mu_key(sg.Database().StorageKey(fluidViscosity)),
      alpha_key(sg.Database().StorageKey(fluidExpansivity)),
      beta_key(sg.Database().StorageKey(fluidCompressibility)),
      s_key((with_steam==true)? sg.Database().StorageKey("saturation steam") : csmp::Index())
{
    if (rho_key.place == ELEMENT_INTEGRATION_POINT ) Visitor<dim>::ApplicationTarget(ELEMENT);
    this->PlacementChecks();
    // allocating the memory for the PROST output structures
    props = newProp('p', 't', 2);

    if ( with_steam ) {
        sprops = newProp('p', 't', 2);
        lprops = newProp('p', 't', 2);
    }
}



/**
     For output of only the fluid density and the fluid viscosity.
     The other variables are computed but not reported.
     
     @attention Variable output is controlled by the variable keys.
     If variables are not initialized (place==UNDEFINED), they will not be output.
     
     @attention the only extra variable that is used if steam is enabled,
     is the 'saturation steam'.
*/
template<uint32_t dim>
IAPWS_H2OPropertiesVisitor<dim>::IAPWS_H2OPropertiesVisitor( Model<dim>& sg,
                                                             const char* temperature,
                                                             const char* fluidPressure,
                                                             const char* fluidDensity,
                                                             const char* fluidViscosity,
                                                             bool with_steam )
    : Visitor<dim>( MODEL, NODE ),
      verbose(false),
      with_steam(with_steam),
      kelvin(273.15),
      bar(1.0e+05),
      dt(1.0e-02),
      dp(1.0e-07), // convergence criterion for PROST
      // getting csmp::Index values for the involved Node and integration point variables
      P_key(sg.Database().StorageKey(fluidPressure)),
      T_key(sg.Database().StorageKey(temperature)),
      rho_key(sg.Database().StorageKey(fluidDensity)),
      mu_key(sg.Database().StorageKey(fluidViscosity)),
      s_key((with_steam==true)? sg.Database().StorageKey("saturation steam") : csmp::Index())
{
    if (rho_key.place == ELEMENT_INTEGRATION_POINT ) Visitor<dim>::ApplicationTarget(ELEMENT);
   // no placement checks here
    // allocating the memory for the PROST output structures
    props = newProp('p', 't', 2);

    if ( with_steam ) {
        sprops = newProp('p', 't', 2);
        lprops = newProp('p', 't', 2);
    }
}






template<uint32_t dim>
void IAPWS_H2OPropertiesVisitor<dim>::PlacementChecks()
{
    // testing that the variables exist on the right place and are of right type
    if ( P_key.place != NODE or P_key.type != SCALAR )
        throw csmp::Exception( ERROR, "IAPWS_H2OPropertiesVisitor (constructor)",
                               "'fluid pressure' must be a scalar node variable" );
    if ( T_key.place != NODE or T_key.type != SCALAR )
        throw csmp::Exception( ERROR, "IAPWS_H2OPropertiesVisitor (constructor)",
                               "'temperature' must be a scalar node variable" );
    if (rho_key.place == NODE)
        cout<<"\nIAPWS_H2OPropertiesVisitor<dim>::fluid density is placed at the NODE!"<<endl;
    if (rho_key.place == ELEMENT_INTEGRATION_POINT)
        cout<<"\nIAPWS_H2OPropertiesVisitor<dim>::fluid density is placed at the INTEGRATION POINT!"<<endl;

    if ( !(rho_key.place == NODE or rho_key.place == ELEMENT_INTEGRATION_POINT) or rho_key.type != SCALAR )
        throw csmp::Exception( ERROR, "IAPWS_H2OPropertiesVisitor (constructor)",
                               "'fluid density' must be a scalar node or integration point variable" );
    if ( !(mu_key.place == NODE or mu_key.place == ELEMENT_INTEGRATION_POINT) or mu_key.type != SCALAR )
        throw csmp::Exception( ERROR, "IAPWS_H2OPropertiesVisitor (constructor)",
                               "'fluid viscosity' must be a scalar node or integration point variable" );
    if ( cp_key.place != NODE or cp_key.type != SCALAR )
        throw csmp::Exception( ERROR, "IAPWS_H2OPropertiesVisitor (constructor)",
                               "'fluid heat capacity' must be a scalar node variable" );
    if ( h_key.place != NODE or h_key.type != SCALAR )
        throw csmp::Exception( ERROR, "IAPWS_H2OPropertiesVisitor (constructor)",
                               "'fluid enthalpy' must be a scalar node variable" );
    if ( beta_key.place != NODE or beta_key.type != SCALAR )
        throw csmp::Exception( ERROR, "IAPWS_H2OPropertiesVisitor (constructor)",
                               "'fluid compressibility' must be a scalar node variable" );
    if ( alpha_key.place != NODE or alpha_key.type != SCALAR )
        throw csmp::Exception( ERROR, "IAPWS_H2OPropertiesVisitor (constructor)",
                               "'fluid expansivity' must be a scalar node variable" );
}



template<uint32_t dim>
void IAPWS_H2OPropertiesVisitor<dim>::Verbose() { verbose = (verbose==false) ? true : false; }



template<uint32_t dim>
IAPWS_H2OPropertiesVisitor<dim>::~IAPWS_H2OPropertiesVisitor() 
{
}



template<uint32_t dim>
void IAPWS_H2OPropertiesVisitor<dim>::Visit( Node<dim>* n ) 
{
    // 1. reading input variables
    n->Read( T_key, Tf );
    n->Read( P_key, Pf );
    
    if ( Tf() < 0.0 || Pf() < bar || Tf() > 1500.0 ) {
        ostringstream t, p, id;
        t << Tf();
        p << Pf();
        id << n->Idx();
        string errormessage;
        errormessage  = "\nErratic temperature and/or pressure at node ";
        errormessage += id.str();
        errormessage += "\nValues are T = ";
        errormessage += t.str();
        errormessage += " oC and P = ";
        errormessage += p.str();
        errormessage += " Pa ";
        throw csmp::Exception( ERROR, "\nIAPWS_H2OPropertiesVisitor<fT>::Visit:",
                               errormessage.c_str() );
    }

    // 2. converting from oC to K for PROST
    temk = Tf() + kelvin; // oC -> K
    
    // 3. computing fluid properties using PROST
    dens = 100.0;
    water_tp(temk, Pf(), dens, dp, props);
    
    rho  = props->d;
    cp   = props->cp;
    h    = props->h;
    beta = 1.0/(props->dp->d_CT*props->d);
    mu = viscosity.ViscosityFromTemperatureAndDensity( Tf(), props->d );

    if ( props->error == 1 ) {
        cout << endl;
        dumpProp( stdout, props );
        cout << endl;
        failed = true;
        throw csmp::Exception( ERROR, "IAPWS_H2OPropertiesVisitor<fT>::Visit:",
                               "Error in PROST, duming output, nothing is written to the nodes...");
    }
    else {
        failed = false;
        // 4. store variables
        n->Store( rho_key,  rho ); // always (nodal) fluid density
        if ( cp_key.place != UNDEFINED ) n->Store( cp_key,   cp );    // heat capacity
        if ( h_key.place  != UNDEFINED ) n->Store( h_key,    h );     // enthalpy
        n->Store( mu_key,   mu ); // always viscosity
        if ( beta_key.place != UNDEFINED ) n->Store( beta_key, beta );  // compressibility


        // 5. if steam fraction shall be calculated, check if fluid is liquid or vapor
        if ( with_steam ) {
            sat_p( Pf(), lprops, sprops );
            // liquid phase
            if      ( temk < lprops->T ) x = 0.0;
            // vapor phase
            else                         x = 1.0;
            n->Store( s_key, x );
        }
        
    }
    // 5. compute thermal expansivity
    dens = 100.0;
    water_tp(temk+dt, Pf(), dens, dp, props );
    if ( props->error == 1 && !failed ) {
        cout << endl;
        dumpProp( stdout, props );
        cout << endl;
        alpha = 0.0;
        throw csmp::Exception( ERROR, "IAPWS_H2OPropertiesVisitor<fT>::Visit:",
                               "Error in PROST, duming output, expansivity is set to zero...");
    }
    
    else {
        alpha  = props->d - rho();
        alpha /= -rho();
        alpha /= dt;
    }
    
    if ( alpha_key.place != UNDEFINED ) n->Store( alpha_key, alpha );  // expansivity
    
    // 6. output if desired
    if ( verbose && props->error != 1 )
    {
        cout <<"\nFluid properties computed with tsteam():"<< endl;
        cout.width(35); cout <<"fluid density (kg m-3): "            << rho()        << endl;
        cout.width(35); cout <<"heat capacity of fluid (J kg K-1): " << cp()         << endl;
        cout.width(35); cout <<"fluid enthalpy (J kg): "             << h()          << endl;
        cout.width(35); cout <<"dynamic viscosity (Pa s-1): "        << mu()         << endl;
        cout.width(35); cout <<"fluid compressibility (Pa-1): "      << beta()         << endl;
        cout.width(35); cout <<"fluid expansivity (oC-1): "          << alpha()         << endl;
        if ( with_steam ) {
            cout.width(35); cout <<"steam fraction (-): "            << x()     << endl;
        }
        cout.flush();
    }
}




template<uint32_t dim>
void IAPWS_H2OPropertiesVisitor<dim>::Visit( Element<dim>* e )
{
    for ( auto i{0U}; i<e->IntegrationPoints(); i++ )
    {
        Tf=e->PropertyValueAtIntegrationPoint(T_key,i);
        Pf=e->PropertyValueAtIntegrationPoint(P_key,i);

        //perform checks
        if ( Tf() < 0.0 || Pf() < bar || Tf() > 1500.0 ) {
            ostringstream t, p, id;
            t << Tf();
            p << Pf();
            id << e->Idx();
            string errormessage;
            errormessage  = "\nErratic temperature and/or pressure at element ";
            errormessage += id.str();
            errormessage += "\nValues are T = ";
            errormessage += t.str();
            errormessage += " oC and P = ";
            errormessage += p.str();
            errormessage += " Pa ";
            throw csmp::Exception( ERROR, "\nIAPWS_H2OPropertiesVisitor<fT>::Visit:",
                                   errormessage.c_str() );
        }

        // 2. converting from oC to K for PROST
        temk = Tf() + kelvin; // oC -> K

        // 3. computing fluid properties using PROST
        dens = 100.0;
        water_tp(temk, Pf(), dens, dp, props);

        rho  = props->d;
        cp   = props->cp;
        h    = props->h;
        beta = 1.0/(props->dp->d_CT*props->d);
        mu = viscosity.ViscosityFromTemperatureAndDensity( Tf(), props->d );


        if ( props->error == 1 ) {
            cout << endl;
            dumpProp( stdout, props );
            cout << endl;
            failed = true;
            throw csmp::Exception( ERROR, "IAPWS_H2OPropertiesVisitor<fT>::Visit:",
                                   "Error in PROST, duming output, nothing is written to the nodes...");
        }
        else {
            failed = false;
            // 4. store variables
            e->Store(i,rho_key,rho);
            e->Store(i,mu_key,mu);

            // if steam fraction shall be calculated, check if fluid is liquid or vapor
            if ( with_steam ) {
                sat_p( Pf(), lprops, sprops );
                // liquid phase
                if      ( temk < lprops->T ) x = 0.0;
                // vapor phase
                else                         x = 1.0;
                //n->Store( s_key, x );
            }
        }


        // 5. compute thermal expansivity
        dens = 100.0;
        water_tp(temk+dt, Pf(), dens, dp, props );
        if ( props->error == 1 && !failed ) {
            cout << endl;
            dumpProp( stdout, props );
            cout << endl;
            alpha = 0.0;
            throw csmp::Exception( ERROR, "IAPWS_H2OPropertiesVisitor<fT>::Visit:",
                                   "Error in PROST, duming output, expansivity is set to zero...");
        }

        else {
            alpha  = props->d - rho();
            alpha /= -rho();
            alpha /= dt;
        }

        if ( alpha_key.place != UNDEFINED ) e->Store( alpha_key, alpha );  // thermal expansivity


        // 6. output if desired
        if ( verbose && props->error != 1 )
        {
            cout <<"\nFluid properties computed with tsteam():"<< endl;
            cout.width(35); cout <<"fluid density (kg m-3): "            << rho()        << endl;
            cout.width(35); cout <<"heat capacity of fluid (J kg K-1): " << cp()         << endl;
            cout.width(35); cout <<"fluid enthalpy (J kg): "             << h()          << endl;
            cout.width(35); cout <<"dynamic viscosity (Pa s-1): "        << mu()         << endl;
            cout.width(35); cout <<"fluid compressibility (Pa-1): "      << beta()         << endl;
            cout.width(35); cout <<"fluid expansivity (oC-1): "          << alpha()         << endl;
            if ( with_steam ) {
                cout.width(35); cout <<"steam fraction (-): "            << x()     << endl;
            }
            cout.flush();
        }
    }
}



template class IAPWS_H2OPropertiesVisitor<1U>;
template class IAPWS_H2OPropertiesVisitor<2U>;
template class IAPWS_H2OPropertiesVisitor<3U>;



} 


































