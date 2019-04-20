#include "PhaseStateFinder_H2O_CO2_NaCl.h"
#include "HaliteLiquidus.h"
#include "EOS_CO2H2ONaCl_Spycher2004.h"
#include "ConvertConcentrationUnitsNaCl.h"
#include "ErrorHandler.h"

// #define PHASESTATEFINDER_H2O_CO2_NACL_DEBUG

#include <limits>

// SKM FIX - for using EOS as in library
namespace csmp {

PhaseStateFinder_H2O_CO2_NaCl::PhaseStateFinder_H2O_CO2_NaCl( const double64& t, const double64& p, const double64& total_mass,
                                                              const double64& massfraction_h2o, const double64& massfraction_co2,
                                                              const double64& massfraction_nacl ) :
    t(t),
    p(p),
    total_mass(total_mass),
    bulk_massfraction_h2o(massfraction_h2o),
    bulk_massfraction_co2(massfraction_co2),
    bulk_massfraction_nacl(massfraction_nacl),
    mystate(SYSTEM_STATE::undefined),
    liquidus(t,p)
{
}





SYSTEM_STATE PhaseStateFinder_H2O_CO2_NaCl::Equilibrate()
{
    // 1. Basic conversions for internal use
    // a) mass to mole fractions of bulk composition and check for consistency
    bulk_xh2o   = bulk_massfraction_h2o/molar_mass_h2o;
    bulk_xco2   = bulk_massfraction_co2/molar_mass_co2;
    bulk_xnacl  = bulk_massfraction_nacl/molar_mass_nacl;
    double64 molenorm = bulk_xh2o + bulk_xco2 + bulk_xnacl;
    bulk_xh2o  /= molenorm;
    bulk_xco2  /= molenorm;
    bulk_xnacl /= molenorm;
    if(fabs(bulk_xh2o+bulk_xco2+bulk_xnacl-1.) > 2.*std::numeric_limits<double64>::epsilon())
    {
        cerr << "Error PhaseStateFinder_H2O_CO2_NaCl::Equilibrate() :\n";
        cerr << "Some of mole fractions of bulk composition differs from 1 beyond numeric precision.\n";
        cerr << "Actual difference is " << bulk_xh2o+bulk_xco2+bulk_xnacl-1. << endl;
        cerr << "Please enter any character to continue or terminate the program : "; char mychar; cin >> mychar;
    }
    // b) Cartesian coordinates of bulk composition: x_bulkcomp, y_bulkcomp
    ConvertMolefractionToCartesian(x_bulkcomp, y_bulkcomp, bulk_xh2o, bulk_xco2, bulk_xnacl);


    // 2. Finding the phase state. In this version we simply check each state one after the other
    //    until the correct one has been found. In that case, fluid properties etc. are computed and
    //    the job is done.
    //
    //    Todo for speed-up?:
    //    For fastest state finding, starting with the previous state as guess and
    //    then using the info obtained in EvaluateFractionsInTriangularRegion()
    //    seems a good candidate for the fastest route. For now I did not put
    //    any priority on this.

    bool equilibrated{false};

    // some convenient variables
    double64 molar_mass_aq, molar_mass_carb;
    double64 massnorm;
    double64 moleFractionOfCarbonicPhase;
    double64 moleFractionOfAqueousPhase;

    // 0. Maybe we should  perform also checks for pure phases H2O, CO2, NaCl for completeness
    // 1. Salt-free system, includes carbonic phase, H2O-CO2, and aqueous phase
//    cout << "Testing salt-free ...\n";
    if( fabs(bulk_xnacl) < 2.*std::numeric_limits<double64>::epsilon() ) // my best guess of how to check for "essentially zero" salt
    {
        double64 xco2  = eos.x_Co2(p,t,0.0);
        double64 yh2o  = eos.y_H2o(p,t,0.0);
        if(bulk_xh2o <= yh2o) // carbonic phase only, all h2o dissolved in there
        {
            mystate   = SYSTEM_STATE::carb;
            massAqueousPhase_  = 0.0;
            massCarbonicPhase_ = total_mass;
            massHalite_        = 0.0;
            rho_carb_  = eos.Rho_CarbonicPhase(p,t);
            mu_carb_   = eos.mu_CarbonicPhase(p,t);
            beta_carb_ = eos.C_CarbonicPhase(p,t);
            Y_h2o_     = bulk_xh2o*molar_mass_h2o;
            Y_co2_     = bulk_xco2*molar_mass_co2;
            massnorm   = Y_h2o_ + Y_co2_;
            Y_h2o_    /= massnorm;
            Y_co2_     = 1.-Y_h2o_;
            // set all aqueous properties to zero
            X_h2o_     = std::numeric_limits<double64>::quiet_NaN();
            X_co2_     = std::numeric_limits<double64>::quiet_NaN();
            X_nacl_    = std::numeric_limits<double64>::quiet_NaN();
            rho_aq_    = std::numeric_limits<double64>::quiet_NaN();
            mu_aq_     = std::numeric_limits<double64>::quiet_NaN();
            beta_aq_   = std::numeric_limits<double64>::quiet_NaN();
            D_Co2_     = std::numeric_limits<double64>::quiet_NaN();
            return mystate;
        }
        else if(bulk_xco2 <= xco2) // aqeuous phase only, all CO2 dissolved
        {
            mystate = SYSTEM_STATE::aq;
            massAqueousPhase_   = total_mass;
            massCarbonicPhase_  = 0.0;
            massHalite_         = 0.0;
            rho_aq_    = eos.densityAqueousPhase(eos.volumePartialMolarCo2(t),eos.densityBrine( p, t, 0. ),xco2);
            mu_aq_     = eos.mu_AqueousPhase( p, t, 0. );
            beta_aq_   = eos.C_AqueousPhase(p,t,0.);
            X_h2o_     = bulk_xh2o*molar_mass_h2o;
            X_co2_     = bulk_xco2*molar_mass_co2;
            massnorm   = X_h2o_ + X_co2_;
            X_h2o_    /= massnorm;
            X_co2_     = 1.-X_h2o_;
            X_nacl_    = 0.;
            D_Co2_     = eos.D_Co2( p, t, 0.0 );
            Y_h2o_     = std::numeric_limits<double64>::quiet_NaN();
            Y_co2_     = std::numeric_limits<double64>::quiet_NaN();
            rho_carb_  = std::numeric_limits<double64>::quiet_NaN();
            mu_carb_   = std::numeric_limits<double64>::quiet_NaN();
            beta_carb_ = std::numeric_limits<double64>::quiet_NaN();
            return mystate;
        }
        else
        {
            mystate = SYSTEM_STATE::aq_carb;
            moleFractionOfCarbonicPhase  = (bulk_xco2-xco2)/(1.-yh2o-xco2);
            moleFractionOfAqueousPhase   = 1.0-moleFractionOfCarbonicPhase;
            molar_mass_carb              = yh2o*molar_mass_h2o+(1.0-yh2o)*molar_mass_co2;
            molar_mass_aq                = xco2*molar_mass_co2 + (1.0-xco2)*molar_mass_h2o;
            massnorm                     = moleFractionOfCarbonicPhase*molar_mass_carb + moleFractionOfAqueousPhase*molar_mass_aq;
            massAqueousPhase_            = moleFractionOfAqueousPhase*molar_mass_aq/massnorm*total_mass;
            massCarbonicPhase_           = total_mass-massAqueousPhase_;
            massHalite_                  = 0.0;
            Y_h2o_                       = yh2o*molar_mass_h2o/molar_mass_carb;
            Y_co2_                       = 1.0-Y_h2o_;
            X_co2_                       = xco2*molar_mass_co2/molar_mass_aq;
            X_h2o_                       = 1.0-X_co2_;
            X_nacl_                      = 0.0;
            rho_aq_                      = eos.densityAqueousPhase(eos.volumePartialMolarCo2(t),eos.densityBrine( p, t, 0. ),xco2);
            mu_aq_                       = eos.mu_AqueousPhase( p, t, 0. );
            beta_aq_                     = eos.C_AqueousPhase(p,t,0.);
            rho_carb_                    = eos.Rho_CarbonicPhase(p,t);
            mu_carb_                     = eos.mu_CarbonicPhase(p,t);
            beta_carb_                   = eos.C_CarbonicPhase(p,t);
            D_Co2_                       = eos.D_Co2( p, t, 0.0 );
            return mystate;
        }
    }

    // Now the four triangular regions; I avoided looping just for my own clarity
    // This code can clearly be simplified
    // NOTICE THAT THE FOURTH (aq_carb) IS DIFFERENT AS THERE ARE ONLY TWO (NOT THREE) PHASES AND
    // THEIR COMPOSITIONS ARE A FUNCTION OF THE BULK COMPOSITION
//    cout << "Testing aq_salt ...\n";
    mystate = SYSTEM_STATE::aq_salt;
    SetCornerPointsTriangularRegions(mystate);
    equilibrated = EvaluateFractionsInTriangularRegion();
    if(equilibrated)
    {
#ifdef PHASESTATEFINDER_H2O_CO2_NACL_DEBUG
        cerr << "PhaseStateFinder_H2O_CO2_NaCl::Equilibrate: Point is in aq_salt\n";
#endif
        // First, compute composition of aqueous phase, which is the weighted sum/average of corners a+b
        molenorm                    = weight_phase_a + weight_phase_b;
        double64 xnacl                = (weight_phase_a*xsat + weight_phase_b*corner_xnacl[1])/molenorm;
        double64 xh2o                 = (weight_phase_a*(1.0-xsat) + weight_phase_b*corner_xh2o[1])/molenorm;
        double64 xco2                 = 1.0-xnacl-xh2o;
        molar_mass_aq               = xco2*molar_mass_co2 + xh2o*molar_mass_h2o + xnacl*molar_mass_nacl;
        X_h2o_                      = xh2o*molar_mass_h2o/molar_mass_aq;
        X_co2_                      = xco2*molar_mass_co2/molar_mass_aq;
        X_nacl_                     = xnacl*molar_mass_nacl/molar_mass_aq;
        const double64 mnacl        = xnacl*55.508/xh2o;
        rho_aq_                     = eos.densityAqueousPhase(eos.volumePartialMolarCo2(t),eos.densityBrine( p, t, mnacl ),xco2);
        mu_aq_                      = eos.mu_AqueousPhase( p, t, mnacl );
        beta_aq_                    = eos.C_AqueousPhase(p,t,mnacl);
        D_Co2_                      = eos.D_Co2( p, t, mnacl );
        massnorm                    = molenorm*molar_mass_aq + weight_phase_c*molar_mass_nacl;
        massAqueousPhase_           = molenorm*molar_mass_aq/massnorm*total_mass;
        massCarbonicPhase_          = 0.0;
        massHalite_                 = total_mass-massAqueousPhase_;
        Y_h2o_     = std::numeric_limits<double64>::quiet_NaN();
        Y_co2_     = std::numeric_limits<double64>::quiet_NaN();
        rho_carb_  = std::numeric_limits<double64>::quiet_NaN();
        mu_carb_   = std::numeric_limits<double64>::quiet_NaN();
        beta_carb_ = std::numeric_limits<double64>::quiet_NaN();
        return mystate;
    }

    else mystate = SYSTEM_STATE::aq_carb_salt;
//    cout << "Testing aq_carb_salt ...\n";
    SetCornerPointsTriangularRegions(mystate);
    equilibrated = EvaluateFractionsInTriangularRegion();
    if(equilibrated)
    {
#ifdef PHASESTATEFINDER_H2O_CO2_NACL_DEBUG
        cerr << "PhaseStateFinder_H2O_CO2_NaCl::Equilibrate: Point is in aq_carb_salt\n";
#endif
        MoleToMassfractionAtTriangularCorners();
        X_h2o_             = corner_massfraction_h2o[a];
        X_co2_             = corner_massfraction_co2[a];
        X_nacl_            = corner_massfraction_nacl[a];
        Y_h2o_             = corner_massfraction_h2o[b];
        Y_co2_             = corner_massfraction_co2[b];
        molar_mass_aq      = corner_xh2o[a]*molar_mass_h2o + corner_xco2[a]*molar_mass_co2 + corner_xnacl[a]*molar_mass_nacl;
        molar_mass_carb    = corner_xh2o[b]*molar_mass_h2o + corner_xco2[b]*molar_mass_co2;
        massnorm           = weight_phase_a*molar_mass_aq + weight_phase_b*molar_mass_carb +weight_phase_c*molar_mass_nacl;
        massAqueousPhase_  = weight_phase_a*molar_mass_aq/massnorm*total_mass;
        massCarbonicPhase_ = weight_phase_b*molar_mass_carb/massnorm*total_mass;
        massHalite_        = total_mass - massAqueousPhase_ - massCarbonicPhase_;
        const double64 mnacl = corner_xnacl[a]*55.508/corner_xh2o[a]; // get molality in aqueous phase
        rho_aq_            = eos.densityAqueousPhase(eos.volumePartialMolarCo2(t),eos.densityBrine( p, t, mnacl ),corner_xco2[a]);
        mu_aq_             = eos.mu_AqueousPhase( p, t, mnacl );
        beta_aq_           = eos.C_AqueousPhase(p,t,mnacl);
        D_Co2_             = eos.D_Co2( p, t, mnacl );
        rho_carb_          = eos.Rho_CarbonicPhase(p,t);
        mu_carb_           = eos.mu_CarbonicPhase(p,t);
        beta_carb_         = eos.C_CarbonicPhase(p,t);
        return mystate;
    }

    else mystate = SYSTEM_STATE::carb_salt;
//    cout << "Testing carb_salt ...\n";
    SetCornerPointsTriangularRegions(mystate);
    equilibrated = EvaluateFractionsInTriangularRegion();
    if(equilibrated)
    {
#ifdef PHASESTATEFINDER_H2O_CO2_NACL_DEBUG
        cerr << "PhaseStateFinder_H2O_CO2_NaCl::Equilibrate: Point is in carb_salt\n";
#endif
        // First, compute composition of carbonic phase, which is the weighted sum/average of corners a+b
        molenorm                    = weight_phase_a + weight_phase_b;
        double64 xh2o                 = (weight_phase_a*corner_xh2o[a] + weight_phase_b*corner_xh2o[b])/molenorm;
        double64 xco2                 =  1.0-xh2o;
        molar_mass_carb             =  xco2*molar_mass_co2 + xh2o*molar_mass_h2o;
        Y_h2o_                      =  xh2o*molar_mass_h2o/molar_mass_carb;
        Y_co2_                      =  xco2*molar_mass_co2/molar_mass_carb;
        massnorm                    = molenorm*molar_mass_carb + weight_phase_c*molar_mass_nacl;
        massAqueousPhase_           = 0.0;
        massCarbonicPhase_          = molenorm*molar_mass_carb/massnorm*total_mass;
        massHalite_                 = total_mass-massCarbonicPhase_;
        rho_carb_                   = eos.Rho_CarbonicPhase(p,t);
        mu_carb_                    = eos.mu_CarbonicPhase(p,t);
        beta_carb_                  = eos.C_CarbonicPhase(p,t);
        X_h2o_     = std::numeric_limits<double64>::quiet_NaN();
        X_co2_     = std::numeric_limits<double64>::quiet_NaN();
        X_nacl_    = std::numeric_limits<double64>::quiet_NaN();
        rho_aq_    = std::numeric_limits<double64>::quiet_NaN();
        mu_aq_     = std::numeric_limits<double64>::quiet_NaN();
        beta_aq_   = std::numeric_limits<double64>::quiet_NaN();
        D_Co2_     = std::numeric_limits<double64>::quiet_NaN();
        return mystate;
    }

    else mystate = SYSTEM_STATE::aq_carb;
//    cout << "Testing aq_carb ...\n";
    SetCornerPointsTriangularRegions(mystate);
    equilibrated = EvaluateFractionsInTriangularRegion(); // This is ok but we are not yet done if equilibrated
    if(equilibrated)
    {
#ifdef PHASESTATEFINDER_H2O_CO2_NACL_DEBUG
        cerr << "PhaseStateFinder_H2O_CO2_NaCl::Equilibrate: Point is in aq_carb\n";
#endif
        EvaluateSaltyAqCarb(); // does all updates
        return mystate;
    }

    else // this is NOT foolproof!!! There should be a true test if we are in aqueous condition  ...
    {
        mystate = SYSTEM_STATE::aq;
        massAqueousPhase_  = total_mass;
        massCarbonicPhase_ = 0.0;
        massHalite_        = 0.0;
        double64 xh2o        = bulk_xh2o;
        double64 xco2        = bulk_xco2;
        double64 xnacl       = bulk_xnacl;
        const double64 mnacl = xnacl*55.508/xh2o;
        rho_aq_            = eos.densityAqueousPhase(eos.volumePartialMolarCo2(t),eos.densityBrine( p, t, mnacl ),xco2);
        mu_aq_             = eos.mu_AqueousPhase( p, t, mnacl );
        beta_aq_           = eos.C_AqueousPhase(p,t,mnacl);
        X_h2o_             = bulk_massfraction_h2o;
        X_co2_             = bulk_massfraction_co2;
        X_nacl_            = bulk_massfraction_nacl;
        D_Co2_     = eos.D_Co2( p, t, mnacl );
        Y_h2o_     = std::numeric_limits<double64>::quiet_NaN();
        Y_co2_     = std::numeric_limits<double64>::quiet_NaN();
        rho_carb_  = std::numeric_limits<double64>::quiet_NaN();
        mu_carb_   = std::numeric_limits<double64>::quiet_NaN();
        beta_carb_ = std::numeric_limits<double64>::quiet_NaN();
        return mystate;
    }
}

void PhaseStateFinder_H2O_CO2_NaCl::EvaluateSaltyAqCarb()
{
    xsat = liquidus.MoleFractionNaCl();
    msat = XNaCl2Molal(xsat);

    // define cartesian coordinates of aqueous saturation point (xaqsat,yaqsat)
    double64 mco2 = eos.m_Co2(p,t,msat);
    double64 xco2 = mco2/(55.508+2.0*msat+mco2); // CAUTION: the factor 2 is a fix for Spycher's convention. I hope it is correct
    double64 xh2o = 55.508/(55.508+msat+mco2); // here that factor is NOT needed ...!
    double64 xnacl = 1.-xh2o-xco2;
    double64 xaqsat,yaqsat;
    ConvertMolefractionToCartesian(xaqsat, yaqsat, xh2o, xco2, xnacl);
    // scaling factor to convert between xnacl and cartesian y coordinate: y = xnacl*scale
    double64 scale = nacl_y; double64 scale2 = scale*scale;
    // define more points
    double64 xco2_carb_0 = 1.0-eos.y_H2o(p,t,0.); double64 xco2_carb_0_2 = xco2_carb_0*xco2_carb_0;
    double64 xco2_aq_0 = eos.x_Co2(p,t,0.);
    const double64 B = -1.73966; // a fit, valid for all t and p
    const double64 A = (xaqsat-xco2_aq_0)/(yaqsat); double64 A2 = A*A; // specific for t and p
    // next one is a Maple evaluation, can certainly be simplified
    double64 y2 = y_bulkcomp; double64 x2 = x_bulkcomp;
    double64 y2_2 = y2*y2; double64 x2_2 =x2*x2;
    double64 B2y2 = B*B*y2*y2;
    double64 twoBscaley2 = 2.*B*scale*y2;
    double64 twoBscaley2xcarb = twoBscaley2*xco2_carb_0;
    double64 Ay2         = A*y2;
    double64 xnacl_aq1 = (-0.5)/(scale*B*(xco2_carb_0-1.));
    xnacl_aq1       *=    ( A*scale*y2 +B*(- y2*xco2_carb_0 + y2) + scale*(- x2 + xco2_carb_0)
                            - sqrt( A2*scale2*y2_2 + Ay2*(-twoBscaley2*xco2_carb_0 + twoBscaley2)
                                    + B2y2*(xco2_carb_0_2-2*xco2_carb_0+1.)
                                    + scale2*(- 2.*Ay2*x2 + 2.*Ay2*xco2_carb_0 + x2_2 - 2.*x2*xco2_carb_0 + xco2_carb_0_2 )
                                    + twoBscaley2xcarb*(x2 - 2.*xco2_aq_0 + xco2_carb_0)
                                    + twoBscaley2*(2.*xco2_aq_0-xco2_carb_0 - x2)
                                    )
                            );

    // Convert result to cartesian coordinates for aqueous phase
    // backup variables: x2, y2 could be used instead, this looks cleaner
    double64 backup_x2 = x_bulkcomp;
    double64 backup_y2 = y_bulkcomp;
    x_bulkcomp = xco2_aq_0+A*(xnacl_aq1*scale);
    y_bulkcomp = xnacl_aq1*scale;
    double64 x1(x_bulkcomp);
    // Now convert to mole fractions of aqueous phase
    SetCornerPointsTriangularRegions( SYSTEM_STATE::full );
    EvaluateFractionsInTriangularRegion();
    corner_xh2o[0]  = weight_phase_a;
    corner_xco2[0]  = weight_phase_b;
    corner_xnacl[0] = weight_phase_c;
    x_bulkcomp = backup_x2;
    y_bulkcomp = backup_y2;

    // Now mole fractions of carbonic phase;
    corner_xco2[1] = xco2_carb_0*(1.+B*xnacl_aq1)-B*xnacl_aq1;
    corner_xh2o[1] = 1.0-corner_xco2[1];

    // Now get weights of phases
    double64 weightaq   = (corner_xco2[1]-x_bulkcomp)/(corner_xco2[1]-x1);
    double64 weightcarb = (1.0-weightaq);

    double64 molar_mass_aq, molar_mass_carb;
    double64 massnorm;

    molar_mass_aq      = corner_xh2o[0]*molar_mass_h2o + corner_xco2[0]*molar_mass_co2 + corner_xnacl[0]*molar_mass_nacl;
    molar_mass_carb    = corner_xh2o[1]*molar_mass_h2o + corner_xco2[1]*molar_mass_co2;
    massnorm           = weightaq*molar_mass_aq + weightcarb*molar_mass_carb;
    massCarbonicPhase_ = weightcarb*molar_mass_carb/massnorm*total_mass;
    massAqueousPhase_  = weightaq*molar_mass_aq/massnorm*total_mass;
    massHalite_        = 0.0;
    const double64 mnacl = corner_xnacl[0]*55.508/corner_xh2o[0];
    X_h2o_             = corner_xh2o[0]*molar_mass_h2o/molar_mass_aq;
    Y_h2o_             = corner_xh2o[1]*molar_mass_h2o/molar_mass_carb;
    X_co2_             = corner_xco2[0]*molar_mass_co2/molar_mass_aq;
    Y_co2_             = corner_xco2[1]*molar_mass_co2/molar_mass_carb;
    X_nacl_            = corner_xnacl[0]*molar_mass_nacl/molar_mass_aq;
    rho_aq_            = eos.densityAqueousPhase(eos.volumePartialMolarCo2(t),eos.densityBrine( p, t, mnacl ),corner_xco2[0]);
    mu_aq_             = eos.mu_AqueousPhase( p, t, mnacl);
    beta_aq_           = eos.C_AqueousPhase(p,t,mnacl);
    rho_carb_          = eos.Rho_CarbonicPhase(p,t);
    mu_carb_           = eos.mu_CarbonicPhase(p,t);
    beta_carb_         = eos.C_CarbonicPhase(p,t);
    D_Co2_     = eos.D_Co2( p, t, mnacl );
    return;
}

void PhaseStateFinder_H2O_CO2_NaCl::MoleToMassfractionAtTriangularCorners()
{
    for(int i = 0; i < 3; ++i)
    {
        corner_massfraction_h2o[i]   = corner_xh2o[i]*molar_mass_h2o;
        corner_massfraction_co2[i]   = corner_xco2[i]*molar_mass_co2;
        corner_massfraction_nacl[i]  = corner_xnacl[i]*molar_mass_nacl;
        const double64 massnorm        = corner_massfraction_h2o[i] + corner_massfraction_co2[i] + corner_massfraction_nacl[i];
        corner_massfraction_h2o[i]  /= massnorm;
        corner_massfraction_co2[i]  /= massnorm;
        corner_massfraction_nacl[i] /= massnorm;
    }
}


void PhaseStateFinder_H2O_CO2_NaCl::ConvertMolefractionTriangularCornersToCartesian(corner mycorner)
{
    ConvertMolefractionToCartesian(x[mycorner], y[mycorner], corner_xh2o[mycorner], corner_xco2[mycorner], corner_xnacl[mycorner]);
}


void PhaseStateFinder_H2O_CO2_NaCl::ConvertMolefractionToCartesian(double64& x, double64& y, double64 xh2o, double64 xco2, double64 xnacl)
{
    x = xh2o*h2o_x + xco2*co2_x + xnacl*nacl_x;
    y = xh2o*h2o_y + xco2*co2_y + xnacl*nacl_y;
}

void PhaseStateFinder_H2O_CO2_NaCl::InputComposition()
{
    cout << "Please enter mole fractions at point of interest (xnacl will be computed from xh2o and xco2):\n";
    cout << "xh2o = "; cin >> bulk_xh2o;
    cout << "xco2 = "; cin >> bulk_xco2;
    bulk_xnacl = 1.0-bulk_xh2o-bulk_xco2;
    ConvertMolefractionToCartesian(x_bulkcomp, y_bulkcomp, bulk_xh2o, bulk_xco2, bulk_xnacl);
}



void PhaseStateFinder_H2O_CO2_NaCl::InputXcompYcomp()
{
    cout << "x = "; cin >> x_bulkcomp;
    cout << "y = "; cin >> y_bulkcomp;
}



bool PhaseStateFinder_H2O_CO2_NaCl::EvaluateFractionsInTriangularRegion()
{
    weight_phase_b = -(x[0]*(y[2]-y_bulkcomp) + x[2]*(-y[0]+y_bulkcomp) + x_bulkcomp*(y[0]-y[2]));
    weight_phase_b /= (x[0]*(y[1]-y[2]) + x[1]*(-y[0]+y[2])+x[2]*(y[0]-y[1]));

    weight_phase_a = -(weight_phase_b*y[1]+y[2]*(-weight_phase_b+1.)-y_bulkcomp)/(y[0]-y[2]);
    weight_phase_c = 1.0 - weight_phase_b - weight_phase_a;
    /* epsilon seems too strict?
    if(weight_phase_b >= -std::numeric_limits<double64>::epsilon() &&
       weight_phase_a >= -std::numeric_limits<double64>::epsilon() &&
       weight_phase_c >= -std::numeric_limits<double64>::epsilon() ) return true;
    */   
    if(weight_phase_b >= -1.0e-6 &&
       weight_phase_a >= -1.0e-6 &&
       weight_phase_c >= -1.0e-6 ) return true;       
    
    // if one of the weights is negative, the point is outside the triangle
    // it should be possible to identify then the next region to search!
    return false;
}



void PhaseStateFinder_H2O_CO2_NaCl::SetCornerPointsTriangularRegions( SYSTEM_STATE mystate ) // tested 1.1.19, TD
{
    xsat = liquidus.MoleFractionNaCl();
    msat = XNaCl2Molal(xsat);
    double64 mco2,yh2o;

    switch(mystate)
    {
    case SYSTEM_STATE::aq_salt :
        corner_xh2o[0]=1.0-xsat; corner_xco2[0]=0.0; corner_xnacl[0]=xsat;
        ConvertMolefractionTriangularCornersToCartesian(a);
        mco2 = eos.m_Co2(p,t,msat);
        corner_xco2[1] = mco2/(55.508+2.0*msat+mco2); // CAUTION: the factor 2 is a fix for Spycher's convention. I hope it is correct
        corner_xh2o[1] = 55.508/(55.508+msat+mco2);
        corner_xnacl[1] = 1.-corner_xh2o[1]-corner_xco2[1];
        ConvertMolefractionTriangularCornersToCartesian(b);
        corner_xh2o[2] = 0.; corner_xco2[2]=0.; corner_xnacl[2]=1.;
        x[2] = nacl_x; y[2] = nacl_y;
        break;
    case SYSTEM_STATE::carb_salt :
        yh2o = corner_xh2o[0]=eos.y_H2o(p,t,msat); corner_xco2[0]=1.0-yh2o; corner_xnacl[0]=0.;
        ConvertMolefractionTriangularCornersToCartesian(a);
        corner_xh2o[1]=0.; corner_xco2[1]=1.; corner_xnacl[1]=0.;
        corner_xh2o[2]=0.; corner_xco2[2]=0.; corner_xnacl[2]=1.;
        x[1] = co2_x; y[1] = co2_y;
        x[2] = nacl_x; y[2] = nacl_y;
        break;
    case SYSTEM_STATE::aq_carb_salt :
        mco2 = eos.m_Co2(p,t,msat);
        corner_xco2[0] = mco2/(55.508+2.0*msat+mco2); // CAUTION: the factor 2 is a fix for Spycher's convention. I hope it is correct
        corner_xh2o[0] = 55.508/(55.508+msat+mco2);
        corner_xnacl[0] = 1.-corner_xco2[0]-corner_xh2o[0];
        ConvertMolefractionTriangularCornersToCartesian(a);
        corner_xh2o[1] = eos.y_H2o(p,t,msat);
        corner_xco2[1] = 1.0-corner_xh2o[1];
        corner_xnacl[1]= 0.0;
        ConvertMolefractionTriangularCornersToCartesian(b);
        corner_xh2o[2]=0.;
        corner_xco2[2]=0.;
        corner_xnacl[2]=1.;
        x[2] = nacl_x; y[2] = nacl_y;
        break;
    case SYSTEM_STATE::aq_carb :
        corner_xco2[0]  = eos.x_Co2(p,t,0.0);
        corner_xh2o[0]  = 1.0 - corner_xco2[0];
        corner_xnacl[0] = 0.0;
        ConvertMolefractionTriangularCornersToCartesian(a);
        corner_xh2o[1] = eos.y_H2o(p,t,msat);
        corner_xco2[1] = 1.0-corner_xh2o[1];
        corner_xnacl[1]= 0.0;
        ConvertMolefractionTriangularCornersToCartesian(b);
        mco2 = eos.m_Co2(p,t,msat);
        corner_xco2[2]  = mco2/(55.508+2.0*msat+mco2); // CAUTION: the factor 2 is a fix for Spycher's convention. I hope it is correct
        corner_xh2o[2]  = 55.508/(55.508+msat+mco2);
        corner_xnacl[2] = 1.-corner_xco2[2]-corner_xh2o[2];
        ConvertMolefractionTriangularCornersToCartesian(c);
        break;
    case SYSTEM_STATE::full:
        corner_xh2o[0]=1.; corner_xco2[0]=0.; corner_xnacl[0]=0.;
        corner_xh2o[1]=0.; corner_xco2[1]=1.; corner_xnacl[1]=0.;
        corner_xh2o[2]=0.; corner_xco2[2]=0.; corner_xnacl[2]=1.;
        x[0] = h2o_x; y[0] = h2o_y;
        x[1] = co2_x; y[1] = co2_y;
        x[2] = nacl_x; y[2] = nacl_y;
        break;
    default:
        cerr << "ERROR: PhaseStateFinder_H2O_CO2_NaCl::SetCornerPointsTriangularRegions(state mystate):\n";
        cerr << "parameter mystate must be aq_salt (=3), carb_salt (=5), or aq_carb_salt (=6)\n";
        cerr << "state provided was " << parseState(mystate) << endl;
    }
    return;
}



void PhaseStateFinder_H2O_CO2_NaCl::DisplayFractions( SYSTEM_STATE mystate )
{
    SetCornerPointsTriangularRegions(mystate);
    double64 xh2o,xco2,xnacl;
    cout << "xh2o : "; cin >> xh2o;
    cout << "xco2 : "; cin >> xco2;
    xnacl = 1.0-xh2o-xco2;
    cout << "xnacl: " << xnacl << endl;
    ConvertMolefractionToCartesian(x_bulkcomp, y_bulkcomp, xh2o, xco2, xnacl);
    EvaluateFractionsInTriangularRegion();
}



/**
    Reports the system state as a string.
    State refers to how many phases co-exist and what they are composed of.
    See phase diagram from Thomas.
 
    @attention the aqueous phase can contain various amounts of dissolved salt.
*/
std::string parseState( SYSTEM_STATE phase_state )
 {
    if ( phase_state == SYSTEM_STATE::aq ) return "aqueous phase only";
    if ( phase_state == SYSTEM_STATE::carb ) return "carbonic phase only";
    if ( phase_state == SYSTEM_STATE::salt ) return "solid salt only";
    if ( phase_state == SYSTEM_STATE::aq_salt ) return "aqueous + salt phase";
    if ( phase_state == SYSTEM_STATE::aq_carb ) return "aqueous + carbonic phase";
    if ( phase_state == SYSTEM_STATE::carb_salt ) return "carbonic + salt phase";
    if ( phase_state == SYSTEM_STATE::aq_carb_salt ) return "aqueous + carbonic + salt phase";
    if ( phase_state == SYSTEM_STATE::full ) return "full system";
    if ( phase_state == SYSTEM_STATE::undefined ) return "undefined";

    return "parseState (H2O-CO2-NaCl system): failed to apply any conversion of system state to string.";
 }






/** Out()

    const double64& t;
    const double64& p;
    const double64& total_mass;
    const double64& bulk_massfraction_h2o;
    const double64& bulk_massfraction_co2;
    const double64& bulk_massfraction_nacl;

    const double64 PI;
    const double64 h2o_x;
    const double64 h2o_y;
    const double64 co2_x;
    const double64 co2_y;
    const double64 nacl_x;
    const double64 nacl_y;
    const double64 molar_mass_h2o;
    const double64 molar_mass_co2;
    const double64 molar_mass_nacl;
    const double64 moles_in_kg_h2o;

    double64 weight_phase_a,weight_phase_b,weight_phase_c; // weights (mole-based) of phases on corner points in triangular regions
    double64 bulk_xh2o,bulk_xco2,bulk_xnacl;               // mole fractions of bulk composition
    double64 bulk_mnacl;                                   // molality of NaCl for the bulk composition
    double64 x_bulkcomp,y_bulkcomp;                        // cartesian coordinates of bulk composition
    double64 x[3],y[3];                                    // corner points of triangular regions within phase diagram
    double64 m;
    double64 xsat, msat, mnacl;

    // SKM FIX - there were 2 versions of total mass. This one is disambiguated by underscore for class member
    //double64 total_mass;

    double64 massCarbonicPhase_;
    double64 rho_carb_;
    double64 mu_carb_;
    double64 beta_carb_;
    double64 Y_h2o_;
    double64 Y_co2_;

    double64 massAqueousPhase_;
    double64 rho_aq_;
    double64 mu_aq_;
    double64 beta_aq_;
    double64 X_h2o_;
    double64 X_co2_;
    double64 X_nacl_;
    double64 D_Co2_;

    double64 massHalite_;

*/
void PhaseStateFinder_H2O_CO2_NaCl::Out( bool print_private_state ) const
 {
    cout <<"\nPhaseStateFinder_H2O_CO2_NaCl::Out: current values:";
     // Access to results, always run Equilibrate before
    cout <<"\n\tmass of carbonic phase (kg): "<<  massCarbonicPhase();
    cout <<"\n\tdensity (kg/m3):             "<<  rho_carb();
    cout <<"\n\tviscosity (Pa.s):            "<<  mu_carb();
    cout <<"\n\tcompressibility (1/Pa):      "<<  beta_carb();
    cout <<"\n\tmass fraction CO2:           "<<  Y_co2();
    cout <<"\n\tmass fraction H2O:           "<<  Y_h2o();

    cout <<"\n\tmass of aqueous phase (kg):  "<<  massAqueousPhase();
    cout <<"\n\tdensity (kg/m3):             "<<  rho_aq();
    cout <<"\n\tviscosity (Pa.s):            "<<  mu_aq();
    cout <<"\n\tcompressibility (1/Pa):      "<<  beta_aq();
    cout <<"\n\tmass fraction H2O:           "<<  X_h2o();
    cout <<"\n\tmass fraction CO2:           "<<  X_co2();
    cout <<"\n\tmass fraction NaCl:          "<<  X_nacl();
    cout <<"\n\tdiffusivity of CO2 (m2/s):   "<<  D_Co2();

    cout <<"\n\tmass of halite (salt) (kg):  "<<  massHalite();
 
    if ( print_private_state )
     cerr <<"\nPhaseStateFinder_H2O_CO2_NaCl::Out: private state cannot be printed yet.\n";
    
} // end Out


} // end csmp
