// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "MagmaModel.h"
#include "Model.h"
#include "Region.h"
#include "PropertyDatabase.h"

using namespace std;
namespace csmp {

template<size_t dim>
MagmaModel<dim>::MagmaModel(Model<dim>& model,
                            double initial_wtpercent_volatile_in_magma,
                            double crystal_density,
                            double crystallization_curve_exponent,
                            double initial_volatile_vol_frac_chamber,
                            double compressibility_rock,
                            double compressibility_magma,
                            double sigma1,
                            int saturation_model,
                            double saturation,
                            double salinity_magmatic_fluids)

    : Visitor<dim>(MODEL, NODE),

    with_pressure_volume_change (true), // allows for pore volume change during degassing

    model_ref_(model),
    prop_ref_(model.Database()),

    liquidus_temperature_key                      (model.Database().StorageKey("liquidus temperature")),
    solidus_temperature_key                       (model.Database().StorageKey("solidus temperature")),
    bulk_volume_key                               (model.Database().StorageKey("bulk volume")),
    minimum_pore_volume_key                       (model.Database().StorageKey("mini pore volume")),
    minimum_mass_fluid_key                        (model.Database().StorageKey("mini mass fluid")),
    volume_deformation_key                        (model.Database().StorageKey("volume deformation")),
    total_volume_deformation_key                  (model.Database().StorageKey("total volume deformation")),
    melt_volume_key                               (model.Database().StorageKey("melt volume")),
    crystal_volume_key                            (model.Database().StorageKey("crystal volume")),
    delta_melt_volume_key                         (model.Database().StorageKey("melt volume change")),
    delta_crystal_volume_key                      (model.Database().StorageKey("crystal volume change")),
    crystallinity_key                             (model.Database().StorageKey("crystallinity")),
    mass_volatile_produced_key                    (model.Database().StorageKey("mass volatile produced")),
    total_mass_volatile_produced_key              (model.Database().StorageKey("total mass volatile produced")),
    melt_mass_key                                 (model.Database().StorageKey("melt mass")),
    crystal_mass_key                              (model.Database().StorageKey("crystal mass")),
    pressure_key                                  (model.Database().StorageKey("fluid pressure")),
    magma_density_key                             (model.Database().StorageKey("magma density")),
    melt_density_key                              (model.Database().StorageKey("melt density")),
    temperature_key                               (model.Database().StorageKey("temperature")),
    previous_temperature_key                      (model.Database().StorageKey("previous temperature")),
    temperature_difference_key                    (model.Database().StorageKey("temperature difference")),
    nodal_compressibility_rock_key                (model.Database().StorageKey("nodal compressibility rock")),
    nodal_porosity_key                            (model.Database().StorageKey("nodal porosity")),
    volatile_dissolved_fraction_key               (model.Database().StorageKey("volatile dissolved fraction")),
    volatile_dissolved_mass_key                   (model.Database().StorageKey("volatile dissolved mass")),
    volatile_dissolved_eq_key                     (model.Database().StorageKey("volatile dissolved eq")),
    pore_volume_key                               (model.Database().StorageKey("pore volume")),
    volatile_saturation_key                       (model.Database().StorageKey("volatile saturation")),
    relative_saturation_key                       (model.Database().StorageKey("relative saturation")),
    volume_change_factor_key                      (model.Database().StorageKey("pore volume change factor")),
    crystal_fraction_key                          (model.Database().StorageKey("crystal fraction")),
    melt_fraction_key                             (model.Database().StorageKey("melt fraction")),
    volatile_fraction_key                         (model.Database().StorageKey("volatile fraction")),
    rock_density_scaling_key                      (model.Database().StorageKey("rock density scaling")),
    nodal_density_rock_key                        (model.Database().StorageKey("nodal density rock")),
    threshold_pressure_key                        (model.Database().StorageKey("threshold pressure")),
    lithostatic_pressure_key                      (model.Database().StorageKey("lithostatic pressure")),
    minimum_nodal_density_rock_key                (model.Database().StorageKey("mini nodal density rock")),
    porous_flag_key                               (model.Database().StorageKey("porous flag node")),
    volumetric_mass_key                           (model.Database().StorageKey("fluid density")),
    Htp_key                                       (model.Database().StorageKey("previous total enthalpy")),
    msp_key                                       (model.Database().StorageKey("previous mass salt")),
    mf_key                                        (model.Database().StorageKey("fluid mass")),
    ml_key                                        (model.Database().StorageKey("fluid mass liquid")),
    mv_key                                        (model.Database().StorageKey("fluid mass vapor")),
    mm_key                                        (model.Database().StorageKey("magmatic fluid mass")),
    mms_key                                       (model.Database().StorageKey("magmatic mass salt")),
    kcl_ms_key                                    (model.Database().StorageKey("mass KCl")),
    ddt_magmatic_ratio_key                        (model.Database().StorageKey("ddt magmatic ratio")),
    mpr_key                                       (model.Database().StorageKey("magmatic production rate")),

    // for consistecy checks:
    m_key                                         (model.Database().StorageKey("mass")),
    H_key                                         (model.Database().StorageKey("enthalpy")),

    hCl_key                                       (model.Database().StorageKey( "enthalpy content liquid" ) ),
    hCv_key                                       (model.Database().StorageKey( "enthalpy content vapor" ) ),

    wt_key                                        (model.Database().StorageKey("salinity")),
    depth_key                                     (model.Database().StorageKey("nodal depth")),
    average_crystallinity_key                     (model.Database().StorageKey("average crystallinity")),
    delta_average_crystallinity_key               (model.Database().StorageKey("delta average crystallinity")),


    compressibility_rock_                         (compressibility_rock),
    compressibility_magma_                        (compressibility_magma),
    crystallization_curve_exponent_               (crystallization_curve_exponent),
    sigma1_                                       (sigma1),   //Parameter for erf function crystallinity
    crystal_density_                              (crystal_density),
    initial_volatile_vol_frac_chamber_            (initial_volatile_vol_frac_chamber),
    initial_wtpercent_volatile_in_magma_          (initial_wtpercent_volatile_in_magma),//in mass %
    saturation_model                              (saturation_model),
    saturation                                    (saturation),
    salinity_magmatic_fluids                      (salinity_magmatic_fluids),

    // initialization; parameters are read for each node
    tem_b(0.),              // temperature of MVP
    pre_b(0.),              // lithostatic pressure
    sal_b(0.),              // NaCl mass fraction
    x_NaCl_b(0.),           // NaCl mole fraction

    //brine(tem_b, pre_b, sal_b),

    // Fluid object from which to fetch thermodynamic properties, can replace "brine" as it is much faster
    fluid ( tem_b, pre_b, x_NaCl_b, 0., 0., 0., 1., false )


{} // end MagmaModel

template<size_t dim>
MagmaModel<dim>::~MagmaModel()
{} // end ~MagmaModel

template<size_t dim>
void MagmaModel<dim>::Visit(Model<dim>* node)
{}

template<size_t dim>
void MagmaModel<dim>::Visit(Region<dim>* node)
{}

template<size_t dim>
void MagmaModel<dim>::Visit(Node<dim>* node)
{

    ReadVariables( node );

    double volumetric_mass_before= volumetric_mass();
    double mass_produced_theory=0.0;

    delta_melt_volume()=0.0;
    delta_crystal_volume()=0.0;
    volume_deformation()=0.0;
    mass_volatile_produced()=0.0;

    nodal_compressibility_rock()=compressibility_magma_; // Base value if not rock
    temperature_difference()=temperature()-previous_temperature();

    sal_b = wt()/100.;                   // salinity of MVP in mass fraction NaCl
    tem_b = temperature();               // temperature of MVP
    pre_b = threshold_pressure();        // threshold pressure
    x_NaCl_b  = Weight2XNaCl(wt());

    double rhothresh=fluid.BulkProperties().rho;
    pre_b = lithostatic_pressure();      // lithostatic pressure
    double thermodynamic_density=fluid.BulkProperties().rho;

    double magma_deformation=0.0;
    double deformation_overpressure=0.0;
    double density_deformation=0.0;
    double minimum_pore_volume_before=minimum_pore_volume();

    double solid_volume_before=crystal_volume_before()+melt_volume_before();
    double bulk_volume_before=bulk_volume();

    double fluidization_limit=0.4;     //Crystal fraction above which the matrix do not behave as a fluid to maintain Plith (low bound)

    bool pause_at_end=false;
    bool display=false;


    bool small_pore_volume=false;

    cerr<<defaultfloat;
    cerr.precision(6);


    ///////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    /// Calculate crystallinity

    // Reversible instantaneous crystallinity calculation
    if (temperature()<=solidus_temperature())
    {
        crystallinity()=1.0;
    }
    else if (temperature()>=liquidus_temperature())
    {
        crystallinity()=0.0;
    }
    else
    {
        crystallinity()=1-pow(((temperature()-solidus_temperature())/(liquidus_temperature()-solidus_temperature())),crystallization_curve_exponent_);
    }

    // Calculate new melt and crystal volumes based on new crystallinity, only if there is melt left
    if (melt_volume_before()>0.0)
    {
        // First we calculate the mass source sink depending on melt under/over saturation
        // after VolatileCalc for a granite
        volatile_dissolved_eq()=6.833e-7*pow(temperature(),2)-0.001356*temperature()+0.7273;
        volatile_dissolved_eq()*=pow(pressure()*1.0e-5,-6.831e-7*pow(temperature(),2)+0.001416*temperature()-0.1255);
        volatile_dissolved_eq()/=100;

        volatile_dissolved_fraction()=volatile_dissolved_eq();
        volatile_dissolved_mass()=melt_mass_before()*volatile_dissolved_eq();

        mass_volatile_produced()=volatile_dissolved_mass_before() - volatile_dissolved_mass();

        mass_produced_theory=mass_volatile_produced();

        volumetric_mass()+= mass_volatile_produced()/pore_volume_before();         //Source-Sink

        // If we have sinks, we need to avoid "drying"
        if(definitelyLessThan(volumetric_mass()*pore_volume_before(),minimum_mass_fluid(),1.e-06))
        {
            mass_volatile_produced()+=(minimum_mass_fluid()-volumetric_mass()*pore_volume_before());
            volumetric_mass() = volumetric_mass_before+ mass_volatile_produced()/pore_volume_before();
        }

        // new amount of dissolved volatile mass after "production"
        volatile_dissolved_mass()=volatile_dissolved_mass_before()-mass_volatile_produced();

        delta_melt_volume() -= crystallinity() * ( melt_volume_before() + crystal_volume_before() - mass_volatile_produced() / crystal_density_ ) ;
        delta_melt_volume() += crystal_volume_before() - mass_volatile_produced() / crystal_density_;
        delta_melt_volume() /= ( crystallinity() * ( 1. - melt_density() / crystal_density_ ) + melt_density() / crystal_density_ );

        delta_crystal_volume()-=delta_melt_volume() * melt_density() / crystal_density_;
        delta_crystal_volume()-=mass_volatile_produced() / crystal_density_;

        if(crystallinity_before()<=0.0 && crystallinity()<=0.0)
        {
            delta_crystal_volume()=0.0;
            delta_melt_volume()=0.0;
        }
    }

    injected_fluid()= mass_volatile_produced()/pore_volume_before(); //kg m-3

    ///////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    /// Calculate melt and crystal volumes

    // update melt and crystal volumes by volume differences just calculated
    melt_volume()=melt_volume_before()+delta_melt_volume();
    crystal_volume()=crystal_volume_before()+delta_crystal_volume();

    // apply threshold: volumes below 1e-6 m³ are treated as zero for numerical stability
    if (melt_volume()<=1.e-6)
    {
        melt_volume()=0.0;
        delta_melt_volume()=-melt_volume_before();
        melt_density()=0.0;
        volatile_dissolved_mass()=0.0;
        volatile_dissolved_fraction()=0.0;
        volatile_dissolved_eq()=0.0;
    }

    if (crystal_volume()<=1.e-6)
    {
        crystal_volume()=0.0;
        delta_crystal_volume()=-crystal_volume_before();
    }

    // Critical line: change in pore volume due to crystallization or melting
    // (density difference between melt and crystals)
    magma_deformation = delta_melt_volume()+delta_crystal_volume();

    // Calculate masses
    melt_mass()     = melt_volume() * melt_density();
    crystal_mass()  = crystal_volume() * crystal_density_;

    // Calculate change in melt volume with density change
    if (melt_volume()>0.0)
    {
        volatile_dissolved_fraction()=volatile_dissolved_mass()/melt_mass();

        // parametrization for silicic melt
        // it is valid for silicic melt, and uses volatile_dissolved_fraction() as input
        // the previous (wrong) one was for dacitic melt and for percentages. it gave wrong values the way it was used in CSMP

        // for Gilgai Granite (leucogranite) as a proxy for the Mole Granite. Best fit compiled from data by Juniper & Kleeman (1979):
        //  c1 = 2.416e3    c2 = 1.982e3    c3 = 1.786e-7
        //  c4 = 1.162e-7   c5 = 9.316e-2   c6 = 1.862
        //  density = c1 - (c2*xh2o) + (c3 + (c4*xh2o))*pressure - (c5 + (c6*xh2o))*temperature

        double c1(2.416e3), c2(1.982e3), c3(1.786e-7), c4(1.162e-7), c5(9.316e-2), c6(1.862);
        double T(temperature()),P(0.),X(volatile_dissolved_fraction());
        if(pressure()<lithostatic_pressure())
        {
            P = lithostatic_pressure();
        }
        else
        {
            P = pressure();
        }

        melt_density() = c1 - (c2*X) + (c3 + (c4*X))*P - (c5+  (c6*X))*T;

        melt_volume()=melt_mass()/melt_density();
        delta_melt_volume()=melt_volume()-melt_volume_before();

        // Critical line: change in pore volume due to crystallization or melting
        // (density difference between melt and crystals)
        magma_deformation=delta_melt_volume()+delta_crystal_volume();

        if((minimum_pore_volume()-magma_deformation)/bulk_volume()<initial_volatile_vol_frac_chamber_)
        {
            magma_deformation=-bulk_volume()*initial_volatile_vol_frac_chamber_+minimum_pore_volume();
            delta_melt_volume()=magma_deformation-delta_crystal_volume();
            melt_volume()=melt_volume_before()+delta_melt_volume();
            melt_density()=melt_mass()/melt_volume();
            small_pore_volume=true;
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Temperature-Pressure-Volume DEFORMATIONS

    // Crystallization related volume change
    pore_volume()=pore_volume_before();

    if( magma_deformation!=0.0  )
    {
        //Adapted pore volume due to crystallization/melting
        pore_volume()-=magma_deformation;
        minimum_pore_volume()-=magma_deformation;
    }

    // Calculate factor by which pore volume has changed;
    // Needed to ensure mass balance of masses of advection variables
    volume_change_factor() = pore_volume()/pore_volume_before();
    volumetric_mass()/=volume_change_factor();

    double MC_PV=pore_volume();
    double MC_rho=volumetric_mass();

    ///////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    /// Calculate pore volume deformation in case of underpressure
    /// Change in volume to maintain Plith low bound in the suspension (below fluidization limit)

    if ( MC_rho<thermodynamic_density && with_pressure_volume_change && crystal_fraction()<fluidization_limit  && !small_pore_volume )
    {
        density_deformation = pore_volume()*(MC_rho/thermodynamic_density-1.0);
        //Collapse porosity by very much if element switch after channels formation because volumetric_mass decreased...

        //if pore volume is smaller than minimum pore volume, clamp pore volume to minimum pore volume
        if(pore_volume()+density_deformation<minimum_pore_volume())
        {
            density_deformation=minimum_pore_volume()-pore_volume();
            small_pore_volume=true;
        }

        bulk_volume()+=density_deformation;
        pore_volume()+=density_deformation;
        volume_deformation()+=density_deformation;
        total_volume_deformation()+=density_deformation;

        // Calculate factor by which pore volume has changed;
        // Needed to ensure mass balance of masses of advection variables
        volume_change_factor() = pore_volume()/MC_PV;
        volumetric_mass()=MC_rho/volume_change_factor();

        if(essentiallyEqual(volumetric_mass(),thermodynamic_density,1.0e-06))
        {
            volumetric_mass()=thermodynamic_density;
        }
        if(pore_volume()<minimum_pore_volume())
        {
            minimum_pore_volume()=pore_volume();
        }

        pause_at_end = FractionSumChecks(pause_at_end);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    /// Calculate pore volume deformation in case of overpressure
    /// Change in volume to force maximum pressure Pthreshold by pore space dilataion
    /// Pore space can contract back to minimum pore volume if pressure is released

    if(MC_rho>rhothresh && with_pressure_volume_change && crystal_fraction() < 0.7 )
    {
        deformation_overpressure=pore_volume()*(MC_rho/rhothresh-1.0);

        // If pore volume is at minimum, negative pressure-volume deformation is cancelled
        if(deformation_overpressure<0.0 && (pore_volume()<=minimum_pore_volume() or small_pore_volume))
        {
            deformation_overpressure=0.0;
        }

        // At high crystallinities, magma will not deform, it can be fractured in the Parmigiani model instead
        if (crystal_fraction()>0.7)
        {
            nodal_compressibility_rock()=compressibility_rock_;
            deformation_overpressure=0.0;
        }

        if(deformation_overpressure !=0.0 )
        {
            //if deformation leads to pore volume smaller than minimum pore volume, clamp deformation to diff between pore volume and the minimum
            if (pore_volume()+deformation_overpressure<minimum_pore_volume())
            {
                deformation_overpressure=minimum_pore_volume()-pore_volume();
            }

            // Update all volumes and volume differences by calculated deformation due to overpressure
            bulk_volume()+=deformation_overpressure;
            pore_volume()+=deformation_overpressure;
            volume_deformation()+=deformation_overpressure;
            total_volume_deformation()+=deformation_overpressure;
        }

        // Calculate factor by which pore volume has changed;
        // Needed to ensure mass balance of masses of advection variables
        volume_change_factor() = pore_volume()/MC_PV;
        volumetric_mass() = MC_rho/volume_change_factor();

        pause_at_end = FractionSumChecks(pause_at_end);
    }

    pause_at_end = ConsistencyChecks1(pause_at_end, density_deformation, deformation_overpressure);

    // Rock density scaling
    double solid_volume=crystal_volume()+melt_volume();

    rock_density_scaling()=bulk_volume()/bulk_volume_before;//bulk volume change scaling (not solid volume)
    //                                                      //this is to conserve rock mass if written rock mass=(1-phi)*rho rock
    nodal_density_rock()*=rock_density_scaling();

    minimum_rock_density_scaling()=(solid_volume+minimum_pore_volume_before)/(solid_volume_before+minimum_pore_volume());
    minimum_nodal_density_rock()*=minimum_rock_density_scaling();

    // Calculate volume fractions, porosity and final pore volume factor
    crystal_fraction()      = crystal_volume()/bulk_volume();
    melt_fraction()         = melt_volume()/bulk_volume();
    volatile_fraction()     = pore_volume()/bulk_volume();
    nodal_porosity()        = volatile_fraction();
    volatile_saturation()   = pore_volume()/(bulk_volume()-crystal_volume());
    volume_change_factor() = pore_volume()/pore_volume_before();

    pause_at_end = FractionSumChecks(pause_at_end);
    pause_at_end = ConsistencyChecks2(pause_at_end);

    // Sources++
    pre_b = pressure(); // in Pa
    mass() += injected_fluid();
    magmatic_mass() += injected_fluid();
    fluid_mass() += injected_fluid();
    magma_density()=(crystal_mass()+melt_mass())/bulk_volume()+volatile_fraction()*volumetric_mass();

    brine_enthalpy=fluid.BulkProperties().h;

    Htp() += mass_volatile_produced()*brine_enthalpy/bulk_volume_before;   //ideally we should modify the enthalpy content "fluid" instead!!!!!
    enthalpy() += mass_volatile_produced()*brine_enthalpy/bulk_volume_before;

    msp() += injected_fluid()*salinity_magmatic_fluids*0.01;
    magmatic_salt() += injected_fluid()*salinity_magmatic_fluids*0.01;

    if (with_halite_trend) mass_KCl() += injected_fluid()*salinity_magmatic_fluids*0.01*fraction_KCl_magma;



    nmf() = mass_volatile_produced();
    mpr() = nmf()/dt;

    model_ref_.InputPropertyValue("magmatic production rate", mpr);

    // Calculate total mass produced at a node
    total_mass_volatile_produced()+=mass_volatile_produced();


    if (volatile_dissolved_eq()!=0.)
    {
        relative_saturation() = std::min(1.,volatile_dissolved_fraction()/volatile_dissolved_eq());
    }
    else
    {
        relative_saturation() = 0.;
    }


    StoreVariables(node);

    if(pause_at_end) {
        cerr<<"This is final pause"<<endl;
        string dummy;
        cerr << " press ENTER to continue... ";
        getline( cin, dummy );}

} // end Visit

template<size_t dim>
void MagmaModel<dim>::InitializeMagma(Model<dim>* node, std::string regionname)
{
    const Region<dim>& ref = node->Region(regionname);

    double minimum_volatile_dissolved(1.);

    // In case of saturation model 2, this loop over all nodes is used to determine the minimal value of thermodynamic equilibrium volatile saturation encountered at any node
    // in the chamber.
    if (saturation_model == 2)
    {
        for( typename std::vector<Node<dim>*>::const_iterator is = ref.NodesBegin(); is != ref.NodesEnd(); ++is )
        {
            (*is)->Read(temperature_key,temperature);
            (*is)->Read(pressure_key,pressure);

            volatile_dissolved_eq()=6.833e-7*pow(temperature(),2)-0.001356*temperature()+0.7273; // this coincides with 'alpha'
            volatile_dissolved_eq()*=pow(pressure()*1.0e-5,-6.831e-7*pow(temperature(),2)+0.001416*temperature()-0.1255); // this value is in percent
            volatile_dissolved_eq()/=100; // now it is a fraction

            if (volatile_dissolved_eq()<minimum_volatile_dissolved)
            {
                minimum_volatile_dissolved = volatile_dissolved_eq();
            }
        }
    }

    for( typename std::vector<Node<dim>*>::const_iterator it = ref.NodesBegin(); it != ref.NodesEnd(); ++it )
    {
        (*it)->Read(bulk_volume_key,bulk_volume);
        (*it)->Read(volumetric_mass_key,volumetric_mass);
        (*it)->Read(pressure_key,pressure);
        (*it)->Read(lithostatic_pressure_key,lithostatic_pressure);
        (*it)->Read(melt_density_key,melt_density);
        (*it)->Read(temperature_key,temperature);

        volatile_dissolved_eq()=6.833e-7*pow(temperature(),2)-0.001356*temperature()+0.7273; // this coincides with 'alpha'
        volatile_dissolved_eq()*=pow(pressure()*1.0e-5,-6.831e-7*pow(temperature(),2)+0.001416*temperature()-0.1255); // this value is in percent
        volatile_dissolved_eq()/=100; // now it is a fraction


        // Saturation model 0: The minimum value of the "initial wtpercent volatile in magma" and the thermodynamical volatile saturation is used.
        if (saturation_model == 0)
        {
            if (initial_wtpercent_volatile_in_magma_/100.<volatile_dissolved_eq()) // the read value is actually a fraction, the name is misleading...
            {
                volatile_dissolved_fraction()=initial_wtpercent_volatile_in_magma_/100.; // this makes sure that the value from the config file is used. if initial_wtpercent_volatile_in_chamber is set to a high value, the eq (saturated) value is used
            }
            else
            {
                volatile_dissolved_fraction()=volatile_dissolved_eq();
            }
        }

        // Saturation model 1: Dissolved volatile fraction is set to a fraction of the thermodynamic volatile saturation
        else if (saturation_model == 1)
        {
            volatile_dissolved_fraction = saturation*volatile_dissolved_eq();
        }

        // Saturation model 2: Dissolved volatile fraction is set to a fraction of the minimal thermodynamic volatile saturation calculated at any node in the chamber
        else if (saturation_model == 2)
        {
            volatile_dissolved_fraction() = saturation*minimum_volatile_dissolved;
        }

        else
        {
            cerr<<"Magma saturation model is not configured correctly."<<endl;
            string dummy;
            cerr << " Press ENTER to continue... ";
            getline( cin, dummy );
        }

        if (volatile_dissolved_eq()!=0.)
        {
            relative_saturation() = std::min(1.,volatile_dissolved_fraction()/volatile_dissolved_eq());
        }
        else
        {
            relative_saturation()=0.;
        }

        // new parametrization for silicic melt (see separate jupyter notebook)
        // it is valid for silicic melt, and uses volatile_dissolved_fraction() as input
        // the previous one was for dacitic melt and for percentages.

        // for Gilgai Granite (leucogranite) as a proxy for the Mole Granite. Best fit compiled from data by Juniper & Kleeman (1979):
        //  c1 = 2.416e3    c2 = 1.982e3    c3 = 1.786e-7
        //  c4 = 1.162e-7   c5 = 9.316e-2   c6 = 1.862
        //  density = c1 - (c2*xh2o) + (c3 + (c4*xh2o))*pressure - (c5 + (c6*xh2o))*temperature

        double c1(2.416e3), c2(1.982e3), c3(1.786e-7), c4(1.162e-7), c5(9.316e-2), c6(1.862);
        double T(temperature()),P(0.),X(volatile_dissolved_fraction());
        if(pressure()<lithostatic_pressure())
        {
            P = lithostatic_pressure();}
        else
        {
            P = pressure();}

        //rhyolite
        melt_density() = c1 - (c2*X) + (c3 + (c4*X))*P - (c5+  (c6*X))*T; // TS


        volatile_fraction()=initial_volatile_vol_frac_chamber_; //the initial porosity

        melt_fraction()=1.-volatile_fraction();

        // Initial volumes calculations assuming X=0 for now
        melt_volume()= melt_fraction()*bulk_volume();
        pore_volume()=volatile_fraction()*bulk_volume();

        minimum_pore_volume() = pore_volume();


        crystal_volume()=0.0;
        delta_melt_volume()=0.0;
        delta_crystal_volume()=0.0;

        // Initial fractions
        crystal_fraction()=crystal_volume()/bulk_volume();
        volatile_saturation()=pore_volume()/(bulk_volume()-crystal_volume());

        // Initial masses
        melt_mass()=melt_volume()*melt_density();
        crystal_mass()=crystal_volume()*crystal_density_;
        volatile_dissolved_mass()=melt_mass()*volatile_dissolved_fraction();

        minimum_mass_fluid()=volumetric_mass()*pore_volume();

        // Initial porosities
        nodal_porosity()=pore_volume()/bulk_volume();

        //Initial salt
        msp() = volumetric_mass()*salinity_magmatic_fluids*0.01;
        magmatic_salt() = volumetric_mass()*salinity_magmatic_fluids*0.01;

        //magmatic_ratio() = 1.;

        magmatic_mass() = volumetric_mass();

        ddt_magmatic_ratio() = 0.;

        //Initial rock compressibility
        nodal_compressibility_rock()=compressibility_magma_;
        if (crystal_fraction()>0.7 /*&& volatile_fraction()>(1.0-crystal_fraction())*0.9*/)//Comment Benoit Nov 2022
        {
            nodal_compressibility_rock()=compressibility_rock_;
        }

        cerr<<"Crystal fraction   = "<<crystal_fraction()<<endl;
        cerr<<"Melt fraction      = "<<melt_fraction()<<endl;
        cerr<<"Fluid fraction     = "<<volatile_fraction()<<endl;
        cerr<<"Total fractions    = "<<crystal_fraction()+melt_fraction()+volatile_fraction()<<endl<<endl;
        cerr<<"Crystal volume     = "<<crystal_volume()<<endl;
        cerr<<"Melt volume        = "<<melt_volume()<<endl;
        cerr<<"Fluid volume       = "<<pore_volume()<<endl;
        cerr<<"Mini mass fluid    = "<<minimum_mass_fluid()<<endl;
        cerr<<"Total volume       = "<<crystal_volume()+melt_volume()+pore_volume()<<endl;
        cerr<<"Bulk volume        = "<<bulk_volume()<<endl<<endl;

        (*it)->Store(minimum_mass_fluid_key,minimum_mass_fluid);
        (*it)->Store(delta_melt_volume_key,delta_melt_volume);
        (*it)->Store(melt_volume_key,melt_volume);
        (*it)->Store(delta_crystal_volume_key,delta_crystal_volume);
        (*it)->Store(crystal_volume_key,crystal_volume);
        (*it)->Store(pore_volume_key,pore_volume);
        (*it)->Store(minimum_pore_volume_key,minimum_pore_volume);
        (*it)->Store(nodal_porosity_key,nodal_porosity);
        (*it)->Store(melt_mass_key,melt_mass);
        (*it)->Store(crystal_mass_key,crystal_mass);
        (*it)->Store(volatile_dissolved_mass_key,volatile_dissolved_mass);
        (*it)->Store(volatile_dissolved_eq_key,volatile_dissolved_eq);
        (*it)->Store(volatile_dissolved_fraction_key,volatile_dissolved_fraction);
        (*it)->Store(volatile_saturation_key,volatile_saturation);
        (*it)->Store(relative_saturation_key,relative_saturation);
        (*it)->Store(crystal_fraction_key,crystal_fraction);
        (*it)->Store(melt_fraction_key,melt_fraction);
        (*it)->Store(volatile_fraction_key,volatile_fraction);
        (*it)->Store(melt_density_key,melt_density);
        (*it)->Store(nodal_compressibility_rock_key,nodal_compressibility_rock);
        (*it)->Store(msp_key, msp);
        (*it)->Store(mm_key,  magmatic_mass );
        (*it)->Store(mms_key,  magmatic_salt );
        (*it)->Store(ddt_magmatic_ratio_key,ddt_magmatic_ratio);
    }
}

template<size_t dim>
void MagmaModel<dim>::CalculateAverageCrystallinity(Model<dim>* mdl, std::string regionname)
{
    const Region<dim>& ref = mdl->Region(regionname);

    mdl->Read(average_crystallinity_key,average_crystallinity);

    double previous_average_crystallinity(average_crystallinity());
    double total_crystal_volume(0.);
    double total_melt_volume(0.);

    for( typename std::vector<Node<dim>*>::const_iterator is = ref.NodesBegin(); is != ref.NodesEnd(); ++is )
    {
        (*is)->Read(crystal_volume_key,crystal_volume);
        (*is)->Read(melt_volume_key,melt_volume);

        total_crystal_volume += crystal_volume();
        total_melt_volume += melt_volume();
    }

    average_crystallinity() = total_crystal_volume/(total_crystal_volume+total_melt_volume);
    delta_average_crystallinity() = average_crystallinity() - previous_average_crystallinity;

    mdl->Store(average_crystallinity_key,average_crystallinity);
    mdl->Store(delta_average_crystallinity_key,delta_average_crystallinity);

}


template<size_t dim>
void MagmaModel<dim>::WithHaliteTrend(double fraction_magma)
{
    with_halite_trend = true;
    fraction_KCl_magma = fraction_magma;
}

template<size_t dim>
void MagmaModel<dim>::SetTimeIncrement( double& time_step )
{
    dt = time_step;
}

template<size_t dim>
void MagmaModel<dim>::ReadVariables(Node<dim>* node)
{
    node->Read(liquidus_temperature_key,liquidus_temperature);
    node->Read(solidus_temperature_key,solidus_temperature);
    node->Read(temperature_key,temperature);
    node->Read(previous_temperature_key,previous_temperature);
    node->Read(total_mass_volatile_produced_key,total_mass_volatile_produced);
    node->Read(bulk_volume_key,bulk_volume);
    node->Read(minimum_pore_volume_key,minimum_pore_volume);
    node->Read(pore_volume_key,pore_volume);
    node->Read(minimum_mass_fluid_key,minimum_mass_fluid);
    node->Read(melt_volume_key,melt_volume_before);
    node->Read(crystal_volume_key,crystal_volume_before);
    node->Read(volatile_dissolved_mass_key,volatile_dissolved_mass_before);
    node->Read(volatile_dissolved_fraction_key,volatile_dissolved_fraction_before);
    node->Read(pore_volume_key,pore_volume_before);
    node->Read(nodal_porosity_key,nodal_porosity_before);
    node->Read(pressure_key,pressure);
    node->Read(melt_mass_key,melt_mass_before);
    node->Read(crystal_mass_key,crystal_mass_before);
    node->Read(crystal_fraction_key,crystal_fraction);
    node->Read(melt_fraction_key,melt_fraction);
    node->Read(volatile_fraction_key,volatile_fraction);
    node->Read(volatile_saturation_key,volatile_saturation);
    node->Read(nodal_density_rock_key,nodal_density_rock);
    node->Read(threshold_pressure_key,threshold_pressure);
    node->Read(lithostatic_pressure_key,lithostatic_pressure);
    node->Read( volumetric_mass_key,  volumetric_mass );
    node->Read( ml_key,  ml );
    node->Read( mv_key,  mv );
    node->Read( Htp_key, Htp );
    node->Read( msp_key, msp );
    node->Read(mm_key,  magmatic_mass );
    node->Read(mf_key,  fluid_mass );
    node->Read(mms_key,  magmatic_salt );
    node->Read(m_key,  mass );
    node->Read(H_key, enthalpy );
    node->Read(hCl_key, hCl );
    node->Read(hCv_key, hCv );
    if (with_halite_trend) node->Read(kcl_ms_key,  mass_KCl );
    node->Read(total_volume_deformation_key,  total_volume_deformation );
    node->Read(melt_density_key,  melt_density );
    node->Read(minimum_nodal_density_rock_key,  minimum_nodal_density_rock );
    node->Read(wt_key,  wt );
    node->Read(depth_key,  depth );
    node->Read(porous_flag_key, porous_flag );
    node->Read(crystallinity_key,crystallinity_before);
}

template<size_t dim>
void MagmaModel<dim>::StoreVariables(Node<dim>* node)
{
    node->Store(volumetric_mass_key, volumetric_mass);
    node->Store(Htp_key, Htp );
    node->Store(msp_key, msp );
    node->Store(m_key,  mass );
    node->Store(H_key, enthalpy );
    node->Store(hCl_key, hCl );
    node->Store(hCv_key, hCv );
    node->Store(mm_key,  magmatic_mass );
    node->Store(mf_key,  fluid_mass );
    node->Store(mms_key,  magmatic_salt );

    if (with_halite_trend) node->Store(kcl_ms_key,  mass_KCl );

    node->Store(temperature_difference_key,temperature_difference);
    node->Store(crystallinity_key,crystallinity);
    node->Store(delta_melt_volume_key,delta_melt_volume);
    node->Store(melt_volume_key,melt_volume);
    node->Store(delta_crystal_volume_key,delta_crystal_volume);
    node->Store(crystal_volume_key,crystal_volume);
    node->Store(melt_mass_key,melt_mass);
    node->Store(crystal_mass_key,crystal_mass);
    node->Store(volatile_dissolved_fraction_key,volatile_dissolved_fraction);
    node->Store(volatile_dissolved_mass_key,volatile_dissolved_mass);
    node->Store(volatile_dissolved_eq_key,volatile_dissolved_eq);
    node->Store(mass_volatile_produced_key,mass_volatile_produced);
    node->Store(total_mass_volatile_produced_key,total_mass_volatile_produced);
    node->Store(pore_volume_key,pore_volume);
    node->Store(volume_change_factor_key,volume_change_factor);
    node->Store(volatile_saturation_key,volatile_saturation);
    node->Store(relative_saturation_key,relative_saturation);
    node->Store(crystal_fraction_key,crystal_fraction);
    node->Store(melt_fraction_key,melt_fraction);
    node->Store(volatile_fraction_key,volatile_fraction);
    node->Store(nodal_porosity_key,nodal_porosity);
    node->Store(bulk_volume_key,bulk_volume);
    node->Store(volume_deformation_key,volume_deformation);
    node->Store(total_volume_deformation_key,total_volume_deformation);
    node->Store(rock_density_scaling_key,rock_density_scaling);
    node->Store(nodal_density_rock_key,nodal_density_rock);
    node->Store(pressure_key,pressure);
    node->Store(magma_density_key,magma_density);
    node->Store(melt_density_key,melt_density);
    node->Store(nodal_compressibility_rock_key,nodal_compressibility_rock);
    node->Store(minimum_nodal_density_rock_key,minimum_nodal_density_rock);
    node->Store(minimum_pore_volume_key,minimum_pore_volume);
}

template<size_t dim>
bool MagmaModel<dim>::ConsistencyChecks1(bool pause_at_end, double density_deformation, double deformation_overpressure)
{
    if (volume_change_factor()>2. or volume_change_factor()<0.5)
    {
        cerr<<"!!!!!!!!! High volume change factor !!!!!!!!!"<<endl;
    }

    if (volumetric_mass()<=0.0)
    {
        cerr<<"!!!!!!!!!volumetric_mass is below zero!!!!!!!!!"<<endl;
        if(pore_volume()<0.0) cerr<<"!!!!!!!!!pore volume is below zero!!!!!!!!! "<<pore_volume()<<endl;
        if(minimum_pore_volume()<0.0) cerr<<"!!!!!!!!!mini pore volume is below zero!!!!!!!!! "<<minimum_pore_volume()<<endl;
        string dummy;
        cerr << " press ENTER to continue... ";
        getline( cin, dummy );
        pause_at_end=true;
    }

    if (minimum_pore_volume()>(pore_volume()+0.1))
    {
        cerr<<"!!!!!!!!!The pore volume is definitely below mini!!!!!!!!!"<<endl;
        cerr<<"the difference is "<<pore_volume()-minimum_pore_volume()<<endl;
        string dummy;
        cerr << " press ENTER to continue... ";
        getline( cin, dummy );
        pause_at_end=true;
    }

    if(density_deformation<0.0 && deformation_overpressure>0.0)
    {
        cerr<<"density_deformation<0.0 && deformation_overpressure>0.0"<<endl;
        cerr<<"Overpressure with underpressure?"<<endl ;
        string dummy;
        cerr << " press ENTER to continue... ";
        getline( cin, dummy );
        pause_at_end=true;
    }

    if(!approximatelyEqual(density_deformation,0.0,1.e-6) && !approximatelyEqual(deformation_overpressure,0.0,1.e-6))
    {
        cerr<<"density_deformation!=0.0 && deformation_overpressure!=0.0"<<endl;
        cerr<<"Overpressure with underpressure?"<<endl ;
        string dummy;
        cerr << " press ENTER to continue... ";
        getline( cin, dummy );
        pause_at_end=true;
    }

    return pause_at_end;
}

template<size_t dim>
bool MagmaModel<dim>::ConsistencyChecks2(bool pause_at_end)
{
    if (volatile_fraction()>(1.0-crystal_fraction())*0.9 && crystal_fraction()<0.4)
    {
        cerr<<"crystal fraction is below 0.4 while volatile saturation is above 0.9"<<endl;
        cerr<<"crystal fraction is: "<<crystal_fraction()<<endl;
        cerr<<"melt fraction is   : " << melt_fraction()<<endl;
        cerr<<"volatile fraction is: "<<volatile_fraction()<<endl;
        cerr<<"volatile_saturation() = "<<volatile_saturation()<<endl;
        cerr<<"volume deformation is: "<<volume_deformation()<<endl;
        //system("pause");pause_at_end=true;
    }

    if(minimum_pore_volume()<0.0)
    {
        cerr<<"Mini pore volume below 0!!!!!"<<endl;
        string dummy;
        cerr << " press ENTER to continue... ";
        getline( cin, dummy );
        pause_at_end=true;
    }

    if(nodal_porosity()<0.0)
    {
        cerr<<"Porosity below 0.0!!!!!"<<endl;
        string dummy;
        cerr << " press ENTER to continue... ";
        getline( cin, dummy );
        pause_at_end=true;
    }

    if(volumetric_mass()<0.0)
    {
        cerr<<"Negative fluid density!!!!!"<<endl;
        string dummy;
        cerr << " press ENTER to continue... ";
        getline( cin, dummy );
        pause_at_end=true;
    }

    if(volume_change_factor()<0.0)
    {
        cerr<<" Negative pore volume change factor!!!!!"<<endl;
        string dummy;
        cerr << " press ENTER to continue... ";
        getline( cin, dummy );
        pause_at_end=true;
    }

    return pause_at_end;
}


template<size_t dim>
bool MagmaModel<dim>::FractionSumChecks(bool pause_at_end)
{
    crystal_fraction()=crystal_volume()/bulk_volume();
    melt_fraction()=melt_volume()/bulk_volume();
    volatile_fraction()=pore_volume()/bulk_volume();

    if (!(essentiallyEqual(crystal_fraction()+melt_fraction()+volatile_fraction(),1.0,1.0e-06)))
    {
        cerr<<"FRACTIONS DO NOT SUM UP TO 1 (D)!!!!!!!!!!!"<<endl ;
        cerr << "crystal fraction     : " << crystal_fraction()<<endl;
        cerr << "melt fraction        : " << melt_fraction()<<endl;
        cerr << "volatile fraction    : " << volatile_fraction() << endl;
        cerr << "total fractions      : " << crystal_fraction()+melt_fraction()+volatile_fraction()<<endl<<endl;
        cerr << "crystallinity        : " << crystallinity()<<endl;
        cerr << "melt volume          : " << melt_volume()<<endl;
        cerr << "crystal volume       : " << crystal_volume()<<endl;
        cerr << "melt volume change   : " << melt_volume()-melt_volume_before()<<endl;
        cerr << "crystal volume change: " << delta_crystal_volume()<<endl;
        cerr << "bulk volume          : " << bulk_volume()<<endl;
        string dummy;
        cerr << " press ENTER to continue... ";
        getline( cin, dummy );
        pause_at_end=true;
    }

    return pause_at_end;
}



template class MagmaModel<1U>;
template class MagmaModel<2U>;
template class MagmaModel<3U>;

} // end namespace csmp
