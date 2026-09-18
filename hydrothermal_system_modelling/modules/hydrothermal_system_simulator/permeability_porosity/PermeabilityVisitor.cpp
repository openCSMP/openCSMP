// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "PermeabilityVisitor.h"

using namespace std;

namespace csmp
{
template<uint32_t dim>
PermeabilityVisitor<dim>::PermeabilityVisitor(Model<dim> &model)
    : Visitor<dim>( MODEL, ELEMENT ),
    model                     ( model ),
    temperature_e_key_        ( model.Database().StorageKey("temperature element") ),
    permeability_key_         ( model.Database().StorageKey("permeability") ),
    vertical_permeability_key_( model.Database().StorageKey("vertical permeability") ),
    horizontal_permeability_key_( model.Database().StorageKey("horizontal permeability") ),
    permeability_tensor_key_ ( model.Database().StorageKey("permeability tensor") ),
    anisotropy_factor_key_   ( model.Database().StorageKey("permeability anisotropy factor") ),
    pore_volume_element_key_ ( model.Database().StorageKey("pore volume element") ),
    bulk_volume_element_key_ ( model.Database().StorageKey("bulk volume element") ),
    porosity_key_            ( model.Database().StorageKey("porosity") ),
    k_start_key_             ( model.Database().StorageKey("BDT start permeability element") ),

    log_k_increase_key_        ( model.Database().StorageKey("log permeability increase") ),
    fluid_pressure_key_        ( model.Database().StorageKey("fluid pressure") ),
    failure_pressure_key_      ( model.Database().StorageKey("failure pressure") ),
    fracturing_ref_key_        ( model.Database().StorageKey("fracturing reference") ),
    fracturing_events_key_     ( model.Database().StorageKey("fracturing events") ),
    gp_sc_key_                 ( model.Database().StorageKey("gradP scaling") ),

    permeability_ID_key_      ( model.Database().StorageKey("permeability ID") ),
    region_ID_key_            ( model.Database().StorageKey("region ID") ),


    depth_dependent_                        ( false ),
    with_reference_depth_                   ( false ),
    temperature_dependent_                  ( false ),
    mineral_dependent_                      ( false ),
    pore_fluid_factor_dependent_            ( false ),
    average_Pff_                            ( false ),
    hydro_fracturing_                       ( false ),
    immediate_closure_                      ( true ),
    avoid_shallow_and_reverse_gradients_    ( true ), // hard-coded !
    overpressured_                          ( false ),
    anisotropic_k_                          ( false ),
    no_frac_below_intrusion_                ( false ),
    with_intrusion_                         ( false ),

    log_k_back_          ( -14. ),
    log_k_max_           ( -10. ),
    log_k_D_             ( 0. ), // this value is defined in the config file
    log_k_T_             ( 0. ),
    log_k_Pff_           ( 0. ),
    log_k_ductile_       ( -17. ),
    log_k_min_           ( -22. ),
    log_k_Hydrofracture_ ( 0. ),
    reference_depth_     ( 0. ),


    hydro_failure_temperature_   (360.),
    ductile_failure_temperature_ (400.),
    lithos_failure_temperature_  (500.)

{

}

template<uint32_t dim>
PermeabilityVisitor<dim>::~PermeabilityVisitor()
{

}

template<uint32_t dim>
void PermeabilityVisitor<dim>::Visit( Model<dim> *m )
{

}

template<uint32_t dim>
void PermeabilityVisitor<dim>::Visit( Region<dim> *r )
{

}


template<uint32_t dim>
void PermeabilityVisitor<dim>::Visit(Element<dim> *e)
{
    // Region ID is used to set permeability of a cooled fractured intrusion
    // CARFEUL: the intrusion region ID set via: SetIntrusionRegionID( uint32_t intrusion_ID )
    ScalarVariable region_ID;
    if (with_intrusion_)
        e->Read(region_ID_key_, region_ID);

    // initialise permeability ID used to visualise which permeability law is dominating and written on the element
    // -1 = minimum permeability because rock is not fracturable;
    //  0 = background permeability; 1 = depth-dependent; 2 = temperature-dependent; 3 = hydro-fracturing; 4 = pore-fluid-factor-dependent; 5 = mineralisation-dependent
    // 11 = depth-dependent after CalculateGradPScaling; 22 = temperature-dependent after CalculateGradPScaling
    k_ID_ = 0.;

    // Define anisotropy factor: values >1 = increased vertical permeability; values<1 = increased horizontal permeability
    anisotropy_factor_ = 1.;

    if (anisotropic_k_)
    {
        e->Read(anisotropy_factor_key_, anisotropy_factor_);
    }

    ScalarVariable fracturing_events, k_start;

    // reset k = 0 for both depth-depenendet and temperature-dependent permeability
    log_k_D_ = log_k_T_ = 0.;

    e->Read(fracturing_ref_key_, fracturing_ref_); // fracturing reference; fracturing_ref<1 --> crystal fraction <0.7 --> no fracturing
    e->Read(fracturing_events_key_, fracturing_events);

    e->Read(permeability_key_, k_); // permeability
    log_k_ = log_k_previous_ = log10(k_());

    e->Read(k_start_key_, k_start); // initial permeability assigned to regions before applying permeability visitor
    log_k_start_ = log10(k_start());


    double depth;
    depth = e->BaryCenter()[1];  // y-coordinate in m; surface has to be at y=0

    // Calculate pore fluid factor used for PFF-dependent permeability and hydrofracturing
    CalculatePoreFluidFactor(e, average_Pff_); // true = average Pff; false = max Pff

    if (definitelyLessThan(fracturing_ref_(), 1., numeric_limits<double>::epsilon()))
    {
        // Element is not fracturable (i.e., crystal fraction < 0.7).
        // log_k is set to minimum value.

        log_k_ = log_k_min_;
        k_ = pow(10., log_k_);
        k_ID_ = -1.;
    }
    else //Element is fracturable (i.e., crystal fraction >= 0.7).
    {
        //! Calculate permeability
        if (depth_dependent_)
        {
            log_k_D_ = CalculateDepthDependentPermeability(e);
            log_k_D_ = std::min( log_k_D_, log_k_start_); // this is to make sure that permeability is not larger than pre-assigned permeability due to depth dependency
        }

        if (mineral_dependent_)
            log_k_mineral_ = CalculatePoroAndPermChangeByMineralPrecipitationAndDissolution(e);

        if (temperature_dependent_)
            log_k_T_ = CalculateTemperatureDependentPermeability(e);

        if (pore_fluid_factor_dependent_)
            log_k_Pff_ = CalculatePoreFluidFactorDependentPermeability(e);

        if (hydro_fracturing_)
            log_k_Hydrofracture_ = CalculateHydrofracturedPermeability(e);

        ////////////////////////////
        if (depth_dependent_)
        {
            log_k_ = log_k_D_;
            k_ID_ = 1.;
        }

        if (mineral_dependent_)
        {
            if (definitelyLessThan(log_k_mineral_, log_k_, numeric_limits<double>::epsilon()))
                k_ID_ = 5.;

            log_k_ = std::min(log_k_, log_k_mineral_);
        }

        if (temperature_dependent_)
        {
            if (definitelyLessThan(log_k_T_, log_k_, numeric_limits<double>::epsilon()))
                k_ID_ = 2.;

            if ( with_intrusion_ && essentiallyEqual(region_ID(), intrusion_ID_, numeric_limits<double>::epsilon()))
                log_k_ = std::min(log_k_start_, log_k_T_); // no depth-dependency in case of magma chamber
            else
                log_k_ = std::min(log_k_, log_k_T_);
        }

        if (hydro_fracturing_ && overpressured_)
        {
            // hydrofracturing (potentially healing) dominates permeability
            if (definitelyGreaterThan(log_k_Hydrofracture_, log_k_, numeric_limits<double>::epsilon()))
            {
                // apply if fracturing is generally allowed or if the element is located above the minimum y coordinate of the intrusion
                if ( !no_frac_below_intrusion_ || depth >= chamber_y_min_)
                {
                    log_k_ = log_k_Hydrofracture_;
                    fracturing_events += 1.;
                    k_ID_ = 3.;
                }
            }

            // hydrofracturing is limited to 2 orders of magnitude higher than depth-dependent
            if (depth_dependent_)
            {
                if (definitelyGreaterThan(log_k_, log_k_D_ + 2., numeric_limits<double>::epsilon()))
                    log_k_ = std::min(log_k_, log_k_D_ + 2.);
            }
        } // end hydro_fracturing

        else if (pore_fluid_factor_dependent_ && !overpressured_) // pore fluid factor
        {
            if (definitelyGreaterThan(log_k_Pff_, log_k_, numeric_limits<double>::epsilon()))
                k_ID_ = 4.;

            log_k_ = std::max(log_k_, log_k_Pff_);
        }

    } // end fracturing_ref >= 0.7

    k_ = pow(10., log_k_);


    //! Define horizontal and vertical permeability
    if (!anisotropic_k_) // isotropic permeability
    {
        if (avoid_shallow_and_reverse_gradients_
            && definitelyGreaterThan(fracturing_ref_(), 0., numeric_limits<double>::epsilon()))
        {
            // perform correction to avoid shallow or reverse gradients
            // gradP_sc is caulculated in PressureGradientVisitor
            CalculateGradPScaling(e);
        }

        kV_ = k_();
        kH_ = k_();
    }

    else // isotropic permeability
    {
        if (k_ID_() == 1 || k_ID_() == 0) // depth-dependent or no dependency
        {
            if (avoid_shallow_and_reverse_gradients_ && definitelyGreaterThan(fracturing_ref_(), 0., numeric_limits<double>::epsilon()))
                CalculateGradPScaling(e);

            if (anisotropy_factor_() >= 1.) // increased vertical permeability
            {
                kH_ = k_();
                kV_ = k_() * anisotropy_factor_;
            }
            else   // increased horizontal permeability
            {
                kH_ = k_() / anisotropy_factor_;
                kV_ = k_();
            }
        }

        else if (k_ID_() == 3 || k_ID_() == 4) // hydrofracturing or PFF: assuming a higher vertical permeability increase than in the horizontal direction
        {
            if (avoid_shallow_and_reverse_gradients_ && definitelyGreaterThan(fracturing_ref_(), 0., numeric_limits<double>::epsilon()))
                CalculateGradPScaling(e);

            kV_ = k_();
            kH_ = k_() / pore_fluid_factor_;
        }

        // temperature-dependent or Pff-dependent or mineral-dependent: isotropic permeability
        // Could be improved: t-dependency could be calculated for different start permeability --> anisotropic k-decrease between T_low and T_high
        else
        {
            if (avoid_shallow_and_reverse_gradients_ && definitelyGreaterThan(fracturing_ref_(), 0., numeric_limits<double>::epsilon()))
                CalculateGradPScaling(e);

            kV_ = k_();
            kH_ = k_();
        }
    } // end isotropic permeability


    //! Store variables
    // store permeability
    e->Store(permeability_key_, k_);
    e->Store(vertical_permeability_key_, kV_);
    e->Store(horizontal_permeability_key_, kH_);

    // Permeability changes
    // TO DO: Add permeability changes for vertical and horizontal permeability
    e->Read(log_k_increase_key_, log_k_increase_);
    log_k_increase_ += log_k_ - log_k_previous_;
    e->Store(log_k_increase_key_, log_k_increase_);

    // Permeability ID
    e->Store(permeability_ID_key_, k_ID_);

    // Fracturing event
    e->Store(fracturing_events_key_, fracturing_events);

} // end visit



//! Functions called within the visitor
template<uint32_t dim>
double PermeabilityVisitor<dim>::CalculateDepthDependentPermeability( Element<dim> *e )
{
    double depth_in_km;

    if (with_reference_depth_)
        depth_in_km = (reference_depth_ - e->BaryCenter()[1]) / 1000;  // surface at y=reference_depth
    else
        depth_in_km = -1. * e->BaryCenter()[1] / 1000;  // surface has to be at y=0

    if (depth_in_km < 0)
        depth_in_km = 1; // avoids NaNs in case of topography

    return -14 - 3.2 * log10( depth_in_km ); // Ingebritsen & Manning (https://doi.org/10.1130/0091-7613(1999)027<1107:GIOAPD>2.3.CO;2)
}


template<uint32_t dim>
double PermeabilityVisitor<dim>::CalculateTemperatureDependentPermeability( Element<dim> *e )
{
    // temperature-dependent permeability following Hayba & Ingebritsen, 1997 (https://doi.org/10.1029/97JB00552)
    double T_onset    = hydro_failure_temperature_;
    double T_ductile  = ductile_failure_temperature_;
    double T_min      = lithos_failure_temperature_;
    double av_k        ( 0. );

    // element temperature (T at element bary center) is used as it is the most accurate
    e->Read(temperature_e_key_, T_);

    if ( T_() < T_onset )
        av_k = log_k_start_;

    else if ( T_() < T_ductile )
        av_k = ( log_k_start_ * ( T_ductile - T_() ) + log_k_ductile_ * ( T_() - T_onset ) ) / ( T_ductile - T_onset );

    else if ( T_() < T_min )
        av_k = ( log_k_ductile_ * ( T_min - T_() ) + log_k_min_ * ( T_() - T_ductile ) ) / ( T_min - T_ductile );

    else
        av_k = log_k_min_;

    return av_k;

}


template<uint32_t dim>
double PermeabilityVisitor<dim>::CalculatePoroAndPermChangeByMineralPrecipitationAndDissolution( Element<dim> *e )
{
    csmp::Index deltaQuartzPrecipitated_element_key ( model.Database().StorageKey("delta quartz precipitated element") );
    csmp::Index deltaGoldPrecipitated_element_key ( model.Database().StorageKey("delta gold precipitated element") );

    e->Read(pore_volume_element_key_, pore_volume_element_);
    e->Read(bulk_volume_element_key_, bulk_volume_element_);
    e->Read(porosity_key_, porosity_);             //kg/m³

    e->Read(deltaGoldPrecipitated_element_key, deltaGoldPrecipitated_element_);             //kg/m³
    e->Read(deltaQuartzPrecipitated_element_key, deltaQuartzPrecipitated_element_);         //kg/m³

    //update porosity_ and permeability according to description in Scott & Driesner, 2018 (https://doi.org/10.1155/2018/6957306)
    double previous_porosity = porosity_();
    double critical_porosity = 0.5 * previous_porosity;

    porosity_() -= deltaQuartzPrecipitated_element_() * bulk_volume_element_() / quartz_density_ / bulk_volume_element_();
    porosity_() -= deltaGoldPrecipitated_element_() * pore_volume_element_() / au_density_ / bulk_volume_element_();
    porosity_() = std::max(porosity_(), 0.005);

    double initial_permeability(k_());
    double permfactor(0.);
    permfactor = (porosity_() - critical_porosity) / (previous_porosity - critical_porosity);
    initial_permeability *= pow(1 - (1 - pow(permfactor, 1.58)), 0.46);
    initial_permeability = std::max(initial_permeability, 1e-22);

    e->Store(porosity_key_, porosity_);

    double log_updated_perm = log10(initial_permeability);

    return log_updated_perm;
}


template<uint32_t dim>
void PermeabilityVisitor<dim>::CalculatePoreFluidFactor( Element<dim> *e, bool average_pff )
{
    // Using maximum pff of element instead of average pff enhances fracturing
    // as it uses the maximum pff at each element (not the average).

    pore_fluid_factor_ = 0.;

    e->NodePropertyVector( fluid_pressure_key_, fluid_pressure_vec_ );
    e->NodePropertyVector( failure_pressure_key_, failure_pressure_vec_ );

    if (average_pff)
    {
        for ( uint32_t i = 0; i < e->Nodes(); i++ )
            pore_fluid_factor_ += fluid_pressure_vec_[i]() / failure_pressure_vec_[i]();

        pore_fluid_factor_ /= static_cast<double>( e->Nodes() );
    }

    else // calculate max pore fluid factor
    {
        for ( uint32_t i = 0; i < e->Nodes(); i++ )
            pore_fluid_factor_ = std::max(pore_fluid_factor_, fluid_pressure_vec_[i]() / failure_pressure_vec_[i]());
    }

    // check for overpressure
    if ( pore_fluid_factor_ > 1.0 )
        overpressured_ = true;
    else
        overpressured_ = false;

}


template<uint32_t dim>
double PermeabilityVisitor<dim>::CalculatePoreFluidFactorDependentPermeability( Element<dim> *e )
{
    // Dependence of permeability on pore-fluid factor assuming a strain of 10% shortening during room temperature deformation
    // after Cox, 2005 (Fig. 9a) doi.org/10.5382/AV100.04

    double log_k_pff(0.);
    pff_ = 0.;

    pff_ = std::min(pore_fluid_factor_, 1.0);
    pff_ = std::max(pff_, 0.3);
    pff_ -= 0.3;
    pff_ /= 0.7;

    log_k_pff = std::min(log_k_D_, log_k_back_);

    if (log_k_T_ < log_k_)
    {
        log_k_pff *= (pff_ * pff_);
        log_k_pff += (1. - pff_ * pff_) * log_k_T_;
    }

    return log_k_pff;
}


template<uint32_t dim>
double PermeabilityVisitor<dim>::CalculateHydrofracturedPermeability( Element<dim> *e )
{
    // Case 1:
    //================================
    // Element is not fracturable
    //================================

    if ( definitelyLessThan( fracturing_ref_(), 1., numeric_limits<double>::epsilon() ) )
        return log_k_min_;

    // end Case 1

    // Case 2:
    //================================
    // No overpressure
    //================================

    if ( !overpressured_ )
    {
        if ( !immediate_closure_ )
            return (log_k_ / pore_fluid_factor_);
        else
            return log_k_min_;
    } // end Case 2


    // Case 3:
    //================================
    // Overpressure
    // Permeability will be scaled
    //================================

    else //hydrofracturing
    {
        double k_temp = k_();
        k_temp *= pore_fluid_factor_;
        k_temp *= pore_fluid_factor_;

        // limit to log_k_max
        if ( log10( k_temp ) > log_k_max_ )
            return log_k_max_;
        else
            return log10( k_temp );

    } // end new hydrofracturing

}


template<uint32_t dim>
void PermeabilityVisitor<dim>::CalculateGradPScaling(Element<dim> *e)
{
    // perform correction to avoid shallow or reverse gradients
    // recheck whether this is necessary
    // gradP_sc is caulculated in PressureGradientVisitor

    // Region ID is used to set permeability of a cooled fractured intrusion
    // CARFEUL: the intrusion region ID set via: SetIntrusionRegionID( uint32_t intrusion_ID )
    ScalarVariable region_ID;
    if (with_intrusion_)
        e->Read(region_ID_key_, region_ID);

    if (!with_intrusion_
        || !essentiallyEqual(region_ID(), intrusion_ID_, numeric_limits<double>::epsilon() )
        || definitelyLessThan(k_(), pow(10., log_k_start_), numeric_limits<double>::epsilon()) )
    {
        e->Read(gp_sc_key_, gradP_sc_);
        // Caution, this is computed in CVFEM_PressureGradientVisitor!
        double grad_factor(1.0);

        if (gradP_sc_() < 0.0)
            grad_factor = 0.0;
        else if (gradP_sc_() < 1.0)
            grad_factor = gradP_sc_();

        k_() *= grad_factor; //grad_factor between 0 and 1

        if (k_() < pow(10., std::min(log_k_T_, log_k_D_)))
        {
            k_() = pow(10., std::min(log_k_T_, log_k_D_));

            if (definitelyLessThan( log_k_T_, log_k_D_, numeric_limits<double>::epsilon() ))
                k_ID_ = 22; // temperature-dependent
            else
                k_ID_ = 11; // depth-dependent
        }
    }


}



//! Setter functions
template<uint32_t dim>
void PermeabilityVisitor<dim>::DepthDependent( bool d_dep )
{
    depth_dependent_ = d_dep;
}

template<uint32_t dim>
void PermeabilityVisitor<dim>::SetReferenceDepthTo( bool with_ref_depth, double ref_depth )
{
    with_reference_depth_ = with_ref_depth;
    reference_depth_      = ref_depth;
}


template<uint32_t dim>
void PermeabilityVisitor<dim>::TemperatureDependent( bool T_dep )
{
    temperature_dependent_ = T_dep;
}


template<uint32_t dim>
void PermeabilityVisitor<dim>::ChangeBrittleDuctileTransitionTemperature(
    bool change_values,
    double transition_start,
    double transition_ductile,
    double transition_end,
    double log_k_start,
    double log_k_duct,
    double log_k_end )
{

    if (change_values)
    {
        hydro_failure_temperature_   = transition_start;
        ductile_failure_temperature_ = transition_ductile;
        lithos_failure_temperature_  = transition_end;

        log_k_back_                  = log_k_start;
        log_k_ductile_               = log_k_duct;
        log_k_min_                   = log_k_end;
    }
}


template<uint32_t dim>
void PermeabilityVisitor<dim>::MineralDependent( bool min_dep )
{
    mineral_dependent_ = min_dep;
}


template<uint32_t dim>
void PermeabilityVisitor<dim>::PoreFluidFactorDependent( bool Pff_dep, bool average_Pff)
{
    pore_fluid_factor_dependent_ = Pff_dep;
    average_Pff_ = average_Pff;
}


template<uint32_t dim>
void PermeabilityVisitor<dim>::Hydrofracturing( bool hydro_frac, double log_perm_max, double log_perm_min)
{
    hydro_fracturing_    = hydro_frac;
    log_k_max_           = log_perm_max;
    log_k_min_           = log_perm_min;
}


template<uint32_t dim>
void PermeabilityVisitor<dim>::SetImmediateClosureTo( bool closure )
{
    immediate_closure_ = closure;
}


template<uint32_t dim>
void PermeabilityVisitor<dim>::SetIntrusionRegionID( uint32_t intrusion_ID )
{
    with_intrusion_ = true;
    intrusion_ID_ = intrusion_ID;
}


template<uint32_t dim>
void PermeabilityVisitor<dim>::AnisotropicPermeabilityTensor( bool anisotropic_k )
{
    anisotropic_k_ = anisotropic_k;
}


template<uint32_t dim>
void PermeabilityVisitor<dim>::NoFracturingBelowIntrusion( bool no_frac_below_intrusion, std::string region_name )
{
    // Prevents fracturing artifacts in lithostatically pressured host rock below the intrusion
    no_frac_below_intrusion_ = no_frac_below_intrusion;

    if (no_frac_below_intrusion_)
    {
        // get maximum depth of intrusion
        double chamber_y_max;
        model.Region(region_name).MinMaxOf( "coordinate y", chamber_y_min_, chamber_y_max);
    }
}



//! Standalone batch operations (not related to Visit function)
template<uint32_t dim>
void PermeabilityVisitor<dim>::AssignPermeabilityTensor(Model<dim> &model)
{
    ScalarVariable vert_perm (ANY, 0.0);
    ScalarVariable horiz_perm (ANY, 0.0);
    TensorVariable<dim> perm_tensor(ANY, 0.);

    const Region<dim> &mref = model.Region("Model");

    typename vector<Element<dim>*>::const_iterator eit;

    for ( eit = mref.CellsBegin(); eit != mref.CellsEnd(); ++eit )
    {
        vert_perm() = (*eit)->Read(vertical_permeability_key_);
        horiz_perm() = (*eit)->Read(horizontal_permeability_key_);

        perm_tensor (0, 0) = horiz_perm();
        perm_tensor (1, 1) = vert_perm();

        if (dim == 3U)
            perm_tensor(2, 2) = horiz_perm();

        (*eit)->Store(permeability_tensor_key_, perm_tensor);
    }
}


template<uint32_t dim>
void PermeabilityVisitor<dim>::ComputeSimpleFracturableFlag()
{
    const Region<dim> &mref = model.Region("Model");
    csmp::Index       fracKey( model.Database().StorageKey("fracturing reference") );
    csmp::Index       TKey( model.Database().StorageKey("temperature element") );
    csmp::Index       TsKey( model.Database().StorageKey("solidus temperature") );

    ScalarVariable    t (ANY, 0.0); // temperature
    ScalarVariable    Ts (ANY, 0.0); // solidus temperature
    ScalarVariable    fref (ANY, 0.0); // fracture reference

    Ts = model.Region("Model").NodeVector().front()->Read(TsKey);

    typename vector<Element<dim>*>::const_iterator eit;
    for ( eit = mref.CellsBegin(); eit != mref.CellsEnd(); ++eit )
    {
        t()    = (*eit)->Read(TKey); // read temperature
        fref() = 1.; // fracturable rock

        if (t() > Ts())
            fref() = 0.; // not fracturable -> minimum permeability applied

        (*eit)->Store(fracKey, fref);
    }
}


template class PermeabilityVisitor<1U>;
template class PermeabilityVisitor<2U>;
template class PermeabilityVisitor<3U>;

}
