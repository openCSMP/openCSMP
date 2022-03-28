#include "TwoPhaseVelocityAndVolumeFlux.h"

#include "MathOperatorLHS.h"
#include "Model.h"
#include "TwoPhaseModel.h"
#include "STL_utilities.h"
#include "Element.h"
#include "Face.h"
#include "variableOperations.h"

using namespace std;

namespace csmp {

/** Constructs the post-processing math operator for velocity/ flux calculations.

@section arguments Input Arguments 

Three character strings for the input variables of the calculation
'conductivity', 'porosity', and 'fluid pressure', a boolean
variable to toggle node averaging on in which case the velocities and
fluxes are extrapolated to the nodes.  

@section messages Messages 

The names of the output variables are hard-wired in this operator. Thus,
checks are performed whether these variables actually exist in the
variable database and errors are reported if they don't or if they
have the wrong placement or type.  

tested:  */


template<uint32_t dim,class CELL>
TwoPhaseVelocityAndVolumeFlux<dim,CELL>::TwoPhaseVelocityAndVolumeFlux(const Model<dim>& sg,
                                                                   TwoPhaseModel<dim>& satfunc,
                                                                   const char* oper,           // conductivity
                                                                   const char* basic,          // porosity
                                                                   const char* test,           // fluid pressure
                                                                   const char* density_nw,
                                                                   const char* density_w,
                                                                   const char* gravity_proj,
                                                                   const char* gravity_term,
                                                                   bool  with_gravity,
                                                                   bool  with_capillary_spreading,
                                                                   bool  phase_velocities,
                                                                   bool  node_averaging,
                                                                   const char* thickness,
                                                                   const char* prop_multiplier,
                                                                   const char* velocity,
                                                                   const char* pore_velocity,
                                                                   const char* volume_flux,
                                                                   const char* velocity_nw,
                                                                   const char* velocity_w,
                                                                   const char* nodal_velocity,
                                                                   const char* nodal_pore_velocity,
                                                                   const char* nodal_volume_flux,
                                                                   const char* nodal_velocity_nw,
                                                                   const char* nodal_velocity_w )

    : MathOperatorLHS<dim>(sg.Database(),oper,basic,test),
      satFunc_ ( satfunc ),
      PF_(3),
      VELOFLUX_(dim+1),
      IVELOFLUX_(dim),
      IPVF_(3*(dim*2+1)),
      NVF_(3*(dim*2+1)),
      veloflux_(dim*2+1),
      temp_veloflux_(sg.Mesh().Nodes()),
      RESULT_(dim*2+1,3),
      IPOL_(3),
      mult_vec_(3),
      rho_w_vec_(3),
      rho_nw_vec_(3),
      components_(dim*2+1),
      node_output_(sg.Mesh().Nodes(),false),
      ac_gravity_(9.80665),  // m s-2
      VERTICAL_AXIS_( (dim==1u) ? 0u : 1u ),
      nodal_averaging_(node_averaging),
      phase_velocities_(phase_velocities),
      with_gravity_(with_gravity),
      with_capillary_(with_capillary_spreading),
      with_multiplier_( ((prop_multiplier==NULL) ? false : true) ),
      multiply_with_cell_thickness_( ((thickness==NULL) ? false : true) ),
      verbose_(false)
{
    MathOperatorLHS<dim>::Name("VelocityAndVolumeFlux", oper, basic, test );

    // getting the necessary csmp::Index keys

    velo_key_			 = sg.Database().StorageKey( velocity );
    ivelo_key_			 = sg.Database().StorageKey( pore_velocity );
    flux_key_			 = sg.Database().StorageKey( volume_flux );

    rho_nw_key_          = sg.Database().StorageKey( density_nw );
    rho_w_key_           = sg.Database().StorageKey( density_w );

    if(with_gravity_){
        grav_key_    = sg.Database().StorageKey( gravity_proj );
        gravity_term_key_    = sg.Database().StorageKey( gravity_term );
    }
    if(with_multiplier_)
        mult_key_			 = sg.Database().StorageKey( prop_multiplier );

    thi_key_       = ( multiply_with_cell_thickness_ ? sg.Database().StorageKey(thickness) : Index());
    
    // Facet Data
    facet_normal_idx_   = sg.Database().StorageKey("facet normal");
    facet_area_idx_     = sg.Database().StorageKey("facet area"  );
    sector_volume_idx_  = sg.Database().StorageKey("facet volume");

    // getting ranges for output variables
    sg.Database().RangeOf(velocity, minmaxV_.first, minmaxV_.second );
    sg.Database().RangeOf(volume_flux, minmaxF_.first, minmaxF_.second );
    
    // testing the Operands
    if ( MathOperatorLHS<dim>::MaterialOperandType() != SCALAR )
        throw csmp::Exception( ERROR, "TwoPhaseVelocityAndVolumeFlux::(constructor)",
                               oper, "Operand must be a scalar property." );

    if ( MathOperatorLHS<dim>::BasicOperandType() != SCALAR )
        throw csmp::Exception( ERROR, "TwoPhaseVelocityAndVolumeFlux::(constructor)",
                               basic, "Basic Operand must be a scalar property." );

    if ( MathOperatorLHS<dim>::TestOperandPlacement() != NODE || MathOperatorLHS<dim>::TestOperandType() != SCALAR )
        throw csmp::Exception( ERROR, "TwoPhaseVelocityAndVolumeFlux::(constructor)",
                               test, "Operand 'fluid pressure' must be a scalar property placed on the nodes." );

    if ( velo_key_.place != ELEMENT || velo_key_.type != VECTOR )
        throw csmp::Exception( ERROR, "TwoPhaseVelocityAndVolumeFlux::(constructor)",
                               "Operand 'velocity' must be a vector property placed on the element." );

    if ( ivelo_key_.place != ELEMENT || ivelo_key_.type != VECTOR )
        throw csmp::Exception( ERROR, "TwoPhaseVelocityAndVolumeFlux::(constructor)",
                               "Operand 'pore velocity' must be a vector property placed on the element." );

    if ( flux_key_.place != ELEMENT || flux_key_.type != SCALAR )
        throw csmp::Exception( ERROR, "TwoPhaseVelocityAndVolumeFlux::(constructor)",
                               "Operand 'volume flux' must be a scalar property placed on the element." );

    if ( (rho_nw_key_.type != SCALAR) || (rho_w_key_.type != SCALAR) )
        throw csmp::Exception( ERROR, "TwoPhaseVelocityAndVolumeFlux::(constructor)",
                               "Fluid density operand must be a scalar property." );

    if( with_gravity_){

        if ( grav_key_.type != VECTOR || (grav_key_.place != ELEMENT && grav_key_.place != FACE) )
            throw csmp::Exception( ERROR, "TwoPhaseVelocityAndVolumeFlux::(constructor)",
                                   "Gravity Projection operand must be a vector property and placed on the element or face." );

        if ( gravity_term_key_.type != VECTOR || (gravity_term_key_.place != ELEMENT && gravity_term_key_.place != FACE) )
            throw csmp::Exception( ERROR, "TwoPhaseVelocityAndVolumeFlux::(constructor)",
                                   "Gravity Term operand must be a vector property and placed on the element or face." );

    }

    if( with_multiplier_){

        if (  mult_key_.place != NODE || mult_key_.type != SCALAR )
            throw csmp::Exception( ERROR, "TwoPhaseVelocityAndVolumeFlux::(constructor)",
                                   "Property multiplier must be a scalar property placed on the node." );

    }

    if( phase_velocities ){

        velo_nw_key_ = sg.Database().StorageKey( velocity_nw );
        velo_w_key_  = sg.Database().StorageKey( velocity_w );

        if ( velo_nw_key_.place != ELEMENT || velo_nw_key_.type != VECTOR )
            throw csmp::Exception( ERROR, "TwoPhaseVelocityAndVolumeFlux::(constructor)",
                                   "Operand 'velocity non-wetting phase' must be a vector property placed on the element." );
        if ( velo_w_key_.place != ELEMENT || velo_w_key_.type != VECTOR )
            throw csmp::Exception( ERROR, "TwoPhaseVelocityAndVolumeFlux::(constructor)",
                                   "Operand 'velocity wetting phase' must be a vector property placed on the element." );
    }

    if ( node_averaging )
    {
        nvelo_key_  = sg.Database().StorageKey(nodal_velocity);
        nivelo_key_ = sg.Database().StorageKey(nodal_pore_velocity);
        nflux_key_  = sg.Database().StorageKey(nodal_volume_flux);

        if ( nvelo_key_.place != NODE || nvelo_key_.type != VECTOR )
            throw csmp::Exception( ERROR, "TwoPhaseVelocityAndVolumeFlux::(constructor)",
                                   "Operand 'nodal velocity' must be a vector property placed on the node." );

        if ( nivelo_key_.place != NODE || ivelo_key_.type != VECTOR )
            throw csmp::Exception( ERROR, "TwoPhaseVelocityAndVolumeFlux::(constructor)",
                                   "Operand 'nodal pore velocity' must be a vector property placed on the node." );

        if ( nflux_key_.place != NODE || nflux_key_.type != SCALAR )
            throw csmp::Exception( ERROR, "TwoPhaseVelocityAndVolumeFlux::(constructor)",
                                   "Operand 'nodal volume flux' must be a scalar property placed on the node." );


        if( phase_velocities ){

            nvelo_nw_key_ = sg.Database().StorageKey( nodal_velocity_nw );
            nvelo_w_key_  = sg.Database().StorageKey( nodal_velocity_w );

            if ( nvelo_nw_key_.place != NODE || nvelo_nw_key_.type != VECTOR )
                throw csmp::Exception( ERROR, "TwoPhaseVelocityAndVolumeFlux::(constructor)",
                                       "Operand 'nodal velocity non-wetting phase' must be a vector property placed on the node." );
            if ( nvelo_w_key_.place != NODE || nvelo_w_key_.type != VECTOR )
                throw csmp::Exception( ERROR, "TwoPhaseVelocityAndVolumeFlux::(constructor)",
                                       "Operand 'nodal velocity wetting phase' must be a vector property placed on the node." );
        }

        // if values are to be averagded on the nodes, the operator needs to be applied twice
        MathOperatorLHS<dim>::ApplicationCycles(2);
    }
    
} // end constructor



/**

Using the physically meaningful range which was specified in the Property-
Database, the output variables are tested for the correctness. Note that this
method assumes that the interstital velocity "pore velocity" should be
in the same range as the Darcy velocity.

@section messages Messages

If the range constraint is violated, a message is printed to 'cout'.
*/
template<uint32_t dim,class CELL>
void TwoPhaseVelocityAndVolumeFlux<dim,CELL>::TestRangeOfOutputVariables() const
{
    // velocity, interstitial velocity
    for ( auto i{0U}; i<dim; i++ )
    {
        if ( vt_[i] < minmaxV_.first || vt_[i] > minmaxV_.second ){
            cerr<<" 'velocity' value= "<< vt_[i]<<endl;
            throw csmp::Exception( ERROR, "TwoPhaseVelocityAndVolumeFlux::TestRangeOfOutputVariables",
                                   "Result variable 'velocity' outside of range specified in database file." );
        }

        if ( ivelo_[i] < minmaxV_.first || ivelo_[i] > minmaxV_.second ){
            cerr<<" 'pore velocity' value= "<< ivelo_[i]<<endl;
            throw csmp::Exception( ERROR, "TwoPhaseVelocityAndVolumeFlux::TestRangeOfOutputVariables",
                                   "Result variable 'pore velocity' outside of range specified in database file." );
        }
    }
    // volume flux
    if ( flux_() < minmaxF_.first || flux_() > minmaxF_.second ){
        cerr<<" 'volume flux' value= "<< flux_() <<endl;
        throw csmp::Exception( ERROR, "TwoPhaseVelocityAndVolumeFlux::TestRangeOfOutputVariables",
                               "Result variable 'volume flux' outside of range specified in database file." );
    }
}

/** Switch to verbose mode (results are reported to stdout).
*/
template<uint32_t dim,class CELL>
void TwoPhaseVelocityAndVolumeFlux<dim,CELL>::Verbose( bool stdoutput ) { verbose_=stdoutput; }

template<uint32_t dim,class CELL>
void TwoPhaseVelocityAndVolumeFlux<dim,CELL>::ExtractVelocity( const DenseMatrix<DM_MIN>&      INP,
                                                                  uint32_t        col,
                                                                  VectorVariable<dim>& vc )
{
    for ( auto i{0U}; i<dim; i++ ) vc(i) = INP(i,col);
}



template<uint32_t dim,class CELL>
void TwoPhaseVelocityAndVolumeFlux<dim,CELL>::ExtractVolumeFlux( const DenseMatrix<DM_MIN>&   INP,
                                                                    uint32_t    col,
                                                                    ScalarVariable& sc )
{
    sc() = INP(dim,col);
}



template<uint32_t dim,class CELL>
void TwoPhaseVelocityAndVolumeFlux<dim,CELL>::ExtractInterstitialVelocity( const DenseMatrix<DM_MIN>& INP,
                                                                              uint32_t         col,
                                                                              VectorVariable<dim>& vc )
        {
    for ( auto i{0U}; i<dim; i++ ) vc(i) = INP(i+dim+1,col);
}



/**

Outputs the computed velocities and fluxes to the Model storage
checking the range of each variable using the ranges associated with each
of the Operands.


@section arguments Input Arguments

A reference to the property memory manager and the Element for which the
variables are output.
*/
template<uint32_t dim,class CELL>
void TwoPhaseVelocityAndVolumeFlux<dim,CELL>::WriteOperands( CELL& e )
{
    if ( MathOperatorLHS<dim>::ApplicationCycle() == 1 )
    {
        // 1. outputting element properties first
        // --------------------------------------
        TestRangeOfOutputVariables();
        e.Store(  velo_key_, vt_ );
        e.Store( ivelo_key_, ivelo_ );
        e.Store(  flux_key_, flux_ );

        if (phase_velocities_)
        {
            e.Store(  velo_nw_key_, velo_nw_ );
            e.Store(  velo_w_key_, velo_w_ );
        }
    }
    if ( MathOperatorLHS<dim>::ApplicationCycle() == 2 )
    {
        if ( nodal_averaging_ )
            for ( auto i{0U}; i<e.Nodes(); i++ )
                // doing this operation only once per node
                if ( !node_output_[ e.N(i)->Idx() ] )
                {
                    // 1. assigning the strain & stress values
                    // ---------------------------------------
                    ExtractVelocity( RESULT_, i, vt_ );
                    e.N(i)->Store( nvelo_key_, vt_ );
                    ExtractVolumeFlux( RESULT_, i, flux_ );
                    e.N(i)->Store( nflux_key_, flux_ );
                    ExtractInterstitialVelocity( RESULT_, i, ivelo_ );
                    e.N(i)->Store( nivelo_key_, ivelo_ );

                    // 3. flagging the node to prevent further computations
                    // ----------------------------------------------------
                    node_output_[ e.N(i)->Idx() ] = true;
                }
    }

} // end WriteOperands


template<uint32_t dim,class CELL>
void TwoPhaseVelocityAndVolumeFlux<dim,CELL>::ComputeTotalMobilityRelativeDensityAndGravityTerm( CELL& e )
{

    // Saturation Functions
    satFunc_.Initialize( e );
    satFunc_.InitializeForBaryCenter( e );
    satFunc_.EffectiveSaturation();

    // storing the total mobility value
    e.Store(  MathOperatorLHS<dim>::MaterialOperandKey(), makeScalar( PLAIN, satFunc_.TotalMobility() ) );

    if (with_gravity_){

        const size_t non_wet_phase(2U),wet_phase(1U);

        if ( rho_nw_key_.place == ELEMENT && rho_w_key_.place == ELEMENT ){

            e.Read( rho_nw_key_, rho_nw_ );
            e.Read( rho_w_key_, rho_w_ );

            rhor_ = satFunc_.Permeability()* satFunc_.G() * (rho_w_ - rho_nw_);
            rhot_ = satFunc_.Permeability()* (satFunc_.MobilityPhase( wet_phase ) * rho_w_ + satFunc_.MobilityPhase( non_wet_phase ) * rho_nw_);

        }
        else if ( rho_nw_key_.place == NODE && rho_w_key_.place == NODE ) {

            e.PropertyValueAtBaryCenter( rho_nw_key_, rho_nw_ );
            e.PropertyValueAtBaryCenter( rho_w_key_, rho_w_ );
            e.NodePropertyVector( rho_nw_key_, rho_nw_vec_ );
            e.NodePropertyVector( rho_w_key_, rho_w_vec_ );

            rhor_ = satFunc_.Permeability() *  satFunc_.G() *(rho_w_ - rho_nw_);
            rhot_ = satFunc_.Permeability() * (satFunc_.MobilityPhase( wet_phase ) * rho_w_ + satFunc_.MobilityPhase( non_wet_phase ) * rho_nw_);

        }else
            throw csmp::Exception( ERROR, "VelocityAndVolumeFlux<dim>::ComputeRelativeDensityAndGravityTerm",
                                   "fluid density is neither a node nor element variable; can't deal with this.");

        // ===============================================================

        // Ascertaining that the velocity vector remains in the plane of the surface element

        e.Read(grav_key_, gproj_);

        // storing the gravity term
        e.Store( gravity_term_key_, gproj_* ( ac_gravity_* rhot_() ) );
    }

} // end ComputeGravityTerm



template<uint32_t dim,class CELL>
void TwoPhaseVelocityAndVolumeFlux<dim,CELL>::GetOperands( const CELL& e )
{
    if ( MathOperatorLHS<dim>::ApplicationCycle() == 1 ) {

        // relative density, total mobility and gravity term
        ComputeTotalMobilityRelativeDensityAndGravityTerm( const_cast<CELL&>(e) );

        // fluid pressure
        e.NodePropertyVector( MathOperatorLHS<dim>::TestOperandKey(), PF_ );

        // porosity
        e.Read( MathOperatorLHS<dim>::BasicOperandKey(), phi_ );

        // nodal multipliers
        if ( with_multiplier_ ) e.NodePropertyVector( mult_key_, mult_vec_ );
        cell_thickness_ = ( multiply_with_cell_thickness_ ? e.Read(thi_key_): 1.0);

        // conductivity
        if ( MathOperatorLHS<dim>::MaterialOperandPlacement() == ELEMENT ){

            if ( MathOperatorLHS<dim>::MaterialOperandType() == SCALAR ) {

                MathOperatorLHS<dim>::MTRL[0].AssignToDiagonal( dim,  e.Read( MathOperatorLHS<dim>::MaterialOperandKey() ) );
            }
            else if ( MathOperatorLHS<dim>::MaterialOperandType() == VECTOR ) {
                VectorVariable<dim>  vc;
                e.Read( MathOperatorLHS<dim>::MaterialOperandKey(), vc );
                MathOperatorLHS<dim>::MTRL[0].AssignToDiagonal( vc );
            }
            else if ( MathOperatorLHS<dim>::MaterialOperandType() == TENSOR ) {
                TensorVariable<dim>  ts;
                e.Read( MathOperatorLHS<dim>::MaterialOperandKey(), ts );
                MathOperatorLHS<dim>::MTRL[0] = ts;
            }

        }else if ( MathOperatorLHS<dim>::MaterialOperandPlacement() == ELEMENT_INTEGRATION_POINT ){

            for ( auto i{0U}; i<e.FE()->IntegrationPoints(); i++ ) {
                if ( MathOperatorLHS<dim>::MaterialOperandType() == SCALAR ) {
                    MathOperatorLHS<dim>::MTRL[i].AssignToDiagonal( dim,
                                                                    e.Read( i, MathOperatorLHS<dim>::MaterialOperandKey() ) );
                }
                else if ( MathOperatorLHS<dim>::MaterialOperandType() == VECTOR ) {
                    VectorVariable<dim>  vc;
                    e.Read( i, MathOperatorLHS<dim>::MaterialOperandKey(), vc );
                    MathOperatorLHS<dim>::MTRL[i].AssignToDiagonal( vc );
                }
                else if ( MathOperatorLHS<dim>::MaterialOperandType() == TENSOR ) {
                    TensorVariable<dim>  ts;
                    e.Read( i, MathOperatorLHS<dim>::MaterialOperandKey(), ts );
                    MathOperatorLHS<dim>::MTRL[i] = ts;
                }
            }
        }
        else // if a nodal variable is dealt with
        {
            if ( e.IntegrationPoints() == 0U )
                throw csmp::Exception( FATAL_ERROR, "VelocityAndVolumeFlux<dim>::GetOperands",
                                       "The current finite element has no integration points",
                                       "Therefore nodal properties cannot be integrated.");

            for ( auto i{0U}; i<e.IntegrationPoints(); i++ )
            {
                MathOperatorLHS<dim>::MTRL[i].Resize(dim,dim);
                MathOperatorLHS<dim>::MTRL[i].Zero();
                MathOperatorLHS<dim>::PropertyAtIntegrationPoint( e, MathOperatorLHS<dim>::MaterialOperandKey(), i, MathOperatorLHS<dim>::MTRL[i] );
            }
        }

    }
    
} // end GetOperands




/**

@section arguments Input Arguments 

A reference to the Element for which the post-processing is done.  
*/
template<uint32_t dim,class CELL>
void TwoPhaseVelocityAndVolumeFlux<dim,CELL>::ComputeContribution( const CELL& e )
{
    typename list<vector<double> >::const_iterator  lit;

    if ( MathOperatorLHS<dim>::ApplicationCycle() == 1 ) {

        if ( verbose_ )
            cout <<"\n\nVelocityAndVolumeFlux::ComputeContribution: Element: "<< e.Idx() << endl;

        // 1. Case of linear elements: Constant velocity
        // ---------------------------------------------
        if ( (e.Interpolation() == 1) && this->MaterialOperandPlacement()==ELEMENT )
        {

            const size_t non_wet_phase(2U),wet_phase(1U);

            // Saturation Functions
            satFunc_.Initialize( e );
            satFunc_.InitializeForBaryCenter( e );
            satFunc_.EffectiveSaturation();

            e.dN_AtBaryCenter(DERIV_);

            vt_ = 0.;

            // interpolation function derivatives at the nodes

            for ( auto i{0U}; i<e.Nodes(); i++ )
                for ( auto j{0U}; j<dim; j++ )
                    vt_(j) += PF_[i]() * -DERIV_(j,i) * MathOperatorLHS<dim>::MTRL[0](j,j);

            if ( with_gravity_ ){

                for( uint32_t xyz = 0; xyz < dim; ++xyz )
                    vt_( xyz ) += rhot_() * ac_gravity_ * gproj_[xyz];
            }

            if( with_capillary_ ){

                for ( auto i{0U}; i<e.Nodes(); i++ ){
                    satFunc_.InitializeForNode(e,i);
                    satFunc_.EffectiveSaturation();
                    const double pc ( satFunc_.pc_Phase( ));
                    for ( auto j{0U}; j<dim; j++ )
                        vt_( j ) += pc  * -DERIV_(j,i) * satFunc_.Permeability() * satFunc_.MobilityPhase( non_wet_phase );
                }

            }

            // interstitial velocity ( pore velocity )
            ivelo_  = vt_;
            ivelo_ /= phi_;

            // flux
            flux_() = vt_.Length();

            // phase velocities
            if ( phase_velocities_ ){

                satFunc_.Initialize( e );
                satFunc_.InitializeForBaryCenter( e );
                satFunc_.EffectiveSaturation();

                velo_nw_ = 0.;
                velo_w_  = 0.;

                for ( auto i{0U}; i<e.Nodes(); i++ )
                    for ( auto j{0U}; j<dim; j++ ){
                        velo_nw_(j) += PF_[i]() * -DERIV_(j,i) * satFunc_.Permeability() * satFunc_.MobilityPhase( non_wet_phase );
                        velo_w_(j)  += PF_[i]() * -DERIV_(j,i) * satFunc_.Permeability() * satFunc_.MobilityPhase( wet_phase );
                    }

                if ( with_gravity_ ){

                    for( auto xyz = 0; xyz < dim; ++xyz ){
                        velo_nw_( xyz ) += satFunc_.Permeability() * satFunc_.MobilityPhase( non_wet_phase )* rho_nw_() * ac_gravity_ * gproj_[xyz];
                        velo_w_( xyz )  += satFunc_.Permeability() * satFunc_.MobilityPhase( wet_phase )    * rho_w_()  * ac_gravity_ * gproj_[xyz];
                    }
                }

                if( with_capillary_ ){

                    for ( auto i{0U}; i<e.Nodes(); i++ ){
                        satFunc_.InitializeForNode(e,i);
                        satFunc_.EffectiveSaturation();
                        const double pc ( satFunc_.pc_Phase());
                        for ( auto j{0U}; j<dim; j++ )
                            velo_nw_( j ) += pc  * -DERIV_(j,i) * satFunc_.Permeability()* satFunc_.MobilityPhase( non_wet_phase);
                    }

                }

            }

            if ( verbose_ ) {
                cout <<"\ncomputed element variable 'velocity':"<< endl;
                vt_.Out();
                cout <<"\ncomputed element variable 'pore velocity':"<< endl;
                ivelo_.Out();
                cout <<"\ncomputed element variable 'volume flux': "<< flux_() << endl;
            }
        }

        // 2. Case numerically integrated elements:
        //    Computing the velocities and the volume flux at the integration points
        // -------------------------------------------------------------------------
        else if ( (this->MaterialOperandPlacement()==ELEMENT_INTEGRATION_POINT) || (this->MaterialOperandPlacement() == NODE) )
        {
            const size_t non_wet_phase(2U),wet_phase(1U);

            // Saturation Functions
            satFunc_.Initialize( e );
            satFunc_.InitializeForBaryCenter( e );
            satFunc_.EffectiveSaturation();

            DERIV_.Resize(dim,e.Nodes());
            IPVF_.resize( e.FE()->IntegrationPoints()*(dim*2+1) );
            NVF_.resize( e.Nodes()*(dim*2+1) );

            // collecting average data for element variables
            vt_    = 0.0;
            ivelo_   = 0.0;
            flux_()  = 0.0;
            velo_nw_ = 0.0;
            velo_w_  = 0.0;

            for ( auto i{0U}; i<e.FE()->IntegrationPoints(); i++ )
            {
                // if density driven flow is computed, calculate rho * g * z
                // compute density at integration points and any multiplier as well

                if ( with_gravity_ ) {

                    // density element variable
                    if ( rho_nw_key_.place == ELEMENT && rho_w_key_.place == ELEMENT ) {

                        rho_fac_      = rhor_();
                        rhot_fac_     = rhot_();
                        rho_nw_fac_   = rho_w_();
                        rho_w_fac_    = rho_nw_();

                        if ( !with_multiplier_ ){

                            mult_fac_ = 1.; // more likely case first

                        }else{

                            e.N_AtIntegrationPoint( i, IPOL_ );
                            mult_fac_=0.;
                            for ( auto j{0U}; j<e.Nodes(); j++ )
                                mult_fac_ += IPOL_[j] * mult_vec_[j]();

                        }
                    }
                    // density nodal variable
                    else if ( rho_nw_key_.place == NODE && rho_w_key_.place == NODE ) {

                        e.N_AtIntegrationPoint( i, IPOL_ );

                        rho_fac_      = 0.0;
                        rhot_fac_     = 0.0 ;
                        mult_fac_     = 0.0;
                        rho_nw_fac_   = 0.0;
                        rho_w_fac_    = 0.0;

                        for ( auto j{0U}; j<e.Nodes(); j++ ) {

                            if ( with_multiplier_ ) mult_fac_ += IPOL_[j] * mult_vec_[j]();
                            else                    mult_fac_ = 1.0;

                            rho_fac_    += IPOL_[j] * satFunc_.Permeability()* satFunc_.G() *( rho_w_vec_[j]() - rho_nw_vec_[j]() );
                            rhot_fac_   += IPOL_[j] * satFunc_.Permeability()* (satFunc_.MobilityPhase( wet_phase ) * rho_w_vec_[j]() + satFunc_.MobilityPhase( non_wet_phase ) * rho_nw_vec_[j]());
                            rho_nw_fac_ += IPOL_[j] * rho_nw_vec_[j]();
                            rho_w_fac_  += IPOL_[j] * rho_w_vec_[j]();
                        }

                    }
                    else {
                        cout <<"\nVelocityAndVolumeFlux<"<<  dim <<">::ComputeContribution: ";
                        cout <<"placement of density variable not recognized."<< endl;
                    }

                }

                // getting DN matrices at the node points
                e.dN_AtIntegrationPoint( DERIV_, i, 1 );

                // compute velocity and flux VELOFLUX
                fill( VELOFLUX_.begin(), VELOFLUX_.end(), 0.0 );

                if( with_capillary_ ){

                    for ( auto n=0; n<e.Nodes(); n++ ){
                        satFunc_.InitializeForNode(e,n);
                        satFunc_.EffectiveSaturation();
                        const double pc ( satFunc_.pc_Phase( ));
                        for ( auto j{0U}; j<dim; j++ )
                            VELOFLUX_[j] += ( PF_[n]() + pc ) * -DERIV_(j,n) * MathOperatorLHS<dim>::MTRL[i](j,j);
                    }

                }else{

                    for ( auto n=0; n<e.Nodes(); n++ )
                        for ( auto j{0U}; j<dim; j++ )
                            VELOFLUX_[j] += PF_[n]() * -DERIV_(j,n) * MathOperatorLHS<dim>::MTRL[i](j,j);

                }

                if ( with_gravity_ ){

                    for( auto xyz = 0; xyz < dim; ++xyz )
                        VELOFLUX_[ xyz ] += mult_fac_ * rhot_fac_ * ac_gravity_ * gproj_[ xyz ];
                }

                // interstitial velocity & volume flux
                for ( auto n=0; n<dim; n++ ) {
                    IVELOFLUX_[n]   = VELOFLUX_[n]/phi_();
                    VELOFLUX_[dim] += VELOFLUX_[n]*VELOFLUX_[n];
                    // summing integration point values for later averaging
                    vt_(n)       += VELOFLUX_[n];
                    ivelo_(n)      += IVELOFLUX_[n];
                }

                VELOFLUX_[dim] = sqrt(VELOFLUX_[dim]);
                flux_() += VELOFLUX_[dim];

                // phase velocities
                if ( phase_velocities_ ){

                    satFunc_.Initialize( e );
                    satFunc_.InitializeForBaryCenter( e );
                    satFunc_.EffectiveSaturation();

                    fill( VELOFLUX_W_.begin(), VELOFLUX_W_.end(), 0.0 );
                    fill( VELOFLUX_NW_.begin(), VELOFLUX_NW_.end(), 0.0 );


                    for ( auto n=0; n<e.Nodes(); n++ )
                        for ( auto j{0U}; j<dim; j++ ){
                            VELOFLUX_NW_[j] += PF_[n]() * -DERIV_(j,n) * satFunc_.Permeability() * satFunc_.MobilityPhase( non_wet_phase );
                            VELOFLUX_W_[j]  += PF_[n]() * -DERIV_(j,n) * satFunc_.Permeability() * satFunc_.MobilityPhase( wet_phase );
                        }

                    if ( with_gravity_ ){

                        for( auto xyz = 0; xyz < dim; ++xyz ){
                            VELOFLUX_NW_[ xyz ] += satFunc_.Permeability() * satFunc_.MobilityPhase( non_wet_phase )* mult_fac_ * rho_nw_fac_ * ac_gravity_ * gproj_[xyz];
                            VELOFLUX_W_ [ xyz ] += satFunc_.Permeability() * satFunc_.MobilityPhase( wet_phase )    * mult_fac_ * rho_w_fac_  * ac_gravity_ * gproj_[xyz];
                        }
                    }

                    if( with_capillary_ ){

                        for ( auto n=0; n<e.Nodes(); n++ ){
                            satFunc_.InitializeForNode(e,n);
                            satFunc_.EffectiveSaturation();
                            const double pc ( satFunc_.pc_Phase( ));
                            for ( auto j{0U}; j<dim; j++ )
                                VELOFLUX_NW_[ j ] += pc  * -DERIV_(j,n) * satFunc_.Permeability()* satFunc_.MobilityPhase( non_wet_phase);
                        }

                    }

                    // interstitial velocity & volume flux
                    for ( auto n=0; n<dim; n++ ) {
                        // summing integration point values for later averaging
                        velo_nw_(n) += VELOFLUX_NW_[n];
                        velo_w_(n)  += VELOFLUX_W_[n];
                    }

                }

                if ( verbose_ ) {
                    cout <<"\nVelocity x, y, (z), magnitude at integration point: "<< i << endl;
                    printVector( VELOFLUX_ );
                    cout <<"\nInterstitial velocity at same point:                "<< endl;
                    printVector( IVELOFLUX_ );
                }

                if ( nodal_averaging_ )
                {
                    // inserting velocity, volume flux, and interstitial velocity
                    // into single STL vector:
                    // velocity & volume flux
                    for ( auto k=0; k<dim; k++ ) {
                        IPVF_[ i*components_ + k ] = VELOFLUX_[k];
                        IPVF_[ i*components_ + dim + 1 + k ] = IVELOFLUX_[k];
                    }
                    IPVF_[ i*components_ + dim ] = VELOFLUX_[dim];
                }
            }
            // Averaging integration point values to get the element variables
            // ---------------------------------------------------------------
            vt_     /= static_cast<double>(e.FE()->IntegrationPoints());
            flux_() /= static_cast<double>(e.FE()->IntegrationPoints());
            ivelo_  /= static_cast<double>(e.FE()->IntegrationPoints());

            if( phase_velocities_ ){
                velo_nw_  /= static_cast<double>(e.FE()->IntegrationPoints());
                velo_w_   /= static_cast<double>(e.FE()->IntegrationPoints());
            }

        } else {
            throw csmp::Exception( ERROR, "ThoPhaseVelocityAndVolumeFlux::(constructor)",
                                   this->MaterialOperandName(),
                                   " must be either an element, node or an integration-point property." );
        }
        
        if ( nodal_averaging_ ){

            // if the computed properties are constant on the element
            if ( e.Interpolation() == 1 )
            {
                for ( auto k=0; k<dim; k++ )
                {
                    veloflux_[k]           = vt_[k];
                    veloflux_[dim + 1 + k] = ivelo_[k];
                }
                veloflux_[dim] = flux_();

                for ( auto i{0U}; i<e.Nodes(); i++ )
                    temp_veloflux_[ e.N(i)->Idx() ].push_back( veloflux_ );
            }
            else
            {
                // 1.2 Extrapolating velocity and volume flux from the
                //     integration points to the nodes.
                // ---------------------------------------------------
                e.ExtrapolateIntegrationPointVariableToNodes( components_, IPVF_, NVF_ );

                for ( auto i{0U}; i<e.Nodes(); i++ )
                {
                    for ( auto k=0; k<components_; k++ ) veloflux_[k] = NVF_[ i*components_ + k ];
                    temp_veloflux_[ e.N(i)->Idx() ].push_back( veloflux_ );
                }
            }
        }

    }  // end of first application cycle








    // 2. During the second visitation, the nodal velocities and fluxes computed for each element
    //    node are averaged and stored in a vector for output.
    //    -------------------------------------------
    if ( (MathOperatorLHS<dim>::ApplicationCycle() == 2) && nodal_averaging_ ) {

        RESULT_.Resize(components_,e.Nodes());

        for ( auto i{0U}; i<e.Nodes(); i++ ) {
            // duplicate calculations are avoided via the boolean vector
            if ( !node_output_[ e.N(i)->Idx() ] )
            {
                for ( auto j{0U}; j<components_; j++ )
                {
                    // averaging velocity/flux components
                    for ( sum_=0.0,lit =temp_veloflux_[ e.N(i)->Idx() ].begin(); lit!=temp_veloflux_[ e.N(i)->Idx() ].end(); lit++ )
                        sum_ += (*lit)[j];

                    sum_ /= static_cast<double>(temp_veloflux_[ e.N(i)->Idx() ].size());

                    RESULT_(j,i) = sum_;
                }
            }

            // now the vector list is no longer needed and therefore erased for the
            // next application of the post-processing operator
            temp_veloflux_[ e.N(i)->Idx() ].erase( temp_veloflux_[ e.N(i)->Idx() ].begin(), temp_veloflux_[ e.N(i)->Idx() ].end() );
        }

    } // end application cycle 2





} // end ComputeContribution

//////////////////////////////////////////////////
// For node-centered calculations
//////////////////////////////////////////////////

template<uint32_t dim,class CELL>
void TwoPhaseVelocityAndVolumeFlux<dim,CELL>::ComputeContribution( const Node<dim>& n_ref )
{
    // initialize output value
    double value = 0.0;

    uint32_t inside_node,outside_node;
    double vtn(0.0),facetArea(1.0);
    size_t advected_phase_n(2U), advected_phase_w(1U);
    const double zero(0.0); //can be num_epsilon or so
    //const double zero( std::numeric_limits<double>::min()); //can be num_epsilon or so
    //const double zero(1.0e-15); //can be num_epsilon or so

    Element<dim>* eptr(nullptr);
    double sign(0.0);

    this->ApplicationCycle(1);

    if( (!with_gravity_) && (!with_capillary_) ){

        double fn(1.0);

        for ( auto t=0U; t<n_ref.Parents(); t++ )
        {
            eptr = n_ref.Parent(t);
            this->GetOperands(*eptr);
            auto pnid(n_ref.ParentNodeNumber(t));
            cell_thickness_ = ( multiply_with_cell_thickness_ ? eptr->Read(this->thi_key_): 1.0);
            eptr->Read(this->velo_key_,vt_);

            satFunc_.Initialize(*eptr);
            // --------------------------------------------------------------------------
            // For all finite-volume facets
            for ( auto k{0U}; k<eptr->FV()->FacetsPerSector(pnid); k++ )
            {
                auto i( eptr->FV()->FacetSurroundingSector(pnid,k) );

                // --------------------------------------------------------------------------
                // Finite Volume Stencil information & total velocity projection
                eptr->FV()->FacetEdgeNodes( i, inside_node, outside_node );
                eptr->Read(i,0,facet_area_idx_,sc_);
                eptr->Read(i,0,facet_normal_idx_,facet_n_);
                facetArea = sc_()*cell_thickness_;
                vtn  = dotProduct( vt_, facet_n_);
                // --------------------------------------------------------------------------
                // Calculate values of saturations for inside & ouside nodes
                const double sn_inside_node  = 1.0 - eptr->N(inside_node)->Read(this->satFunc_.WettingPhaseSaturationKey() );
                const double sn_outside_node = 1.0 - eptr->N(outside_node)->Read(this->satFunc_.WettingPhaseSaturationKey() );

                // --------------------------------------------------------------------------
                // Initilize saturation functions for Facet Integration Point
                satFunc_.InitializeForFacetIntegrationPoint(i,0U,*eptr);

                // --------------------------------------------------------------------------
                // Define Upstream Saturation based on Phase Velocity
                if ( vtn > zero )
                    satFunc_.SaturationWettingPhase( 1.0 - sn_inside_node );
                else
                    satFunc_.SaturationWettingPhase( 1.0 - sn_outside_node );

                satFunc_.EffectiveSaturation();

                // --------------------------------------------------------------------------
                // Finally based on the Upstream Saturation Calculate Phase Velocity of Advected Phase
                fn = satFunc_.f_Phase( advected_phase_n );

                flux_() = vtn*fn*facetArea;

                sign = ( (pnid == inside_node) ? 1.0 : -1.0);    // +: inside node, -: outside node

                value  += sign*flux_();
            }
        }

    }else if( !with_gravity_ ){

        double sn_inside_node (0.0), ln_inside_node (0.0),  lw_inside_node(0.0);
        double sn_outside_node(0.0), ln_outside_node (0.0), lw_outside_node(0.0);

        //double upstream_sn(0.0), upstream_sw(0.0);
        double upstream_mobility_n(0.0),upstream_mobility_w(0.0),total_mobility(0.0);
        double upstream_fn(0.0);
        //double upstream_lambda_overbar(0.0);

        double dsdn(1.0),dpcdn(1.0);
        //double viscous_velocity_component(1.0), capillary_velocity_component(1.0);
        double vn_capillary_component_of_velocity (1.0), vw_capillary_component_of_velocity (1.0);
        double vn_at_facet_int_point (1.0), vw_at_facet_int_point (1.0);


        for ( auto t=0U; t<n_ref.Parents(); t++ )
        {
            eptr = n_ref.Parent(t);
            this->GetOperands(*eptr);
            auto pnid(n_ref.ParentNodeNumber(t));
//            size_t eid(eptr->Idx());

            cell_thickness_ = ( multiply_with_cell_thickness_ ? eptr->Read(this->thi_key_): 1.0);
            eptr->Read(this->velo_key_,vt_);

            satFunc_.Initialize(*eptr);

            // --------------------------------------------------------------------------
            // Compute Capillary Pressure Gradient

            eptr->dN_AtBaryCenter( DN_ );
            dsdn_ = 0.0;
            for ( auto j{0U}; j<eptr->Nodes(); j++ ) {
                 const double sn = eptr->N(j)->Read( this->TestOperandKey() );
                 for ( auto k{0U}; k<dim; k++ ) dsdn_[k] += DN_(k,j) * sn;
            }

            // --------------------------------------------------------------------------
            // For all finite-volume facets
            for ( auto k{0U}; k<eptr->FV()->FacetsPerSector(pnid); k++ )
            {
                auto i( eptr->FV()->FacetSurroundingSector(pnid,k) );

                // --------------------------------------------------------------------------
                // Finite Volume Stencil information & total velocity projection
                eptr->FV()->FacetEdgeNodes( i, inside_node, outside_node );
                eptr->Read(i,0,facet_area_idx_,sc_);
                eptr->Read(i,0,facet_normal_idx_,facet_n_);
                facetArea = sc_()*cell_thickness_;
                vtn  = dotProduct( vt_, facet_n_);

                // --------------------------------------------------------------------------
                // Calculate values of saturations & mobilities for inside node
                satFunc_.InitializeForNode( *eptr, inside_node );
                satFunc_.EffectiveSaturation();
                sn_inside_node = eptr->N(inside_node)->Read( this->TestOperandKey()  );
                ln_inside_node = satFunc_.MobilityPhase( advected_phase_n );
                lw_inside_node = satFunc_.MobilityPhase( advected_phase_w );

                // --------------------------------------------------------------------------
                // Calculate values of saturations & mobilities for outside node
                satFunc_.InitializeForNode( *eptr, outside_node );
                satFunc_.EffectiveSaturation();
                sn_outside_node = eptr->N(outside_node)->Read( this->TestOperandKey() );
                ln_outside_node = satFunc_.MobilityPhase( advected_phase_n );
                lw_outside_node = satFunc_.MobilityPhase( advected_phase_w );

                // --------------------------------------------------------------------------
                // Initilize saturation functions for Facet Integration Point
                satFunc_.InitializeForFacetIntegrationPoint( i, 0U, *eptr );
                satFunc_.EffectiveSaturation();

                // --------------------------------------------------------------------------
                // Calculate Phase Velocities at the facet integration points, define Capillary Component of the Phase Velocities if necessary

                dsdn = facet_n_.DotProduct(dsdn_);
                dpcdn = -satFunc_.dpcds_Phase( )*dsdn;

                vn_capillary_component_of_velocity = satFunc_.MobilityPhase(advected_phase_w) * satFunc_.Permeability() * dpcdn;
                vw_capillary_component_of_velocity = satFunc_.MobilityPhase(advected_phase_n) * satFunc_.Permeability() * dpcdn;

                vn_at_facet_int_point = vtn - vn_capillary_component_of_velocity;
                vw_at_facet_int_point = vtn + vw_capillary_component_of_velocity;


                // --------------------------------------------------------------------------
                // Define Upstream Mobilities based on Phase Velocities

                // mixture moving inside the CV        mixture moving outside the CV
                //
                //       |  vt                                 /|\ vt
                //       |                                      |
                //      \|/              /|\                    |
                //     -----              |  N_up             -----
                //   /       \                              /       \
                //   \       /                              \       /
                //     -----              |                   -----
                //      /|\              \|/  N_down            |
                //       |                                      |
                //       |  vt                                 \|/ vt

                if( (vn_at_facet_int_point >zero ) && ( vw_at_facet_int_point >zero ) ){

                    upstream_mobility_n=ln_inside_node;
                    upstream_mobility_w=lw_inside_node;

                    //upstream_sn = sn_inside_node;
                    //upstream_sw = 1.0 - sn_inside_node;

                }else if ((vn_at_facet_int_point>zero)&&(vw_at_facet_int_point<zero) ){

                    upstream_mobility_n=ln_inside_node;
                    upstream_mobility_w=lw_outside_node;

                    //upstream_sn = sn_inside_node;
                    //upstream_sw = 1.0 - sn_outside_node;

                }else if ( (vn_at_facet_int_point<zero) && (vw_at_facet_int_point>zero) ){

                    upstream_mobility_n=ln_outside_node;
                    upstream_mobility_w=lw_inside_node;

                    //upstream_sn = sn_outside_node;
                    //upstream_sw = 1.0 - sn_inside_node;

                }else if ( (vn_at_facet_int_point<zero) && (vw_at_facet_int_point<zero) ){

                    upstream_mobility_n=ln_outside_node;
                    upstream_mobility_w=lw_outside_node;

                    //upstream_sn = sn_outside_node;
                    //upstream_sw = 1.0 - sn_outside_node;

                }else{

                    upstream_mobility_n = 0.5*(ln_inside_node + ln_outside_node);
                    upstream_mobility_w = 0.5*(lw_inside_node + lw_outside_node);

                    //upstream_sn = 0.5 * ( sn_inside_node + sn_outside_node );
                    //upstream_sw = 1.0 - upstream_sn;

                }

                // --------------------------------------------------------------------------
                // Finally based on the Upstream Mobilities Calculate Phase Velocity of Advected Phase, include Capillary Component if necessary

                total_mobility=upstream_mobility_n+upstream_mobility_w;
                upstream_fn=(total_mobility!=0.0? upstream_mobility_n/total_mobility : 0.0);

                flux_() = upstream_fn * vn_at_facet_int_point * facetArea;

                sign = ( (pnid == inside_node) ? 1.0 : -1.0);    // +: inside node, -: outside node

                value += sign*flux_();
            }
        }

    }else{

        double sn_inside_node (0.0), ln_inside_node (0.0),  lw_inside_node(0.0);
        double sn_outside_node(0.0), ln_outside_node (0.0), lw_outside_node(0.0);

        //double upstream_sn(0.0), upstream_sw(0.0);
        double upstream_mobility_n(0.0),upstream_mobility_w(0.0),total_mobility(0.0);
        double upstream_fn(0.0), upstream_lambda_overbar(0.0);

        double dsdn(1.0),dpcdn(1.0),gproj_n(0.0);
        double viscous_velocity_component(1.0), capillary_velocity_component(1.0), gravity_velocity_component(1.0);
        double vn_capillary_component_of_velocity (1.0), vw_capillary_component_of_velocity (1.0);
        double vn_gravity_component_of_velocity (1.0), vw_gravity_component_of_velocity (1.0);
        double vn_at_facet_int_point (1.0), vw_at_facet_int_point (1.0);

        for ( auto t=0U; t<n_ref.Parents(); t++ )
        {
            eptr = n_ref.Parent(t);
            this->GetOperands(*eptr);
            auto pnid(n_ref.ParentNodeNumber(t));
//            size_t eid(eptr->Idx());

            cell_thickness_ = ( multiply_with_cell_thickness_ ? eptr->Read(this->thi_key_): 1.0);
            eptr->Read(this->velo_key_,vt_);
            eptr->Read( grav_key_, gproj_);

            satFunc_.Initialize(*eptr);

            // --------------------------------------------------------------------------
            // Compute Capillary Pressure Gradient if necessary
            if( with_capillary_ ){

                eptr->dN_AtBaryCenter( DN_ );

                dsdn_ = 0.0;
                for ( auto j{0U}; j<eptr->Nodes(); j++ ) {
                     const double sn = eptr->N(j)->Read( this->TestOperandKey() );
                     for ( auto k{0U}; k<dim; k++ ) dsdn_[k] += DN_(k,j) * sn;
                }

            }

            // --------------------------------------------------------------------------
            // For all finite-volume facets
            for ( auto k{0U}; k<eptr->FV()->FacetsPerSector(pnid); k++ )
            {
                auto i( eptr->FV()->FacetSurroundingSector(pnid,k) );

                // --------------------------------------------------------------------------
                // Finite Volume Stencil information & total velocity projection
                eptr->FV()->FacetEdgeNodes( i, inside_node, outside_node );
                eptr->Read(i,0,facet_area_idx_,sc_);
                eptr->Read(i,0,facet_normal_idx_,facet_n_);
                facetArea = sc_()*cell_thickness_;
                vtn  = dotProduct(vt_, facet_n_);

                // --------------------------------------------------------------------------
                // Calculate values of saturations & mobilities for inside node
                satFunc_.InitializeForNode( *eptr, inside_node );
                satFunc_.EffectiveSaturation();
                sn_inside_node = eptr->N(inside_node)->Read( this->TestOperandKey()  );
                ln_inside_node = satFunc_.MobilityPhase( advected_phase_n );
                lw_inside_node = satFunc_.MobilityPhase( advected_phase_w );

                // --------------------------------------------------------------------------
                // Calculate values of saturations & mobilities for outside node
                satFunc_.InitializeForNode( *eptr, outside_node );
                satFunc_.EffectiveSaturation();
                sn_outside_node = eptr->N(outside_node)->Read( this->TestOperandKey() );
                ln_outside_node = satFunc_.MobilityPhase( advected_phase_n );
                lw_outside_node = satFunc_.MobilityPhase( advected_phase_w );


                // --------------------------------------------------------------------------
                // Initilize saturation functions for Facet Integration Point
                satFunc_.InitializeForFacetIntegrationPoint( i, 0U, *eptr );
                satFunc_.EffectiveSaturation();

                // --------------------------------------------------------------------------
                // Define Gravity Component of the Phase Velocities
                gproj_n  =  -dotProduct(gproj_, facet_n_);

                vn_gravity_component_of_velocity = satFunc_.MobilityPhase(advected_phase_w) * satFunc_.GravityTerm() * gproj_n;
                vw_gravity_component_of_velocity = satFunc_.MobilityPhase(advected_phase_n) * satFunc_.GravityTerm() * gproj_n;

                // --------------------------------------------------------------------------
                // Calculate Phase Velocities at the facet integration points, define Capillary Component of the Phase Velocities if necessary

                if( with_capillary_ ){

                    dsdn = facet_n_.DotProduct( dsdn_);
                    dpcdn = -satFunc_.dpcds_Phase( )*dsdn;

                    vn_capillary_component_of_velocity = satFunc_.MobilityPhase(advected_phase_w) * satFunc_.Permeability() * dpcdn;
                    vw_capillary_component_of_velocity = satFunc_.MobilityPhase(advected_phase_n) * satFunc_.Permeability() * dpcdn;
                    vn_at_facet_int_point = vtn - vn_gravity_component_of_velocity - vn_capillary_component_of_velocity;
                    vw_at_facet_int_point = vtn + vw_gravity_component_of_velocity + vw_capillary_component_of_velocity;

                }else{

                    vn_at_facet_int_point = vtn - vn_gravity_component_of_velocity;
                    vw_at_facet_int_point = vtn + vw_gravity_component_of_velocity;

                }

                // --------------------------------------------------------------------------
                // Define Upstream Mobilities based on Phase Velocities

                // mixture moving inside the CV        mixture moving outside the CV
                //
                //       |  vt                                 /|\ vt
                //       |                                      |
                //      \|/              /|\                    |
                //     -----              |  N_up             -----
                //   /       \                              /       \
                //   \       /                              \       /
                //     -----              |                   -----
                //      /|\              \|/  N_down            |
                //       |                                      |
                //       |  vt                                 \|/ vt

                if( (vn_at_facet_int_point >zero ) && ( vw_at_facet_int_point >zero ) ){

                    upstream_mobility_n=ln_inside_node;
                    upstream_mobility_w=lw_inside_node;

                    //upstream_sn = sn_inside_node;
                    //upstream_sw = 1.0 - sn_inside_node;

                }else if ((vn_at_facet_int_point>zero)&&(vw_at_facet_int_point<zero) ){

                    upstream_mobility_n=ln_inside_node;
                    upstream_mobility_w=lw_outside_node;

                    //upstream_sn = sn_inside_node;
                    //upstream_sw = 1.0 - sn_outside_node;

                }else if ( (vn_at_facet_int_point<zero) && (vw_at_facet_int_point>zero) ){

                    upstream_mobility_n=ln_outside_node;
                    upstream_mobility_w=lw_inside_node;

                    //upstream_sn = sn_outside_node;
                    //upstream_sw = 1.0 - sn_inside_node;

                }else if ( (vn_at_facet_int_point<zero) && (vw_at_facet_int_point<zero) ){

                    upstream_mobility_n=ln_outside_node;
                    upstream_mobility_w=lw_outside_node;

                    //upstream_sn = sn_outside_node;
                    //upstream_sw = 1.0 - sn_outside_node;

                }else{

                    upstream_mobility_n = 0.5*(ln_inside_node + ln_outside_node);
                    upstream_mobility_w = 0.5*(lw_inside_node + lw_outside_node);

                    //upstream_sn = 0.5 * ( sn_inside_node + sn_outside_node );
                    //upstream_sw = 1.0 - upstream_sn;

                }

                total_mobility=upstream_mobility_n+upstream_mobility_w;
                upstream_fn=(total_mobility!=0.0? upstream_mobility_n/total_mobility : 0.0);
                upstream_lambda_overbar=(total_mobility!=0.0? (upstream_mobility_n*upstream_mobility_w)/total_mobility : 0.0);

                // --------------------------------------------------------------------------
                // Finally based on the Upstream Mobilities Calculate Phase Velocity of Advected Phase, include Capillary Component if necessary

                viscous_velocity_component = upstream_fn * vtn ;
                gravity_velocity_component = upstream_lambda_overbar * satFunc_.GravityTerm() * gproj_n;

                if ( with_capillary_ ) {

                    capillary_velocity_component = upstream_fn * vn_capillary_component_of_velocity ;

                    flux_() = (viscous_velocity_component - gravity_velocity_component - capillary_velocity_component ) * facetArea;

                } else flux_() = (viscous_velocity_component - gravity_velocity_component ) * facetArea;

                sign = ( (pnid == inside_node) ? 1.0 : -1.0);    // +: inside node, -: outside node

                value += sign*flux_();
            }
        }
    }

} // end ComputeContribution


template class TwoPhaseVelocityAndVolumeFlux<1U,Element<1U> >;
template class TwoPhaseVelocityAndVolumeFlux<2U,Element<2U> >;
template class TwoPhaseVelocityAndVolumeFlux<3U,Element<3U> >;

//template class TwoPhaseVelocityAndVolumeFlux<1U,Face<1U> >;
//template class TwoPhaseVelocityAndVolumeFlux<2U,Face<2U> >;
//template class TwoPhaseVelocityAndVolumeFlux<3U,Face<3U> >;

} 
















