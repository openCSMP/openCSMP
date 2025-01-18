#include "VelocityAndVolumeFlux.h"
#include "Model.h"
#include "Region.h"
#include "MathOperatorLHS.h"
#include "STL_utilities.h"
#include "Element.h"
#include "Face.h"
#include "ErrorHandler.h"
#include "CSMP_physical_constants.h"

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
template<uint32_t dim, template<uint32_t> class CELL>
VelocityAndVolumeFlux<dim,CELL>::VelocityAndVolumeFlux( const Model<dim>& sg,
                                                        const char* oper,   // conductivity
                                                        const char* basic,  // porosity
                                                        const char* test,   // fluid pressure
                                                        bool  node_averaging,
                                                        const char* velocity,
                                                        const char* pore_velocity,
                                                        const char* volume_flux,
                                                        const char* nodal_velocity,
                                                        const char* nodal_pore_velocity,
                                                        const char* nodal_volume_flux )
  : MathOperatorLHS<dim,CELL>(sg.Database(),oper,basic,test),
    // getting the necessary csmp::Index keys
    velo_key_(sg.Database().StorageKey(velocity)),
    ivelo_key_(sg.Database().StorageKey(pore_velocity)),
    flux_key_(sg.Database().StorageKey(volume_flux)),
    PF_(3),
    VELOFLUX_(dim+1),
    IVELOFLUX_(dim),
    IPVF_(3*(dim*2+1)),
    NVF_(3*(dim*2+1)),
    veloflux_(dim*2+1),
    temp_veloflux_(sg.Mesh().Nodes()),
    RESULT_(dim*2+1,3),
    components_(dim*2+1),
    node_output_(sg.Mesh().Nodes(),false),
    nodal_averaging_(node_averaging),
    verbose_(false),
    with_gravity_(false),
    with_multiplier_(false),
    ac_gravity_(ACC_GRAVITY),  // m s-2
    VERTICAL_AXIS_( (dim==1u) ? 0u : 1u )
 {
    MathOperatorLHS<dim,CELL>::Name("VelocityAndVolumeFlux", oper, basic, test );

    // getting ranges for output variables
    sg.Database().RangeOf(velocity, minmaxV_.first, minmaxV_.second );
    sg.Database().RangeOf(volume_flux, minmaxF_.first, minmaxF_.second );
    
    // testing the Operands 
    if ( this->MaterialOperandPlacement() != ELEMENT and this->MaterialOperandPlacement() != ELEMENT_INTEGRATION_POINT  )
      throw csmp::Exception( ERROR, "VelocityAndVolumeFlux::(constructor)", 
                   oper, " must be either an element or a constraint-point property." );

     if ( MathOperatorLHS<dim,CELL>::BasicOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "VelocityAndVolumeFlux::(constructor)", 
                   basic, "Basic Operand must be a scalar property." );

    if ( MathOperatorLHS<dim,CELL>::TestOperandPlacement() != NODE ||
         MathOperatorLHS<dim,CELL>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "VelocityAndVolumeFlux::(constructor)", 
                   test, "Operand 'fluid pressure' must be a scalar property placed on the nodes." );

    if ( velo_key_.place != ELEMENT || velo_key_.type != VECTOR )
     throw csmp::Exception( ERROR, "VelocityAndVolumeFlux::(constructor)", 
                   "Operand 'velocity' must be a vector property placed on the element." );

    if ( ivelo_key_.place != ELEMENT || ivelo_key_.type != VECTOR )
     throw csmp::Exception( ERROR, "VelocityAndVolumeFlux::(constructor)", 
                   "Operand 'pore velocity' must be a vector property placed on the element." );

    if ( flux_key_.place != ELEMENT || flux_key_.type != SCALAR )
      throw csmp::Exception( ERROR, "VelocityAndVolumeFlux::(constructor)", 
                   "Operand 'volume flux' must be a scalar property placed on the element." );

    if ( node_averaging ) 
      {
         nvelo_key_  = sg.Database().StorageKey(nodal_velocity);
         nivelo_key_ = sg.Database().StorageKey(nodal_pore_velocity);
         nflux_key_  = sg.Database().StorageKey(nodal_volume_flux);
        
         if ( nvelo_key_.place != NODE || nvelo_key_.type != VECTOR )
           throw csmp::Exception( ERROR, "VelocityAndVolumeFlux::(constructor)", 
                   "Operand 'nodal velocity' must be a vector property placed on the node." );

         if ( nivelo_key_.place != NODE || ivelo_key_.type != VECTOR )
           throw csmp::Exception( ERROR, "VelocityAndVolumeFlux::(constructor)", 
                   "Operand 'nodal pore velocity' must be a vector property placed on the node." );

         if ( nflux_key_.place != NODE || nflux_key_.type != SCALAR )
           throw csmp::Exception( ERROR, "VelocityAndVolumeFlux::(constructor)", 
                   "Operand 'nodal volume flux' must be a scalar property placed on the node." );

         // if values are to be averagded on the nodes, the operator needs to be applied twice
         MathOperatorLHS<dim,CELL>::ApplicationCycles(2);
      }
   
   csmp::ErrorHandler& csmp_error( ErrorHandler::Instance() );
   
   if ( WithLowerDimensionalElements(sg) ) {
        csmp_error.Note( WARNING, "VelocityAndVolumeFlux<dim,CELL>::VelocityAndVolumeFlux:",
                          "your model contains lower dimensional elements, have you taken care of a thickness attribute?" );
     }
    
 } // end constructor




/**
    Like in first case, but considering buoyancy forces.
    
    @param relative_density 'relative density' is the fluid density - a reference density, for instance 1000 kg/m3
*/
template<uint32_t dim, template<uint32_t> class CELL>
VelocityAndVolumeFlux<dim,CELL>::VelocityAndVolumeFlux( const Model<dim>& sg,
                                                        const char* oper,   // conductivity
                                                        const char* basic,  // porosity
                                                        const char* test,   // fluid pressure
                                                        const char* relative_density, // relative fluid density
                                                        bool  node_averaging,
                                                        const char* velocity,
                                                        const char* pore_velocity,
                                                        const char* volume_flux,
                                                        const char* nodal_velocity,
                                                        const char* nodal_pore_velocity,
                                                        const char* nodal_volume_flux)

  : MathOperatorLHS<dim,CELL>(sg.Database(),oper,basic,test),
    PF_(3),
    VELOFLUX_(dim+1),
    IVELOFLUX_(dim),
    IPVF_(3*(dim*2+1)),
    NVF_(3*(dim*2+1)),
    veloflux_(dim*2+1),
    temp_veloflux_(sg.Mesh().Nodes()),
    RESULT_(dim*2+1,3),
    components_(dim*2+1),
    node_output_(sg.Mesh().Nodes(),false),
    nodal_averaging_(node_averaging),
    verbose_(false),
    with_gravity_(true),
    with_multiplier_(false),
    ac_gravity_(ACC_GRAVITY),  // m s-2
    VERTICAL_AXIS_( (dim==1u) ? 0u : 1u )
 {
    MathOperatorLHS<dim,CELL>::Name("VelocityAndVolumeFlux", oper, basic, test );

    // getting the necessary csmp::Index keys
    velo_key_   = sg.Database().StorageKey( velocity );
    ivelo_key_  = sg.Database().StorageKey( pore_velocity );
    flux_key_   = sg.Database().StorageKey( volume_flux );
    rhor_key_   = sg.Database().StorageKey( relative_density );
    
    // getting ranges for output variables
    sg.Database().RangeOf( velocity, minmaxV_.first, minmaxV_.second );
    sg.Database().RangeOf( volume_flux, minmaxF_.first, minmaxF_.second );
    
    // testing the Operands 
    if ( this->MaterialOperandPlacement() != ELEMENT and this->MaterialOperandPlacement() != ELEMENT_INTEGRATION_POINT  )
      throw csmp::Exception( ERROR, "VelocityAndVolumeFlux::(constructor)", 
                   oper, " must be either an element or a constraint-point property." );

    if ( MathOperatorLHS<dim,CELL>::BasicOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "VelocityAndVolumeFlux::(constructor)", 
                   basic, "Basic Operand must be a scalar property." );

    if ( MathOperatorLHS<dim,CELL>::TestOperandPlacement() != NODE ||
         MathOperatorLHS<dim,CELL>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "VelocityAndVolumeFlux::(constructor)", 
                   test, "Operand 'fluid pressure' must be a scalar property placed on the nodes." );

    if ( velo_key_.place != ELEMENT || velo_key_.type != VECTOR )
     throw csmp::Exception( ERROR, "VelocityAndVolumeFlux::(constructor)", 
                   "Operand 'velocity' must be a vector property placed on the element." );

    if ( ivelo_key_.place != ELEMENT || ivelo_key_.type != VECTOR )
     throw csmp::Exception( ERROR, "VelocityAndVolumeFlux::(constructor)", 
                   "Operand 'pore velocity' must be a vector property placed on the element." );

    if ( flux_key_.place != ELEMENT || flux_key_.type != SCALAR )
      throw csmp::Exception( ERROR, "VelocityAndVolumeFlux::(constructor)", 
                   "Operand 'volume flux' must be a scalar property placed on the element." );

    if ( rhor_key_.type != SCALAR )
      throw csmp::Exception( ERROR, "VelocityAndVolumeFlux::(constructor)", 
                   "Fluid density operand must be a scalar property." );

    if ( node_averaging ) 
      {
         nvelo_key_  = sg.Database().StorageKey(nodal_velocity);
         nivelo_key_ = sg.Database().StorageKey(nodal_pore_velocity);
         nflux_key_  = sg.Database().StorageKey(nodal_volume_flux);
        
         if ( nvelo_key_.place != NODE || nvelo_key_.type != VECTOR )
           throw csmp::Exception( ERROR, "VelocityAndVolumeFlux::(constructor)", 
                   "Operand 'nodal velocity' must be a vector property placed on the node." );

         if ( nivelo_key_.place != NODE || ivelo_key_.type != VECTOR )
           throw csmp::Exception( ERROR, "VelocityAndVolumeFlux::(constructor)", 
                   "Operand 'nodal pore velocity' must be a vector property placed on the node." );

         if ( nflux_key_.place != NODE || nflux_key_.type != SCALAR )
           throw csmp::Exception( ERROR, "VelocityAndVolumeFlux::(constructor)", 
                   "Operand 'nodal volume flux' must be a scalar property placed on the node." );

         // if values are to be averagded on the nodes, the operator needs to be applied twice
         MathOperatorLHS<dim,CELL>::ApplicationCycles(2);
      }
   
   if ( WithLowerDimensionalElements(sg) )
     throw csmp::Exception( ERROR, "VelocityAndVolumeFlux<dim,CELL>::VelocityAndVolumeFlux:",
                "post-processing of gravity-influenced flows in the lower-dimensional elements will not work with this visitor." );
   
 } // end constructor





template<uint32_t dim, template<uint32_t> class CELL>
VelocityAndVolumeFlux<dim,CELL>::VelocityAndVolumeFlux( const Model<dim>& sg,
                                                        const char* oper,   // conductivity
                                                        const char* basic,  // porosity
                                                        const char* test,   // fluid pressure
                                                        const char* relative_density, // relative fluid density
                                                        const char* prop_multiplier, // operand multiplier
                                                        bool  node_averaging,
                                                        const char* velocity,
                                                        const char* pore_velocity,
                                                        const char* volume_flux,
                                                        const char* nodal_velocity,
                                                        const char* nodal_pore_velocity,
                                                        const char* nodal_volume_flux)
  : MathOperatorLHS<dim,CELL>(sg.Database(),oper,basic,test),
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
    rho_vec_(3),
    components_(dim*2+1),
    node_output_(sg.Mesh().Nodes(),false),
    nodal_averaging_(node_averaging),
    verbose_(false),
    with_gravity_(true),
    with_multiplier_(true),
    ac_gravity_(ACC_GRAVITY),  // m s-2
    VERTICAL_AXIS_( (dim==1u) ? 0u : 1u )
 {
    MathOperatorLHS<dim,CELL>::Name("VelocityAndVolumeFlux", oper, basic, test );

    // getting the necessary csmp::Index keys
    velo_key_   = sg.Database().StorageKey( velocity );
    ivelo_key_  = sg.Database().StorageKey( pore_velocity );
    flux_key_   = sg.Database().StorageKey( volume_flux );
    rhor_key_   = sg.Database().StorageKey( relative_density );
    mult_key_   = sg.Database().StorageKey( prop_multiplier );
    
    // getting ranges for output variables
    sg.Database().RangeOf(velocity, minmaxV_.first, minmaxV_.second );
    sg.Database().RangeOf(volume_flux, minmaxF_.first, minmaxF_.second );
    
    // testing the Operands 
    if ( MathOperatorLHS<dim,CELL>::MaterialOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "VelocityAndVolumeFlux::(constructor)", 
                     oper, "Operand must be a scalar property." );

    if ( MathOperatorLHS<dim,CELL>::BasicOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "VelocityAndVolumeFlux::(constructor)", 
                   basic, "Basic Operand must be a scalar property." );

    if ( MathOperatorLHS<dim,CELL>::TestOperandPlacement() != NODE ||
         MathOperatorLHS<dim,CELL>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "VelocityAndVolumeFlux::(constructor)", 
                   test, "Operand 'fluid pressure' must be a scalar property placed on the nodes." );

    if ( velo_key_.place != ELEMENT || velo_key_.type != VECTOR )
     throw csmp::Exception( ERROR, "VelocityAndVolumeFlux::(constructor)", 
                   "Operand 'velocity' must be a vector property placed on the element." );

    if ( ivelo_key_.place != ELEMENT || ivelo_key_.type != VECTOR )
     throw csmp::Exception( ERROR, "VelocityAndVolumeFlux::(constructor)", 
                   "Operand 'pore velocity' must be a vector property placed on the element." );

    if ( flux_key_.place != ELEMENT || flux_key_.type != SCALAR )
      throw csmp::Exception( ERROR, "VelocityAndVolumeFlux::(constructor)", 
                   "Operand 'volume flux' must be a scalar property placed on the element." );

    if ( rhor_key_.type != SCALAR )
      throw csmp::Exception( ERROR, "VelocityAndVolumeFlux::(constructor)", 
                   "Fluid density operand must be a scalar property." );

    if (  mult_key_.place != NODE || mult_key_.type != SCALAR )
      throw csmp::Exception( ERROR, "VelocityAndVolumeFlux::(constructor)", 
                   "Property multiplier must be a scalar property placed on the node." );

    if ( node_averaging ) 
      {
         nvelo_key_  = sg.Database().StorageKey(nodal_velocity);
         nivelo_key_ = sg.Database().StorageKey(nodal_pore_velocity);
         nflux_key_  = sg.Database().StorageKey(nodal_volume_flux);
        
         if ( nvelo_key_.place != NODE || nvelo_key_.type != VECTOR )
           throw csmp::Exception( ERROR, "VelocityAndVolumeFlux::(constructor)", 
                   "Operand 'nodal velocity' must be a vector property placed on the node." );

         if ( nivelo_key_.place != NODE || ivelo_key_.type != VECTOR )
           throw csmp::Exception( ERROR, "VelocityAndVolumeFlux::(constructor)", 
                   "Operand 'nodal pore velocity' must be a vector property placed on the node." );

         if ( nflux_key_.place != NODE || nflux_key_.type != SCALAR )
           throw csmp::Exception( ERROR, "VelocityAndVolumeFlux::(constructor)", 
                   "Operand 'nodal volume flux' must be a scalar property placed on the node." );

         // if values are to be averagded on the nodes, the operator needs to be applied twice
         MathOperatorLHS<dim,CELL>::ApplicationCycles(2);
      }                                  

   if ( WithLowerDimensionalElements(sg) )
     throw csmp::Exception( ERROR, "VelocityAndVolumeFlux<dim,CELL>::VelocityAndVolumeFlux:",
                "post-processing of gravity-influenced flows in the lower-dimensional elements will not work with this visitor." );
   
 } // end constructor




template<uint32_t dim, template<uint32_t> class CELL>
VelocityAndVolumeFlux<dim,CELL>::~VelocityAndVolumeFlux()
 {
    // no dynamically allocated variables
 }



/** Switch to verbose mode (results are reported to stdout).
*/
template<uint32_t dim, template<uint32_t> class CELL>
void VelocityAndVolumeFlux<dim,CELL>::Verbose( bool stdoutput ) { verbose_=stdoutput; }


/**

Using the physically meaningful range which was specified in the Property-
Database, the output variables are tested for the correctness. Note that this
method assumes that the interstital velocity "pore velocity" should be
in the same range as the Darcy velocity.

@section messages Messages

If the range constraint is violated, a message is printed to 'cout'.
*/
template<uint32_t dim, template<uint32_t> class CELL>
void VelocityAndVolumeFlux<dim,CELL>::TestRangeOfOutputVariables() const
  {
     // velocity, interstitial velocity
     for ( uint32_t i{0U}; i<dim; i++ )
       {
          if ( velo_[i] < minmaxV_.first || velo_[i] > minmaxV_.second )
            throw csmp::Exception( ERROR, "VelocityAndVolumeFlux::TestRangeOfOutputVariables:",
                          "Result variable 'velocity' outside of range specified in database file." );

          if ( ivelo_[i] < minmaxV_.first || ivelo_[i] > minmaxV_.second )
            throw csmp::Exception( ERROR, "VelocityAndVolumeFlux::TestRangeOfOutputVariables:",
                          "Result variable 'pore velocity' outside of range specified in database file." );
       }
     // volume flux
     if ( flux_() < minmaxF_.first || flux_() > minmaxF_.second )
            throw csmp::Exception( ERROR, "VelocityAndVolumeFlux::TestRangeOfOutputVariables:",
                          "Result variable 'volume flux' outside of range specified in database file." );
  }





/**
    If there are elements from more than one spatial dimension in the model, we have a problem.
*/
template<uint32_t dim, template<uint32_t> class CELL>
bool VelocityAndVolumeFlux<dim,CELL>::WithLowerDimensionalElements( const Model<dim>& m ) const
 {
    std::pair<int32_t,int32_t>  dimensionality = m.Region("Model").ElementSpatialDimensions();
    if ( dimensionality.first > 1 ) return true;
   
    return false;
   
 } // end WithLowerDimensionalElements




template<uint32_t dim, template<uint32_t> class CELL>
void VelocityAndVolumeFlux<dim,CELL>::GetOperands( const CELL<dim>& e )
{
    if ( MathOperatorLHS<dim,CELL>::ApplicationCycle() == 1 ) {
    
    // fluid pressure
    e.NodePropertyVector( MathOperatorLHS<dim,CELL>::TestOperandKey(), PF_ );
 
    // relative density
    if ( with_gravity_ ) {
        if ( rhor_key_.place == ELEMENT )
          e.Read( rhor_key_, rhor_ );
        else if ( rhor_key_.place == NODE ) {
             e.PropertyValueAtBaryCenter( rhor_key_, rhor_ );
             e.NodePropertyVector( rhor_key_, rho_vec_ );
          }      
        else
        throw csmp::Exception( ERROR, "VelocityAndVolumeFlux<dim>::GetOperands", 
                       "fluid density is neither a node nor element variable; can't deal with this.");
      }
      
    // nodal multipliers
    if ( with_multiplier_ ) e.NodePropertyVector( mult_key_, mult_vec_ );
    
    // conductivity
    if ( MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT )
      {
          if ( MathOperatorLHS<dim,CELL>::MaterialOperandType() == SCALAR ) {
               MathOperatorLHS<dim,CELL>::MTRL[0].AssignToDiagonalAndZeroOffDiagonal( dim,
                                                 e.Read( MathOperatorLHS<dim,CELL>::MaterialOperandKey() ) );
            }
          else if ( MathOperatorLHS<dim,CELL>::MaterialOperandType() == VECTOR ) {
               VectorVariable<dim>  vc;
               e.Read( MathOperatorLHS<dim,CELL>::MaterialOperandKey(), vc );
               MathOperatorLHS<dim,CELL>::MTRL[0].AssignToDiagonal( vc );
            }
          else if ( MathOperatorLHS<dim,CELL>::MaterialOperandType() == TENSOR ) {
               TensorVariable<dim>  ts;
               e.Read( MathOperatorLHS<dim,CELL>::MaterialOperandKey(), ts );
               MathOperatorLHS<dim,CELL>::MTRL[0] = ts;
           }
      }
    else if ( MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT_INTEGRATION_POINT )
      {
         for ( auto i{0U}; i<e.FE()->IntegrationPoints(); i++ ) {
	          if ( MathOperatorLHS<dim,CELL>::MaterialOperandType() == SCALAR ) {
	               MathOperatorLHS<dim,CELL>::MTRL[i].AssignToDiagonalAndZeroOffDiagonal( dim,
	                                                    e.Read( i, MathOperatorLHS<dim,CELL>::MaterialOperandKey() ) );
	            }
	          else if ( MathOperatorLHS<dim,CELL>::MaterialOperandType() == VECTOR ) {
	               VectorVariable<dim>  vc;
	               e.Read( i, MathOperatorLHS<dim,CELL>::MaterialOperandKey(), vc );
	               MathOperatorLHS<dim,CELL>::MTRL[i].AssignToDiagonal( vc );
	            }
	          else if ( MathOperatorLHS<dim,CELL>::MaterialOperandType() == TENSOR ) {
	               TensorVariable<dim>  ts;
	               e.Read( i, MathOperatorLHS<dim,CELL>::MaterialOperandKey(), ts );
	               MathOperatorLHS<dim,CELL>::MTRL[i] = ts;
	           }
          }
      }
    else // if a nodal variable is dealt with
      {
         if ( e.FE()->IntegrationPoints() == 0U ) 
           throw csmp::Exception( FATAL_ERROR, "VelocityAndVolumeFlux<dim>::GetOperands", 
                                        "The current finite element has no integration points",
                                        "Therefore nodal properties cannot be integrated.");
      
         for ( auto i=0U; i<e.FE()->IntegrationPoints(); i++ )
           {   
              MathOperatorLHS<dim,CELL>::MTRL[i].Resize(dim,dim);
              MathOperatorLHS<dim,CELL>::MTRL[i].Zero();
              MathOperatorLHS<dim,CELL>::PropertyAtIntegrationPoint( e,
                                  MathOperatorLHS<dim,CELL>::MaterialOperandKey(), i, MathOperatorLHS<dim,CELL>::MTRL[i] );
           }
      }
         // porosity
         e.PropertyValueAtBaryCenter( MathOperatorLHS<dim,CELL>::BasicOperandKey(), phi_ );
     } 
    
} // end GetOperands






/**
 
@section arguments Input Arguments 

A reference to the Element for which the post-processing is done.  
*/
template<uint32_t dim, template<uint32_t> class CELL>
void VelocityAndVolumeFlux<dim,CELL>::ComputeContribution( const CELL<dim>& e )
{
   if ( MathOperatorLHS<dim,CELL>::ApplicationCycle() == 1 ) {
        if ( verbose_ )
          cout <<"\n\nVelocityAndVolumeFlux::ComputeContribution: Element: "<< e.Idx() << endl; 
      
        // 1. Case of linear (analytically integrated) elements: constant velocity
        // -----------------------------------------------------------------------
        if ( e.Interpolation() == 1 && 
             e.IntegrationPoints() == 0 )
          {
              // interpolation function derivatives at the nodes
              e.dN_AtBaryCenter( DERIV_ );
              velo_ = 0.;
              for ( uint32_t i{0U}; i<e.Nodes(); i++ )
                for ( uint32_t j{0U}; j<dim; j++ )
                  velo_(j) += PF_[i]() * -DERIV_(j,i) * MathOperatorLHS<dim,CELL>::MTRL[0](j,j);
              
              if ( with_gravity_ )
                velo_(VERTICAL_AXIS_) -= ac_gravity_ * rhor_() * MathOperatorLHS<dim,CELL>::MTRL[0](VERTICAL_AXIS_,VERTICAL_AXIS_);

              // interstitial velocity    
              ivelo_  = velo_;
              ivelo_ /= phi_;
              
              // flux
              flux_ = velo_.Length();
              
              if ( verbose_ ) {
                   cout <<"\ncomputed element variable 'velocity':"<< endl;
                   velo_.Out();
                   cout <<"\ncomputed element variable 'pore velocity':"<< endl;
                   ivelo_.Out();
                   cout <<"\ncomputed element variable 'volume flux': "<< flux_() << endl;
                }
          }

        // 2. Case of numerically integrated elements:
        //    Computing the velocities and the volume flux at the integration points
        // -------------------------------------------------------------------------
        else 
          {
             DERIV_.Resize(dim,static_cast<uint32_t>(e.Nodes()));
             IPVF_.resize( e.FE()->IntegrationPoints()*(dim*2+1) );
             NVF_.resize( e.Nodes()*(dim*2+1) );
             // collecting averadge data for element variables
             velo_ = ivelo_ = flux_ = 0.;
             
             for ( uint32_t i{0U}; i<e.IntegrationPoints(); i++ )
               {
                  // if density driven flow is computed, calculate rho * g * z
                  // comput density at integration points and any multiplier as well
                  if ( with_gravity_ ) {
                      // density element variable
                      if ( rhor_key_.place == ELEMENT ) {
                          rho_fac_ = rhor_();
                          if ( !with_multiplier_ ) mult_fac_ = 1.; // more likely case first
                          else {
                               e.N_AtIntegrationPoint( i, IPOL_ );
                               mult_fac_=0.;
                               for ( uint32_t j{0U}; j<e.Nodes(); j++ )
                                 mult_fac_ += IPOL_[j] * mult_vec_[j]();
                            }
                        }                      
                      // density nodal variable  
                      else if ( rhor_key_.place == NODE ) {
                          e.N_AtIntegrationPoint( i, IPOL_ );
                          rho_fac_=mult_fac_=0.0;
                          for ( auto j{0U}; j<e.Nodes(); j++ ) {
                              if ( with_multiplier_ ) mult_fac_ += IPOL_[j] * mult_vec_[j]();
                              else                   mult_fac_ = 1.0;
                              rho_fac_  += IPOL_[j] * rho_vec_[j]();
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
                  
                  if ( MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT )
                    for ( uint32_t n=0; n<e.Nodes(); n++ )
                      for ( uint32_t j{0U}; j<dim; j++ )
                        // -DERIV because fluid flows down pressure
                        VELOFLUX_[j] += PF_[n]() * -DERIV_(j,n) * MathOperatorLHS<dim,CELL>::MTRL[0](j,j);
                  else
                     for ( uint32_t n=0u; n<e.Nodes(); n++ )
                      for ( uint32_t j{0U}; j<dim; j++ )
                        VELOFLUX_[j] += PF_[n]() * -DERIV_(j,n) * MathOperatorLHS<dim,CELL>::MTRL[i](j,j);
            
                  if ( with_gravity_ ) VELOFLUX_[VERTICAL_AXIS_] -= ac_gravity_ *
                                                                 rho_fac_ * MathOperatorLHS<dim,CELL>::MTRL[0](VERTICAL_AXIS_,VERTICAL_AXIS_) *
                                                                 mult_fac_;

                  // interstitial velocity & volume flux
                  for ( uint32_t n=0u; n<dim; n++ ) {
                       IVELOFLUX_[n]   = VELOFLUX_[n]/phi_();
                       VELOFLUX_[dim] += VELOFLUX_[n]*VELOFLUX_[n];
                       // summing integration point values for later averaging
                       velo_(n)              += VELOFLUX_[n];
                       ivelo_(n)             += IVELOFLUX_[n];
                    }
                  flux_ += VELOFLUX_[dim] = sqrt(VELOFLUX_[dim]);

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
                       for ( uint32_t k=0u; k<dim; k++ ) {
                            IPVF_[ i*components_ + k ] = VELOFLUX_[k];
                            IPVF_[ i*components_ + dim + 1 + k ] = IVELOFLUX_[k];
                         }
                       IPVF_[ i*components_ + dim ] = VELOFLUX_[dim];
                    }
               }              
             // Averaging integration point values to get the element variables
             // ---------------------------------------------------------------
             velo_  /= static_cast<double>(e.FE()->IntegrationPoints());
             flux_  /= static_cast<double>(e.FE()->IntegrationPoints());
             ivelo_ /= static_cast<double>(e.FE()->IntegrationPoints());
         }
        
       if ( nodal_averaging_ )
         {
             // if the computed properties are constant on the element
             if ( e.Interpolation() == 1 and
                  e.IntegrationPoints() == 0 )
               {
                  for ( uint32_t k=0; k<dim; k++ )
                    {
                       veloflux_[k]           = velo_[k];
                       veloflux_[dim + 1 + k] = ivelo_[k];
                    }
                  veloflux_[dim] = flux_();
                    
                  for ( uint32_t i{0U}; i<e.Nodes(); i++ )
                    temp_veloflux_[ e.N(i)->Idx() ].push_back( veloflux_ );
               }
             else
               {
                  // 1.2 Extrapolating velocity and volume flux from the 
                  //     integration points to the nodes.
                  // ---------------------------------------------------
                  e.ExtrapolateIntegrationPointVariableToNodes( components_, IPVF_, NVF_ );
     
                  for ( uint32_t i{0U}; i<e.Nodes(); i++ )
                    {
                       for ( uint32_t k=0; k<components_; k++ ) veloflux_[k] = NVF_[ i*components_ + k ];
                       temp_veloflux_[ e.N(i)->Idx() ].push_back( veloflux_ );
                    }
               }
          }

     }  // end of first application cycle  


   // 2. During the second visitation, the nodal velocities and fluxes computed for each element
   //    node are averaged and stored in a vector for output.
   //    -------------------------------------------
   if ( MathOperatorLHS<dim,CELL>::ApplicationCycle() == 2 && nodal_averaging_ ) {
        RESULT_.Resize(components_,static_cast<uint32_t>(e.Nodes()));
        for ( auto i{0U}; i<e.Nodes(); i++ ) {
          // duplicate calculations are avoided via the boolean vector
          if ( !node_output_[ e.N(i)->Idx() ] )
            {
               for ( uint32_t j{0U}; j<components_; j++ )
                 {  
                    // averaging velocity/flux components
                    sum_=0.0;
                    for ( auto lit =temp_veloflux_[ e.N(i)->Idx() ].begin();
                          lit!=temp_veloflux_[ e.N(i)->Idx() ].end(); lit++ )
                          sum_ += (*lit)[j];
                    sum_ /= static_cast<double>(temp_veloflux_[ e.N(i)->Idx() ].size());
                    RESULT_(j,i) = sum_;
                 }
             }

          // now the vector list is no longer needed and therefore erased for the 
          // next application of the post-processing operator
          temp_veloflux_[ e.N(i)->Idx() ].erase(
                         temp_veloflux_[ e.N(i)->Idx() ].begin(),
                         temp_veloflux_[ e.N(i)->Idx() ].end() );
       }

   } // end application cycle 2

} // end ComputeContribution






/**

Outputs the computed velocities and fluxes to the Model storage
checking the range of each variable using the ranges associated with each
of the Operands.


@section arguments Input Arguments

A reference to the property memory manager and the Element for which the
variables are output.
*/
template<uint32_t dim, template<uint32_t> class CELL>
void VelocityAndVolumeFlux<dim,CELL>::WriteOperands( CELL<dim>& e )
 {
    if ( MathOperatorLHS<dim,CELL>::ApplicationCycle() == 1 )
      {
         // 1. outputting element properties first
         // --------------------------------------
         TestRangeOfOutputVariables();
         e.Store(  velo_key_, velo_ );
         e.Store( ivelo_key_, ivelo_ );
         e.Store(  flux_key_, flux_ );
      }
    if ( MathOperatorLHS<dim,CELL>::ApplicationCycle() == 2 )
      {
         if ( nodal_averaging_ )
           for ( uint32_t i{0U}; i<e.Nodes(); i++ )
             // doing this operation only once per node
             if ( !node_output_[ e.N(i)->Idx() ] )
               {
                  // 1. assigning the strain & stress values
                  // ---------------------------------------
                  ExtractVelocity( RESULT_, i, velo_ );
                  e.N(i)->Store( nvelo_key_, velo_ );
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



 





template<uint32_t dim, template<uint32_t> class CELL>
void VelocityAndVolumeFlux<dim,CELL>::ExtractVelocity( const DenseMatrix<DM_MIN>& INP,
                                                       uint32_t        col,
                                                       VectorVariable<dim>& vc )
 {
    for ( auto i{0U}; i<dim; i++ ) vc(i) = INP(i,col);
 }







template<uint32_t dim, template<uint32_t> class CELL>
void VelocityAndVolumeFlux<dim,CELL>::ExtractVolumeFlux( const DenseMatrix<DM_MIN>&   INP,
                                                            uint32_t    col,
                                                            ScalarVariable& sc )
 {
    sc() = INP(dim,col);
 }

 

template<uint32_t dim, template<uint32_t> class CELL>
void VelocityAndVolumeFlux<dim,CELL>::ExtractInterstitialVelocity( const DenseMatrix<DM_MIN>& INP,
                                                                   uint32_t col,
                                                                   VectorVariable<dim>& vc )
 {
    for ( auto i{0U}; i<dim; i++ ) vc(i) = INP(i+dim+1,col);
 }
 


template class VelocityAndVolumeFlux<1U>;
template class VelocityAndVolumeFlux<2U>;
template class VelocityAndVolumeFlux<3U>;

template class VelocityAndVolumeFlux<1U,Face>;
template class VelocityAndVolumeFlux<2U,Face>;
template class VelocityAndVolumeFlux<3U,Face>;

}
 
 
 
 
 
 
 
 
 



