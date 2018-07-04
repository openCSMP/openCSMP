#include "FractureComputeVelocityVisitor.h"
#include "Model.h"

namespace csmp{

  template<size_t dim>
  FractureComputeVelocityVisitor<dim>::FractureComputeVelocityVisitor( Model<dim>& model, TwoPhaseModel<dim>& saturationFunctions,
                                                                       const char* vt_Tag,
                                                                       const char* vn_Tag,
                                                                       const char* fluidPressureTag,
                                                                       const char* fractureCapillaryPressureTag,
                                                                       const char* volFluxTag,
                                                                       const char* prevVolFluxTag,
                                                                       const char* gravityVectorTag,
                                                                       const char* permeabilityTag,
                                                                       bool withCapillaryGradient,
                                                                       bool withGravity )
      : Visitor<dim>( MODEL, ELEMENT ), model_( model ), saturationFunctions_( &saturationFunctions ),
        vt_Key_( model.Database().StorageKey( vt_Tag ) ),
        vn_Key_( model.Database().StorageKey( vn_Tag ) ),
        fluidPressureKey_( model.Database().StorageKey( fluidPressureTag ) ),
        fracCapillaryPressureKey_( model.Database().StorageKey( fractureCapillaryPressureTag ) ),
        volumeFluxKey_( model.Database().StorageKey( volFluxTag ) ),
        previousVolumeFluxKey_( model.Database().StorageKey( prevVolFluxTag ) ),
        gravityVectorKey_ (model.Database().StorageKey( gravityVectorTag ) ),
        permeabilityKey_( model.Database().StorageKey( permeabilityTag ) ),
        withCapillaryGradient_( withCapillaryGradient ),
        withGravity_( withGravity ),
        gravityAcc_( 9.80665 ),
        DERIV_( dim, dim )
    {
    }


  template<size_t dim>
  void FractureComputeVelocityVisitor<dim>::Visit( Element<dim>* element )
  {
    velo_ = 0.;

    saturationFunctions_->Initialize( *element );
    saturationFunctions_->InitializeForBaryCenter( *element);
    saturationFunctions_->EffectiveSaturation();

    if( withCapillaryGradient_ )
    {
      // Computing the total velocity: vt = -k (lt ( grad p + grad pc) - (s_w rho_w + s_o rho_o) g)
      // first without gravity vt = -k (lt grad p)
      element->dN_AtBaryCenter( DERIV_, 1U );

      for ( size_t i = 0; i < element->Nodes(); ++i )
      {
        double pf = element->N(i)->Read( fluidPressureKey_ );
        double nodalFracPc = element->N(i)->Read( fracCapillaryPressureKey_ );
        for( size_t xyz = 0; xyz < dim; ++xyz )
          velo_( xyz ) += ( ( pf + nodalFracPc ) * -DERIV_( xyz, i ) ) * saturationFunctions_->TotalMobility(); // total mobility already includes multiplication with k
      }
    }
    else
    {
        // Computing the total velocity: vt = -k (lt grad p - (s_w rho_w + s_o rho_o) g)
        // first without gravity vt = -k (lt grad p)
        element->dN_AtBaryCenter( DERIV_, 1U );
        for ( size_t i = 0; i < element->Nodes(); ++i )
        {
          double64 pf = element->N(i)->Read( fluidPressureKey_ );
          for( size_t xyz = 0; xyz < dim; ++xyz )
            velo_( xyz ) += pf * -DERIV_( xyz, i ) * saturationFunctions_->TotalMobility(); // total mobility already includes multiplication with k
        }
    }

    if( withGravity_ )
    {
      // Taking into account the gravitational flow component for calculation of velocity of non-wetting phase
      // k ( Lambda Water * rho_w + Lambda Oil * Rho_o) g
      element->Read( gravityVectorKey_, gravityVector_ );

      //adding the gravitational flow component to the velocity vector
      for( size_t xyz = 0; xyz < dim; ++xyz )
        velo_( xyz ) += gravityVector_( xyz );
    }

    // storing the computed velocity
    element->Store( vt_Key_, velo_ );

    // backing up the previous volume flux
    element->Read( volumeFluxKey_, flux_ );
    element->Store( previousVolumeFluxKey_, flux_ );

    // volume flux
    flux_ = velo_.Length();
    element->Store( volumeFluxKey_, flux_ );

    // computing vn from vt
    // vn = dfds vt + k  d l_overbar / ds (rhow - rhor) g
    velo_ *= saturationFunctions_->dfds();
    velo_ -= element->Read( permeabilityKey_ ) * saturationFunctions_->dGds()
           * ( saturationFunctions_->DensityNonWettingPhase() - saturationFunctions_->DensityWettingPhase() ) * gravityAcc_;

    element->Store( vn_Key_, velo_);
  }



  template class FractureComputeVelocityVisitor<3>;

} //csmp
