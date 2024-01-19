#include "ComputeGravityTermVisitor.h"
#include "Model.h"
#include "Region.h"
#include "TwoPhaseModel.h"
#include "ScalarVariable.h"

using namespace std;

namespace csmp {

template<uint32_t dim>
ComputeGravityTermVisitor<dim>::ComputeGravityTermVisitor( Model<dim>& model,
                                                           TwoPhaseModel<dim>& saturationFunctions,
                                                           const char* gravityVectorTag,
                                                           const char* permeabilityTag )
    : Visitor<dim>( MODEL, ELEMENT ), model_( model ),
      saturationFunctions_( &saturationFunctions ),
      gravityVectorKey_ (model.Database().StorageKey( gravityVectorTag ) ),
      permeabilityKey_( model.Database().StorageKey( permeabilityTag ) ),
      singlePhaseViscosityKey_ ( ),
      singlePhaseDensityKey_( ),
      gravityAcc_( 9.80665 )
{
}

template<uint32_t dim>
ComputeGravityTermVisitor<dim>::ComputeGravityTermVisitor( Model<dim>& model,
                                                           const char* gravityVectorTag,
                                                           const char* permeabilityTag, 
                                                           const char* viscosityTag,
                                                           const char* densityTag )
    : Visitor<dim>( MODEL, ELEMENT ), model_( model ),
      saturationFunctions_( NULL ),
      gravityVectorKey_ (model.Database().StorageKey( gravityVectorTag ) ),
      permeabilityKey_( model.Database().StorageKey( permeabilityTag ) ),
      singlePhaseViscosityKey_ ( model.Database().StorageKey( viscosityTag ) ),
      singlePhaseDensityKey_( model.Database().StorageKey( densityTag ) ),
      gravityAcc_( 9.80665 )
{
}



template<uint32_t dim>
ComputeGravityTermVisitor<dim>::ComputeGravityTermVisitor(Model<dim>& model,
                                                           Index gravityVector,
                                                           Index permeability,
                                                           Index viscosity,
                                                           Index density )
    : Visitor<dim>( MODEL, ELEMENT ), model_( model ),
      saturationFunctions_( NULL ),
      gravityVectorKey_ (gravityVector ),
      permeabilityKey_( permeability ),
      singlePhaseViscosityKey_ ( viscosity ),
      singlePhaseDensityKey_( density ),
      gravityAcc_( 9.80665 )
{
}


template<uint32_t dim>
void ComputeGravityTermVisitor<dim>::Visit( Element<dim>* element )
{
    //! two cases are currently supported by the code
    //! 1) twophase with gravity on the ELEMENT
    //! 2) singlephase with gravity on the ELEMENT_INTEGRATION_POINT
    
    ScalarVariable      mu, rho;
    double              gravityTerm = std::numeric_limits<double>::quiet_NaN();
    VectorVariable<dim> gravityVector;
    
    if (saturationFunctions_ != NULL)
    {
      saturationFunctions_->Initialize( *element );
      saturationFunctions_->InitializeForBaryCenter( *element);
      saturationFunctions_->EffectiveSaturation();
    
      // Taking into account the gravitational flow component for calculation of velocity of non-wetting phase
      // k ( Lambda Water * rho_w + Lambda Oil * Rho_o) g
      gravityTerm = element->Read( permeabilityKey_ ) * gravityAcc_ *
                   ( saturationFunctions_->MobilityPhase(1) * saturationFunctions_->DensityWettingPhase() +
                     saturationFunctions_->MobilityPhase(2) * saturationFunctions_->DensityNonWettingPhase() );
    }

    if( element->FE()->IsLine() )
    {
      //! projection of a vector g(0,-1,0) onto a vector r(r[0],r[1],r[2]) 
      //! (where |r|=1) is equal to
      //! (r, g)r = -r[1] r

      Point<dim> line_vector( element->N(0)->Coordinate() - element->N(1)->Coordinate() );

      line_vector.NormalizeLengthTo(1.);
      
      line_vector *= -line_vector[1];

      //! r*=(r[1]>0.0?-1:1) Note that the direction of the gravity force should be oriented downwards
      if(line_vector[1]>0) line_vector*=-1;

      gravityVector(0) = line_vector[0];
      gravityVector(1) = line_vector[1];
      gravityVector(2) = line_vector[2];

    }
    else if( element->FE()->IsSurface())
    {
      //! projection of vector g(0,-1, 0) onto a surface with a normal vector n(n[0], n[1], n[2]) 
      //! (where |n|=1) is equal to
      //! n x (n x g) = (n, g) n - g = -n[1] n - g 

      Point<dim> line_vector1( element->N(0)->Coordinate() - element->N(1)->Coordinate() );
      Point<dim> line_vector2( element->N(1)->Coordinate() - element->N(2)->Coordinate() );

      Point<dim> surface_normal( crossProduct( line_vector1,line_vector2) );
      
      surface_normal.NormalizeLengthTo (1.);
      
      surface_normal *= -surface_normal[1];
      surface_normal[1]+=1.0;

      //! n*=(n[1]>0.0?-1:1) Note that the direction of the gravity force should be oriented downwards
      if( surface_normal[1] > 0. ) surface_normal *= -1.;

      gravityVector(0) = surface_normal[0];
      gravityVector(1) = surface_normal[1];
      gravityVector(2) = surface_normal[2];
    }
    else
    {
      gravityVector(0) = 0.;
      gravityVector(1) = -1.;
      gravityVector(2) = 0.;
    }

    if (saturationFunctions_ == NULL) // single phase -> INTEGRATION_POINT
    {
      // k/mu * rho * g
      gravityTerm = element->Read( permeabilityKey_ ) * gravityAcc_;
      // loop over element integration points
      for ( auto ip=0;ip<element->IntegrationPoints(); ++ip)
      {
        element->PropertyValueAtIntegrationPoint( singlePhaseDensityKey_, ip, rho );
        element->PropertyValueAtIntegrationPoint( singlePhaseViscosityKey_, ip, mu );
        
        element->Store( ip, gravityVectorKey_, gravityVector*gravityTerm/mu()*rho() );
      }
    }
    else // two phase -> ELEMENT
      element->Store( gravityVectorKey_, gravityVector * gravityTerm );

}


template class ComputeGravityTermVisitor<3>;
template class ComputeGravityTermVisitor<2>;
template class ComputeGravityTermVisitor<1>;

} //csmp
