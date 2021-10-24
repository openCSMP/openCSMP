#include "ComputeGravityTermVisitor.h"
#include "Model.h"
#include "Region.h"
#include "TwoPhaseModel.h"
#include "ScalarVariable.h"
#if defined(_OPENMP )
#include "omp.h"
#endif

using namespace std;

namespace csmp{

template<size_t dim>
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
#if defined(_OPENMP )
    this->femgrs_.resize(omp_get_max_threads());
    for (long int tid = 0 ; tid < omp_get_max_threads();tid ++)
        this->femgrs_[tid].InitializeElements(dim,model.FE_Manager().InterpolationOrder(),true);
#endif
}

template<size_t dim>
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
#if defined(_OPENMP )
    this->femgrs_.resize(omp_get_max_threads());
    for (long int tid = 0 ; tid < omp_get_max_threads();tid ++)
        this->femgrs_[tid].InitializeElements(dim,model.FE_Manager().InterpolationOrder(),true);
#endif
}

template<size_t dim>
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
#if defined(_OPENMP )
    this->femgrs_.resize(omp_get_max_threads());
    for (long int tid = 0 ; tid < omp_get_max_threads();tid ++)
        this->femgrs_[tid].InitializeElements(dim,model.FE_Manager().InterpolationOrder(),true);
#endif
}

template<size_t dim>
void ComputeGravityTermVisitor<dim>::Visit( Model<dim>* m ){
    if (this->Verbose()) cout <<" ComputeGravityTermVisitor<dim>::Visit(Model<dim>*)"<<endl;
#if defined(_OPENMP )
    this->Visit(&(m->Region("Model")));
#endif
}

template<size_t dim>
void ComputeGravityTermVisitor<dim>::Visit(Region<dim>* region ){
    if (this->Verbose()) cout <<" ComputeGravityTermVisitor<dim>::Visit(Region<dim>*) : "<<region->Name()<<endl;

#if defined(_OPENMP )
#pragma omp parallel
    {
        Element<dim>* ep;
        FiniteElement* fe_tmp;
        size_t tid=omp_get_thread_num();
#pragma omp for
        for ( long int e= 0 ; e < region->Elements(); e++ ){
            ep = region->E(e);
            fe_tmp=ep->FE();
            // change pointer here
            ep->Assign(femgrs_[tid].E(ep->FE_Type()));
            this->ComputeContribution(ep);
            // put it back here
            ep->Assign(fe_tmp);
        }
    }
#endif
}

template<size_t dim>
void ComputeGravityTermVisitor<dim>::Visit( Element<dim>* element )
{
#if !defined(_OPENMP)
    this->ComputeContribution(element);
#endif
}

template<size_t dim>
void ComputeGravityTermVisitor<dim>::ComputeContribution( Element<dim>* element )
{

    //! two cases are currently supported by the code
    //! 1) twophase with gravity on the ELEMENT
    //! 2) singlephase with gravity on the ELEMENT_INTEGRATION_POINT
    
    ScalarVariable mu, rho;
    double gravityTerm;
    VectorVariable<dim> gravityVector;
    
    if (saturationFunctions_ != NULL)
    {
#if defined(_OPENMP )
        ErrorHandler& error_handler (ErrorHandler::Instance());
        string errmsg="Unfortunately two-phase functions are not yet Openmp-ized.";
        error_handler.notice(FATAL_ERROR,"ComputeGravityTermVisitor<dim>::ComputeContribution(element)",errmsg);
#endif
      saturationFunctions_->Initialize( *element );
      saturationFunctions_->InitializeForBaryCenter( *element);
      saturationFunctions_->EffectiveSaturation();
    
      // Taking into account the gravitational flow component for calculation of velocity of non-wetting phase
      // k ( Lambda Water * rho_w + Lambda Oil * Rho_o) g
      gravityTerm = element->Read( permeabilityKey_ ) * gravityAcc_ *
                   ( saturationFunctions_->MobilityPhase(1) * saturationFunctions_->DensityWettingPhase() +
                     saturationFunctions_->MobilityPhase(2) * saturationFunctions_->DensityNonWettingPhase() );
    }

    if( element->FE()->IsLineElement() )
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
    else if( element->FE()->IsSurfaceElement())
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
      for (size_t ip=0;ip<element->IntegrationPoints (); ++ip)
      {
        element->PropertyValueAtIntegrationPoint( singlePhaseDensityKey_, ip, rho );
        element->PropertyValueAtIntegrationPoint( singlePhaseViscosityKey_, ip, mu );
        
        element->Store( ip, gravityVectorKey_, gravityVector*gravityTerm/mu()*rho() );
      }
    }
    else // two phase -> ELEMENT
      element->Store( gravityVectorKey_, gravityVector * gravityTerm);

}

template class ComputeGravityTermVisitor<3>;
template class ComputeGravityTermVisitor<2>;
template class ComputeGravityTermVisitor<1>;

} //csmp
