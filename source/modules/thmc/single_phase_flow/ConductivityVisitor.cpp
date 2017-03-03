#include "ConductivityVisitor.h"

#include "CapillaryNumberVisitor.h"
#include "Model.h"
#include "ScalarVariable.h"
#include "VectorVariable.h"
#if defined(_OPENMP )
#include "omp.h"
#endif

namespace csmp{


template<size_t dim>
ConductivityVisitor<dim>::ConductivityVisitor(Model<dim>& model,
                                              const char* specific_saturated_hydraulic_conductivity, // this is without mult. by density
                                              const char* permeability,
                                              const char* viscosity,
                                              const char* density,
                                              const char* saturated_hydraulic_conductivity,
                                              const char* compressibility,
                                              const char* porosity,
                                              const char* specific_saturated_hydraulic_diffusivity)
    : Visitor<dim>( MODEL, ELEMENT ), model_(model),
      sshcKey_( model.Database().StorageKey(specific_saturated_hydraulic_conductivity) ),
      kKey_   ( model.Database().StorageKey(permeability) ),
      muKey_  ( model.Database().StorageKey(viscosity) ),
      rhoKey_ ( ((density==NULL) ? csmp::Index() : model.Database().StorageKey(density))),
      shcKey_ ( ((saturated_hydraulic_conductivity==NULL) ? csmp::Index() : model.Database().StorageKey(saturated_hydraulic_conductivity))),
      ctKey_  ( ((compressibility==NULL) ? csmp::Index() : model.Database().StorageKey(compressibility))),
      phiKey_ ( ((porosity==NULL) ? csmp::Index() : model.Database().StorageKey(porosity))),
      sshdKey_( ((specific_saturated_hydraulic_diffusivity==NULL) ? csmp::Index() : model.Database().StorageKey(specific_saturated_hydraulic_diffusivity)))
{
    if ( sshcKey_.place != ELEMENT )
        throw csmp::Exception( ERROR, "ConductivityVisitor", specific_saturated_hydraulic_conductivity, " must be an element property." );
    if ( kKey_.place != ELEMENT )
        throw csmp::Exception( ERROR, "ConductivityVisitor", permeability, " must be an element property." );
    if ( muKey_.place != NODE )
        throw csmp::Exception( ERROR, "ConductivityVisitor", viscosity, " must be a node property." );
    if ( rhoKey_ != csmp::Index() && rhoKey_.place != ELEMENT )
        throw csmp::Exception( ERROR, "ConductivityVisitor", density, " must be an element property." );
    if ( shcKey_ != csmp::Index() && shcKey_.place != ELEMENT )
        throw csmp::Exception( ERROR, "ConductivityVisitor", saturated_hydraulic_conductivity, " must be an element property." );
    if ( ctKey_ != csmp::Index()  && ctKey_.place != ELEMENT )
        throw csmp::Exception( ERROR, "ConductivityVisitor", compressibility, " must be an element property." );
    if ( phiKey_ != csmp::Index() && phiKey_.place != ELEMENT )
        throw csmp::Exception( ERROR, "ConductivityVisitor", porosity, " must be an element property." );
    if ( sshdKey_ != csmp::Index() && sshdKey_.place != ELEMENT )
        throw csmp::Exception( ERROR, "ConductivityVisitor", specific_saturated_hydraulic_diffusivity, " must be an element property." );

#if defined(_OPENMP )
    thread_result_.resize(omp_get_max_threads());
    femgrs_.resize(omp_get_max_threads());
#pragma omp parallel
    {
        size_t tid = omp_get_thread_num();
        femgrs_[tid].InitializeElements(dim,model.FE_Manager().InterpolationOrder(),true);
    }
#endif
}

template<size_t dim>
void ConductivityVisitor<dim>::Visit( Model<dim>* m ){
#if defined(_OPENMP )
    this->Visit(&m->Region("Model"));
#endif
}

template<size_t dim>
void ConductivityVisitor<dim>::Visit(Region<dim>* region ){
#if defined(_OPENMP )
#pragma omp parallel
    {
        size_t tid=omp_get_thread_num();
        Element<dim>* ep;
#pragma omp for
        for ( long int e= 0 ; e < region->Elements(); e++ ){
            ep = region->E(e);
            FiniteElement* fe_tmp=ep->FE();
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
void ConductivityVisitor<dim>::Visit( Element<dim>* element )
{
#if !defined(_OPENMP)
    // first we calculate sshc
    element->PropertyValueAtBaryCenter( muKey_, result );
    result = element->Read( kKey_ ) / result();
    element->Store( sshcKey_, result );

    // if key exists, calculate sshd
    if (sshdKey_!=csmp::Index()){
        resultd = result;
        resultd /= element->Read(phiKey_);
        resultd /= element->Read(ctKey_);
        element->Store( sshdKey_, resultd);
    }

    // if key exists, calculate shc
    if (shcKey_!=csmp::Index()){
        result=result* element->Read(rhoKey_);
        element->Store( shcKey_, result );
    }
#endif
}

template<size_t dim>
void ConductivityVisitor<dim>::ComputeContribution( Element<dim>* element )
{
#if defined(_OPENMP )
    // first we calculate sshc
    element->PropertyValueAtBaryCenter( muKey_, thread_result_[omp_get_thread_num()] );
    thread_result_[omp_get_thread_num()] = element->Read( kKey_ ) / thread_result_[omp_get_thread_num()]();
    element->Store( sshcKey_, thread_result_[omp_get_thread_num()] );

    // if key exists, calculate sshd
    if (sshdKey_!=csmp::Index()){
        result = thread_result_[omp_get_thread_num()];
        result /= element->Read(phiKey_);
        result /= element->Read(ctKey_);
        element->Store( sshdKey_, result );
    }

    // if key exists, calculate shc
    if (shcKey_!=csmp::Index()){
        thread_result_[omp_get_thread_num()]=thread_result_[omp_get_thread_num()]* element->Read(rhoKey_);
        element->Store( shcKey_, thread_result_[omp_get_thread_num()] );
    }
#endif
}


template class ConductivityVisitor<1>;
template class ConductivityVisitor<2>;
template class ConductivityVisitor<3>;

}

