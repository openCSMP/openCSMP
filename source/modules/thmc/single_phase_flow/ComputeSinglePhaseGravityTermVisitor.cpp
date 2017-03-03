#include "ComputeSinglePhaseGravityTermVisitor.h"
#include "Model.h"
#include "TwoPhaseModel.h"
#include "ScalarVariable.h"
#include "CSMP_physical_constants.h"
#if defined(_OPENMP )
#include "omp.h"
#endif

using namespace std;

namespace csmp{

template<size_t dim>
ComputeSinglePhaseGravityTermVisitor<dim>::ComputeSinglePhaseGravityTermVisitor(Model<dim>& model,
                                                                                const char* permeabilityTag,
                                                                                const char* viscosityTag,
                                                                                const char* densityTag ,
                                                                                const char* gravityVectorTag,
                                                                                const char* model_gravity_vector,
                                                                                const char* densityTag2) // this is meant to post multiply the gravity term if provided.
    : Visitor<dim>( MODEL, ELEMENT ), model_( model ),
      gravityVectorKey_ (model.Database().StorageKey( gravityVectorTag ) ),
      permeabilityKey_( model.Database().StorageKey( permeabilityTag ) ),
      viscosityKey_ ( model.Database().StorageKey( viscosityTag ) ),
      densityKey_( model.Database().StorageKey( densityTag ) ),
      densityKey2_((densityTag2==NULL) ? csmp::Index() : model.Database().StorageKey(densityTag2)),
      gravitational_acceleration_(ACC_GRAVITY)
{
    model.Read(model.Database().StorageKey(model_gravity_vector),gravity_unit_vector_);

#if defined(_OPENMP )
    this->femgrs_.resize(omp_get_max_threads());
    for (long int tid = 0 ; tid < omp_get_max_threads();tid ++)
        this->femgrs_[tid].InitializeElements(dim,model.FE_Manager().InterpolationOrder(),true);
#endif
}

template<size_t dim>
void ComputeSinglePhaseGravityTermVisitor<dim>::Visit( Model<dim>* m ){
    if (this->Verbose()) cout <<" ComputeSinglePhaseGravityTermVisitor<dim>::Visit(Model<dim>*)"<<endl;
#if defined(_OPENMP )
    this->Visit(&(m->Region("Model")));
#endif
}

template<size_t dim>
void ComputeSinglePhaseGravityTermVisitor<dim>::Visit(Region<dim>* region ){
    if (this->Verbose()) cout <<" ComputeSinglePhaseGravityTermVisitor<dim>::Visit(Region<dim>*) : "<<region->Name()<<endl;

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
void ComputeSinglePhaseGravityTermVisitor<dim>::Visit( Element<dim>* element )
{
#if !defined(_OPENMP)
    this->ComputeContribution(element);
#endif
}

template<size_t dim>
void ComputeSinglePhaseGravityTermVisitor<dim>::ComputeContribution( Element<dim>* element )
{

    ScalarVariable mu, rho;
    double gravityTerm;
    VectorVariable<dim> gravityVector;
    
    if( element->FE()->IsLineElement() && (dim == 2 || dim ==3))
    {
        //! projection of a gravity vector (e.g. (0,-1,0) ) onto a vector r(r[0],r[1],r[2])
        //! r is a vector from one node to the other of the line element.
        //! (where |r|=1)
        //! We then perform the projection of the gravity vector on the line element

        Point<dim> line_vector( element->N(0)->Coordinate() - element->N(1)->Coordinate() );
        line_vector.NormalizeLengthTo(1.);
        VectorVariable<dim> projection = gravity_unit_vector_.ProjectOnto(line_vector.Coordinates());

        for (size_t j = 0 ; j <dim; j++)
            gravityVector(j) = projection[j];

    }
    else if( element->FE()->IsSurfaceElement() && dim ==3)
    {
        //! projection of vector g(0,-1, 0) onto a surface with a normal vector n(n[0], n[1], n[2])
        //! (where |n|=1) is equal to
        //! n x (n x g) = (n, g) n - g = -n[1] n - g
        //! This code obviously assumes all nodes of the element are co-planar
        Point<dim> line_vector1( element->N(0)->Coordinate() - element->N(1)->Coordinate() );
        Point<dim> line_vector2( element->N(1)->Coordinate() - element->N(2)->Coordinate() );

        Point<dim> surface_normal( crossProduct( line_vector1,line_vector2) );
        surface_normal.NormalizeLengthTo (1.);

        Point<dim> firstcross(crossProduct(surface_normal,gravity_unit_vector_.P()));
        Point<dim> secondcross(crossProduct(surface_normal,firstcross));

        for (size_t j = 0 ; j <dim; j++)
            gravityVector(j) = secondcross[j];
    }
    else
    {
        for (size_t j = 0 ; j <dim; j++)
            gravityVector(j) = gravity_unit_vector_(j);
    }
    //-----------
    // Calculation of:
    // k/mu * rho * g
    // if densityTag2 is supplied, then that variable will be multiplied (e.g. rho*k/mu * rho * g )

    gravityTerm =  element->Read( permeabilityKey_ ) * gravitational_acceleration_;
    // loop over element integration points
    if (gravityVectorKey_.place==ELEMENT_INTEGRATION_POINT){
        for (size_t ip=0;ip<element->IntegrationPoints (); ++ip)
        {
            element->PropertyValueAtIntegrationPoint( densityKey_, ip, rho );
            element->PropertyValueAtIntegrationPoint( viscosityKey_, ip, mu );
            element->Store( ip, gravityVectorKey_, gravityVector*rho()*gravityTerm/mu() );
        }
        if (densityKey2_!=csmp::Index())
            for (size_t ip=0;ip<element->IntegrationPoints (); ++ip)
            {
                element->Read( ip, gravityVectorKey_, gravityVector);
                element->PropertyValueAtIntegrationPoint( densityKey2_, ip, rho );
                element->Store( ip, gravityVectorKey_, gravityVector*rho());
            }
    }
    else if (gravityVectorKey_.place == ELEMENT){
        element->PropertyValueAtBaryCenter( densityKey_, rho );
        element->PropertyValueAtBaryCenter( viscosityKey_, mu );
        element->Store( gravityVectorKey_, gravityVector*rho()*gravityTerm/mu() );

        if (densityKey2_!=csmp::Index()){
            element->Read( gravityVectorKey_, gravityVector);
            element->PropertyValueAtBaryCenter( densityKey2_, rho );
            element->Store( gravityVectorKey_, gravityVector*rho());
        }
    }
    else {
        throw csmp::Exception(CSMP_FATAL_ERROR,"ComputeSinglePhaseGravityTermVisitor<dim>::ComputeContribution",
                              " Placement of the gravity term vector is not supported",parsePlacement(gravityVectorKey_.place).c_str());
    }
//    cin.get();
}

template class ComputeSinglePhaseGravityTermVisitor<3>;
template class ComputeSinglePhaseGravityTermVisitor<2>;
template class ComputeSinglePhaseGravityTermVisitor<1>;

} //csmp
