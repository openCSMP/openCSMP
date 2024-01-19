#include "ComputeSinglePhaseGravityTermVisitor.h"
#include "Model.h"
#include "Region.h"
#include "TwoPhaseModel.h"
#include "ScalarVariable.h"
#include "CSMP_physical_constants.h"

using namespace std;

namespace csmp {

template<uint32_t dim, template<uint32_t> class CELL>
ComputeSinglePhaseGravityTermVisitor<dim,CELL>::ComputeSinglePhaseGravityTermVisitor(Model<dim>& model,
                                                                                const char* permeabilityTag,
                                                                                const char* viscosityTag,
                                                                                const char* densityTag ,
                                                                                const char* gravityVectorTag,
                                                                                const char* model_gravity_vector,
                                                                                const char* densityTag2) // this is meant to post multiply the gravity term if provided.
    : Visitor<dim>( MODEL, ELEMENT ),
      gravityVectorKey_ (model.Database().StorageKey( gravityVectorTag ) ),
      permeabilityKey_( model.Database().StorageKey( permeabilityTag ) ),
      viscosityKey_ ( model.Database().StorageKey( viscosityTag ) ),
      densityKey_( model.Database().StorageKey( densityTag ) ),
      densityKey2_((densityTag2==NULL) ? csmp::Index() : model.Database().StorageKey(densityTag2)),
      gravitational_acceleration_(ACC_GRAVITY)
{
    model.Read(model.Database().StorageKey(model_gravity_vector),gravity_unit_vector_);
    
    // adjusting application target if necessary
    if constexpr ( TypeMatchesVariablePlacement<CELL,FACE>::value )
      this->ApplicationTarget( FACE );
    else if constexpr ( TypeMatchesVariablePlacement<CELL,INTER_FACE>::value )
      this->ApplicationTarget( INTER_FACE );
}


template<uint32_t dim, template<uint32_t> class CELL>
void ComputeSinglePhaseGravityTermVisitor<dim,CELL>::Visit( CELL<dim>* cell )
{
    ScalarVariable mu, rho;
    double gravityTerm;
    VectorVariable<dim> gravityVector;
    
    if ( cell->IsLine() && dim != 2u )
    {
        //! projection of a gravity vector (e.g. (0,-1,0) ) onto a vector r(r[0],r[1],r[2])
        //! r is a vector from one node to the other of the line element.
        //! (where |r|=1)
        //! We then perform the projection of the gravity vector on the line element

        Point<dim> line_vector( cell->N(0)->Coordinate() - cell->N(1)->Coordinate() );
        line_vector.NormalizeLengthTo(1.);
        VectorVariable<dim> projection = gravity_unit_vector_.ProjectOnto(line_vector.Coordinates());

        for ( uint32_t j{0U} ; j <dim; j++)
            gravityVector(j) = projection[j];

    }
    else if( cell->IsSurface() && dim == 3u )
    {
        //! projection of vector g(0,-1, 0) onto a surface with a normal vector n(n[0], n[1], n[2])
        //! (where |n|=1) is equal to
        //! n x (n x g) = (n, g) n - g = -n[1] n - g
        //! This code obviously assumes all nodes of the element are co-planar
        Point<dim> line_vector1( cell->N(0)->Coordinate() - cell->N(1)->Coordinate() );
        Point<dim> line_vector2( cell->N(1)->Coordinate() - cell->N(2)->Coordinate() );

        Point<dim> surface_normal( crossProduct( line_vector1,line_vector2) );
        surface_normal.NormalizeLengthTo (1.);

        Point<dim> firstcross(crossProduct(surface_normal,gravity_unit_vector_.P()));
        Point<dim> secondcross(crossProduct(surface_normal,firstcross));

        for (uint32_t j = 0U; j <dim; j++)
            gravityVector(j) = secondcross[j];
    }
    else
    {
        for (uint32_t j = 0U; j <dim; j++)
            gravityVector(j) = gravity_unit_vector_(j);
    }
    //-----------
    // Calculation of:
    // k/mu * rho * g
    // if densityTag2 is supplied, then that variable will be multiplied (e.g. rho*k/mu * rho * g )

    gravityTerm =  cell->Read( permeabilityKey_ ) * gravitational_acceleration_;
    // loop over element integration points
    if (gravityVectorKey_.place==ELEMENT_INTEGRATION_POINT){
        for (uint32_t ip=0U;ip<cell->IntegrationPoints (); ++ip)
        {
            cell->PropertyValueAtIntegrationPoint( densityKey_, ip, rho );
            cell->PropertyValueAtIntegrationPoint( viscosityKey_, ip, mu );
            cell->Store( ip, gravityVectorKey_, gravityVector*rho()*gravityTerm/mu() );
        }
        if (densityKey2_!=csmp::Index())
            for (uint32_t ip=0U;ip<cell->IntegrationPoints (); ++ip)
            {
                cell->Read( ip, gravityVectorKey_, gravityVector);
                cell->PropertyValueAtIntegrationPoint( densityKey2_, ip, rho );
                cell->Store( ip, gravityVectorKey_, gravityVector*rho());
            }
    }
    else if (gravityVectorKey_.place == ELEMENT){
        cell->PropertyValueAtBaryCenter( densityKey_, rho );
        cell->PropertyValueAtBaryCenter( viscosityKey_, mu );
        cell->Store( gravityVectorKey_, gravityVector*rho()*gravityTerm/mu() );

        if (densityKey2_!=csmp::Index()){
            cell->Read( gravityVectorKey_, gravityVector);
            cell->PropertyValueAtBaryCenter( densityKey2_, rho );
            cell->Store( gravityVectorKey_, gravityVector*rho());
        }
    }
    else {
        throw csmp::Exception(FATAL_ERROR,"ComputeSinglePhaseGravityTermVisitor<dim>::ComputeContribution",
                              " Placement of the gravity term vector is not supported",parsePlacement(gravityVectorKey_.place).c_str());
    }
}


template class ComputeSinglePhaseGravityTermVisitor<3,Element>;
template class ComputeSinglePhaseGravityTermVisitor<2,Element>;
template class ComputeSinglePhaseGravityTermVisitor<1,Element>;

template class ComputeSinglePhaseGravityTermVisitor<3,Face>;
template class ComputeSinglePhaseGravityTermVisitor<2,Face>;
template class ComputeSinglePhaseGravityTermVisitor<1,Face>;

} //csmp
