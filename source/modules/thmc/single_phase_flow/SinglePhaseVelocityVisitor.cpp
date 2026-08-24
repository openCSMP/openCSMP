#include "SinglePhaseVelocityVisitor.h"
#include "Model.h"
#include "Region.h"
#include "CSMP_mathUtilities.h"
#include "IsoparametricLinearTetrahedron.h"
#include "CSMP_physical_constants.h"
#include "Exception.h"
#include "ErrorHandler.h"

using namespace std;

namespace csmp {

template<uint32_t dim, template<uint32_t> class CELL>
SinglePhaseVelocityVisitor<dim,CELL>::SinglePhaseVelocityVisitor(Model<dim>& model,
                                                            const char* porosity,
                                                            const char* conductivity,
                                                            const char* fluid_density,
                                                            const char* fluid_pressure,
                                                            const char* velocity,                      // results:
                                                            const char* model_gravity_vector,          // input from model: example (0,-1.0,0) or (0,1.0,0)
                                                            const char* pore_velocity,                 // results: True velocity (v/phi)
                                                            const char* volume_flux,                   // results:
                                                            const char* nodal_velocity,                // results:
                                                            const char* nodal_pore_velocity,           // results:
                                                            const char* nodal_volume_flux):            //

    components_(dim*2+1),
    conductivity_key_  (model.Database().StorageKey(conductivity)),
    phi_key_           (model.Database().StorageKey(porosity)),
    rho_key_           (model.Database().StorageKey(fluid_density)),
    fpres_key_         (model.Database().StorageKey(fluid_pressure)),
    velo_key_          (model.Database().StorageKey(velocity)),
    ivelo_key_         ((pore_velocity==NULL) ? csmp::Index() : model.Database().StorageKey(pore_velocity)),
    flux_key_          ((volume_flux==NULL) ? csmp::Index() : model.Database().StorageKey(volume_flux)),
    nvelo_key_         ((nodal_velocity==NULL) ? csmp::Index() : model.Database().StorageKey(nodal_velocity)),
    nivelo_key_        ((nodal_pore_velocity==NULL) ? csmp::Index() : model.Database().StorageKey(nodal_pore_velocity)),
    nflux_key_         ((nodal_volume_flux==NULL) ? csmp::Index() : model.Database().StorageKey(nodal_volume_flux)),
    with_gravity_      (false),
    gravitational_acceleration_(ACC_GRAVITY),
    application_cycle_(1)
{
    if (model_gravity_vector == NULL) gravity_unit_vector_=0.;
    else {
        model.Read(model.Database().StorageKey(model_gravity_vector),gravity_unit_vector_);
        if (gravity_unit_vector_.Length()>0.0)
            with_gravity_=true;
     }

    this->ApplicationLevel(MODEL);
    this->ApplicationTarget(ELEMENT);

    // getting ranges for output variables

    if ( conductivity_key_.type != SCALAR || conductivity_key_.place != ELEMENT)
        ErrorHandler::Instance().Note( FATAL_ERROR, "SinglePhaseVelocityVisitor::(constructor)",
                             conductivity, " must be a scalar property placed on the nodes. (anything else needs implementation)" );

    if ( phi_key_.type != SCALAR )
        ErrorHandler::Instance().Note( FATAL_ERROR, "SinglePhaseVelocityVisitor::(constructor)",
                             porosity, " must be a scalar property." );
    if ( rho_key_.type != SCALAR )
        ErrorHandler::Instance().Note( FATAL_ERROR, "SinglePhaseVelocityVisitor::(constructor)",
                             fluid_density, " must be a scalar property on the element." );

    if ( fpres_key_.place != NODE || fpres_key_.type != SCALAR )
        ErrorHandler::Instance().Note( FATAL_ERROR, "SinglePhaseVelocityVisitor::(constructor)",
                             fluid_pressure, " must be a scalar property placed on the nodes." );

    if ( velo_key_.place != ELEMENT || velo_key_.type != VECTOR )
        ErrorHandler::Instance().Note( FATAL_ERROR, "SinglePhaseVelocityVisitor::(constructor)",
                             velocity, " must be a vector property placed on the element." );
    model.Database().RangeOf(velocity, minmaxV_.first, minmaxV_.second );

    if ( nvelo_key_!=csmp::Index()){
        if ( ivelo_key_.place != ELEMENT || ivelo_key_.type != VECTOR ){
            ErrorHandler::Instance().Note( FATAL_ERROR, "SinglePhaseVelocityVisitor::(constructor)",
                                 pore_velocity, " must be a vector property placed on the element." );
        }
        else{
            model.Database().RangeOf(pore_velocity, minmaxF_.first, minmaxF_.second );
        }
    }
    if ( flux_key_!=csmp::Index()){
        if ( flux_key_.place != ELEMENT || flux_key_.type != SCALAR ){
            ErrorHandler::Instance().Note( FATAL_ERROR, "SinglePhaseVelocityVisitor::(constructor)",
                                 volume_flux, " must be a scalar property placed on the element." );
        }
        else{
            model.Database().RangeOf(volume_flux, minmaxF_.first, minmaxF_.second );
        }
    }
    if ( nvelo_key_!=csmp::Index()){
        if ( nvelo_key_.place != NODE || nvelo_key_.type != VECTOR ){
            ErrorHandler::Instance().Note( FATAL_ERROR, "SinglePhaseVelocityVisitor::(constructor)",
                                 nodal_velocity, " must be a vector property placed on the node." );
        }
        else {
            model.Database().RangeOf(nodal_velocity, minmaxF_.first, minmaxF_.second );
        }
    }

    if ( nivelo_key_!=csmp::Index()){
        if ( nivelo_key_.place != NODE || ivelo_key_.type != VECTOR ){
            ErrorHandler::Instance().Note( FATAL_ERROR, "SinglePhaseVelocityVisitor::(constructor)",
                                 nodal_pore_velocity, " must be a vector property placed on the node." );
        }
        else{
            model.Database().RangeOf(nodal_pore_velocity, minmaxF_.first, minmaxF_.second );
        }
    }

    if ( nflux_key_!=csmp::Index()){
        if ( nflux_key_.place != NODE || nflux_key_.type != SCALAR ){
            ErrorHandler::Instance().Note( FATAL_ERROR, "SinglePhaseVelocityVisitor::(constructor)",
                                 nodal_volume_flux, " must be a scalar property placed on the node." );
        }
        else{
            model.Database().RangeOf(nodal_volume_flux, minmaxF_.first, minmaxF_.second );
        }
    }

    const Region<dim>& domain = model.Region("Model");

    // Check if there are LDE's present
    bool volume_elements_present(false);
    for ( auto e = domain.CellsBegin(); e != domain.CellsEnd(); e++ )
      if ( (*e)->IsVolume() )
        volume_elements_present=true;

    for ( auto e = domain.CellsBegin() ; e != domain.CellsEnd(); e++ ){

        if ((dim == 2 && (*e)->IsLine()) || (dim == 3 && (*e)->IsSurface() && volume_elements_present))
        {
            ErrorHandler::Instance().Note(FATAL_ERROR,"SinglePhaseVelocityVisitor, constructor","You have lower dimensional elements in your model",
                                " Either implement code to support them in this visitor or use another that does.");
        }
    }

} // end constructor








template<uint32_t dim, template<uint32_t> class CELL>
void SinglePhaseVelocityVisitor<dim,CELL>::Visit( CELL<dim>* e )
{
    std::vector<ScalarVariable >    PF, RHO;
    std::vector<DenseMatrix<dim> >  MTRL(1);// fluid pressure
    ScalarVariable                  phi, flux, rho;
    VectorVariable<dim>             velo,ivelo;
    DenseMatrix<DM_MIN>             DERIV;
    e->NodePropertyVector(fpres_key_, PF );

    // fluid density if needed
    if ( with_gravity_ ) {

        if ( rho_key_.place == ELEMENT )
            e->Read( rho_key_, rho );
        else if ( rho_key_.place == NODE ) {
            if ( e->Interpolation() == 1U )
                e->PropertyValueAtBaryCenter( rho_key_, rho );
            else
                e->NodePropertyVector( rho_key_, RHO );
        }
        else
            ErrorHandler::Instance().Note( FATAL_ERROR, "SinglePhaseVelocityVisitor<dim>::GetOperands",
                                 "fluid density is neither a node nor element variable; can't deal with this.");
    }

    // conductivity
    if ( conductivity_key_.place == ELEMENT )
    {
        if ( conductivity_key_.type == SCALAR )
            MTRL[0].AssignToDiagonalAndZeroOffDiagonal( dim, e->Read( conductivity_key_ ) );
        else if ( conductivity_key_.type == VECTOR ) {
            VectorVariable<dim>  vc;
            e->Read( conductivity_key_, vc );
            MTRL[0].AssignToDiagonal( vc );
        }
        else if ( conductivity_key_.type == TENSOR ) {
            TensorVariable<dim>  ts;
            e->Read( conductivity_key_, ts );
            MTRL[0] = ts;
        }
    }
    else if ( conductivity_key_.place == ELEMENT_INTEGRATION_POINT )
    {
        for ( uint32_t i{0U}; i<e->IntegrationPoints(); i++ )
        {
            if ( conductivity_key_.type == SCALAR )
                MTRL[i].AssignToDiagonalAndZeroOffDiagonal( dim, e->Read(i, conductivity_key_));
            else if ( conductivity_key_.type == VECTOR ) {
                VectorVariable<dim>  vc;
                e->Read(i, conductivity_key_, vc );
                MTRL[i].AssignToDiagonal( vc );
            }
            else if ( conductivity_key_.type == TENSOR ) {
                TensorVariable<dim>  ts;
                e->Read(i, conductivity_key_, ts );
                MTRL[i] = ts;
            }
        }
    }
    else if ( conductivity_key_.place == NODE )
    {
        for ( uint32_t i{0U}; i<e->IntegrationPoints(); i++ )
        {
            if ( conductivity_key_.type == SCALAR )
                MTRL[i].AssignToDiagonalAndZeroOffDiagonal( dim, e->PropertyValueAtIntegrationPoint( conductivity_key_, i ) );
            else if ( conductivity_key_.type == VECTOR ) {
                VectorVariable<dim>  vc;
                e->PropertyValueAtIntegrationPoint( conductivity_key_, i, vc );
                MTRL[i].AssignToDiagonal( vc );
            }
            else if ( conductivity_key_.type == TENSOR ) {
                TensorVariable<dim>  ts;
                e->PropertyValueAtIntegrationPoint( conductivity_key_, i, ts );
                MTRL[i] = ts;
            }
        }
    }

    // porosity
    e->Read( phi_key_, phi );

    if ( e->Interpolation() == 1U && velo_key_.place == ELEMENT){
        // ---------------------------------------------
        // 1. Case of linear elements: Constant velocity
        // ---------------------------------------------
        // vel = -(k/mu)*(grad p - rho*g*{g_unit})
        velo = 0.;
        for ( uint32_t i{0U}; i<e->Nodes(); i++ )
            for ( uint32_t j{0U}; j<dim; j++ )
                velo(j) += PF[i]() * DERIV(j,i) ;

        if (with_gravity_) {
            // 1. Case of linear elements: Constant velocity
            // ---------------------------------------------

            if ( e->Interpolation() == 1U && velo_key_.place == ELEMENT){
                for (uint32_t i{0U} ; i < dim ; i++)
                    velo(i) -= gravity_unit_vector_(i) *gravitational_acceleration_* rho();

            }
            else{
                ErrorHandler::Instance().Note(FATAL_ERROR,"SinglePhaseVelocityVisitor<dim>::AddGravitationalComponentToVelocity( Element<dim>)",
                                    " Interpolation for this element is higher than 1, or velocity is not placed on the element",
                                    " Whatever functionality you are attempting to use probably needs implementation.");
            }
        }

        for ( uint32_t j{0U}; j<dim; j++ )
            velo(j) *= -MTRL[0](j,j);

        if (ivelo_key_!=csmp::Index()) {
            ivelo  = velo;
            ivelo /= phi;
            e->Store( ivelo_key_, ivelo );
        }
        if (flux_key_!=csmp::Index()) {
            flux = velo.Length();
            e->Store( flux_key_, flux );
        }

        e->Store( velo_key_, velo );

    }
    else{
        ErrorHandler::Instance().Note(FATAL_ERROR,"SinglePhaseVelocityVisitor<dim>::ComputeVelocity( Element<dim>)",
                            " Interpolation for this element is higher than 1, or velocity is not placed on the element",
                            " Whatever functionality you are attempting to use probably needs implementation.");
    }
}


template class SinglePhaseVelocityVisitor<1U,Element>;
template class SinglePhaseVelocityVisitor<2U,Element>;
template class SinglePhaseVelocityVisitor<3U,Element>;

template class SinglePhaseVelocityVisitor<1U,Face>;
template class SinglePhaseVelocityVisitor<2U,Face>;
template class SinglePhaseVelocityVisitor<3U,Face>;

} // csmp






















