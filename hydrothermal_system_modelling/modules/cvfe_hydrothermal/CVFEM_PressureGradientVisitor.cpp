// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "CVFEM_PressureGradientVisitor.h"
#include "Exception.h"

using namespace std;

namespace csmp {

/** custom constructor */
template<uint32_t dim>
CVFEM_PressureGradientVisitor<dim>::CVFEM_PressureGradientVisitor( Model<dim>& model,
                                                                  const char* fluid_pressure,
                                                                  const char* lithostatic_pressure,
                                                                  const char* permeability,
                                                                  const char* KgradP,
                                                                  const char* gradP_scaling)
    : gradient_scaling(false)
{

    this->ApplicationLevel(REGION);
    this->ApplicationTarget(ELEMENT);

    k_key_          = model.Database().StorageKey(permeability);
    pres_key_       = model.Database().StorageKey(fluid_pressure);
    lith_pres_key_  = model.Database().StorageKey(lithostatic_pressure);
    KgradP_key_     = model.Database().StorageKey(KgradP);
    scale_key_      = model.Database().StorageKey(gradP_scaling);

    if ((k_key_.type != TENSOR && k_key_.type != SCALAR) || k_key_.place != ELEMENT)
        throw csmp::Exception( ERROR, "CVFEM_PressureGradientVisitor::(constructor)",
                              permeability, " must be a tensor or scalar property placed on the elements." );

    if ( pres_key_.type != SCALAR || pres_key_.place != NODE )
        throw csmp::Exception( ERROR, "CVFEM_PressureGradientVisitor::(constructor)",
                              fluid_pressure, " must be a scalar property placed on the nodes." );

    if ( lith_pres_key_.type != SCALAR || lith_pres_key_.place != NODE )
        throw csmp::Exception( ERROR, "CVFEM_PressureGradientVisitor::(constructor)",
                              lithostatic_pressure, " must be a scalar property placed on the nodes." );

    if ( KgradP_key_.type != VECTOR || KgradP_key_.place != ELEMENT )
        throw csmp::Exception( ERROR, "CVFEM_PressureGradientVisitor::(constructor)",
                              KgradP, " must be a vector property placed on the element." );

    if ( scale_key_.type != SCALAR || scale_key_.place != ELEMENT )
        throw csmp::Exception( ERROR, "CVFEM_PressureGradientVisitor::(constructor)",
                              gradP_scaling, " must be a scalar property placed on the element." );

}


/** default destructor */
template<uint32_t dim>
CVFEM_PressureGradientVisitor<dim>::~CVFEM_PressureGradientVisitor()
{}

/** visit function for Model */
template<uint32_t dim>
void CVFEM_PressureGradientVisitor<dim>::Visit(Model<dim>* n)
{}

/** visit function for Region */
template<uint32_t dim>
void CVFEM_PressureGradientVisitor<dim>::Visit(Region<dim>* n)
{}

/** visit function for Element */
template<uint32_t dim>
void CVFEM_PressureGradientVisitor<dim>::Visit(Element<dim>* n)
{

    n->NodePropertyVector( pres_key_, p_ );
    n->NodePropertyVector( lith_pres_key_, lp_ );

    if ( k_key_.type == SCALAR )
        n->Read( k_key_, k_ );
    else
        n->Read( k_key_, KT_ );

    if (gradient_scaling) ComputeGradient( *n );
    else ComputeGradient2( *n );

    n->Store( KgradP_key_, KgradP_ );
    n->Store( scale_key_, gradP_factor_ );

}

/** compute gradient with gradient scaling */
template<uint32_t dim>
void CVFEM_PressureGradientVisitor<dim>::ComputeGradient( Element<dim>& e )
{
    //current implementation sometimes lead to gradP_factor_= NaN, not sure if it matters
    DERIV_.Resize(dim,e.Nodes());
    KgradP_ = 0.;

    for ( uint32_t i{0}; i<e.IntegrationPoints(); i++ )
    {
        e.dN_AtIntegrationPoint( DERIV_, i, 1 );//Benoit: this is slow in OPEN-CSMP

        for ( uint32_t in{0u}; in<e.Nodes(); in++ )
        {
            for ( uint32_t j{0U}; j<dim; j++ )
            {
                KgradP_(j) += lp_[in]() * -DERIV_(j,in);
            }
        }
    }

    KgradP_  /= static_cast<double>(e.IntegrationPoints());
    gradP_factor_() = 1./KgradP_(dim-1);

    KgradP_ = 0.;
    for ( uint32_t i{0}; i<e.IntegrationPoints(); i++ )
    {
        e.dN_AtIntegrationPoint( DERIV_, i, 1 );//Benoit: this is slow in OPEN-CSMP

        for ( uint32_t in{0u}; in<e.Nodes(); in++ )
        {
            for ( uint32_t j{0U}; j<dim; j++ )
            {
                KgradP_(j) += p_[in]() * -DERIV_(j,in);
            }
        }
    }

    KgradP_ /= static_cast<double>(e.IntegrationPoints());
    gradP_factor_() *= KgradP_(dim-1);

    //------------------
    // Multiply by permeability

    if ( k_key_.type == SCALAR )
    {
        KgradP_ *= k_(); //multiply by scalar permeability
    }

    else
        for ( uint32_t j=0; j<dim; j++ )
        {
            KgradP_(j) *= KT_(j, j);//multiply by diagonal tensor permeability (term by term)
        }

} // end ComputeContribution


/** compute gradient without gradient scaling */
template<uint32_t dim>
void CVFEM_PressureGradientVisitor<dim>::ComputeGradient2( Element<dim>& e )
{

    DERIV_.Resize(dim,e.Nodes());

    gradP_factor_() = 0.;
    for ( uint32_t i{0U}; i<e.Nodes(); i++ )
    {
        if (p_[i]()>lp_[i]())
            gradP_factor_() += p_[i]()/lp_[i]();
        else
            gradP_factor_() += 1.;
    }
    gradP_factor_() /= static_cast<double>(e.Nodes());

    if (definitelyLessThan( gradP_factor_(), 1.0, numeric_limits<double>::epsilon()))
        gradP_factor_() = 1.0;

    KgradP_ = 0.;
    for ( uint32_t i{0}; i<e.IntegrationPoints(); i++ )
    {
        e.dN_AtIntegrationPoint( DERIV_, i, 1 );

        for ( uint32_t in=0; in<e.Nodes(); in++ )
            for ( uint32_t j{0U}; j<dim; j++ )
                KgradP_(j) += p_[in]() * -DERIV_(j,in);
    }

    KgradP_  /= static_cast<double>(e.IntegrationPoints());

    if ( k_key_.type == SCALAR )
        KgradP_  *= k_(); //multiply by scalar permeability

    else
        for ( uint32_t j=0; j<dim; j++ )
        {
            KgradP_(j) *=KT_(j,j);//multiply by diagonal tensor permeability (term by term)
        }

} // end ComputeContribution2


template class CVFEM_PressureGradientVisitor<1U>;
template class CVFEM_PressureGradientVisitor<2U>;
template class CVFEM_PressureGradientVisitor<3U>;

} // csmp
