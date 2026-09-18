// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV.h"
#include "UpwindControlVisitor.h"

using namespace std;

namespace csmp {

/** default destructor */
template<uint32_t dim, template<uint32_t> class CELL>
CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV<dim,CELL>::~CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV() {}


/** custom constructor */
template<uint32_t dim, template<uint32_t> class CELL>
CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV<dim,CELL>::CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV( const PropertyDatabase<dim>& pref,
                                                                                            UpwindControlVisitor<dim>& upwind_visitor,
                                                                                            ExplicitFiniteVolumeTransportPHX<dim>& fv_transport,
                                                                                            const char* oper,
                                                                                            const char* test,
                                                                                            const char* upwind,
                                                                                            const char* grav_trigger,
                                                                                            const char* thickness)//Benoit 2025 add
    : CVFEM_MathOperatorRHS<dim,CELL>(pref,oper,test),
      UpwindVisitor(upwind_visitor),
      finite_volume(fv_transport),
      B(3,3),
      gravity(-9.80665), // scalar acts to increase the pressure
      xyz(dim-1),
      upwind_(pref.StorageKey(upwind)),
      rho_key(pref.StorageKey(grav_trigger)),
      thickness_(pref.StorageKey(thickness))//Benoit 2025 add
{
    if (dim==3) xyz = 1;

    string name = "CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV ";
    name += upwind;
    char * cname = new char[name.length()+1];
    strcpy(cname, name.c_str());
    MathOperatorRHS<dim,CELL>::Name(cname, oper, test);
    
    // testing the Operands
    if ( MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() != ELEMENT ||
         MathOperatorRHS<dim,CELL>::MaterialOperandType() != SCALAR )
        throw csmp::Exception( ERROR, "Upwind_Integral_dNT_rhsop_dN_dV_2<dim>::(constructor)",
                               oper, "Operand must be a scalar property placed on the element." );

    if ( MathOperatorRHS<dim,CELL>::TestOperandPlacement() != NODE ||
         MathOperatorRHS<dim,CELL>::TestOperandType() != SCALAR )
        throw csmp::Exception( ERROR, "Upwind_Integral_dNT_rhsop_dN_dV_2<dim>::(constructor)",
                               test, "Dependent variable must be a scalar property placed on the nodes." );

    if ( upwind_.Placement() != NODE || upwind_.Type() != SCALAR )
        throw csmp::Exception( ERROR, "Upwind_Integral_dNT_rhsop_dN_dV_2<dim>::(constructor)",
                               upwind, "Upwind variable must be a scalar property placed on the nodes." );

    if ( rho_key.place != NODE || rho_key.type != SCALAR )
        throw csmp::Exception( ERROR, "Upwind_Integral_dNT_rhsop_dN_dV_2<dim>::(constructor)",
                               grav_trigger, "Upwind variable must be a scalar property placed on the nodes." );
}


template<uint32_t dim, template<uint32_t> class CELL>
void CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV<dim,CELL>::GetOperands( const CELL<dim>& e )
{
    MathOperatorRHS<dim,CELL>::GetOperands(e);

    e.NodePropertyVector(upwind_.Key(), upwind_var_);
}


template<uint32_t dim, template<uint32_t> class CELL>
void CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV<dim,CELL>::GetOperandsCVFEM( const CELL<dim>& e, csmp::Index upwind_var_key )
{
    MathOperatorRHS<dim,CELL>::GetOperands(e);

    e.NodePropertyVector(upwind_.Key(), upwind_var_);
    e.NodePropertyVector(upwind_var_key, upwind_var_multiplier);

    for (auto i{0U}; i < upwind_var_multiplier.size(); i++)
        upwind_var_[i]() *= upwind_var_multiplier[i]();
}


template<uint32_t dim, template<uint32_t> class CELL>
void CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
{

    GetUpwindMatrix( e );
    
    // initialize output matrix
    MathOperatorRHS<dim,CELL>::RHS.resize(e.Nodes());
    fill(MathOperatorRHS<dim,CELL>::RHS.begin(), MathOperatorRHS<dim,CELL>::RHS.end(), 0.0);

    //    e.dN( B );
    operand = MathOperatorRHS<dim,CELL>::MTRL[0](0,0);
    //    B *= operand;
    
    // loop over facets of the element
    for ( uint32_t i{0}; i<e.FV()->Facets(); i++ )
    {

        //e.dN_AtIntegrationPoint( B, i, 1 );
        e.dN(B);
        //       e.dN( B, e.ConnectedFiniteVolumeStencil()->FacetEdgeMidPoint(i ).Coordinates() );
        B *= operand;
        e.FV()->FacetEdgeNodes( i, inside_node_, outside_node_ );
        area = finite_volume.GetFacetArea( e.Idx(), i );
        area*=e.Read(thickness_);//Benoit 2025 add

        normal_component = finite_volume.GetFacetNormalComponent( e.Idx(), i, xyz );
        grav = gravity * normal_component;
        contribution = area*MathOperatorRHS<dim,CELL>::MTRL[0](0,0);
        if (upwind(inside_node_, outside_node_) == 1) contribution *= upwind_var_[inside_node_]();
        else if (upwind(inside_node_, outside_node_) == 2) contribution *= upwind_var_[outside_node_]();
        else contribution = 0.0;
        contribution *= grav;
        MathOperatorRHS<dim,CELL>::RHS[inside_node_]  -= contribution;
        MathOperatorRHS<dim,CELL>::RHS[outside_node_] += contribution;
    }
}    


template<uint32_t dim, template<uint32_t> class CELL>
void CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV<dim,CELL>::GetUpwindMatrix( const CELL<dim>& e )
{
    upwind = UpwindVisitor.UpwindMatrix( rho_key, e.Idx() );

} // GetUpwindMatrix


template class CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV<1,Element>;
template class CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV<2,Element>;
template class CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV<3,Element>;

} // csmp
