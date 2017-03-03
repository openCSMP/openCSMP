#include "CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV.h"

using namespace std;

namespace csmp {

/** default destructor */
template<size_t dim>
CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV<dim>::~CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV() {}

/** custom constructor */
template<size_t dim>
CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV<dim>::CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV(const PropertyDatabase<dim>& pref, 
                                                                       UpwindControlVisitor<dim>& upwind_visitor,
                                                                       ExplicitFiniteVolumeTransportPHX<dim>& fv_transport,
                                                                       const char* oper,
                                                                       const char* test,
                                                                       const char* upwind,
                                                                       const char* grav_trigger)
  : CVFEM_MathOperatorRHS<dim>(pref,oper,test),
    UpwindVisitor(upwind_visitor),
    finite_volume(fv_transport),
    B(3,3),
    gravity(-9.80665), // scalar acts to increase the pressure
    xyz(dim-1),
    upwind_(pref.StorageKey(upwind)),
    rho_key( pref.StorageKey(grav_trigger) )
 {
    if (dim==3) xyz = 1;

    string name = "CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV ";
    name += upwind;
    char * cname = new char[name.length()+1];
    strcpy(cname, name.c_str());
    MathOperatorRHS<dim>::Name(cname, oper, test);
    
	// testing the Operands 
    if ( MathOperatorRHS<dim>::MaterialOperandPlacement() != ELEMENT || 
         MathOperatorRHS<dim>::MaterialOperandType() != SCALAR )
    throw csmp::Exception( CSMP_ERROR, "Upwind_Integral_dNT_rhsop_dN_dV_2<dim>::(constructor)", 
                   oper, "Operand must be a scalar property placed on the element." );

    if ( MathOperatorRHS<dim>::TestOperandPlacement() != NODE || 
         MathOperatorRHS<dim>::TestOperandType() != SCALAR )
    throw csmp::Exception( CSMP_ERROR, "Upwind_Integral_dNT_rhsop_dN_dV_2<dim>::(constructor)", 
                   test, "Dependent variable must be a scalar property placed on the nodes." );
     
    if ( upwind_.Placement() != NODE || upwind_.Type() != SCALAR )
    throw csmp::Exception( CSMP_ERROR, "Upwind_Integral_dNT_rhsop_dN_dV_2<dim>::(constructor)", 
                   upwind, "Upwind variable must be a scalar property placed on the nodes." );

    if ( rho_key.place != NODE || rho_key.type != SCALAR )
    throw csmp::Exception( CSMP_ERROR, "Upwind_Integral_dNT_rhsop_dN_dV_2<dim>::(constructor)", 
                   grav_trigger, "Upwind variable must be a scalar property placed on the nodes." );
                   
 }

template<size_t dim>
void CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV<dim>::GetOperands( Element<dim>& e )
{

  MathOperatorRHS<dim>::GetOperands(e);

  e.NodePropertyVector(upwind_.Key(), upwind_var_);
  
}

template<size_t dim>
void CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV<dim>::GetOperandsCVFEM( Element<dim>& e, csmp::Index upwind_var_key )
{

  MathOperatorRHS<dim>::GetOperands(e);
   
  e.NodePropertyVector(upwind_.Key(), upwind_var_);
  e.NodePropertyVector(upwind_var_key, upwind_var_multiplier);

  for (size_t i = 0; i < upwind_var_multiplier.size(); i++)
     upwind_var_[i]() *= upwind_var_multiplier[i]();
  
}


template<size_t dim>
void CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV<dim>::ComputeContribution( Element<dim>& e )
{
 
    GetUpwindMatrix( e );
    
    // initialize output matrix
    MathOperatorRHS<dim>::RHS.resize(e.Nodes());
    fill(MathOperatorRHS<dim>::RHS.begin(), MathOperatorRHS<dim>::RHS.end(), 0.0);

//    e.dN( B );
    operand = MathOperatorRHS<dim>::MTRL[0](0,0);
//    B *= operand;
    
    // loop over facets of the element
    for ( size_t i=0U; i<e.FV()->Facets(); i++ )
       {

       e.dN_AtIntegrationPoint( B, i, 1 );
//       e.dN( B, e.ConnectedFiniteVolumeStencil()->FacetEdgeMidPoint(i ).Coordinates() );
    B *= operand;
          e.FV()->FacetEdgeNodes( i, inside_node_, outside_node_ );
          area = finite_volume.GetFacetArea( e.Idx(), i );
          
          normal_component = finite_volume.GetFacetNormalComponent( e.Idx(), i, xyz );
          grav = gravity * normal_component;
          contribution = area*MathOperatorRHS<dim>::MTRL[0](0,0);
    	  if (upwind(inside_node_, outside_node_) == 1) contribution *= upwind_var_[inside_node_]();
    	  else if (upwind(inside_node_, outside_node_) == 2) contribution *= upwind_var_[outside_node_]();
    	  else contribution = 0.0;
    	  contribution *= grav;
          MathOperatorRHS<dim>::RHS[inside_node_]  -= contribution;
          MathOperatorRHS<dim>::RHS[outside_node_] += contribution;
       }
}    

template<size_t dim>
void CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV<dim>::GetUpwindMatrix( Element<dim>& e )
{

  upwind = UpwindVisitor.UpwindMatrix( rho_key, e.Idx() );    

} // GetUpwindMatrix


template class CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV<1>;
template class CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV<2>;
template class CVFEM_Upwind_NumIntegral_dNT_rhsop_g_dV<3>;

} // csmp
