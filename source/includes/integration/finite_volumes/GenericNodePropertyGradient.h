#ifndef GENERIC_NODE_PROPERTY_GRADIENT_H
#define GENERIC_NODE_PROPERTY_GRADIENT_H

#include "PropertyHandle.h"
#include "ErrorHandler.h"
#include "FiniteVolumeStencilManager.h"
#include "FiniteVolumeStencil.h"
#include "Model.h"
 
namespace csmp {

template<size_t> class Region;
/**

@brief Computes the gradient of a node variable at the node using a least-squares approach.

The NodePropertyGradient object can be used to compute the nodal gradient of
a property (solute concentration, fluid phase volume fraction, etc...) such
that the property can be advected with second order accuarcy. It is used by the
NodePropertyGradientLimiter as well in second order accurate advection schemes
to avoid spurious oscillations.

The resulting gradient is stored at the node.

*/
template<size_t dim>
class GenericNodePropertyGradient {

public:
    GenericNodePropertyGradient( Model<dim>&, const char* region, const char* prop="node property gradient");
    GenericNodePropertyGradient( Model<dim>&, const char* region, std::vector<char*> advected_props );
    
    ~GenericNodePropertyGradient(); 
    
    /// provide public access to center of mass of fv:
    void GenericCenterOfMass( size_t global_node_id, VectorVariable<dim>& mass_center ) const;
  
    /// provide public access to distance between facet and mass centre:
    void GenericDistanceFacetFVBarycenter( size_t global_el_id,
                                           size_t local_facet_id,
                                           size_t local_node_id,
                                           VectorVariable<dim>& distance ) const;

    /// computation of property gradient
    void CalculateGenericLeastSquareSums();

    /// computation of property gradient
    void CalculateGenericNodalGradient();

    /// initializing the property key
    void SetPropertyKey( csmp::Index& key );
    
    // storage requirements
    double SizeOf() const;
    
    // TODO: why are these not private?
    Region<dim>&                        gref_;
    double                            tolerance;
    csmp::Index                         u_key;
    const csmp::Index                   grad_key, mctr_key;
    std::vector<PropertyHandle<dim>* >  gradient;
    PropertyHandle<dim>                 mass_center;
    std::vector<double>               det, sum_x2, sum_y2, sum_z2, sum_xy,sum_xz, sum_yz;
    std::vector<size_t>                 zero_grad_index_;
    std::vector<std::vector<size_t> >   neighbors_; ///< global indices of parent nodes; [node_id] -> vector with id's
    
    // 
    /// [EL_Id] [local_facet_id] [local_node_id], distance between facet with ID local_facet_id and node with ID local_node_id
    std::vector<std::vector<std::vector<VectorVariable<dim> > > > distance_facet_FVBarycenter_;
    
 private:
    /// computation of centers of masses of all FV's
    void CalculateCenterOfMass();
     
    void CalculateDistanceFacetFVBary( std::vector<std::vector<std::vector<VectorVariable<dim> > > >& d_facet_FVBarycenter );
     
    void PushBackAvoidDuplicate( std::vector<size_t>& old_vector, std::vector<size_t>& possible_new_entries );
     
    void ConvertToGlobalCoordinates( Element<dim>& el, const Point<dim>& local_c, std::vector<double>& global_c );
     
    std::vector<std::vector<VectorVariable<dim> > >  distance;
    std::vector<std::pair<VectorVariable<dim>, double> > center_of_mass_;
    std::vector<VectorVariable<dim> >                 inner_;
    std::vector< std::vector<VectorVariable<dim> > >  middle_;
};


}
#endif
