#include "PoreVolumeVisitor.h"
#include "Exception.h"

using namespace std;

namespace csmp {

  /** custom constructor */
  template<size_t dim>
  PoreVolumeVisitor<dim>::PoreVolumeVisitor( Model<dim>& model, 
                                             const char* porosity,
                                             const char* volume,
                                             const char* pore_volume )
  { 

    this->ApplicationLevel(REGION);
    this->ApplicationTarget(NODE);
    // pore_volume_key_ !!
    property_key_ = model.Database().StorageKey(pore_volume);
    phi_key_      = model.Database().StorageKey(porosity);
    vol_key_      = model.Database().StorageKey(volume);

    if ( property_key_.type != SCALAR || property_key_.place != NODE )
      throw csmp::Exception( CSMP_ERROR, "PoreVolumeVisitor::(constructor)", 
                             pore_volume, " must be a nodal scalar property." );

    if ( phi_key_.type != SCALAR || phi_key_.place != NODE )
      throw csmp::Exception( CSMP_ERROR, "PoreVolumeVisitor::(constructor)", 
                             porosity, " must be a scalar property." );

    if ( vol_key_.type != SCALAR || vol_key_.place != NODE )
      throw csmp::Exception( CSMP_ERROR, "PoreVolumeVisitor::(constructor)", 
                             volume, " must be a scalar property." );

  }

  /** deconstructor */
  template<size_t dim>
  PoreVolumeVisitor<dim>::~PoreVolumeVisitor() 
  {}

  /** visit function for Node */
  template<size_t dim>
  void PoreVolumeVisitor<dim>::Visit(Node<dim>* n) 
  { 

    // add the functionality of porosity as an element variable?

    n->Read( phi_key_, phi );

    vol() = 0.;
    pore_vol() = 0.;

    for (unsigned int i = 0; i < n->Parents(); i++)
      {
        volume  = n->Parent(i)->FE()->Volume()/n->Parent(i)->FE()->Nodes();
        pore_vol() += volume;
        vol() += volume;    
      }
    pore_vol() *= phi();
	
    n->Store( property_key_, pore_vol );
    n->Store( vol_key_, vol );

  }

  /** visit function for Region */
  template<size_t dim>
  void PoreVolumeVisitor<dim>::Visit(Region<dim>* n)
  {
	  // no calcution for the region
  }

  template class PoreVolumeVisitor<1U>;
  template class PoreVolumeVisitor<2U>;
  template class PoreVolumeVisitor<3U>;

} // csmp
