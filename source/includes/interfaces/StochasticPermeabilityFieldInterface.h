#ifndef STOCHASTIC_PERMEABILITY_FIELD_INTERFACE_H
#define STOCHASTIC_PERMEABILITY_FIELD_INTERFACE_H

#include "CSMP_definitions.h"
#include "PropertyHandle.h"

namespace csmp {

template<uint32_t> class Model;
template<uint32_t> class PropertyHandle;

template<uint32_t dim>
class StochasticPermeabilityFieldInterface {
    
  public:
    StochasticPermeabilityFieldInterface( Model<dim>& );
    ~StochasticPermeabilityFieldInterface();
    bool Read2DStochasticPermeabilityFieldOnRegularGrid( Model<dim>&, const char* fname );
    bool Read2DStochasticPermeabilityField( Model<dim>&, const char* fname );
    
  private:
    void Tokenize( const std::string& str, std::vector<std::string>& tkns );
    PropertyHandle<dim> nodal_perm, log_perm;
 };

}

/** 
class StochasticPermeabilityFieldInterface  StochasticPermeabilityFieldInterface "interfaces\StochasticPermeabilityFieldInterface.h"
@author S. Geiger
@date 2003

 
@section motivation Motivation

To map HYDROGEN (Bellin & Rubin, 1996) generated 2D stochastic permeability 
fields to CSMP for structured and unstructured finite element meshes.


@section participants Participants

Employs the FiniteDifferenceGrid and FEMToGridVisitor to map the permeability
values to the finite element mesh. PropertyHandles are used to create new variables
'nodal permeability, 'log nodal permeability', and a temporary 'backup permeability'.  
'log nodal permeability' can be used from the main() file to visualize the
stochastic permeability field.
 

@section collaboration Collaboration

Object interacts with the Model. 


@section examples Application Examples
 
The object should be used from the main() file before the comptuation of the
fluid pressure but potentially after creating new Regions, as the initial permeability
field is completely overwritten and cannot be used to identify Regions any longer.

@code
  // Build the 2D Model
  Model< 2>  my_sg( my_mesh_container );
  
  // Setup the interface to read a HYDROGEN generated stochastic permeability field
  StochasticPermeabilityFieldInterface<2U> new_stochastic_k_field( my_sg );
  
  // Read in k-field from file 'k_field.perm' and compute the element permeability values on a regular grid
  new_stochastic_k_field.Read2DStochasticPermeabilityFieldOnRegularGrid( my_sg, "k_field" );
@endcode
 */


#endif
