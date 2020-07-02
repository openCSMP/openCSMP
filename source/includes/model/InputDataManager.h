#ifndef CSMP_INPUT_DATA_MANAGER_H
#define CSMP_INPUT_DATA_MANAGER_H

#include "CSMP_definitions.h"

namespace csmp {

class ComputationalSettings;
template<size_t> class Model;

/** @brief input of parameter values from block-structured text file.
 
@author S.K. Matthai
@author S. Geiger
@author S.G. Roberts
@date 2001
 
@section motivation Motivation
 
To have an ASCII-file based interface to CSMP which allows to assign
material properties and initial conditions to models and their
subregions as identified by unique names (of groups).  

Boundary and essential conditions should also be assignable via condition
flags and this approach needs to work in 1-3D and both for box- and 
for freeform models.  

Since many meshing tools do not allow to assign names to subregions of
models / meshes, but subregions and model internal and external boundaries
are identified by integer or floating point (material property) values,
this interface should also provide the possibility to create named
model subregions from property values and ranges. The only possible
properties for this purpose inside of CSP are element properties. This
is because each subregion must at least consist of one finite element.
 
 
@section design Design Intent

A simple design is achieved by dividing the input data into five 
different data blocks, The user decides which blocks to read,
by setting boolean flag arguments to ConfigureFromFile():

Block 1 - to define named regions of elements with characteristic values
of element properties.  

Block 2 - to set default property values for the entire model.
These can get overwritten in individual subregions by the values 
and flags supplied in the data blocks 3-5.

Block 3 - propert values for  model subregions. The user can determine
whether these values are applied to the entire subregion (option=COMPLETE),
its boundary (option=PERIMETER), or only the interior of the region 
(option=INTERIIR). The CSMP user's guide explain how these parts 
of a region are defined and how Region obkects work in general. 
The perimeter elements of regions are only those which
have at least one of their faces on the region boundary. 

@attention Since perimeter nodes are shared between adjacent regions
the last value assignment instruction will determine the values 
that the target variables at these boundary nodes will have
after initialisation is complete. This is not the case for 
other property placements.

Block 4 - (box-shaped models only): assigns boundary conditions 
using the boundary identifiers defined by the enumeriation BOX_BOUNDARY
in 'Box.h'. Box-shaped models offer a few extra configuration 
possibilities. For instance,
linear boundary property variations can be assigned,
which is not possible for free-form models for which only
uniform boundary values can be assigned.
For the latter, computations performed on the boundaries can 
achieve this goal.

Block 5 - Essential conditions that are assigned by setting property
condition flags for specifc regions. Again the discriminators 'complete, interior,
perimeter' can be used. To assign boundary condition values to irregular
models, one first sets the desired property value in the target region
and then assigns the intended condition flag to it, i.e. DIRICH
or NEUMANN etc. in Block 5.

Any of the data blocks are started and terminated by an empty line.  

Apart from the block-structured approach to the configuration file the 
InputDataManager provides several methods to initialize models created with
specific meshing tools, such as Shewchuk's triangle mesher and ICEM's
suite of meshing tools. 

Computational Settings - is an additional block that can be appended,
allowing the user to define the time stepping strategy, during and
times when simulation results shall be output to file.

Again, a blank line indicates the end of this block and users can 
write any comments, references, and observations made on models
as kind of a documentation, following the configuration.
 
@section applicability Applicability

The InputDataManager can be used to configure any CSMP model, but its
functionality is restricted to that of public interfaces of the
Model.
 
@section participants Participants
 
The data input manager requires a Model to act on- and
interacts with the ComputationalSettings object that is used to 
store runtime information. 
 
@section consequences Consequences
 
The decision that CSMP variable names may contain white space forced us
to restrict the use of token (string item) separators to the tab keys.
Thus every item in any configuration line must be separated by a tab
key and any other white space must be avoided (also at the end of lines).
Comments can be inserted anywhere into configuration files.
To write a comment, start the line with a # sign.

The data input manager allows the user to configure CSMP models without
the need for compilation of a program.
 
@section examples Application Examples

An example configuration file as would be used to configure a model
created with the Triangle mesher and including comments is shown
in the following (the CSMP Example suite contains many others):

 @code
'mymodel-configuration.txt' file created X/X/X - this is its title
@endcode

blank line (thereafter region identifications in terms of their permeability
and in alphabetical order)

@code
# Block 1: creating a region called joint containing all finite elements whose
# permeability ranges between 1.0e-8 and 1.0e-8 m2.
joint				tab		permeability	tab		1.0e-8 1.0e-8
left fault zone		tab		permeability	tab		1.0e-12 1.0e-11
@endcode

blank line (thereafter default properties and initial conditions assigned to whole model)

@code
# Block 2: default properties for the entire the model
permeability   	tab   1.0e-13
storativity    	tab   1.0e-9
@endcode

blank line (thereafter region descriptions and assigned properties)
discriminators here are: 'complete', 'interior', or 'boundary'

@code
# Block 3: specific properties of certain model subregions
left fault zone		tab  	interior 		tab		porosity	tab   0.25
left fault zone		tab  	boundary 		tab		porosity	tab   0.25
@endcode

blank line (thereafter essential conditions, 'TOP' or 'top' will both work,
so will 'Dirichlet' or 'DIRICH', or 'DIRICHLET'). Conditions only apply to boundaries flagged using the mechanism for box-shaped models

@code
# Block 4: Dirichlet boundary conditions assigned to box-shaped models via BOX_BOUNDARY flags
top		tab		Dirichlet	tab		fluid pressure		tab		1.0 1.0 4.3 4.7
bottom	tab		Dirichlet	tab		fluid pressure		tab		1.0 1.0 1.0 1.0
bottom	tab		Dirichlet	tab		displacement		tab		1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 DIRICH PLAIN
# (here the model will know whether to read either 2(2D) or 4(3D) values
@endcode

blank line (thereafter regional condition flags: 'Dirichlet' or 'Neumann')

@code
# Block 5: essential conditions assigned to regions
well 	tab		interior	tab 	fluid pressure 	tab 	Dirichlet
@endcode

blank line (thereafter boundary condition flags for free-form boundaries)

@code
# Block 6: boundart conditions applied to CSMP Boundary objects (that can have arbitrary shape)
# Boundary name       part of       flag           variable name      uniform value on boundary
BOUNDARY_TOP     tab  complete tab  Dirichlet  tab fluid pressure tab 16495146.
@endcode

The configuration file could be read configuring a Model object either with the
specific method:  

 @code
    bool ConfigureFromFile( Model<dim>& sg, const char* fname );
 @endcode

Or with the generic method, but the specific settings: 

@code
    bool ConfigureFromFile( Model<dim>& sg, const char* fname,
                            bool block1, 
                            bool block2,
                            bool block3,
                            bool block4,
                            bool block5,
                            bool block6 );
@endcode

*/
template<size_t dim>
class  InputDataManager {
  public:
    /// writes configuration options to screen
    void Help() const;
  
    /// reads data for blocks 1-5, assigning them to box-shaped models
    bool ConfigureFromFile( Model<dim>&, const char* fname );

    /// key generic configuration method for time-dependent models that can contain Boundary objects
    bool ConfigureFromFile( Model<dim>&, const char* fname,
                            bool block1,            ///< region name from parameter range
                            bool block2,            ///< default property values
                            bool block3,            ///< regional property values
                            bool block4,            ///< boundary conditions for box-shaped model
                            bool block5,            ///< essential conditions for regions
                            bool block6,            ///< boundary conditions for arbitrary-shaped model
                            ComputationalSettings& settings );

    /// costumized configuration of CSMP models including a particualr configuration file
    bool ConfigureFromSpecificFile( Model<dim>&, const char* configuration_fname,
                            bool block1,            ///< region name from parameter range
                            bool block2,            ///< default property values
                            bool block3,            ///< regional property values
                            bool block4,            ///< boundary conditions for box-shaped model
                            bool block5,            ///< essential conditions for regions
                            bool block6,            ///< boundary conditions for arbitrary-shaped model
                            ComputationalSettings& settings );

    /// as above but without runtime information
    bool ConfigureFromFile( Model<dim>&, const char* fname,
                            bool block1,            ///< region name from parameter range
                            bool block2,            ///< default property values
                            bool block3,            ///< regional property values
                            bool block4,            ///< boundary conditions for box-shaped model
                            bool block5,            ///< essential conditions for regions
                            bool block6 = false );  ///< boundary conditions for arbitrary-shaped model

    /// configuration restricted to unique regions in the input model
    bool ConfigureRegionsFromFile( Model<dim>&, const char* fname, std::set<std::string>& regions );

    /// FRED was a FRACMAN consortium of Golder Associates, this interface is used in the Fracman GUI
    bool ConfigureFRED_ModelFromFile( Model<dim>&, const char* fname,
                                      std::map<std::string,std::vector<double64> >& well_data,
                                      ComputationalSettings& settings );

    /// configures box-shaped ANSYS models output using ANSYS' csp input interface
    bool Configure_ANSYS_ModelFromFile( Model<dim>&, const char* fname );

    /// configures arbitrarily shaped ANSYS models output to CSMP, recognising boundaries if their name contains BOUNDARY
    bool ConfigureIrregular_ANSYS_ModelFromFile( Model<dim>&, const char* fname );

private:
    bool ReadBlocks( Model<dim>&,
                     std::ifstream&,
                     std::set<std::string>& regions,
                     std::map<std::string,std::vector<double64> >& well_data,
                     ComputationalSettings& settings,
                     bool region_specifications,      // region name from parameter range
                     bool default_property_values,    // default property values
                     bool region_property_values,     // regional property values
                     bool region_property_conditions, // regional property conditions
                     bool box_boundary_conditions,    // boundary conditions for box-shaped model
                     bool boundary_conditions,        // boundary conditions for arbitrary-shaped model
                     bool well_settings,              // well names, locations, rates, ratios
                     bool computational_settings      // computational settings
                   );
};


}// end namespace csmp

#endif
 
 
 
 
 
 
 
 
 
 
 
 
 
 

