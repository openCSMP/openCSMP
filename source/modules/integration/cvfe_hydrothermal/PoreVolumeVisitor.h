#ifndef PORE_VOLUME_VISITOR_H
#define PORE_VOLUME_VISITOR_H

/*   Changelog
     November 2014, Philipp Weis:
	 - initial port to CSMP++ and strong simplification as compared to the csp5-version.
*/

#include "Visitor.h"
//#include "VariableEnums.h"
#include "Model.h"

namespace csmp {

 /// Calculating the volume and pore volume of a node-centered control volume.

template<size_t dim>
class PoreVolumeVisitor : public Visitor<dim> {
  public:
    PoreVolumeVisitor( Model<dim>&, 
	                   const char* porosity,      // nodal variable for porosity [0 to 1]
	                   const char* volume,        // nodal variable to be used in further calculations and post-processing
	                   const char* pore_volume ); // nodal variable to be used in further calculations and post-processing
    ~PoreVolumeVisitor();
    
    virtual void Visit(Node<dim>* n); // Application level and target are set to NODE in the constructor.
	virtual void Visit(Region<dim>* n); // Application level and target are set to NODE in the constructor.

  private:

    csmp::Index property_key_, phi_key_, vol_key_;               
    ScalarVariable phi, pore_vol, vol;
    double64 volume;

};

  /**
     @class PoreVolumeVisitor PoreVolumeVisitor.h

     @author Philipp Weis, ETH Zuerich
     @section contact Contact
     philipp.weis@erdw.ethz.ch

     @changes changes Latest Changes                                                                                  
  
     @section motivation Motivation
      Intitially designed to avoid the usage of the FiniteVolumeManager in csp5.
	  Output as nodal field variables allows to use them for post-processing (e.g. with ParaView or MatLab).

     @section usage Usage
      Used within the CVFEM scheme (Weis et al., Geofluids, 2014).

     @code
	 The visitor reads in a nodal variable for porosity, which has to be initialised before visitation
	 It further requires existing nodal variables for volume and pore volume.
	 The three variable names are provided as constructor arguments.

     The visitor calculates the node-centered volume and pore volume, assuming that all elements contribute equally to the control volumes of their respecetive nodes.

	 The values for volume and pore volume are stored as nodal variables.
          
     @endcode
     
     @section dependencies Dependencies
     
     @section issues Known issues
	 There may be new funcionality in CSMP++ that could make this visitor redundant.
	 Only Visit(Node<dim>* n) is overwritten. However, application target and level are always set to NODE.
	 How to limit to Region (former Group)?

	 Does the Application Level have to be changed?
	 How can it be restricted to a Region?
	 Include the possibility to use porosity as an element variable?
	 Change name to NodalPoreVolumeVisitor?
     
     @section testing Testing
     testing was done in the period before publication in 2014.

  */
 
} // csmp

#endif
