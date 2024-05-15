#ifndef ANSYS_MODEL3D_H
#define ANSYS_MODEL3D_H

#include "Model.h"

namespace csmp {

/**
    Model for  the construction from ANSYS ICEM 3D box-shaped and irregular meshes (including fracture-only DFN models).
    Additional "-regions.txt" file is required to select the mesh domains that shall be incorporated into the CSMP model.
 
    @attention if the variables "element number" and / or "node number" are in the variables file,
    they will be initialised with the data from the original VSet.
    This is important because CSMP's internal element and node numbering varies from run to run
    since its mesh connectivity is based on pointers.
    
    @note boundaries will always be created because they belong to intact CSMP models. This means that correspondingly
    named (see BOX_BOUNDARY or include string 'BOUNDARY' in their names ) surfaces must be present in the mesh. 
*/
class ANSYS_Model3D : public Model<3U> {
  public:
    /// input from ANSYS *.asc, *.dat, and *-variable.txt files.
    /// also loads -regions file using prefix from the configuration file
    ANSYS_Model3D( const char* icem_file_set,
                   const char* regions_file_prefix,
                   const char* variable_file,
                   bool binary_file = true );   /* true = reduce regions according to regions file, false = does not redure regions */

    /// (non-)isoparametric input from ANSYS *.asc, *.dat and *-variable.txt files
    ANSYS_Model3D( bool isoparametric,
                   const char* icem_file_set,
                   const char* variable_file,
                   bool binary_file = true  );   /* true = reduce regions according to regions file, false = does not redure regions */

    /// input from ANSYS *.asc, *.dat and *-variable.txt files
    ANSYS_Model3D( const char* icem_file_set,
                   const char* variable_file,
                   bool binary_file = true );    /* true = reduce regions according to regions file, false = does not redure regions */

    /// input from ANSYS *.asc, *.dat files
    /// creates empty property database
    ANSYS_Model3D( const char* icem_file_set,
                   bool binary_file = true );    /* true = reduce regions according to regions file, false = does not redure regions */

    // To rebuild model from CSMP native binary file do not use an ANSYS model
    virtual ~ANSYS_Model3D();

    /// renumbers the nodes (0..n) as in the original ANSYS model; returns true if changes were made
    bool RestoreOriginalNodeNumbering( bool verbose );

    /// access to the node points of the model in the original order output from ANSYS
    std::vector<Point<3U> >::const_iterator VerticesBegin() const;
    std::vector<Point<3U> >::const_iterator VerticesEnd() const;

private:
    void InitializeANSYS( const char* mesh_file_set,
                           const char* regions_file_prefix,
                           bool binary_input_file );

    void InitializeANSYS( bool isoparametric,
                           const char* mesh_file_set,
                           const char* regions_file_prefix,
                           bool binary_file );

    std::vector<Point<3U> > node_coords_;        ///< node coordinates in VSet order to re-establish original node numbering if necessary
};

} // end csmp

#endif



