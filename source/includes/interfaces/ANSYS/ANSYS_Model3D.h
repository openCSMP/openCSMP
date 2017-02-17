#ifndef ANSYS_MODEL3D_H
#define ANSYS_MODEL3D_H

#include "Model.h"

namespace csmp {


/**
   ANSYS 3D box-shaped and irregular meshes (including fracture-only DFN models)
*/

class ANSYS_Model3D : public Model<3U> {
  public:

    /// input from ANSYS *.asc, *.dat, and *-variable.txt files.
    /// also loads -regions file using prefix from the configuration file
    ANSYS_Model3D( const char* icem_file_set,
                   const char* regions_file_prefix,
                   const char* variable_file,
                   bool irregular_mesh    = false,  /* true = free-form model, but box boundaries will still be picked up; false = only box boundaries */
                   bool binary_file       = true,   /* true = binary, false = ascii */
                   bool use_regions_file  = true,   /* true = reduce regions according to regions file, false = does not redure regions */
                   bool create_boundaries = true);  /* true = creates boundaries around model, false = does not create boundaries */

    /// input from ANSYS *.asc, *.dat and *-variable.txt files
    ANSYS_Model3D( const char* icem_file_set, 
                   const char* variable_file,
                   bool irregular_mesh    = false,  /* true = free-form model, but box boundaries will still be picked up; false = only box boundaries */
                   bool binary_file       = true,   /* true = binary, false = ascii */
                   bool use_regions_file  = true,   /* true = reduce regions according to regions file, false = does not redure regions */
                   bool create_boundaries = true);  /* true = creates boundaries around model, false = does not create boundaries */

    /// input from ANSYS *.asc, *.dat files
    /// creates empty property database
    ANSYS_Model3D( const char* icem_file_set,
                   bool irregular_mesh    = false,  /* true = free-form model, but box boundaries will still be picked up; false = only box boundaries */
                   bool binary_file       = true,   /* true = binary, false = ascii */
                   bool use_regions_file  = true,   /* true = reduce regions according to regions file, false = does not redure regions */
                   bool create_boundaries = true);  /* true = creates boundaries around model, false = does not create boundaries */

    // To rebuild model from CSMP native binary file do not use an ANSYS model

    virtual ~ANSYS_Model3D();
  
    /// renumbers the nodes (0..n) as in the original ANSYS model; returns true if changes were made
    bool RestoreOriginalNodeNumbering( bool verbose );
  
    /// access to the node points of the model in the original order output from ANSYS
    std::vector<Point<3U> >::const_iterator VerticesBegin() const;
    std::vector<Point<3U> >::const_iterator VerticesEnd() const;
  
  private:

    void Initialize( const char* icem_file_set,
                     const char* regions_file_prefix,
                     bool irregular_mesh,
                     bool binary_file,
                     bool use_regions_file,
                     bool create_boundaries );
  
    std::vector<Point<3U> > node_coords_; ///< node coordinates in VSet order to re-establish original node numbering if necessary
};

} // end csmp

#endif



