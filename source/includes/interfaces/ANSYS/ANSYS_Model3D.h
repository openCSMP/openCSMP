#ifndef ANSYS_MODEL3D_H
#define ANSYS_MODEL3D_H

#include "Model.h"

namespace csmp {


/**
    ANSYS ICEM 3D box-shaped and irregular meshes (including fracture-only DFN models).
 
    @attention if the variables "element number" and / or "node number" are in the variables file,
    they will be initialised with the data from the original VSet.
    This is important because CSMP's internal element and node numbering varies from run to run
    since its mesh connectivity is based on pointers.
 
    @note if elements are elimitated from the model in the build process, corresponding ID values dissappear.
    However, the sequence of the element is not changed which means that "element number" can be collapsed
    to achieve a consecutive range again.
    
    @attention boudaries will always be created because they belong to intact CSMP models.
*/
class ANSYS_Model3D : public Model<3U> {
public:

  /// input from ANSYS *.asc, *.dat, and *-variable.txt files.
  /// also loads -regions file using prefix from the configuration file
  ANSYS_Model3D( const char* icem_file_set,
                 const char* regions_file_prefix,
                 const char* variable_file,
                 bool irregular_mesh = false,           /* true = free-form model, but box boundaries will still be picked up; false = only box boundaries */
                 bool binary_file = true,               /* true = binary, false = ascii */
                 bool use_regions_file = true );        /* true = reduce regions according to regions file, false = does not redure regions */

  /// (non-)isoparametric input from ANSYS *.asc, *.dat and *-variable.txt files
  ANSYS_Model3D( bool isoparametric,
                 const char* icem_file_set,
                 const char* variable_file,
                 bool irregular_mesh = false,           /* true = free-form model, but box boundaries will still be picked up; false = only box boundaries */
                 bool binary_file = true,               /* true = binary, false = ascii */
                 bool use_regions_file = true  );       /* true = reduce regions according to regions file, false = does not redure regions */

  /// input from ANSYS *.asc, *.dat and *-variable.txt files
  ANSYS_Model3D( const char* icem_file_set,
                 const char* variable_file,
                 bool irregular_mesh = false,           /* true = free-form model, but box boundaries will still be picked up; false = only box boundaries */
                 bool binary_file = true,               /* true = binary, false = ascii */
                 bool use_regions_file = true );        /* true = reduce regions according to regions file, false = does not redure regions */

  /// input from ANSYS *.asc, *.dat files
  /// creates empty property database
  ANSYS_Model3D( const char* icem_file_set,
                 bool irregular_mesh = false,           /* true = free-form model, but box boundaries will still be picked up; false = only box boundaries */
                 bool binary_file = true,               /* true = binary, false = ascii */
                 bool use_regions_file = true );        /* true = reduce regions according to regions file, false = does not redure regions */

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
                         bool irregular_mesh,
                         bool binary_input_file,
                         bool use_regions_file );

  void InitializeANSYS( bool isoparametric,
                         const char* mesh_file_set,
                         const char* regions_file_prefix,
                         bool irregular_mesh,
                         bool binary_file,
                         bool use_regions_file );

  std::vector<Point<3U> > node_coords_;        ///< node coordinates in VSet order to re-establish original node numbering if necessary
};

} // end csmp

#endif



