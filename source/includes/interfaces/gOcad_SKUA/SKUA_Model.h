#ifndef CSMP_SKUA_MODEL_H
#define CSMP_SKUA_MODEL_H

#include "Model.h"

namespace csmp {


/**
    Imports 3D SKUA  models meshed with SKUA's FiniteElement mesher plugin and using an output procedure
    created by Thomas Jerome (GMDK Inc.).
 
    @attention if the variables "element number" and / or "node number" are in the variables file,
    they will be initialised with the data from the original VSet.
    This is important because CSMP's internal element and node numbering varies from run to run
    since its mesh connectivity is based on pointers.
 
    @note if elements are elimitated from the model in the build process, corresponding ID values dissappear.
    However, the sequence of the element is not changed which means that "element number" can be collapsed
    to achieve a consecutive range again.
    
    @attention To rebuild model from CSMP native binary file do not use this SKUA model class, 
    the standard CSMP model. All functionality should be there.
*/
class SKUA_Model : public Model<3U> {
public:

  /// input from SKUA *.asc, *.dat, and *-variable.txt files.
  /// also loads -regions file using prefix from the configuration file
  SKUA_Model( const char* input_file_set,
              const char* variable_file,
              bool binary_file = true,             ///< true = binary, false = ascii */
              bool use_regions_file = true,        ///< true = reduce regions according to regions file, false = does not redure regions
              bool create_boundaries = true,       ///< true = creates boundaries around model, false = does not create boundaries 
              bool create_splitboundaries = false, ///< TODO: true = creates splitboundaries around model, false = does not create splitboundaries 
              bool isoparametric = true );         ///< only option if model does not only consist of simplex elements

  virtual ~SKUA_Model();

  /// renumbers the nodes (0..n) as in the original SKUA model; returns true if changes were made
  bool RestoreOriginalNodeNumbering( bool verbose );

  /// access to the node points of the model in the original order output from SKUA
  std::vector<Point<3U> >::const_iterator VerticesBegin() const;
  std::vector<Point<3U> >::const_iterator VerticesEnd() const;

private:

  void Initialize( const char* mesh_file_set,
                   bool binary_file,
                   bool use_regions_file,
                   bool create_boundaries,
                   bool create_splitboundaries,
                   bool isoparametric );

  std::vector<Point<3U> > node_coords_;        ///< node coordinates in VSet order to re-establish original node numbering if necessary
};

} // end csmp

#endif /* CSMP_SKUA_MODEL_H */



