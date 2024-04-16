#ifndef ANSYS_MODEL_2D_H
#define ANSYS_MODEL_2D_H

#include "Model.h"

namespace csmp {
/**
       2D specialisation of Model for input files from ANSYS interface.
       
       @attention Boundary objects will be created automatically if there are the  lower-dimensional regions with Box boudary names or BOUNDARY in their name strings.
       
       @attention SplitBoundary objects will be created automatically if the model is not contiguous.
 */
class ANSYS_Model2D : public Model<2U> {
  public:

  /// input from ANSYS *.asc, *.dat, and *-variable.txt files, see source for more extensive documentation
  ANSYS_Model2D( const char* icem_file_set,
                 const char* regions_file_prefix,
                 const char* variable_file,
                 bool irregular_mesh = false,
                 bool binary_file = true,
                 bool use_regions_file = true );

  /// input from ANSYS *.asc, *.dat and *-variable.txt files, as above but allows to choose whether analytically integrated of isoparametric numercially integrated elements are used
  ANSYS_Model2D( bool isoparametric,
                const char* icem_file_set,
                const char* variable_file,
                bool irregular_mesh = false,           /* true = non-box shaped model, false = box shaped model */
                bool binary_file = true,               /* true = binary, false = ascii */
                bool use_regions_file = true );        /* true = reduce regions according to regions file, false = does not redure regions */

  /// input from ANSYS *.asc, *.dat and *-variable.txt files
  ANSYS_Model2D(
    const char* icem_file_set,
    const char* variable_file,
    bool irregular_mesh = false,           /* true = non-box shaped model, false = box shaped model */
    bool binary_file = true,               /* true = binary, false = ascii */
    bool use_regions_file = true );        /* true = reduce regions according to regions file, false = does not redure regions */

  /// input from ANSYS *.asc, *.dat files; creates empty property database
  ANSYS_Model2D(
    const char* icem_file_set,
    bool irregular_mesh = false,           /* true = non-box shaped model, false = box shaped model */
    bool binary_file = true,               /* true = binary, false = ascii */
    bool use_regions_file = true );        /* true = reduce regions according to regions file, false = does not redure regions */

  virtual ~ANSYS_Model2D();

private:

  void InitializeANSYS( bool isoparametric,
                        const char* mesh_file_set,
                        const char* regions_file_prefix,
                        bool irregular_mesh,
                        bool binary_input_file,
                        bool use_regions_file );

  void InitializeANSYS( const char* mesh_file_set,
                        const char* regions_file_prefix,
                        bool irregular_mesh,
                        bool binary_file,
                        bool use_regions_file );

};

} // end csmp

#endif


