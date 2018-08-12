#ifndef ANSYS_MODEL_H
#define ANSYS_MODEL_H

#include "ANSYS_Interface.h"
#include "Model.h"

namespace csmp 
{
  /// Construct a CSMP model from ANSYS *.asc, *.dat, *-variable.txt and *-regions.txt files.
  template<size_t dim>
  class ANSYS_Model : public Model<dim>  {
    public:
      /// input from ANSYS *.asc, *.dat, *-variable.txt and *-regions.txt files.
      /// or from *-1D-mesh.txt
      ANSYS_Model( const std::string& mesh_file_set,
                   const std::string& regions_file_prefix,
                   const std::string& variable_file,
                   bool irregular_mesh,             /* true = non-box shaped model, false = box shaped model */
                   bool binary_file       = true,   /* true = binary, false = ascii */
                   bool use_regions_file  = true,   /* true = reduce regions according to regions file, false = does not redure regions */
                   bool create_boundaries = true,   /* true = creates boundaries around model, false = does not create boundaries */
                   bool ansys_interface   = true ); /* true = use ansys interface, false = use 1D mesher */

      /// input from ANSYS *.asc, *.dat, *-variable.txt and *-regions.txt files.
      /// or from *-1D-mesh.txt
      ANSYS_Model( const std::string& mesh_file_set,
                   const std::string& variable_file,
                   bool irregular_mesh,             /* true = non-box shaped model, false = box shaped model */
                   bool binary_file       = true,   /* true = binary, false = ascii */
                   bool use_regions_file  = true,   /* true = reduce regions according to regions file, false = does not redure regions */
                   bool create_boundaries = true,   /* true = creates boundaries around model, false = does not create boundaries */
                   bool ansys_interface   = true ); /* true = use ansys interface, false = use 1D mesher */

      /// input from ANSYS *.asc, *.dat and *-regions.txt files.
      /// or from *-1D-mesh.txt
      /// creates empty property database
      ANSYS_Model( const std::string& mesh_file_set,
                   bool irregular_mesh,             /* true = non-box shaped model, false = box shaped model */
                   bool binary_file       = true,   /* true = binary, false = ascii */
                   bool use_regions_file  = true,   /* true = reduce regions according to regions file, false = does not redure regions */
                   bool create_boundaries = true,   /* true = creates boundaries around model, false = does not create boundaries */
                   bool ansys_interface   = true ); /* true = use ansys interface, false = use 1D mesher */

      virtual ~ANSYS_Model();

    private:
      void Initialize( const std::string& mesh_file_set,
                       const std::string& regions_file_prefix,
                       bool irregular_mesh,
                       bool binary_file,
                       bool use_regions_file,
                       bool create_boundaries,
                       bool ansys_interface );
  };

} // csmp

#endif
