#ifndef CGNS_MODEL_H
#define CGNS_MODEL_H

#include "CGNS_Interface.h"
#include "Model.h"

#include "ModelTopology.h"
#include "ModelTime.h"

#include "Exception.h"

#include <string>

namespace csmp 
{

  template<uint32_t dim>
  class CGNS_Model : public Model<dim>
    {

    public:

      /// input from CGNS *.cgns, *.dat, *-variable.txt and *-regions.txt files.
      CGNS_Model( const std::string& mesh_file_set,
                  const std::string& regions_file_prefix,
                  const std::string& variable_file,
                  bool use_regions_file  = true,   /* true = reduce regions according to regions file, false = does not redure regions */
                  bool create_boundaries = true );  /* true = creates boundaries around model, false = does not create boundaries */

      /// input from CGNS *.cgns, *.dat, *-variable.txt and *-regions.txt files.
      CGNS_Model( const std::string& mesh_file_set,
                  const std::string& variable_file,
                  bool use_regions_file  = true,   /* true = reduce regions according to regions file, false = does not redure regions */
                  bool create_boundaries = true );  /* true = creates boundaries around model, false = does not create boundaries */

      /// input from CGNS *.cgns, *.dat and *-regions.txt files.
      /// creates empty property database
      CGNS_Model( const std::string& mesh_file_set,
                  bool use_regions_file  = true,   /* true = reduce regions according to regions file, false = does not redure regions */
                  bool create_boundaries = true );  /* true = creates boundaries around model, false = does not create boundaries */

      /// input from CSMP binary *.vset *-regions.dat, *-boundaries.dat and *-variables.txt files
      CGNS_Model( const std::string& vset_dat_files,
                  const std::string& variable_file );

      /// input from CSMP binary *.vset *-regions.dat, *-boundaries.dat files
      CGNS_Model( const std::string& vset_dat_files );

      virtual ~CGNS_Model();

      void Write_CGNS_Mesh( const std::string& mesh_file_set );

    private:

      void Initialize( const std::string& mesh_file_set,
                       const std::string& regions_file_prefix,
                       bool use_regions_file,
                       bool create_boundaries );
  };

  /**
  @class CGNS_Model CGNS_Model "interfaces/CGNS_Model.h"

  @author R. Manasipov
  @date 2015

  */

} // csmp

#endif
