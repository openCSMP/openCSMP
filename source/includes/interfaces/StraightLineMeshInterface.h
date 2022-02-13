#ifndef STRAIGHT_LINE_MESH_INTERFACE_H
#define STRAIGHT_LINE_MESH_INTERFACE_H

#include "VSet.h"

namespace csmp {

/**
@brief Creates mesh based on the parameters provided in text file "*-mesh.txt".

@author R. Manasipov
@date 2014

*/
template< uint32_t dim >
class StraightLineMeshInterface {

  public:
    /// creates mesh based on the parameters provided in text file "*-mesh.txt"
    StraightLineMeshInterface();

    /// creates mesh based on the parameters provided in text file
    StraightLineMeshInterface( csmp::VSet<dim>& vset,
                               const std::string& mesh_file );
    void WriteInfoFile();

  private:


    void WriteSampleFile( std::ofstream& ofs );

    void Initialize( csmp::VSet<dim>& vset,
                     const std::string& mesh_file );

    bool ReadMesh( std::ifstream& ifs, char* text_line, size_t line_length,
                   std::set<std::string>& keywords, std::string& keyword,
                   std::vector<std::string>& keyword_params );

    csmp::VSet<dim>*  vset_;
};


} // end namespace csmp

#endif
