#ifndef CSMP_BINARY_FILE_INTERFACE_H
#define CSMP_BINARY_FILE_INTERFACE_H

#include "Model.h"

namespace csmp {

/// read/write connectivity and property data of a model (class kept in library for backward compatibility)
template<size_t dim>
class BinaryFileInterface {
  public:
    BinaryFileInterface();
    ~BinaryFileInterface();

    /// using the functionality of the VSet, stores the mesh of the current model to disk
    bool WriteConnectivityFile( const Model<dim>&, const char* file_name ) const;

    /// reads the mesh that is stored in a binary file into the VSet
    bool ReadConnectivityFile( const char* file_name, VSet<dim>&, double64& time ) const;
  
    /// using FEM_Data containers, writes the variables associated with a model to disk; @note this storage is inefficient as it is padded
    bool WriteDataTo( const Model<dim>&, 
                      const char* file_name, const char* var_name,
                      long timestep=-1 ) const; ///< timestep is printed using 8 digits and padding with zeros; negative number suppresses output
  
    /// reads stand-alone property data from the binary file into the FEM_Data container
    template<typename Var>
    std::string  ReadDataFrom( const char* file_name, 
                               const PropertyDatabase<dim>&, FEM_Data<Var>& ) const;

    ///
    std::string ReadVariableName( const char* file_name ) const;
 };


} // end namespace csmp

#endif
