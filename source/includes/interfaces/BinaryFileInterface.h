#ifndef CSMP_BINARY_FILE_INTERFACE_H
#define CSMP_BINARY_FILE_INTERFACE_H

#include "Model.h"

namespace csmp {

/// read/write connectivity and data from binary files
template<size_t dim>
class BinaryFileInterface {
  public:
    BinaryFileInterface();
    ~BinaryFileInterface();

    bool WriteConnectivityFile( const Model<dim>&, const char* file_name ) const;

    bool ReadConnectivityFile( const char* file_name, VSet<dim>&, double64& time ) const;
  
    bool WriteDataTo( const Model<dim>&, 
                      const char* file_name, const char* var_name, long timestep ) const;
                      
    std::string ReadVariableName( const char* file_name ) const;                  

    template<typename Var>
    std::string  ReadDataFrom( const char* file_name, 
                               const PropertyDatabase<dim>&, FEM_Data<Var>& ) const;
 };


} // end namespace csmp

#endif
