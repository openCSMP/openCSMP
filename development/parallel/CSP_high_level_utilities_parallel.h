#ifndef CSP_HIGH_LEVEL_UTILITIES_PARALLEL_H
#define CSP_HIGH_LEVEL_UTILITIES_PARALLEL_H

#include "SuperGroup.h"
#include "mpi.h"

namespace csp {

template<stl_index  dim>
csp_float  printModelDimensionsParallel( const SuperGroup<csp_float,dim>& sg, bool intermed_or_max ); 

template<stl_index  dim>
csp_float  printRangeOfVariableParallel( const SuperGroup<csp_float,dim>& sg, 
                                         const char* var, bool max_or_min=true );
template<stl_index  dim>
csp_float  printRangeOfVariableParallel( const SuperGroup<csp_float,dim>& sg, 
	                                     Standard_IO_Handler& io, const char* var,
	                                     bool max_instead_of_min=true );
template<stl_index  dim>
csp_float  printRangeOfVariableParallel( const SuperGroup<csp_float,dim>& sg, 
                                         const char* group, const char* var, bool max_or_min=true );
                                         
template<stl_index  dim>
csp_float  printRangeOfVariableParallel( const SuperGroup<csp_float,dim>& sg, 
                                         Standard_IO_Handler& io, const char* group,
                                         const char* var, bool max_or_min=true );
                                         
}	                                     


#endif
