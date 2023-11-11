#ifndef MPICH_IGNORE_CXX_SEEK
#define MPICH_IGNORE_CXX_SEEK
#endif

#include "CSP_high_level_utilities_parallel.h"

using namespace std;

namespace csp {

template<stl_index  dim>
csp_float  printModelDimensionsParallel( const SuperGroup<csp_float,dim>& sg, bool intermed_or_max )
 {
    vector<csp_float>  xyz(6);
    sg.Dimensions( xyz );
    
    csp_float xmin, ymin, zmin, xmax, ymax, zmax;
    
    MPI::COMM_WORLD.Allreduce( &xyz[0], &xmin, 1, MPI::DOUBLE, MPI::MIN );
    MPI::COMM_WORLD.Allreduce( &xyz[1], &xmax, 1, MPI::DOUBLE, MPI::MAX );
    if ( dim != 1U ) {
        MPI::COMM_WORLD.Allreduce( &xyz[2], &ymin, 1, MPI::DOUBLE, MPI::MIN );
        MPI::COMM_WORLD.Allreduce( &xyz[3], &ymax, 1, MPI::DOUBLE, MPI::MAX );
      }
    if ( dim == 3 ) {  
        MPI::COMM_WORLD.Allreduce( &xyz[4], &zmin, 1, MPI::DOUBLE, MPI::MIN );
        MPI::COMM_WORLD.Allreduce( &xyz[5], &zmax, 1, MPI::DOUBLE, MPI::MAX );
      }

    cout <<"\nprintModelDimensions: Dimensions of model (meters): "<< endl;
    cout <<"xmin, xmax (horizontal right):    "<< xmin <<" "<< xmax << endl;
    if ( dim != 1U ) cout <<"ymin, ymax (vertical upward):     "<< ymin <<" "<< ymax << endl;
    if ( dim == 3U ) cout <<"zmin, zmax (horizontal to front): "<< zmin <<" "<< zmax << endl << endl;

    set<csp_float>  axis;
    axis.insert( xmax - xmin );
    if ( dim != 1U ) axis.insert( ymax - ymin );
    if ( dim == 3U ) axis.insert( zmax - zmin );
    
    set<csp_float>::const_reverse_iterator  it = axis.rbegin();
    
    if ( !intermed_or_max ) return *it;
    
    if ( axis.size() >= 2U ) it++;
    
    return *it;

 } // end printModelDimensions



template<stl_index  dim>
csp_float printRangeOfVariableParallel( const SuperGroup<csp_float,dim>& sg, 
                                        const char* var, bool max_or_min )
 {
     extern csp_float  global_time;
     csp_float         lmin, lmax, gmin, gmax;
     const PropertyDatabase& p_ref = sg.ReferencePropertyDatabase();

     sg.MinMaxOf( var, lmin, lmax );
     
     MPI::COMM_WORLD.Allreduce( &lmin, &gmin, 1, MPI::DOUBLE, MPI::MIN );
     MPI::COMM_WORLD.Allreduce( &lmax, &gmax, 1, MPI::DOUBLE, MPI::MAX );

     cout <<"\nRange of variable ["<< p_ref.Unit(var) <<"]: '";
     cout << var <<"': "<< gmin <<" to "<< gmax << endl;
          
     if ( !max_or_min ) return gmin;
     return gmax;
 }



template<stl_index  dim>
csp_float printRangeOfVariableParallel( const SuperGroup<csp_float,dim>& sg, 
                                        Standard_IO_Handler& io, 
                                        const char* var, bool max_or_min )
 {
     extern csp_float  global_time;
     csp_float         lmin, lmax, gmin, gmax;
     const PropertyDatabase& p_ref = sg.ReferencePropertyDatabase();

     sg.MinMaxOf( var, lmin, lmax );

     MPI::COMM_WORLD.Allreduce( &lmin, &gmin, 1, MPI::DOUBLE, MPI::MIN );
     MPI::COMM_WORLD.Allreduce( &lmax, &gmax, 1, MPI::DOUBLE, MPI::MAX );

     cout <<"\nRange of variable ["<< p_ref.Unit(var) <<"]: '";
     cout << var <<"': "<< gmin <<" to "<< gmax << endl;
     
     // recording the measured variable value range at given timestep
     string var_info("\nt = ");
     char   info[100];
     sprintf( info, "%lf", global_time );
     var_info += info;
     var_info += " secs, range of'";
     var_info += var;
     var_info += "' [";
     var_info += p_ref.Unit(var);
     var_info += "]: ";
     sprintf( info, "%lf", gmin );
     var_info += info;
     var_info += " to ";
     sprintf( info, "%lf", gmax );
     var_info += info;
     var_info += "\n";
     
     io.RecordInformation( var_info ); 
     
     if ( !max_or_min ) return gmin;
     return gmax;
 }




// as above, but for groups
template<stl_index  dim>
csp_float printRangeOfVariableParallel( const SuperGroup<csp_float,dim>& sg, 
                                        const char* group, const char* var, bool max_or_min )
 {
     extern csp_float  global_time;
     csp_float         lmin, lmax, gmin, gmax;
     const PropertyDatabase& p_ref = sg.ReferencePropertyDatabase();
     
     sg.GroupMinMaxOf( group, var, lmin, lmax );
 
     MPI::COMM_WORLD.Allreduce( &lmin, &gmin, 1, MPI::DOUBLE, MPI::MIN );
     MPI::COMM_WORLD.Allreduce( &lmax, &gmax, 1, MPI::DOUBLE, MPI::MAX );

     cout <<"\nRange of variable ["<< p_ref.Unit(var) <<"]: '";
     cout << var <<"' in group '"<< group <<"': "<< gmin <<" to "<< gmax << endl;
          
     if ( !max_or_min ) return lmin;
     return lmax;
 }




template<stl_index  dim>
csp_float printRangeOfVariableParallel( const SuperGroup<csp_float,dim>& sg, 
                                        Standard_IO_Handler& io, 
                                        const char* group,
                                        const char* var, bool max_or_min )
 {
     extern csp_float  global_time;
     csp_float         lmin, lmax, gmin, gmax;
     const PropertyDatabase& p_ref = sg.ReferencePropertyDatabase();

     sg.MinMaxOf( var, lmin, lmax );

     MPI::COMM_WORLD.Allreduce( &lmin, &gmin, 1, MPI::DOUBLE, MPI::MIN );
     MPI::COMM_WORLD.Allreduce( &lmax, &gmax, 1, MPI::DOUBLE, MPI::MAX );

     cout <<"\nRange of variable ["<< p_ref.Unit(var) <<"]: '";
     cout << var <<"': "<< gmin <<" to "<< gmax <<" in group '"<< group <<"'"<< endl;
     
     // recording the measured variable value range at given timestep
     string var_info("\nt = ");
     char   info[100];
     sprintf( info, "%lf", global_time );
     var_info += info;
     var_info += " , group: ";
     var_info += group;
     var_info += ", secs, range of'";
     var_info += var;
     var_info += "' [";
     var_info += p_ref.Unit(var);
     var_info += "]: ";
     sprintf( info, "%lf", gmin );
     var_info += info;
     var_info += " to ";
     sprintf( info, "%lf", gmax );
     var_info += info;
     var_info += "\n";
     
     io.RecordInformation( var_info ); 
     
     if ( !max_or_min ) return gmin;
     return gmax;
 }
 

template
csp_float  printModelDimensionsParallel( const SuperGroup<csp_float,2>& sg, bool intermed_or_max ); 

template
csp_float  printRangeOfVariableParallel( const SuperGroup<csp_float,2>& sg, 
                                         const char* var, bool max_or_min );
template
csp_float  printRangeOfVariableParallel( const SuperGroup<csp_float,2>& sg, 
	                                     Standard_IO_Handler& io, const char* var,
	                                     bool max_instead_of_min );
template
csp_float  printRangeOfVariableParallel( const SuperGroup<csp_float,2>& sg, 
                                         const char* group, const char* var, bool max_or_min );
template
csp_float  printRangeOfVariableParallel( const SuperGroup<csp_float,2>& sg, 
	                                     Standard_IO_Handler& io, 
	                                     const char* group, const char* var,
	                                     bool max_instead_of_min );
  
 
template
csp_float  printModelDimensionsParallel( const SuperGroup<csp_float,3>& sg, bool intermed_or_max ); 

template
csp_float  printRangeOfVariableParallel( const SuperGroup<csp_float,3>& sg, 
                                         const char* var, bool max_or_min );
template
csp_float  printRangeOfVariableParallel( const SuperGroup<csp_float,3>& sg, 
	                                     Standard_IO_Handler& io, const char* var,
	                                     bool max_instead_of_min );
template
csp_float  printRangeOfVariableParallel( const SuperGroup<csp_float,3>& sg, 
                                         const char* group, const char* var, bool max_or_min );
template
csp_float  printRangeOfVariableParallel( const SuperGroup<csp_float,3>& sg, 
	                                     Standard_IO_Handler& io, 
	                                     const char* group, const char* var,
	                                     bool max_instead_of_min );


}
