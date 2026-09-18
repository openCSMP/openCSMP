// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "Model1D.h"
#include "Exception.h"
#include "Region.h"
#include "Model.h"
#include "ModelTopology.h"

using namespace std;

namespace csmp {


/**

Default constructor of Model is called.
Then if binary file option is set to true: the model get build from "bin_file".vset and "bin_file".dat files.

*/
template<uint32_t dim>
Model1D<dim>::Model1D( const std::string& input_file)
 : Model<dim>( input_file.c_str() )
 {
   this->Name( input_file.c_str() );
 }



/**

Default constructor of Model is called.
Then the model is build using the method Initialize.
*/
template<uint32_t dim>
Model1D<dim>::Model1D( const std::string& name,
                       const std::string& variable_file,
                       double length,
                       size_t   elements,
                       const std::vector<double>& splitnode_coordinates,
                       const Point<dim>& origin,
                       const Point<dim>& destination )
 : Model<dim>( variable_file.c_str() )
 {
    this->Name(name.c_str());
    Initialize( length, elements,
                splitnode_coordinates,
                origin, destination );
 }



template<uint32_t dim>
Model1D<dim>::Model1D( const std::string& name,
                       double length,
                       size_t   elements,
                       const std::vector<double>& splitnode_coordinates,
                       const Point<dim>& origin,
                       const Point<dim>& destination )
 {
    this->Name(name.c_str());
    Initialize( length, elements,
                splitnode_coordinates,
                origin, destination );
 }



template<uint32_t dim>
Model1D<dim>::Model1D(const std::string& name,
                       const std::string& variable_file,
                       double length,
                       double min_size, double max_size,
                       double width_of_transition_zone,
                       MeshDensity* density,
                       const std::vector<double>& splitnode_coordinates,
                       const Point<dim>& origin,
                       const Point<dim>& destination )
 : Model<dim>( variable_file.c_str() )
 {
    this->Name(name.c_str());
    Initialize( length, min_size, max_size, width_of_transition_zone, density,
                splitnode_coordinates,
                origin, destination );
 }

template<uint32_t dim>
Model1D<dim>::Model1D(const std::string& name,
                      double length,
                      double dx_min, double dx_max,
                      double width_of_transition_zone,
                      MeshDensity* density,
                      const std::vector<double>& splitnode_coordinates,
                      const Point<dim>& origin,
                      const Point<dim>& destination )
 {
    this->Name(name.c_str());
    Initialize( length, dx_min, dx_max, width_of_transition_zone, density,
                splitnode_coordinates,
                origin, destination );
 }
 
 

template<uint32_t dim>
Model1D<dim>::Model1D( const std::string& name,
                       const std::string& variable_file,
                       const std::vector<double>& node_coordinates,
                       const std::vector<double>& splitnode_coordinates,
                       const Point<dim>& origin,
                       const Point<dim>& destination )
 : Model<dim>( variable_file.c_str() )
{
   this->Name(name.c_str());
   Initialize( node_coordinates,
               splitnode_coordinates,
               origin, destination );
}



template<uint32_t dim>
Model1D<dim>::Model1D( const std::string& name,
                       const std::vector<double>& node_coordinates,
                       const std::vector<double>& splitnode_coordinates,
                       const Point<dim>& origin,
                       const Point<dim>& destination )
{
   this->Name(name.c_str());
   Initialize( node_coordinates,
               splitnode_coordinates,
               origin, destination );
}



template<uint32_t dim>
Model1D<dim>::Model1D( const std::string& name,
                       const std::string& variable_file,
                       const std::vector<Point<dim> >& node_coordinates,
                       const std::vector<Point<dim> >& splitnode_coordinates,
                       const Point<dim>& origin,
                       const Point<dim>& destination )
 : Model<dim>( variable_file.c_str() )
{
   this->Name(name.c_str());
   Initialize( node_coordinates,
               splitnode_coordinates,
               origin, destination );
}



template<uint32_t dim>
Model1D<dim>::Model1D( const std::string& name,
                       const std::vector<Point<dim> >& node_coordinates,
                       const std::vector<Point<dim> >& splitnode_coordinates,
                       const Point<dim>& origin,
                       const Point<dim>& destination )
{
   this->Name(name.c_str());
   Initialize( node_coordinates,
               splitnode_coordinates,
               origin, destination );
}


template<uint32_t dim>
Model1D<dim>::~Model1D()
 {
 }




/// creates mesh from VSet
template<uint32_t dim>
void Model1D<dim>::Initialize( VSet<dim>& vset )
{
    // Create Model from Vset data and conectivity information
//    const bool create_boundaries                    ( false );
//    const bool non_box_shaped_model                 ( false );
//    const bool isoparametric_elements               ( true  );

    Model<dim>::Initialize( vset );
}

/// creates an uniform mesh

template<uint32_t dim>
void Model1D<dim>::Initialize( double length, size_t elements,
                               const std::vector<double>& splitnode_coordinates,
                               const Point<dim>& origin,
                               const Point<dim>& destination )
{
  VSet<dim>      vset;
  LineElementMesher<dim>  mesher;
  mesher.BuildUniformMesh( vset,
                           length, elements,
                           splitnode_coordinates,
                           origin, destination );

  Initialize( vset );

 } // end Initialize




/// creates an exponentially refined mesh

template<uint32_t dim>
void Model1D<dim>::Initialize(double length, double dx_min, double dx_max, double width_of_transition_zone, MeshDensity* density,
                              const std::vector<double>& splitnode_coordinates,
                              const Point<dim>& origin,const Point<dim>& destination )
{
  VSet<dim>      vset;
  LineElementMesher<dim>  mesher;
  mesher.BuildRefinedMesh( vset,
                           length, dx_min, dx_max, width_of_transition_zone, density,
                           splitnode_coordinates,
                           origin, destination );

  Initialize( vset );

 } // end Initialize


/// creates a custom mesh with vector( 1D node coordinates )

template<uint32_t dim>
void Model1D<dim>::Initialize( const std::vector<double>& node_coordinates,
                               const std::vector<double>& splitnode_coordinates,
                               const Point<dim>& origin,
                               const Point<dim>& destination )
{
    VSet<dim>      vset;
    LineElementMesher<dim>  mesher;
    mesher.BuildCustomMesh( vset,
                            node_coordinates,
                            splitnode_coordinates,
                            origin, destination );

    Initialize( vset );

} // end Initialize



/// creates a custom mesh with vector( 1D, 2D or 3D node coordinates )

template<uint32_t dim>
void Model1D<dim>::Initialize( const std::vector<Point<dim> >& node_coordinates,
                               const std::vector<Point<dim> >& splitnode_coordinates,
                               const Point<dim>& origin,
                               const Point<dim>& destination )
{
    VSet<dim>      vset;
    LineElementMesher<dim>  mesher;
    mesher.BuildCustomMesh( vset,
                            node_coordinates,
                            splitnode_coordinates,
                            origin, destination );

    Initialize( vset );

} // end Initialize


template class Model1D<1U>;
template class Model1D<2U>;
template class Model1D<3U>;

} // end csmp






