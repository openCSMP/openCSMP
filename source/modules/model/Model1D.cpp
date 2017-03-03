#include "Model1D.h"
#include "Exception.h"
#include "Region.h"
#include "Model.h"
#include "ModelTopology.h"
#include "CSMP_highLevelUtilities.h"

using namespace std;

namespace csmp {


/**

Default constructor of Model is called.
Then if binary file option is set to true: the model get build from "bin_file".vset and "bin_file".dat files.

*/
template<size_t dim>
Model1D<dim>::Model1D( const std::string& input_file)
 : Model<dim>( input_file.c_str() )
 {
   this->Name( input_file.c_str() );
 }

template<size_t dim>
Model1D<dim>::Model1D( const std::string& input_file,
                       const std::string& variable_file )
 : Model<dim>( input_file.c_str(), variable_file.c_str() )
 {
   this->Name( input_file.c_str() );
 }


/**

Default constructor of Model is called.
Then the model is build using the method Initialize.
*/
template<size_t dim>
Model1D<dim>::Model1D( const std::string& name,
                       const std::string& variable_file,
                       double64 length,
                       size_t   elements,
                       const std::vector<double64>& splitnode_coordinates,
                       const Point<dim>& origin,
                       const Point<dim>& destination )
 : Model<dim>( variable_file.c_str(), false )
 {
    this->Name(name.c_str());
    Initialize( length, elements,
                splitnode_coordinates,
                origin, destination );
 }

template<size_t dim>
Model1D<dim>::Model1D( const std::string& name,
                       double64 length,
                       size_t   elements,
                       const std::vector<double64>& splitnode_coordinates,
                       const Point<dim>& origin,
                       const Point<dim>& destination )
 {
    this->Name(name.c_str());
    Initialize( length, elements,
                splitnode_coordinates,
                origin, destination );
 }


template<size_t dim>
Model1D<dim>::Model1D(const std::string& name,
                       const std::string& variable_file,
                       double64 length,
                       double64 min_size, double64 max_size,
                       double64 width_of_transition_zone,
                       MeshDensity* density,
                       const std::vector<double64>& splitnode_coordinates,
                       const Point<dim>& origin,
                       const Point<dim>& destination )
 : Model<dim>( variable_file.c_str(), false )
 {
    this->Name(name.c_str());
    Initialize( length, min_size, max_size, width_of_transition_zone, density,
                splitnode_coordinates,
                origin, destination );
 }

template<size_t dim>
Model1D<dim>::Model1D(const std::string& name,
                      double64 length,
                      double64 dx_min, double64 dx_max,
                      double64 width_of_transition_zone,
                      MeshDensity* density,
                      const std::vector<double64>& splitnode_coordinates,
                      const Point<dim>& origin,
                      const Point<dim>& destination )
 {
    this->Name(name.c_str());
    Initialize( length, dx_min, dx_max, width_of_transition_zone, density,
                splitnode_coordinates,
                origin, destination );
 }

template<size_t dim>
Model1D<dim>::Model1D( const std::string& name,
                       const std::string& variable_file,
                       const std::vector<double64>& node_coordinates,
                       const std::vector<double64>& splitnode_coordinates,
                       const Point<dim>& origin,
                       const Point<dim>& destination )
 : Model<dim>( variable_file.c_str(), false )
{
   this->Name(name.c_str());
   Initialize( node_coordinates,
               splitnode_coordinates,
               origin, destination );
}

template<size_t dim>
Model1D<dim>::Model1D( const std::string& name,
                       const std::vector<double64>& node_coordinates,
                       const std::vector<double64>& splitnode_coordinates,
                       const Point<dim>& origin,
                       const Point<dim>& destination )
{
   this->Name(name.c_str());
   Initialize( node_coordinates,
               splitnode_coordinates,
               origin, destination );
}

template<size_t dim>
Model1D<dim>::Model1D( const std::string& name,
                       const std::string& variable_file,
                       const std::vector<Point<dim> >& node_coordinates,
                       const std::vector<Point<dim> >& splitnode_coordinates,
                       const Point<dim>& origin,
                       const Point<dim>& destination )
 : Model<dim>( variable_file.c_str(), false )
{
   this->Name(name.c_str());
   Initialize( node_coordinates,
               splitnode_coordinates,
               origin, destination );
}

template<size_t dim>
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


template<size_t dim>
Model1D<dim>::~Model1D()
 {
 }




/// creates mesh from VSet
template<size_t dim>
void Model1D<dim>::Initialize( VSet<dim>& vset )
{
    // Create Model from Vset data and conectivity information
    const bool create_boundaries                    ( false );
    const bool non_box_shaped_model                 ( false );
    const bool isoparametric_elements               ( true  );

    Model<dim>::Initialize( isoparametric_elements,
                            vset,
                            create_boundaries,
                            non_box_shaped_model
                          );
}

/// creates an uniform mesh

template<size_t dim>
void Model1D<dim>::Initialize( double64 length, size_t elements,
                               const std::vector<double64>& splitnode_coordinates,
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

template<size_t dim>
void Model1D<dim>::Initialize(double64 length, double64 dx_min, double64 dx_max, double64 width_of_transition_zone, MeshDensity* density,
                              const std::vector<double64>& splitnode_coordinates,
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

template<size_t dim>
void Model1D<dim>::Initialize( const std::vector<double64>& node_coordinates,
                               const std::vector<double64>& splitnode_coordinates,
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

template<size_t dim>
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






