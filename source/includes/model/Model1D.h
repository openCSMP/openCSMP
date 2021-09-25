#ifndef MODEL_1D_H
#define MODEL_1D_H

#include "Model.h"
#include "LineElementMesher.h"

namespace csmp {

/// Line model that is templatized so that it can also be placed into three-dimensional space
template<size_t dim=1U>
class Model1D : public Model<dim> {

public:

    /// input from csmp native fileset
    Model1D( const std::string& input_file);

    /// creates an uniform mesh
    Model1D( const std::string& model_name,
             double64 length,
             size_t   elements,
             const std::vector<double64>& splitnode_coordinates = std::vector<double64>(),
             const Point<dim>& origin = Point<dim>(),
             const Point<dim>& destination = Point<dim>() );

    /// creates an uniform mesh
    Model1D( const std::string& model_name,
             const std::string& variable_file,
             double64 length,
             size_t   elements,
             const std::vector<double64>& splitnode_coordinates = std::vector<double64>(),
             const Point<dim>& origin = Point<dim>(),
             const Point<dim>& destination = Point<dim>() );

    /// creates a refined mesh
    Model1D( const std::string& model_name,
             double64 length,
             double64 dx_min,
             double64 dx_max,
             double64 width_of_transition_zone,
             MeshDensity* density,
             const std::vector<double64>& splitnode_coordinates = std::vector<double64>(),
             const Point<dim>& origin = Point<dim>(),
             const Point<dim>& destination = Point<dim>() );

    /// creates a refined mesh
    Model1D( const std::string& model_name,
             const std::string& variable_file,
             double64 length,
             double64 min_size,
             double64 max_size,
             double64 width_of_transition_zone,
             MeshDensity* density,
             const std::vector<double64>& splitnode_coordinates = std::vector<double64>(),
             const Point<dim>& origin = Point<dim>(),
             const Point<dim>& destination = Point<dim>() );

    /// creates a custom mesh with vector( 1D node coordinates )
    Model1D( const std::string& model_name,
             const std::vector<double64>& node_coordinates,
             const std::vector<double64>& splitnode_coordinates = std::vector<double64>(),
             const Point<dim>& origin = Point<dim>(),
             const Point<dim>& destination = Point<dim>() );

    /// creates a custom mesh with vector( 1D node coordinates )
    Model1D( const std::string& model_name,
             const std::string& variable_file,
             const std::vector<double64>& node_coordinates,
             const std::vector<double64>& splitnode_coordinates = std::vector<double64>(),
             const Point<dim>& origin = Point<dim>(),
             const Point<dim>& destination = Point<dim>() );

    /// creates a custom mesh with vector( 1D, 2D or 3D node coordinates )
    Model1D( const std::string& model_name,
             const std::vector<Point<dim> >& node_coordinates,
             const std::vector<Point<dim> >& splitnode_coordinates = std::vector<Point<dim> >(),
             const Point<dim>& origin = Point<dim>(),
             const Point<dim>& destination = Point<dim>() );

    /// creates a custom mesh with vector( 1D, 2D or 3D node coordinates )
    Model1D( const std::string& model_name,
             const std::string& variable_file,
             const std::vector<Point<dim> >& node_coordinates,
             const std::vector<Point<dim> >& splitnode_coordinates = std::vector<Point<dim> >(),
             const Point<dim>& origin = Point<dim>(),
             const Point<dim>& destination = Point<dim>() );

    ~Model1D();

  protected:

    void Initialize( VSet<dim>& );

    /// creates an uniform mesh
    void Initialize( double64 length, size_t elements,
                     const std::vector<double64>&,
                     const Point<dim>&, const Point<dim>& );

    /// creates an exponentially refined mesh
    void Initialize( double64 length, double64 dx_min, double64 dx_max,
                     double64 width_of_transition_zone,
                     MeshDensity* density,
                     const std::vector<double64>&,
                     const Point<dim>&, const Point<dim>& );

    /// creates a custom mesh with vector( 1D node coordinates )
    void Initialize( const std::vector<double64>&,
                     const std::vector<double64>&,
                     const Point<dim>&, const Point<dim>& );

    /// creates a custom mesh with vector( 1D, 2D or 3D node coordinates )
    void Initialize( const std::vector<Point<dim> >&,
                     const std::vector<Point<dim> >&,
                     const Point<dim>&, const Point<dim>& );
};

}

#endif


