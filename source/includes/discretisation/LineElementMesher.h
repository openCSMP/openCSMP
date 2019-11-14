#ifndef CSMP_LINE_ELEMENT_MESHER_H
#define CSMP_LINE_ELEMENT_MESHER_H

#include "VSet.h"

namespace csmp {

/// mesh refinement helper base class
class MeshDensity {
public:
    MeshDensity( );
    virtual ~MeshDensity();
    virtual double64 operator()( double64 x );
};

/// helper subclass
class LinDensity: virtual public MeshDensity {
public:
    LinDensity( double64 a = 1.0, double64 b = 0.0 );
    virtual ~LinDensity();
    virtual double64 operator()( double64 x );
    void SetA( double64 );
    void SetB( double64 );
private:
    /// y = a_ * x + b_
    double64 a_;
    double64 b_;
};

/// exponential element size variation
class ExpDensity: virtual public MeshDensity {
public:
    ExpDensity( double64 a = 1.0, double64 b = 1.0  );
    virtual ~ExpDensity();
    virtual double64 operator()(double64 x );
    void SetA( double64 );
    void SetB( double64 );
private:
    /// y = b_ * exp( a_ * x )
    double64 a_;
    double64 b_;
};

/// error function based mesh size variation
class ErfDensity: virtual public MeshDensity {
public:
    ErfDensity( double64 a = 3.0, double64 b = 1.0 );
    ~ErfDensity();
    virtual double64 operator()(double64 x );
    void SetA( double64 );
    void SetB( double64 );
private:
    /// y = b_ * erf( a * x )
    double64 a_;
    double64 b_;
};


/**
     Line element mesher.
     
     @author Stephan Matthai
     @author extended by Roman Manasipov
     @date 2002, 2012
*/
template<size_t dim>
class LineElementMesher {
  public:
    LineElementMesher( );

    /// creates an uniform mesh
    void BuildUniformMesh( VSet<dim>& vset,
                           double64 length,
                           size_t   n_elements,
                           const std::vector<double64>& splitnode_coordinates = std::vector<double64>(),
                           const Point<dim>& origin = Point<dim>(),
                           const Point<dim>& destination = Point<dim>() );
                  
    /// creates a refined mesh
    void BuildRefinedMesh( VSet<dim>& vset,
                           double64 length,
                           double64 dx_min,
                           double64 dx_max,
                           double64 width_of_transition_zone,
                           MeshDensity* density,
                           const std::vector<double64>& splitnode_coordinates = std::vector<double64>(),
                           const Point<dim>& origin = Point<dim>(),
                           const Point<dim>& destination = Point<dim>() );

    /// creates a custom mesh with vector( 1D node coordinates )
    void BuildCustomMesh( VSet<dim>& vset,
                          const std::vector<double64>& node_coordinates,
                          const std::vector<double64>& splitnode_coordinates = std::vector<double64>(),
                          const Point<dim>& origin = Point<dim>(),
                          const Point<dim>& destination = Point<dim>() );

    /// creates a custom mesh with vector( 1D, 2D or 3D node coordinates )
    void BuildCustomMesh( VSet<dim>& vset,
                          const std::vector<Point<dim> >& node_coordinates,
                          const std::vector<Point<dim> >& splitnode_coordinates = std::vector<Point<dim> >(),
                          const Point<dim>& origin = Point<dim>(),
                          const Point<dim>& destination = Point<dim>() );

    /// Method which modifies already existing mesh

    void AssignCornerPoints( VSet<dim>& vset,
                             const Point<dim>& origin,
                             const Point<dim>& destination );

    void InsertSplitNodes( VSet<dim>& vset,
                           const std::vector<double64>& splitnodes );

    void InsertSplitNodes( VSet<dim>& vset,
                           std::set<Point<dim> >& splitnodes );

private:

    /// Methods which include all extra functionality of mesh modificiation

    void CompleteMesh( VSet<dim>& vset,
                       const std::vector<double64>& splitnodes,
                       const Point<dim>& origin = Point<dim>(),
                       const Point<dim>& destination = Point<dim>() );

    void CompleteMesh( VSet<dim>& vset,
                       const std::vector<Point<dim> >& splitnodes,
                       const Point<dim>& origin = Point<dim>(),
                       const Point<dim>& destination = Point<dim>() );

    /// Building blocks
    
    /// creates an uniform mesh
    void UniformMesh( VSet<dim>& vset,
                      double64 length,
                      size_t n_elements );

    /// creates a refined mesh based on errf
    void RefinedMesh( VSet<dim>& vset,
                      double64 length,
                      double64 dx_min,
                      double64 dx_max,
                      double64 width_of_transition_zone,
                      MeshDensity* density );

    /// creates a custom mesh with vector( 1D node coordinates )
    void CustomMesh( VSet<dim>& vset,
                     const std::vector<double64>& node_coordinates );

    /// creates a custom mesh with vector( 1D, 2D or 3D node coordinates )
    void CustomMesh( VSet<dim>& vset,
                     const std::vector<Point<dim> >& node_coordinates );


    void EstablishConnectivity( VSet<dim>& vset );

    void EstablishConnectivity( VSet<dim>& vset,
                                std::set<Point<dim> >& splitnodes );

    VSet<dim>*          vset_;
    std::set<double64>  splitnodes_;

};

} // end namespace csmp

#endif
