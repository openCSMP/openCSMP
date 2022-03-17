#ifndef CSMP_LINE_ELEMENT_MESHER_H
#define CSMP_LINE_ELEMENT_MESHER_H

#include "VSet.h"
#include "Point.h"

namespace csmp {

/// mesh refinement helper base class
class MeshDensity {
public:
    MeshDensity( );
    virtual ~MeshDensity();
    virtual double operator()( double x );
};

/// helper subclass
class LinDensity: virtual public MeshDensity {
public:
    LinDensity( double a = 1.0, double b = 0.0 );
    virtual ~LinDensity();
    virtual double operator()( double x );
    void SetA( double );
    void SetB( double );
private:
    /// y = a_ * x + b_
    double a_;
    double b_;
};

/// exponential element size variation
class ExpDensity: virtual public MeshDensity {
public:
    ExpDensity( double a = 1.0, double b = 1.0  );
    virtual ~ExpDensity();
    virtual double operator()(double x );
    void SetA( double );
    void SetB( double );
private:
    /// y = b_ * exp( a_ * x )
    double a_;
    double b_;
};

/// error function based mesh size variation
class ErfDensity: virtual public MeshDensity {
public:
    ErfDensity( double a = 3.0, double b = 1.0 );
    ~ErfDensity();
    virtual double operator()(double x );
    void SetA( double );
    void SetB( double );
private:
    /// y = b_ * erf( a * x )
    double a_;
    double b_;
};


/**
     Line element mesher.
     
     @author Stephan Matthai
     @author extended by Roman Manasipov
     @date 2002, 2012
*/
template<uint32_t dim>
class LineElementMesher {
  public:
    LineElementMesher( );

    /// creates an uniform mesh
    void BuildUniformMesh( VSet<dim>& vset,
                           double length,
                           size_t   n_elements,
                           const std::vector<double>& splitnode_coordinates = std::vector<double>(),
                           const Point<dim>& origin = Point<dim>(),
                           const Point<dim>& destination = Point<dim>() );
                  
    /// creates a refined mesh
    void BuildRefinedMesh( VSet<dim>& vset,
                           double length,
                           double dx_min,
                           double dx_max,
                           double width_of_transition_zone,
                           MeshDensity* density,
                           const std::vector<double>& splitnode_coordinates = std::vector<double>(),
                           const Point<dim>& origin = Point<dim>(),
                           const Point<dim>& destination = Point<dim>() );

    /// creates a custom mesh with vector( 1D node coordinates )
    void BuildCustomMesh( VSet<dim>& vset,
                          const std::vector<double>& node_coordinates,
                          const std::vector<double>& splitnode_coordinates = std::vector<double>(),
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
                           const std::vector<double>& splitnodes );

    void InsertSplitNodes( VSet<dim>& vset,
                           std::set<Point<dim> >& splitnodes );

private:

    /// Methods which include all extra functionality of mesh modificiation

    void CompleteMesh( VSet<dim>& vset,
                       const std::vector<double>& splitnodes,
                       const Point<dim>& origin = Point<dim>(),
                       const Point<dim>& destination = Point<dim>() );

    void CompleteMesh( VSet<dim>& vset,
                       const std::vector<Point<dim> >& splitnodes,
                       const Point<dim>& origin = Point<dim>(),
                       const Point<dim>& destination = Point<dim>() );

    /// Building blocks
    
    /// creates an uniform mesh
    void UniformMesh( VSet<dim>& vset,
                      double length,
                      size_t n_elements );

    /// creates a refined mesh based on errf
    void RefinedMesh( VSet<dim>& vset,
                      double length,
                      double dx_min,
                      double dx_max,
                      double width_of_transition_zone,
                      MeshDensity* density );

    /// creates a custom mesh with vector( 1D node coordinates )
    void CustomMesh( VSet<dim>& vset,
                     const std::vector<double>& node_coordinates );

    /// creates a custom mesh with vector( 1D, 2D or 3D node coordinates )
    void CustomMesh( VSet<dim>& vset,
                     const std::vector<Point<dim> >& node_coordinates );


    void EstablishConnectivity( VSet<dim>& vset );

    void EstablishConnectivity( VSet<dim>& vset,
                                std::set<Point<dim> >& splitnodes );

    VSet<dim>*          vset_;
    std::set<double>  splitnodes_;

};

} // end namespace csmp

#endif
