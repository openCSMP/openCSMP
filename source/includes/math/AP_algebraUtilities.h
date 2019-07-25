#ifndef AP_ALGEBRA_UTILITIES_H
#define AP_ALGEBRA_UTILITIES_H

#include "CSMP_definitions.h"
#include "FiniteVolumeStencil.h"

namespace csmp {
  
  /** Calculates the dot product between two vectors.
   
   @section arguments Input Arguments
   
   The two vectors to be computed.
   
   @return The result value of the dot product.
   
   @section implementation Implementation
   
   The std::inner_product between the two vectors is returned.
   
   @section application Application
   
   This may be used to calculate the projection of one vector onto another, as
   well as one of the step in the calculation of the vector norm.
   
   */
  inline double64 dotProduct(const std::vector<double64>& vec1, const std::vector<double64>& vec2)
  {
    assert(vec1.size() == vec2.size());
    return std::inner_product(vec1.begin(), vec1.end(), vec2.begin(), 0.);
  }
  
  
  
  
  std::vector<double64> crossProduct( const std::vector<double64>&, const std::vector<double64>&);
  
  
  std::vector<double64> crossProduct3by3(/*1st std::vector*/const double64&, const double64&, const double64&,
                                         /*2nd std::vector*/const double64&, const double64&, const double64&);
  
  
  template<size_t dim>
  void areaCenterOfMass(const std::vector<Point<dim> >&, Point<dim>& );
  
  
  bool normalOfPolygon(const std::vector< Point<3U> >& vecPolygon, const Point<3U>& vecNormalAt, Point<3U>& vecNormal);
  
  template<size_t dim>
  bool areaOfPolygon(const std::vector< Point<dim> >& vecPolygon, const size_t& iNrOfFacetPoints, double64 & fArea);
  
  template<size_t dim>
  void localSurfaceNormal(const std::vector< Point< dim> >& vecPoints, const Point<dim>& vecNormalAt, const size_t& iLevelOfRefinement, Point<dim>& vecNormal);
  
  template<size_t dim>
  bool intersection(/*1st line*/const Point< dim>&, const Point< dim>&,
                    /*2nd line*/const Point< dim>&, const Point< dim>&,
                    Point< dim>&);
  
  
  double64 distanceBetweenPoints(const std::vector<double64>&, const std::vector<double64>&);
  
  
  /** Calculates the norm of a given vector.
   
   @section arguments Input Arguments
   
   The vector to be normed.
   
   @return The normed vector. However, the normed vector is returned in the same input
   variable, which is passed by reference.
   
   @section implementation Implementation
   
   Calculates the square root of the dot product of the vector by itself.
   
   @section application Application
   
   The main objective is to calculate the length / norm of a given vector in
   order to evaluate its magnitude.
   */
  inline double64 normOfVector(const std::vector<double64>& vector1)
  {
    return std::sqrt(dotProduct(vector1, vector1));
  }
  
  
  /**
   
   
   Normalizes the vector using a Euclidean norm, thus dividing each term of
   the vector by the total length of the vector.
   
   @section arguments Input Arguments
   
   @param vectorToNormalize The vector to be normalized.
   
   @section implementation Implementation
   
   At first the norm of the vector is calculated, then each component of the
   vector is divided by this quantity.
   
   @section application Application
   
   This method is used to make the input vector a unitary vector, thus to make
   its length equal to one, preserving its original direction and orientation.
   
   @note SKM - replaced - std::bind2nd(std::multiplies<double64>(), fNorm - by lambda function
   
   */
  inline void euclideanNormalize( std::vector<double64>& vectorToNormalize )
  {
    const double64 fNorm(1.0 / std::sqrt(dotProduct(vectorToNormalize, vectorToNormalize)));
    if(fNorm != 0)
      std::transform( vectorToNormalize.begin(), vectorToNormalize.end(), vectorToNormalize.begin(),
                      [fNorm](const double64 val){ return val * fNorm; } );
  }
  
  /** Binary predicate used to compare elements of a set to be order by magnitude.
   
   original set: -10, -4, 1, 5, 9
   
   This set: -10, 9, 5, -4, 1
   
   @section implementation Implementation
   This is just an operator that takes two input parameters: x and y, and returns
   true if the absolute value of x is greater than the absolute value of y.
   
   @section application Application
   An example of usage is the following:
   
   @code
   std::set< greater_fabs > vals;
   vals.insert(9);
   vals.insert(-10);
   vals.insert(11);
   @endcode
   
   The order of the elements will be: 11, -10, 9
   This is used in TensorVariable.cpp to order the eigenvalues.
   
   */
  //binary predicate for the larger magnitude of vectors
  template <class fT>
  struct greater_fabs : std::binary_function<fT,double64,bool>
  {
    bool operator()(const fT& x, const fT& y) const {return static_cast<bool>(std::fabs(x) > std::fabs(y));}
  };
  
  
  /**
   
   Calculates the center of equal masses at the vertices of a planar, closed,
   convex surface defined by vecPoints.
   
   @section arguments Input Arguments
   
   @param vecPoints The points representing the vertices.
   
   @param vecCentroid the centroid of the point cloud.
   
   @section implementation Implementation
   
   At first all the points are summed, and then each of these summed components
   is divided by the total number of vertices.
   
   @section application Application
   
   The main idea is to calculate the center of mass delimited by the given
   vertices, supposing that the body is of constant density.
   */
  template<size_t dim>
  inline void vertexCenterOfMass( const std::vector<Point<dim> >& vecPoints, Point<dim>& vecCentroid )
  {
    for(size_t i = 0U; i < vecPoints.size(); i++)
      vecCentroid += vecPoints[i];
    
    vecCentroid /= vecPoints.size();
  }
  
  
  
  
  /**
                  
   Calculates the center of equal masses at the vertices of a planar, closed,
   convex surface defined by three points.
   
   @section arguments Input Arguments
   
   @param vecCentroid The points representing the vertices.
   
   @section implementation Implementation
   
   At first all the points are summed, and then each of these summed components
   is divided by the total number of vertices.
   
   @section application Application
   
   The main idea is to calculate the center of mass delimited by the three given
   vertices, supposing that the body is of constant density.
   */
  template<size_t dim>
  inline void vertexCenterOfMass3Vertices( const Point<dim>& pt1, const Point<dim>& pt2,
                                          const Point<dim>& pt3, Point<dim>& vecCentroid )
  {
    vecCentroid = (pt1 + pt2 + pt3) / 3.;
  }
  
} // end csmp

#endif
