#ifndef LINEAR_TETRAHEDRON_H
#define LINEAR_TETRAHEDRON_H

#include "FiniteElement.h"

namespace csmp {

class LinearTetrahedron : public FiniteElement {
  public:

    LinearTetrahedron();
    virtual ~LinearTetrahedron();

    virtual double       Volume();
    virtual void           CounterClockwiseNodes( std::vector<size_t>& ids ) const;
    virtual size_t         CornerNodes() const  { return 4U; }
    virtual void           CornerNodes( std::vector<size_t>& ids ) const;
    virtual void           NodesOfSegment( size_t segm_id, std::vector<size_t>& snids ) const;
    virtual void           NodesOfFace( size_t face_id, std::vector<size_t>& fnids ) const;
    virtual CSMP_FEM_TYPE  ElementTypeOfFace( size_t face ) const;
    virtual CSMP_FEM_TYPE  ElementTypeOfSegment( size_t /* segment */ ) const { return LINEAR_BAR; };
  
    virtual void           UnitNormalToFace( size_t face, std::vector<double>& unrml ) const;

    virtual   void         N( std::vector<double>& N, const std::vector<double>& xyz );
    virtual   void         dN( DenseMatrix<DM_MIN>& M );
    virtual   void         IntegralNN( DenseMatrix<DM_MIN>& M );

  private:
    size_t                 n( size_t i, size_t a );
    void                   UpdateFor();
};

/**

@class LinearTetrahedron  LinearTetrahedron "finite_elements/LinearTetrahedron.h"
@date 1998
@author S.K. Matthaei
@author Stephen G. Roberts */


} // end namespace csmp

#endif


