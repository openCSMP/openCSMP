// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LINEAR_TETRAHEDRON_H
#define LINEAR_TETRAHEDRON_H

#include "FiniteElement.h"

namespace csmp {

class LinearTetrahedron final : public FiniteElement {
  public:

    LinearTetrahedron();

    virtual double  Volume();
    virtual void    CounterClockwiseNodes( std::vector<uint32_t>& ids ) const;
    virtual uint32_t  CornerNodes() const  { return 4U; }
    virtual void    CornerNodes( std::vector<uint32_t>& ids ) const;
    virtual void    NodesOfSegment( uint32_t segm_id, std::vector<uint32_t>& snids ) const;
    
    virtual std::vector<uint32_t>  NodesOfFace( uint32_t face_id ) const;
    virtual std::vector<uint32_t>  CornerNodesOfFace( uint32_t face_id ) const;
    virtual std::vector<uint32_t>  NodesConnectedTo( uint32_t node_id ) const;
    
    virtual CSMP_FEM_TYPE  ElementTypeOfFace( uint32_t face ) const;
    virtual CSMP_FEM_TYPE  ElementTypeOfSegment( uint32_t /* segment */ ) const { return LINEAR_BAR; };
  
    virtual void UnitNormalToFace( uint32_t face, std::vector<double>& unrml ) const;
    
    virtual void N( std::vector<double>& N, const std::vector<double>& xyz );
    virtual void dN( DenseMatrix<DM_MIN>& M );
    virtual void IntegralNN( DenseMatrix<DM_MIN>& M );

  private:
    uint32_t n( uint32_t i, uint32_t a );
    void   UpdateFor();
};

/**

@class LinearTetrahedron  LinearTetrahedron "finite_elements/LinearTetrahedron.h"
@date 1998
@author S.K. Matthaei
@author Stephen G. Roberts */


} // end namespace csmp

#endif


