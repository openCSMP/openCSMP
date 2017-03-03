#ifndef BOUNDARY_CONNECTOR_H
#define BOUNDARY_CONNECTOR_H

#include <map>
#include <vector>

#include "Box.h"

namespace csmp {

template<size_t> class Element;
template<size_t> class Face;
template<size_t> class InterFace;
template<size_t> class Node;
template<size_t> class Model;
template<size_t> class Region;
template<size_t> class Boundary;

/// Establishes connectivity for csmp::SplitBoundary interfaces
template<size_t dim>
class BoundaryConnector {
  public:
    enum ParentOrientation{ OUTER = 1, INNER = -1 };

  public:
    static size_t FaceParentsFromElementEquivalent( Face<dim>* face, Element<dim>* const element, std::map<size_t,int>& elementUnitNormalOrientation );
    static bool   OnOutside( InterFace<dim>&, Face<dim>&, Node<dim>& );
    static bool   OnOutside( InterFace<dim>&, InterFace<dim>&, Node<dim>& );
    static bool   OnOutside( InterFace<dim>&, Element<dim>&, Node<dim>& );
    static bool   BoundaryFaceNodes( const Element<dim>& element, BOX_BOUNDARY boxBoundary, std::vector<Node<dim>*>& boundaryFaceNodes, size_t& bFace );
    static void   OuterAndInnerParent( Element<dim>* const element, std::pair<Element<dim>*,Element<dim>*>& parents, std::map<size_t,int>& elementUnitNormalOrientation );
    static void   GlobalOrientation( const Region<dim>& lowDimRegion, std::map<size_t,int>& elementUnitNormalOrientation );
    static bool   BoundaryNode( Node<dim>&, Model<dim>& );
    static bool   SplitBoundaryNode( Node<dim>&, Model<dim>& );
    static void   RemoveElementNeighborConnectivity( Element<dim>&, Element<dim>& );

  private:
    BoundaryConnector() = delete;
    ~BoundaryConnector() = delete;
    BoundaryConnector( const BoundaryConnector& ) = delete;

    static size_t AddPotentialFaceParentElement( Element<dim>*, std::map<size_t,std::pair<size_t,Element<dim>*> >& );
    
    static size_t EligibleFaceParentElements( size_t faceNodes,
                                              std::pair<Element<dim>*,Element<dim>*>& ,
                                              const std::map<size_t,std::pair<size_t,Element<dim>*> >& );
  };


/**
@class BoundaryConnector BoundaryConnector "main_library/BoundaryConnector.h"

@author P. Lang
@author Alina Yapparova
@date 2011-2012

A class dispatching methods to facilitate the creation of csmp::SplitBoundary, 
focusing on connectivity aspects. This helps in reducing size and coupling of csmp::Face and csmp::Boundary.
*/

} // csmp

#endif 
