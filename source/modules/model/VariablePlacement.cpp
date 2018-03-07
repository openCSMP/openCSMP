#include "VariablePlacement.h"
#include "Element.h"
#include "Node.h"
#include "Exception.h"


namespace csmp {

  template<>
  void calculateN<1u>(const Element<1u>& e, const Point<1u>& p, double64* coeff)
  {
    auto fe = e.FE();
    switch (e.FV()->Geometry()) {
      case LINE:
        fe->Nr( p[0], coeff );
        break;

      default:
        throw csmp::Exception(ERROR,
                          "calculateN",
                          "Element dimension must be 1, 2, or 3");
    }
  }

  template<>
  void calculateN<2u>(const Element<2u>& e, const Point<2u>& p, double64* coeff)
  {
    auto fe = e.FE();
    switch (e.FV()->Geometry()) {
      case LINE:
        fe->Nr( p[0], coeff );
        break;

      case SURFACE:
        fe->Nrs( p[0], p[1], coeff );
        break;

      default:
        throw csmp::Exception(ERROR,
                          "calculateN",
                          "Element dimension must be 1, 2, or 3");
    }
  }

  template<>
  void calculateN<3u>(const Element<3u>& e, const Point<3u>& p, double64* coeff)
  {
    auto fe = e.FE();
    switch (e.FV()->Geometry()) {
      case VOLUME:
        fe->Nrst( p[0], p[1], p[2], coeff );
        break;

      case SURFACE:
        fe->Nrs( p[0], p[1], coeff );
        break;

      case LINE:
        fe->Nr( p[0], coeff );
        break;

      default:
        throw csmp::Exception(ERROR,
                          "calculateN",
                          "Element dimension must be 1, 2, or 3");
    }
  }

  template<>
  void
  calculateDN(const Element<1u>& e, const Point<1u>& p, std::vector<double64>* DN)
  {
    auto fe = e.FE();
    switch (e.FV()->Geometry()) {
      case LINE:
        fe->dNr(p[0], DN[0]);
        fe->Jacobian( DN[0] );
        break;

      default:
        throw csmp::Exception(ERROR,
                          "calculateDN",
                          "Element dimension must be 1, 2, or 3");
    }
  }


  template<>
  void
  calculateDN(const Element<2u>& e, const Point<2u>& p, std::vector<double64>* DN)
  {
    auto fe = e.FE();
    switch (e.FV()->Geometry()) {
      case LINE:
        fe->dNr(p[0], DN[0]);
        fe->Jacobian( DN[0] );
        break;

      case SURFACE:
        fe->dNr(p[0], p[1], DN[0]);
        fe->dNs(p[0], p[1], DN[1]);
        fe->Jacobian( DN[0], DN[1] );
        break;

      default:
        throw csmp::Exception(ERROR,
                          "calculateDN",
                          "Element dimension must be 1, 2, or 3");
    }
  }


  template<>
  void
  calculateDN(const Element<3u>& e, const Point<3u>& p, std::vector<double64>* DN)
  {
    auto fe = e.FE();
    switch (e.FV()->Geometry()) {
      case VOLUME:
        fe->dNr(p[0], p[1], p[2], DN[0]);
        fe->dNs(p[0], p[1], p[2], DN[1]);
        fe->dNt(p[0], p[1], p[2], DN[2]);
        fe->Jacobian( DN[0], DN[1], DN[2] );
        break;
      case SURFACE:
        fe->dNs(p[0], p[1], DN[1]);
        fe->Jacobian( DN[0], DN[1] );
        break;
      case LINE:
        fe->dNr(p[0], DN[0]);
        fe->Jacobian( DN[0] );
        break;
      default:
        throw csmp::Exception(ERROR,
                          "calculateDN",
                          "Element dimension must be 1, 2, or 3");
    }
  }

  template<size_t dim>
  Point<dim>
  directedAreaOfFacet(const Element<dim>& e, size_t iFacet)
  {
    auto fv = e.FV();
    switch (fv->Geometry()) {
      case LINE:
      {
        Point<dim> normal;
        const size_t iNrNodes(e.Nodes());
        for (size_t iNode = 0U; iNode < iNrNodes; ++iNode) {
          const Point<dim> n(e.N(iNode)->Coordinate());
          auto weights = fv->FacetNormalTransformationNodeWeights(iFacet, iNode);
          normal += weights.first * n;
        }
        return normal;
      }

      case SURFACE:
      {
        Point<dim> tangent;
        Point<dim> bitangent;
        const size_t iNrNodes(e.Nodes());
        for (size_t iNode = 0U; iNode < iNrNodes; ++iNode) {
          const Point<dim> n(e.N(iNode)->Coordinate());
          auto weights = fv->FacetNormalTransformationNodeWeights(iFacet, iNode);
          tangent += weights.first * n;
          bitangent += weights.second * n;
        }
        double64 length = exteriorProductLength(tangent, bitangent);
        tangent.NormalizeLengthTo(1.0);
        Point<dim> normal = bitangent - dotProduct(tangent,bitangent) * tangent;
        normal.NormalizeLengthTo(length);
        return normal;
      }

      case VOLUME:
      {
        Point<dim> v0(0.0);
        Point<dim> v1(0.0);
        const size_t iNrNodes(e.Nodes());
        for (size_t iNode = 0; iNode < iNrNodes; ++iNode) {
          auto xform_weights = fv->FacetNormalTransformationNodeWeights(iFacet, iNode);
          const Point<dim> n(e.N(iNode)->Coordinate());
          v0 += xform_weights.first * n;
          v1 += xform_weights.second * n;
        }
        return crossProduct(v1, v0);
      }

      default:
        throw csmp::Exception(ERROR,
                          "directedAreaOfFacet",
                          "Element dimension must be 1, 2, or 3");
    }
  }

  template Point<1u> directedAreaOfFacet(const Element<1u>& e, size_t iFacet);
  template Point<2u> directedAreaOfFacet(const Element<2u>& e, size_t iFacet);
  template Point<3u> directedAreaOfFacet(const Element<3u>& e, size_t iFacet);

}
